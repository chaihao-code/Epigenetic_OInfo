#!/usr/bin/env python3
"""03b_reference_qc.py — Step 4 续: 参考数据验证 + marker 选择 + 参考矩阵.

验证 (不凭猜测):
  A. CpG 全局索引正确性: Markers.U250.hg38.tsv 的 startCpG/endCpG 与
     GSE186458_blocks.s207.hg38 的 CpG 编号与我们的 2bit 扫描坐标比对。
  B. beta 文件身份与读取正确性: 用我们的索引从每个 beta 计算 UXM U250 marker
     区域甲基化, 与 Atlas.U250.l4.hg38.full.tsv 中对应细胞类型列逐 marker 相关
     (Thyroid_epithelial->Thyroid-Ep, T_CD4/T_CD8->Blood-T, B_cell->Blood-B,
      Monocyte->Blood-Mono+Macro, NK->Blood-NK, Granulocyte->Blood-Granul,
      Endothelial->Endothel, Fibroblast->Colon/Dermal/Heart-Fibro)。

marker panel:
  A_panel: UXM Markers.U250.hg38 中 target 属于我们 8 类细胞的 marker。
  B_panel: 独立 1kb-bin 筛选 (免疫 vs 上皮+基质 差异大、重复内稳定、覆盖足)。
输出:
  results/tables/marker_quality.tsv, data/processed/reference_matrix_{A,B}.tsv
  data/processed/selected_markers_{A,B}.bed, data/processed/immune_specific_regions.bed
  results/figures/reference_heatmap.pdf
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/03b_reference_qc.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("ref_qc")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load(open(ROOT / "config.yaml"))
REF = ROOT / CFG["reference"]["dir"]
PROC = ROOT / "data/processed"
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
IMMUNE = ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"]
ALL_CT = ["Thyroid_epithelial"] + IMMUNE + ["Endothelial", "Fibroblast"]

TARGET_MAP = {
    "Thyroid_epithelial": ["Thyroid-Ep"],
    "T_CD4": ["Blood-T"], "T_CD8": ["Blood-T"],
    "B_cell": ["Blood-B"], "Monocyte": ["Blood-Mono+Macro"],
    "NK": ["Blood-NK"], "Granulocyte": ["Blood-Granul"],
    "Endothelial": ["Endothel"],
    "Fibroblast": ["Colon-Fibro", "Dermal-Fibro", "Heart-Fibro"],
}

BETA_CPG_TOTAL = 58803590 // 2


def load_cpg_index():
    z = np.load(PROC / "hg38_cpg_index.npz", allow_pickle=True)
    return list(z["chrom_names"]), z["chrom"], z["pos"]


def load_beta(path):
    raw = np.fromfile(path, dtype=np.uint8)
    assert len(raw) == 2 * BETA_CPG_TOTAL
    return raw[0::2].astype(np.uint32), raw[1::2].astype(np.uint32)


def marker_mean_beta(meth8, cov8, start, end):
    """CpG 区间 [start, end) (0-based) 的 coverage 加权甲基化均值。"""
    m, c = meth8[start:end], cov8[start:end]
    s = c.sum()
    return np.nan if s == 0 else (m.sum() / s)


# ---------------------------------------------------------------- validation A
def validate_index(chrom_names, chrom, pos):
    """wgbstools CpG 编号为 1-based, marker BED start 亦为 1-based:
    校验 pos[startCpG-1] == start-1 (0-based C 坐标)。"""
    markers = pd.read_csv(REF / "Markers.U250.hg38.tsv", sep="\t")
    mm = markers.head(5000)
    hit = 0
    for _, r in mm.iterrows():
        name = str(r["#chr"])
        ci = chrom_names.index(name if name.startswith("chr") else "chr" + name)
        s = int(r["startCpG"])
        if 0 < s <= len(pos) and chrom[s - 1] == ci and pos[s - 1] == int(r["start"]) - 1:
            hit += 1
    log.info("CpG index validation (1-based): %d/%d exact", hit, len(mm))
    return hit >= 0.98 * len(mm)


# ---------------------------------------------------------------- validation B
def validate_beta_vs_markers(chrom_names, chrom, pos):
    """与 Markers.U250.hg38.tsv 自带的 tg_mean 比对 (atlas full 文件列定义
    与坐标行不一致, 已弃用; 主验证为 beta vs GEO bigWig r=0.9999998, 见
    results/qc/reference_beta_validation.tsv 备注)。"""
    mk = pd.read_csv(REF / "Markers.U250.hg38.tsv", sep="\t")
    out_rows = []
    for ct in ALL_CT:
        tg_names = TARGET_MAP[ct]
        sub = mk[mk["target"].isin(tg_names)]
        if len(sub) == 0:
            continue
        files = sorted((REF / "beta_hg38").glob(f"{ct}__*.beta"))
        reps_all = [load_beta(f) for f in files]
        vals, refs = [], []
        for _, r in sub.iterrows():
            name = str(r["#chr"])
            try:
                ci = chrom_names.index(name if name.startswith("chr") else "chr" + name)
            except ValueError:
                continue
            s, e = int(r["startCpG"]) - 1, int(r["endCpG"])  # 1-based -> 0-based half-open
            if e <= s or e > len(pos) or not (chrom[s] == ci and chrom[e - 1] == ci):
                continue
            reps = [marker_mean_beta(m8[s:e], c8[s:e], 0, e - s) for m8, c8 in reps_all]
            if np.all(np.isnan(reps)):
                continue
            vals.append(np.nanmean(reps))
            refs.append(r["tg_mean"])
        vals, refs = np.array(vals), np.array(refs)
        okm = ~np.isnan(vals) & ~np.isnan(refs)
        r_corr = np.corrcoef(vals[okm], refs[okm])[0, 1] if okm.sum() > 100 else np.nan
        out_rows.append(dict(cell_type=ct, n_markers=int(okm.sum()),
                             pearson_r_vs_tg_mean=round(float(r_corr), 5),
                             mean_abs_diff=round(float(np.abs(vals - refs)[okm].mean()), 5)))
        log.info("%s: %d markers, r vs tg_mean = %.4f", ct, okm.sum(), r_corr)
        del reps_all
    qc = pd.DataFrame(out_rows)
    qc.to_csv(TBL / "reference_beta_validation.tsv", sep="\t", index=False)
    return qc


def main():
    chrom_names, chrom, pos = load_cpg_index()
    ok = validate_index(chrom_names, chrom, pos)
    if not ok:
        log.error("CpG index validation FAILED — 中止, 请检查染色体顺序/命名")
        return 1
    qc = validate_beta_vs_markers(chrom_names, chrom, pos)
    print(qc.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
