#!/usr/bin/env python3
"""03_build_reference.py — Step 4: 参考矩阵构建与 marker 质控.

流程:
  1. 用 hg38.2bit 扫描全部 CpG, 构建 wgbstools 风格的 CpG 全局索引
     (染色体顺序 chr1..22,X,Y [+chrM 候选], 以 beta 文件字节数反推的总数校验)。
  2. 读全部 .hg38.beta (#meth,#coverage uint8 对) -> per-CpG beta + coverage。
  3. 用 GSE186458_blocks.s207.hg38.bed.gz(若含逐样本甲基化列)验证索引与文件身份。
  4. 聚合到 1kb 常染色体 bin: 每参考样本的 coverage 加权甲基化均值
     -> data/processed/ref_bin_methylation.npz (bins x 27)。
  5. QC: 生物学重复一致性 (基因组范围 beta 相关), 细胞类型可辨识性热图。
  6. marker panel A: UXM Markers.U250.hg38 (atlas 现成 marker),
     panel B: 独立 1kb-bin 筛选 (免疫 vs 上皮对比大、重复内稳定)。
输出 results/tables/{reference_sample_metadata.tsv, reference_bins_qc.tsv,
marker_quality.tsv}, data/processed/{ref_bin_methylation.npz, selected_markers.bed},
results/figures/reference_heatmap.pdf
"""

# Portable repository paths; no analysis is executed by this helper.
from pathlib import Path as _ReleasePath
_RELEASE_ROOT = next(p for p in _ReleasePath(__file__).resolve().parents
                     if (p / ".epigenetic_oinfo_root").is_file())
def _release_path(value):
    return value.replace("@ROOT@", str(_RELEASE_ROOT))
def _release_config(value):
    if isinstance(value, dict):
        return {k: _release_config(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_config(v) for v in value]
    return _release_path(value) if isinstance(value, str) else value
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import py2bit
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/03_build_reference.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("build_reference")

ROOT = Path(__file__).resolve().parents[1]
CFG = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
REF = ROOT / CFG["reference"]["dir"]
PROC = ROOT / "data/processed"
PROC.mkdir(parents=True, exist_ok=True)

MARKS = CFG["marks"]
IMMUNE_TYPES = ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"]
STROMAL_TYPES = ["Endothelial", "Fibroblast"]
ALL_CELLTYPES = ["Thyroid_epithelial"] + IMMUNE_TYPES + STROMAL_TYPES

BETA_CPG_TOTAL = 58803590 // 2  # 29,401,795 (下载文件实测)

_st = Path(CFG["reference"].get("staging_dir", ""))
STAGE = _st if (_st / "hg38.2bit").exists() and (_st / "beta_hg38").exists() else REF


# ---------------------------------------------------------------- CpG index
def build_cpg_index():
    """扫描 hg38.2bit 建 CpG 索引; 总数必须与 BETA_CPG_TOTAL 一致。"""
    cache = PROC / "hg38_cpg_index.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return z["chrom"], z["pos"]

    tb = py2bit.open(str(STAGE / "hg38.2bit"))
    chroms_all = tb.chroms()
    order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    order = [c for c in order if c in chroms_all]

    chrom_out, pos_out = [], []
    for c in order:
        seq = tb.sequence(c).upper()
        arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        is_c = arr == ord("C")
        is_g_next = np.empty_like(is_c)
        is_g_next[:-1] = arr[1:] == ord("G")
        is_g_next[-1] = False
        pos = np.where(is_c & is_g_next)[0]
        chrom_out.append(np.full(len(pos), order.index(c), dtype=np.int8))
        pos_out.append(pos.astype(np.int64))
        log.info("%s: %d CpGs", c, len(pos))
    tb.close()
    chrom = np.concatenate(chrom_out)
    pos = np.concatenate(pos_out)
    log.info("total CpGs: %d (expect %d)", len(pos), BETA_CPG_TOTAL)
    if len(pos) != BETA_CPG_TOTAL:
        log.warning("CpG count mismatch! 尝试仅常染色体+X/Y 以外的组合请人工检查")
    np.savez_compressed(cache, chrom=chrom, pos=pos,
                        chrom_names=np.array(order))
    return chrom, pos


def load_beta(path):
    """beta 文件 -> (meth uint8, cov uint8) 数组, 长度 = BETA_CPG_TOTAL。"""
    raw = np.fromfile(path, dtype=np.uint8)
    assert len(raw) == 2 * BETA_CPG_TOTAL, f"{path}: {len(raw)}"
    return raw[0::2].astype(np.uint16), raw[1::2].astype(np.uint16)


# ---------------------------------------------------------------- main
def main():
    chrom_idx, cpg_pos = build_cpg_index()
    znames = np.load(PROC / "hg38_cpg_index.npz", allow_pickle=True)["chrom_names"]
    chrom_names = list(znames)

    # --- 1kb 常染色体 bin 边界 ---
    bins = pd.read_csv(CFG["data"]["bins_1kb"], sep="\t", header=None,
                       names=["chrom", "start", "end"])
    bins = bins[bins["chrom"].str.match(r"^chr\d+$")].reset_index(drop=True)
    n_bins = len(bins)
    log.info("autosomal 1kb bins: %d", n_bins)
    bin_chrom_code = bins["chrom"].map({c: i for i, c in enumerate(chrom_names)}).values
    bin_start = bins["start"].values
    bin_end = bins["end"].values

    # --- 加载全部参考样本 beta -> 1kb bin 加权甲基化 ---
    beta_files = sorted((STAGE / "beta_hg38").glob("*.hg38.beta"))
    log.info("loading %d beta files", len(beta_files))
    meth_mat = np.zeros((n_bins, len(beta_files)), dtype=np.float64)
    cov_mat = np.zeros((n_bins, len(beta_files)), dtype=np.float32)
    sample_meta = []
    for j, bf in enumerate(beta_files):
        ct, gsm = bf.name.split("__")[0], bf.name.split("__")[1].split(".")[0]
        meth, cov = load_beta(bf)
        beta = np.where(cov > 0, meth / np.maximum(cov, 1), np.nan)
        w = cov.astype(np.float32)
        # 逐染色体聚合到 bin (用 bincount; np.add.at 对列视图不写回, 已踩坑)
        for ci, cname in enumerate(chrom_names[:22]):
            sel = chrom_idx == ci
            p = cpg_pos[sel]
            v = beta[sel]
            ww = w[sel]
            bsel = bin_chrom_code == ci
            bs, be = bin_start[bsel], bin_end[bsel]
            bin_of_cpg = np.searchsorted(bs, p, side="right") - 1
            valid = (bin_of_cpg >= 0) & (p < be[-1])
            gidx = np.where(bsel)[0][bin_of_cpg[valid]]
            cov_mat[:, j] += np.bincount(gidx, weights=ww[valid], minlength=n_bins)
            wv = np.where(np.isnan(v), 0, ww * v)[valid]
            meth_mat[:, j] += np.bincount(gidx, weights=wv, minlength=n_bins)
        sample_meta.append(dict(cell_type=ct, gsm=gsm, file=bf.name))
        log.info("loaded %s (%s)", bf.name, ct)
    meth_bin = np.where(cov_mat > 0, meth_mat / np.maximum(cov_mat, 1e-9), np.nan)

    smeta = pd.DataFrame(sample_meta)
    smeta.to_csv(ROOT / "results/tables/reference_sample_metadata.tsv",
                 sep="\t", index=False)
    np.savez_compressed(PROC / "ref_bin_methylation.npz",
                        meth=meth_bin.astype(np.float16), cov=cov_mat.astype(np.float16),
                        samples=np.array(smeta["cell_type"] + "__" + smeta["gsm"]))

    # --- QC: 重复一致性 (基因组范围逐bin beta 相关) ---
    qc_rows = []
    for ct in ALL_CELLTYPES:
        cols = [j for j, c in enumerate(smeta["cell_type"]) if c == ct]
        if len(cols) < 2:
            continue
        for a in range(len(cols)):
            for b in range(a + 1, len(cols)):
                x, y = meth_bin[:, cols[a]], meth_bin[:, cols[b]]
                ok = ~np.isnan(x) & ~np.isnan(y)
                r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 1000 else np.nan
                qc_rows.append(dict(cell_type=ct,
                                    sample_a=smeta.iloc[cols[a]]["gsm"],
                                    sample_b=smeta.iloc[cols[b]]["gsm"],
                                    n_shared=int(ok.sum()), pearson_r=round(float(r), 5)))
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(ROOT / "results/tables/reference_replicate_qc.tsv", sep="\t", index=False)
    log.info("replicate r median:\n%s", qc.groupby("cell_type")["pearson_r"].median())

    # --- 细胞类型均值矩阵 (每 cell type 3 重复均值) ---
    ct_mean = np.zeros((n_bins, len(ALL_CELLTYPES)), dtype=np.float32)
    for k, ct in enumerate(ALL_CELLTYPES):
        cols = [j for j, c in enumerate(smeta["cell_type"]) if c == ct]
        ct_mean[:, k] = np.nanmean(meth_bin[:, cols], axis=1)
    np.savez_compressed(PROC / "ref_celltype_bin_mean.npz",
                        mean=ct_mean, celltypes=np.array(ALL_CELLTYPES))

    log.info("DONE step4-part1 (CpG index + bin matrix + replicate QC)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
