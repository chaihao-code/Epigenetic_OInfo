#!/usr/bin/env python3
"""04_select_markers.py — Step 4 续: marker panel 选择 + 参考矩阵 + 免疫特异区域.

Panel B (atlas UXM U250): Markers.U250.hg38.tsv 中 target 属于我们 8 类细胞的
  marker (Thyroid-Ep, Blood-T, Blood-B, Blood-Mono+Macro, Blood-NK, Blood-Granul,
  Endothel, Colon/Dermal/Heart-Fibro)。
  注意: UXM marker 原设计用于 read-level U/X/M 统计; 这里用区域 mean beta
  作为参考值, 属于近似, 已在 methods/limitations 说明 (核对定义适用性)。

Panel A (独立 1kb-bin): 从参考 bin 甲基化矩阵筛选免疫 vs 上皮+基质差异大、
  重复内稳定、覆盖充足的 bin, 方向 (免疫高/低甲基化) 均衡。

R 矩阵: R[g, k] = marker g 在 cell type k 的 mean beta (3 重复均值)。
输出:
  data/processed/reference_matrix_{A,B}.tsv   (长表: marker x cell_type)
  data/processed/selected_markers_{A,B}.bed
  data/processed/immune_specific_regions.bed  (Step 9 用)
  data/processed/thyroid_specific_regions.bed
  results/tables/marker_quality.tsv
  results/figures/reference_heatmap.pdf
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/04_select_markers.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("select_markers")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load(open(ROOT / "config.yaml"))
REF = ROOT / CFG["reference"]["dir"]
PROC = ROOT / "data/processed"
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

IMMUNE = ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"]
STROMAL = ["Endothelial", "Fibroblast"]
ALL_CT = ["Thyroid_epithelial"] + IMMUNE + STROMAL
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
    return raw[0::2].astype(np.uint32), raw[1::2].astype(np.uint32)


# ---------------------------------------------------------------- Panel B (UXM)
def build_panel_B(chrom_names, chrom, pos):
    markers = pd.read_csv(REF / "Markers.U250.hg38.tsv", sep="\t")
    keep = markers[markers["target"].isin({t for v in TARGET_MAP.values() for t in v})].copy()
    log.info("Panel B (UXM): %d markers (of %d)", len(keep), len(markers))
    # 目标细胞类型归属 (一个 marker 的 target 即其特异细胞类型)
    tgt2ct = {t: ct for ct, ts in TARGET_MAP.items() for t in ts}

    # 每 cell type 的 3 个 beta
    ct_betas = {}
    for ct in ALL_CT:
        ct_betas[ct] = [load_beta(f) for f in sorted((REF / "beta_hg38").glob(f"{ct}__*.beta"))]

    rows, bed = [], []
    for _, r in keep.iterrows():
        name = str(r["#chr"])
        ci = chrom_names.index(name if name.startswith("chr") else "chr" + name)
        s, e = int(r["startCpG"]) - 1, int(r["endCpG"])
        if e <= s or e > len(pos) or chrom[s] != ci or chrom[e - 1] != ci:
            continue
        means, sds, covs = {}, {}, {}
        for ct in ALL_CT:
            reps, rc = [], []
            for m8, c8 in ct_betas[ct]:
                m, c = m8[s:e], c8[s:e]
                cv = int(c.sum())
                rc.append(cv)
                reps.append(np.nan if cv == 0 else m.sum() / cv)
            means[ct] = float(np.nanmean(reps))
            sds[ct] = float(np.nanstd(reps))
            covs[ct] = float(np.mean(rc))
        vals = np.array([means[c] for c in ALL_CT])
        row = dict(chrom=chrom_names[ci], start=int(r["start"]), end=int(r["end"]),
                   startCpG=int(r["startCpG"]), endCpG=int(r["endCpG"]),
                   target_cell_type=tgt2ct.get(r["target"], r["target"]),
                   atlas_target=r["target"], direction=r.get("direction", ""),
                   between_contrast=float(np.nanmax(vals) - np.nanmin(vals)),
                   **{f"mean_{c}": means[c] for c in ALL_CT},
                   **{f"sd_{c}": sds[c] for c in ALL_CT},
                   **{f"cov_{c}": covs[c] for c in ALL_CT})
        rows.append(row)
        bed.append((chrom_names[ci], int(r["start"]) - 1, int(r["end"]),
                    f"B{len(bed):05d}"))
    df = pd.DataFrame(rows)
    df.insert(0, "marker_id", [f"B{i:05d}" for i in range(len(df))])
    df.to_csv(TBL / "marker_quality.tsv", sep="\t", index=False)
    pd.DataFrame(bed, columns=["chrom", "start", "end", "name"]).to_csv(
        PROC / "selected_markers_B.bed", sep="\t", index=False, header=False)
    # 长表 R 矩阵
    long = df.melt(id_vars="marker_id", value_vars=[f"mean_{c}" for c in ALL_CT],
                   var_name="cell_type", value_name="methylation")
    long["cell_type"] = long["cell_type"].str.replace("mean_", "", regex=False)
    long.to_csv(PROC / "reference_matrix_B.tsv", sep="\t", index=False)
    log.info("Panel B (UXM) R matrix: %d markers x %d cell types", df["marker_id"].nunique(), len(ALL_CT))
    return df


# ---------------------------------------------------------------- Panel A (independent) + regions
def build_panel_A():
    z = np.load(PROC / "ref_bin_methylation.npz", allow_pickle=True)
    meth, cov, samples = z["meth"].astype(np.float32), z["cov"].astype(np.float32), list(z["samples"])
    zc = np.load(PROC / "ref_celltype_bin_mean.npz", allow_pickle=True)
    ct_mean, cts = zc["mean"], list(zc["celltypes"])
    assert cts == ALL_CT, cts

    bins = pd.read_csv(CFG["data"]["bins_1kb"], sep="\t", header=None,
                       names=["chrom", "start", "end"])
    bins = bins[bins["chrom"].str.match(r"^chr\d+$")].reset_index(drop=True)

    immune_idx = [cts.index(c) for c in IMMUNE]
    other_idx = [cts.index(c) for c in ["Thyroid_epithelial"] + STROMAL]
    epi_idx = cts.index("Thyroid_epithelial")

    immune_mean = ct_mean[:, immune_idx].mean(axis=1)
    other_mean = ct_mean[:, other_idx].mean(axis=1)
    epi_mean = ct_mean[:, epi_idx]
    d = immune_mean - other_mean          # >0: 免疫高甲基化; <0: 免疫低甲基化
    d_thy = epi_mean - np.delete(ct_mean, epi_idx, axis=1).mean(axis=1)

    # 重复内稳定性: 每个 cell type 3 重复的 bin 值 sd, 取 9 类最大值
    rep_sd = np.zeros(meth.shape[0], dtype=np.float32)
    for ct in ALL_CT:
        cols = [samples.index(f"{ct}__{g}") for g in
                [p.name.split("__")[1].split(".")[0]
                 for p in sorted((REF / "beta_hg38").glob(f"{ct}__*.beta"))]]
        rep_sd = np.maximum(rep_sd, np.nanstd(meth[:, cols], axis=1))
    # 覆盖: 9 类平均覆盖度
    mean_cov = np.zeros(meth.shape[0], dtype=np.float32)
    for ct in ALL_CT:
        cols = [samples.index(s) for s in samples if s.startswith(ct + "__")]
        mean_cov += np.nanmean(cov[:, cols], axis=1)
    mean_cov /= len(ALL_CT)

    ok = (mean_cov >= 10) & (rep_sd <= 0.08) & ~np.isnan(d)
    log.info("bins passing filters: %d / %d", ok.sum(), len(ok))

    # --- 免疫特异区域 (Step 9): 强差异 ---
    strong = ok & (np.abs(d) >= 0.4)
    reg = bins.copy()
    reg["d_immune_minus_other"] = d
    reg["mean_coverage"] = mean_cov
    hyper = reg[ok & (d >= 0.4)]
    hypo = reg[ok & (d <= -0.4)]
    out = pd.concat([
        hyper.assign(type="immune_hyper"),
        hypo.assign(type="immune_hypo"),
    ])
    out.to_csv(PROC / "immune_specific_regions.bed", sep="\t", index=False)
    thy = reg[ok & (np.abs(d_thy) >= 0.4)]
    thy["d_thyroid_minus_other"] = d_thy[ok & (np.abs(d_thy) >= 0.4)]
    thy.to_csv(PROC / "thyroid_specific_regions.bed", sep="\t", index=False)
    log.info("immune-specific bins: hyper=%d hypo=%d; thyroid-specific=%d",
             len(hyper), len(hypo), len(thy))

    # --- Panel A markers: 按 |d| 排序取 top, 方向均衡 ---
    cand = np.where(ok)[0]
    order = cand[np.argsort(-np.abs(d[cand]))]
    hyper_sel = [i for i in order if d[i] > 0][:1250]
    hypo_sel = [i for i in order if d[i] < 0][:1250]
    sel = np.array(hyper_sel + hypo_sel)
    mq = bins.iloc[sel].copy()
    mq.insert(0, "marker_id", [f"A{i:05d}" for k, i in enumerate(sel)])
    for k, ct in enumerate(ALL_CT):
        mq[f"mean_{ct}"] = ct_mean[sel, k]
    mq["d_immune_minus_other"] = d[sel]
    mq["between_contrast"] = np.nanmax(ct_mean[sel], axis=1) - np.nanmin(ct_mean[sel], axis=1)
    mq["rep_sd_max"] = rep_sd[sel]
    mq["mean_coverage"] = mean_cov[sel]
    mq["direction"] = np.where(d[sel] > 0, "immune_hyper", "immune_hypo")
    mq.to_csv(PROC / "selected_markers_A.bed", sep="\t", index=False)
    long = mq.melt(id_vars="marker_id", value_vars=[f"mean_{ct}" for ct in ALL_CT],
                   var_name="cell_type", value_name="methylation")
    long["cell_type"] = long["cell_type"].str.replace("mean_", "", regex=False)
    long.to_csv(PROC / "reference_matrix_A.tsv", sep="\t", index=False)
    log.info("Panel A (independent): %d markers (1250 hyper + 1250 hypo)", len(mq))
    return mq


# ---------------------------------------------------------------- heatmap
def heatmap(panelA, panelB):
    ra = panelA.set_index("marker_id")[[f"mean_{c}" for c in ALL_CT]]
    ra.columns = ALL_CT
    rb = panelB.set_index("marker_id")[[f"mean_{c}" for c in ALL_CT]]
    rb.columns = ALL_CT
    fig, axes = plt.subplots(2, 1, figsize=(7, 9))
    for ax, mat, title in ((axes[0], ra, f"Panel A: independent 1kb-bin markers (n={len(ra)})"),
                           (axes[1], rb, f"Panel B: UXM U250 markers (n={len(rb)})")):
        sns.heatmap(mat, ax=ax, cmap="viridis", vmin=0, vmax=1,
                    cbar_kws=dict(label="mean methylation"))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_yticks([])
    plt.tight_layout()
    fig.savefig(FIG / "reference_heatmap.pdf")
    plt.close(fig)
    log.info("saved reference_heatmap.pdf")


def main():
    chrom_names, chrom, pos = load_cpg_index()
    panelB = build_panel_B(chrom_names, chrom, pos)
    panelA = build_panel_A()
    heatmap(panelA, panelB)
    return 0


if __name__ == "__main__":
    sys.exit(main())
