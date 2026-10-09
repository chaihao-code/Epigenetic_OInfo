#!/usr/bin/env python3
"""06_deconvolve_samples.py — Step 6: 甲状腺样本甲基化去卷积.

目标样本提取:
  从每样本 WGB-Seq methylation_profile bigWig (0-100 百分制) 按 marker 坐标
  取 CpG 值均值 -> 0-1 beta。无 coverage 信息, 缺失即该范围无 CpG 记录。
拟合:
  Panel A/B × Model A/B 共 4 组合; 染色体 block bootstrap (200 次) 给 95% CI;
  拟合优度: 残差 RMSE / R^2, marker 缺失率。
输出:
  results/tables/thyroid_cell_proportions.tsv (每行 sample x 组合)
  results/figures/fig1_cell_composition.{pdf,png} (1A-1D 四面板)
  data/processed/target_marker_methylation_{A,B}.tsv
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from deconv_utils import load_R, fit_nnls, fit_bootstrap, MODEL_B, MODEL_A  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/06_deconvolve_samples.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("deconvolve")

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load(open(ROOT / "config.yaml"))
PROC = ROOT / "data/processed"
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
BW_DIR = Path(CFG["data"]["bigwig_dir"])


def methylation_bigwig(sample):
    hits = sorted(BW_DIR.glob(f"*CEEHRC.{sample}.WGB-Seq.methylation_profile.bigWig"))
    hits = [h for h in hits if not h.name.endswith(".1")]
    return hits[0]


def chrom_arrays(bw, chroms):
    """每染色体 (pos, val/100) 数组 (仅 CpG 有值)。"""
    out = {}
    for c in chroms:
        if c not in bw.chroms():
            continue
        iv = bw.intervals(c, 0, bw.chroms()[c])
        if not iv:
            out[c] = (np.array([], dtype=np.int64), np.array([], dtype=np.float32))
            continue
        arr = np.array(iv)
        out[c] = (arr[:, 0].astype(np.int64), arr[:, 2].astype(np.float32) / 100.0)
    return out


def extract_panel_A(sample, cpg_pos_by_chrom, chrom_names):
    """Panel B (UXM) marker 坐标来自 selected_markers_B.bed (无表头, 0-based start)。

    (函数名保留历史命名; 2026-10-09 起 UXM 属 Panel B)
    """
    bed = pd.read_csv(PROC / "selected_markers_B.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "name"])
    bw = pyBigWig.open(str(methylation_bigwig(sample)))
    carr = chrom_arrays(bw, set(bed["chrom"]))
    vals = np.full(len(bed), np.nan)
    for c, grp in bed.groupby("chrom"):
        pos, v = carr.get(c, (np.array([]), np.array([])))
        if len(pos) == 0:
            continue
        for i, r in grp.iterrows():
            m = (pos >= r["start"]) & (pos < r["end"])
            if m.any():
                vals[i] = v[m].mean()
    bw.close()
    return pd.DataFrame({"marker_id": bed["name"], sample: vals})


def extract_panel_B(sample):
    """Panel A (独立筛选) marker 即 1kb bin (selected_markers_A.bed 含表头)。"""
    bed = pd.read_csv(PROC / "selected_markers_A.bed", sep="\t")
    bw = pyBigWig.open(str(methylation_bigwig(sample)))
    carr = chrom_arrays(bw, set(bed["chrom"]))
    vals = np.full(len(bed), np.nan)
    for c, grp in bed.groupby("chrom"):
        pos, v = carr.get(c, (np.array([]), np.array([])))
        if len(pos) == 0:
            continue
        bs = grp["start"].values
        idx = np.searchsorted(pos, bs)
        for k, i in enumerate(grp.index):
            j0, j1 = idx[k], np.searchsorted(pos, grp.loc[i, "end"])
            if j1 > j0:
                vals[i] = v[j0:j1].mean()
    bw.close()
    return pd.DataFrame({"marker_id": bed["marker_id"], sample: vals})


def extract_targets():
    for panel in ("A", "B"):
        out = None
        for s in CFG["samples"]["included"]:
            df = extract_panel_B(s) if panel == "A" else extract_panel_A(s, None, None)
            out = df if out is None else out.merge(df, on="marker_id")
            log.info("extracted sample %s panel %s", s, panel)
        out.to_csv(PROC / f"target_marker_methylation_{panel}.tsv", sep="\t", index=False)


def deconvolve():
    rows = []
    for panel in ("A", "B"):
        targ = pd.read_csv(PROC / f"target_marker_methylation_{panel}.tsv", sep="\t")
        for model in ("A", "B"):
            R, ids, cts, ch = load_R(panel, model)
            lut = dict(zip(ids, ch))
            for s in CFG["samples"]["included"]:
                y_map = dict(zip(targ["marker_id"], targ[s]))
                y = np.array([y_map.get(m, np.nan) for m in ids])
                missing = np.isnan(y).mean()
                p = fit_nnls(y, R)
                if p is None:
                    log.warning("fit failed %s panel%s model%s", s, panel, model)
                    continue
                ok = ~np.isnan(y)
                resid = y[ok] - R[ok] @ p
                ss_res = float((resid ** 2).sum())
                ss_tot = float(((y[ok] - y[ok].mean()) ** 2).sum())
                mean_, ci = fit_bootstrap(y, R, ch, n_boot=200, seed=42)
                for k, ct in enumerate(cts):
                    rows.append(dict(
                        sample_id=s, panel=panel, model=model,
                        cell_type=ct, fraction=float(p[k]),
                        ci_lo=float(ci[0, k]), ci_hi=float(ci[1, k]),
                        marker_missing_rate=float(missing),
                        fit_rmse=float(np.sqrt((resid ** 2).mean())),
                        fit_r2=1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
                    ))
                log.info("%s panel%s model%s: fit ok, missing=%.3f, r2=%.3f",
                         s, panel, model, missing,
                         1 - ss_res / ss_tot if ss_tot > 0 else np.nan)
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t", index=False)
    return df


def figures(df):
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")
    d = df[(df.panel == "A") & (df.model == "A")].merge(meta, on="sample_id")
    order = meta.sort_values("age")["sample_id"].tolist()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # 1A stacked bar
    ax = axes[0, 0]
    pivot = d.pivot(index="sample_id", columns="cell_type", values="fraction").loc[order]
    bottoms = np.zeros(len(pivot))
    colors = plt.cm.tab10.colors
    for i, ct in enumerate(pivot.columns):
        ax.bar(range(len(pivot)), pivot[ct], bottom=bottoms,
               label=ct.replace("_", " "), color=colors[i % 10], edgecolor="white", lw=0.3)
        bottoms += pivot[ct].values
    ax.set_xticks(range(len(pivot)))
    ax.set_xticklabels([f"{s}\n{a:.0f}y" for s, a in
                        zip(pivot.index, meta.set_index("sample_id").loc[pivot.index, "age"])],
                       fontsize=7)
    ax.set_ylabel("estimated fraction")
    ax.set_title("1A: cell composition (Panel A / Model A)", fontsize=10)
    ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_ylim(0, 1)

    # 1B immune vs age
    ax = axes[0, 1]
    imm = d[d.cell_type.isin(["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"])]
    imm_frac = imm.groupby("sample_id")["fraction"].sum().reset_index().merge(meta, on="sample_id")
    ax.scatter(imm_frac.age, imm_frac.fraction, s=60, color="steelblue", edgecolor="k")
    for _, r in imm_frac.iterrows():
        ax.annotate(r["sample_id"], (r.age, r.fraction), fontsize=6, alpha=0.7)
    ax.set_xlabel("age (y)")
    ax.set_ylabel("total immune fraction")
    ax.set_title("1B: immune fraction vs age", fontsize=10)

    # 1C per-celltype vs age
    ax = axes[1, 0]
    for ct in ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"]:
        sub = d[d.cell_type == ct]
        ax.scatter(sub.age, sub.fraction, s=30, label=ct)
    ax.set_xlabel("age (y)")
    ax.set_ylabel("fraction")
    ax.set_title("1C: cell type fractions vs age", fontsize=10)
    ax.legend(fontsize=7)

    # 1D fit quality
    ax = axes[1, 1]
    fq = df.groupby(["panel", "model", "sample_id"])[["fit_r2", "marker_missing_rate"]].first().reset_index()
    sns.boxplot(data=fq, x="model", y="fit_r2", hue="panel", ax=ax)
    ax.set_title("1D: deconvolution fit quality (R²)", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_cell_composition.pdf")
    fig.savefig(FIG / "fig1_cell_composition.png", dpi=200)
    plt.close(fig)
    log.info("saved fig1_cell_composition")


def main():
    for panel in ("A", "B"):
        if not (PROC / f"target_marker_methylation_{panel}.tsv").exists():
            extract_targets()
            break
    df = deconvolve()
    figures(df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
