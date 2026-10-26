#!/usr/bin/env python3
"""09_marker_go_analysis.py — Step 10a: 免疫 marker / GO 富集检查.

三部分:
A. RNA-seq marker 检查 (IHEC mRNA-Seq RPKM bigWig, 6 样本):
   - 免疫标志: PTPRC CD3D CD3E TRAC MS4A1 CD79A LST1 TYROBP
   - 甲状腺上皮标志: TG TPO SLC5A5 TSHR (+PAX8)
   - 基因体平均 RPKM (forward+reverse 求和), 与 age / 估计免疫比例相关
B. 免疫特异区域 overlap 富集:
   - 各区 sig_fdr_0.10 且 slope<0 (losing-redundancy) bins 与
     immune_specific_regions.bed 的 overlap
   - 背景: 同区域全部 bins; 另做区域构成匹配的随机置换 (1000 次)
C. GO before/after filtering:
   - 复现原 GSEA 流程 (run_gsea_tss.py): ranking = -log10(p)*sign(r),
     每基因取最小 p; gp.prerank GO_Biological_Process_2021, seed 固定
   - before: 原始 annotated_tss 表; after: 排除免疫特异 bins 后重算
     per-bin r/p 再排名 (n=6, 与原分析一致)
输出:
  results/tables/immune_marker_analysis.tsv
  results/tables/immune_region_enrichment.tsv
  results/tables/go_before_after_terms.tsv
  results/figures/fig4_marker_go.{pdf,png}
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import statsmodels.api as sm
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/10_marker_go.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("marker_go")

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"

IMMUNE_GENES = ["PTPRC", "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "LST1", "TYROBP"]
THYROID_GENES = ["TG", "TPO", "SLC5A5", "TSHR", "PAX8"]
MARKER_GENES = IMMUNE_GENES + THYROID_GENES
AGE_COLS = ["CEMT0098_30y", "CEMT0130_37y", "CEMT0132_37y",
            "CEMT0042_46y", "CEMT0044_55y", "CEMT0040_67y"]
AGES = [30, 37, 37, 46, 55, 67]
SAMPLES = ["CEMT0098", "CEMT0130", "CEMT0132", "CEMT0042", "CEMT0044", "CEMT0040"]
RNA_DIR = Path(_release_path('@ROOT@/data/aging/thyroid/rna_bigwig'))
IMMUNE_KEYWORDS = re.compile(
    r"t cell|t-cell|leukocyte|lymphocyte|mhc|antigen|immun|interferon|"
    r"myeloid|monocyte|granulocyte|neutrophil|nk cell|natural killer|"
    r"inflammat|chemokine|cytokine|adhesion", re.I)
N_PERM = 1000


# ---------- A. RNA-seq marker ----------
def load_gene_coords(gtf):
    """gencode GTF -> {gene_name: (chrom, start, end, strand)} 蛋白编码优先"""
    coords = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if p[2] != "gene":
                continue
            attrs = dict(re.findall(r'(\S+) "([^"]*)"', p[8]))
            name, gtype = attrs.get("gene_name"), attrs.get("gene_type", "")
            if name in MARKER_GENES and name not in coords and gtype in ("protein_coding", "TR_C_gene"):
                coords[name] = (p[0], int(p[3]) - 1, int(p[4]), p[6])
    return coords


def rna_marker_check(cfg):
    coords = load_gene_coords(cfg["data"]["gencode_gtf"])
    missing = [g for g in MARKER_GENES if g not in coords]
    log.info("gene coords: %d/%d found; missing=%s", len(coords), len(MARKER_GENES), missing)
    rows = []
    bw_cache = {}
    for gene, (chrom, start, end, strand) in coords.items():
        for s in SAMPLES:
            if s not in bw_cache:
                fw = sorted(RNA_DIR.glob(f"*.{s}.mRNA-Seq.rpkm_forward.bigWig"))[0]
                rv = sorted(RNA_DIR.glob(f"*.{s}.mRNA-Seq.rpkm_reverse.bigWig"))[0]
                bw_cache[s] = (pyBigWig.open(str(fw)), pyBigWig.open(str(rv)))
            bwf, bwr = bw_cache[s]
            # IHEC rpkm bigWig 为链特异性有符号值: 正链基因为正,
            # 负链基因存于 reverse 文件且为负值 -> 取绝对值;
            # 无覆盖区间 (未表达) 为 NaN -> 计 0
            vf = np.array([np.nan if v is None else abs(v)
                           for v in (bwf.values(chrom, start, end) or [])], dtype=float)
            vr = np.array([np.nan if v is None else abs(v)
                           for v in (bwr.values(chrom, start, end) or [])], dtype=float)
            rpkm = (np.nanmean(vf) if vf.size and np.isfinite(vf).any() else 0.0) \
                 + (np.nanmean(vr) if vr.size and np.isfinite(vr).any() else 0.0)
            detected = bool((vf.size and np.isfinite(vf).any()) or (vr.size and np.isfinite(vr).any()))
            rows.append(dict(gene=gene, group=("immune" if gene in IMMUNE_GENES else "thyroid"),
                             sample_id=s, age=dict(zip(SAMPLES, AGES))[s],
                             gene_body_rpkm=rpkm, detected=detected))
    df = pd.DataFrame(rows)
    df.to_csv(TBL / "immune_marker_analysis.tsv", sep="\t", index=False)

    # 与 age / 免疫比例相关
    prop = pd.read_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t")
    imm = (prop[(prop.panel == "A") & (prop.model == "A")
                & prop.cell_type.isin(["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"])]
           .groupby("sample_id")["fraction"].sum())
    corr_rows = []
    for gene, g in df.groupby("gene"):
        x_age = g.age.values.astype(float)
        y = g.gene_body_rpkm.values.astype(float)
        ra, pa = stats_spearman(x_age, y)
        ri, pi = stats_spearman(imm.loc[SAMPLES].values, y)
        # 免疫 marker 均值作为 immune score
        corr_rows.append(dict(gene=gene, group=g["group"].iloc[0],
                              spearman_r_age=ra, spearman_p_age=pa,
                              spearman_r_immunefrac=ri, spearman_p_immunefrac=pi,
                              mean_rpkm=np.mean(y)))
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(TBL / "immune_marker_correlation.tsv", sep="\t", index=False)
    log.info("marker correlation:\n%s", corr.to_string(index=False))
    return df, corr


def stats_spearman(x, y):
    from scipy import stats
    with np.errstate(invalid="ignore"):
        r, p = stats.spearmanr(x, y)
    return float(r), float(p)


# ---------- B. overlap 富集 ----------
def overlap_enrichment(cfg):
    imm = pd.read_csv(ROOT / "data/processed/immune_specific_regions.bed", sep="\t")
    imm_bins = set(f"{c}:{s}-{e}" for c, s, e in zip(imm.chrom, imm.start, imm.end))
    ldir = Path(cfg["data"]["local_oinfo_dir"])
    rows = []
    rng = np.random.default_rng(42)
    for d in sorted(ldir.iterdir()):
        f = d / "oinfo_per_sample_fdr.bed"
        if not f.exists():
            continue
        df = pd.read_csv(f, sep="\t")
        n_all = len(df)
        if n_all < 50:
            continue
        in_imm = df.bin_id.isin(imm_bins).values
        # 与原分析候选定义一致: 未校正 p<0.05 且 slope<0 (losing redundancy);
        # 记录: per-bin FDR<0.10 位点过少 (<5/区域), 无法以 FDR 集合作前景
        sig_down = (df.p < 0.05) & (df.slope < 0)
        n_sig = int(sig_down.sum())
        if n_sig < 5:
            continue
        obs_frac = (df.bin_id[sig_down].isin(imm_bins)).mean()
        # 背景: 全部 bins 中免疫比例
        bg_frac = in_imm.mean()
        # 区域构成匹配置换: 从全部 bins 随机抽 n_sig 个
        perm = []
        idx_all = np.arange(n_all)
        for _ in range(N_PERM):
            sub = rng.choice(idx_all, size=n_sig, replace=False)
            perm.append(in_imm[sub].mean())
        perm = np.array(perm)
        emp_p = float((perm >= obs_frac).mean())
        rows.append(dict(region=d.name, n_bins=n_all, n_sig_down=n_sig,
                         obs_immune_frac=round(obs_frac, 4),
                         bg_immune_frac=round(bg_frac, 4),
                         perm_mean=round(perm.mean(), 4),
                         perm_q95=round(np.percentile(perm, 95), 4),
                         empirical_p=max(emp_p, 1 / N_PERM)))
    out = pd.DataFrame(rows)
    out.to_csv(TBL / "immune_region_enrichment.tsv", sep="\t", index=False)
    log.info("overlap enrichment:\n%s", out.to_string(index=False))
    return out


# ---------- C. GO before/after ----------
def gene_ranking(df):
    """复现 run_gsea_tss.py: ranking = -log10(p)*sign(r), 每基因最小 p"""
    d = df[df.gene_name_tss.notna() & (df.gene_name_tss != "")].copy()
    d["p"] = d["p"].clip(lower=1e-300)
    d = d.sort_values("p").drop_duplicates("gene_name_tss")
    return pd.DataFrame({
        "gene": d.gene_name_tss,
        "rank": -np.log10(d.p) * np.sign(d.r),
    }).sort_values("rank", ascending=False)


def recompute_rp(df, drop_bins):
    """排除免疫特异 bins 后, 用每样本 local O-info 重算 per-bin r/p (n=6)"""
    from scipy import stats
    d = df[~df.bin_id.isin(drop_bins)].copy()
    if len(d) < 30:
        return None
    X = d[AGE_COLS].values.astype(float)
    rs, ps = [], []
    for i in range(len(d)):
        r, p = stats.pearsonr(AGES, X[i])
        rs.append(r)
        ps.append(p)
    d["r"] = rs
    d["p"] = ps
    return d


def go_before_after(cfg):
    import gseapy as gp
    imm = pd.read_csv(ROOT / "data/processed/immune_specific_regions.bed", sep="\t")
    imm_bins = set(f"{c}:{s}-{e}" for c, s, e in zip(imm.chrom, imm.start, imm.end))
    ldir = Path(cfg["data"]["local_oinfo_dir"])
    rows = []
    for d in sorted(ldir.iterdir()):
        f = d / "oinfo_per_sample_fdr_annotated_tss.tsv"
        if not f.exists():
            continue
        df = pd.read_csv(f, sep="\t", low_memory=False)
        if "gene_name_tss" not in df.columns or len(df) < 100:
            continue
        for tag, dd in (("before", df), ("after", recompute_rp(df, imm_bins))):
            if dd is None:
                continue
            rk = gene_ranking(dd)
            if len(rk) < 50:
                continue
            pre = gp.prerank(rnk=rk, gene_sets="GO_Biological_Process_2021",
                             permutation_num=1000, seed=42, no_plot=True,
                             outdir=None, verbose=False)
            res = pre.res2d.copy()
            res["region"] = d.name
            res["filter"] = tag
            res["immune_term"] = res.Term.str.contains(IMMUNE_KEYWORDS)
            rows.append(res)
            log.info("GSEA %s %s: terms=%d immune_terms_fdr05=%d", d.name, tag,
                     len(res), ((res["FDR q-val"] < 0.05) & res.immune_term).sum())
    out = pd.concat(rows, ignore_index=True)
    keep = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "region", "filter", "immune_term"]
    out[keep].to_csv(TBL / "go_before_after_terms.tsv", sep="\t", index=False)
    return out


def main():
    import sys as _s
    cfg = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
    if "--overlap-only" in _s.argv:
        df = pd.read_csv(TBL / "immune_marker_analysis.tsv", sep="\t")
        ov = overlap_enrichment(cfg)
        go = None
        go = pd.read_csv(TBL / "go_before_after_terms.tsv", sep="\t")
    else:
        df, corr = rna_marker_check(cfg)
        ov = overlap_enrichment(cfg)
        go = go_before_after(cfg)

    # 汇总: 每 region before/after 的 immune term 计数 (FDR<0.05 与 FDR<0.25)
    go["FDR q-val"] = pd.to_numeric(go["FDR q-val"], errors="coerce")
    summ = (go.groupby(["region", "filter", "immune_term"])["FDR q-val"]
            .agg(n_terms="count", n_fdr05=lambda s: (s < 0.05).sum(),
                 n_fdr25=lambda s: (s < 0.25).sum())
            .reset_index())
    summ.to_csv(TBL / "go_immune_term_summary.tsv", sep="\t", index=False)

    # ---------------- Figure 4 ----------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    # 4A: marker RPKM heatmap (log10)
    ax = axes[0]
    piv = df.pivot_table(index="gene", columns="sample_id", values="gene_body_rpkm")
    piv = piv.loc[MARKER_GENES]
    im = ax.imshow(np.log10(piv[SAMPLES].values + 0.01), aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels([f"{s[-4:]}\n{a}y" for s, a in zip(SAMPLES, AGES)], fontsize=7)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels(piv.index, fontsize=7)
    for i, g in enumerate(piv.index):
        ax.get_yticklabels()[i].set_color(
            "crimson" if g in IMMUNE_GENES else "darkgreen")
    plt.colorbar(im, ax=ax, label="log10 RPKM")
    ax.set_title("4A: immune / thyroid marker expression", fontsize=9)

    # 4B: immune score (mean immune marker RPKM) vs age & immune fraction
    ax = axes[1]
    score = df[df.group == "immune"].groupby("sample_id")["gene_body_rpkm"].mean()
    ax.scatter(AGES, score.loc[SAMPLES].values, s=60, color="crimson", edgecolor="k")
    from scipy import stats as st
    r, p = st.spearmanr(AGES, score.loc[SAMPLES].values)
    ax.set_xlabel("age (y)")
    ax.set_ylabel("mean immune-marker RPKM")
    ax.set_title(f"4B: immune score vs age\nSpearman r={r:.2f}, p={p:.3g}", fontsize=9)

    # 4C: overlap enrichment obs vs perm null
    ax = axes[2]
    if len(ov):
        y = np.arange(len(ov))
        ax.errorbar(ov.obs_immune_frac, y, xerr=[np.maximum(ov.obs_immune_frac - ov.perm_q95, 0)],
                    fmt="o", color="crimson", ms=5, label="observed (q95 perm error)")
        ax.scatter(ov.bg_immune_frac, y, color="gray", s=15, label="region background")
        ax.set_yticks(y)
        ax.set_yticklabels(ov.region, fontsize=6)
        ax.set_xlabel("fraction of losing-redundancy bins\nthat are immune-specific")
        ax.set_title("4C: immune-region overlap", fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_marker_go.pdf")
    fig.savefig(FIG / "fig4_marker_go.png", dpi=200)
    plt.close(fig)
    log.info("saved fig4")
    return 0


if __name__ == "__main__":
    sys.exit(main())
