#!/usr/bin/env python3
"""02_gsea_pathway.py

Pathway-level preranked GSEA on the immune-adjusted age effect.

For each region:
1. Map 1-kb bins to genes via TSS overlap (+/- TSS_WINDOW, Gencode hg38 GTF,
   same mapping as the original unadjusted pipeline).
2. Gene-level ranking metric (mirroring the original analysis, but from the
   immune-adjusted model):
       ranking = -log10(p_age) * sign(beta_age)
   where p_age / beta_age come from 01_regress_oinfo_immune.py. A gene hit by
   several bins keeps the bin with the smallest p_age.
3. gseapy.prerank against GO_Biological_Process_2021 and MSigDB_Hallmark_2020.

Inputs  : results/regression/<region>.tsv
Outputs : results/gsea_pathway/<region>/<library>/  (gseapy report + plots)
          results/gsea_pathway/<region>/<region>.ranked_genes.tsv
          logs/02_gsea_pathway.log
"""

import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C


def build_tss_bed(path):
    """Extract [TSS-WINDOW, TSS+WINDOW] windows per gene from the GTF."""
    rows = []
    with open(C.GTF_PATH) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "transcript":
                continue
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            attr = dict(
                item.strip().split(" ", 1) for item in f[8].rstrip(";").split("; ")
            )
            gene = attr.get("gene_name", "").strip('"')
            if not gene:
                continue
            tss = start if strand == "+" else end
            rows.append((chrom, max(0, tss - C.TSS_WINDOW), tss + C.TSS_WINDOW, gene))
    tss_df = pd.DataFrame(rows, columns=["chr", "start", "end", "gene_name"])
    tss_df.to_csv(path, sep="\t", index=False, header=False)
    logging.info("TSS windows: %d (from %s)", len(tss_df), C.GTF_PATH)


def annotate_bins(df, tss_bed):
    """Return df plus gene_name_tss via bedtools intersect."""
    tmp = Path(tempfile.mkdtemp(prefix="tss_"))
    bins_bed = tmp / "bins.bed"
    df[["chr", "start", "end"]].to_csv(bins_bed, sep="\t", index=False, header=False)
    out = tmp / "intersect.tsv"
    subprocess.run(
        f"bedtools intersect -a {bins_bed} -b {tss_bed} -wa -wb > {out}",
        shell=True, check=True,
    )
    if out.stat().st_size == 0:
        df["gene_name_tss"] = np.nan
        return df
    inter = pd.read_csv(out, sep="\t", header=None,
                        names=["chr", "start", "end", "b_chr", "b_start", "b_end", "gene_name"])
    gene_map = (
        inter.groupby(["chr", "start", "end"])["gene_name"]
        .apply(lambda x: ",".join(sorted(set(x))))
        .reset_index()
        .rename(columns={"gene_name": "gene_name_tss"})
    )
    return df.merge(gene_map, on=["chr", "start", "end"], how="left")


def gene_ranking(reg_df):
    """Gene-level ranking from the immune-adjusted regression.

    ranking = -log10(p_age) * sign(beta_age); multi-hit genes keep min p_age.
    """
    df = reg_df.dropna(subset=["p_age", "beta_age", "gene_name_tss"]).copy()
    df["ranking"] = -np.log10(df["p_age"].clip(lower=1e-300)) * np.sign(df["beta_age"])
    # a bin may map to several TSS genes -> one (gene, bin) row each
    df = df.assign(gene_name_tss=df["gene_name_tss"].str.split(",")).explode("gene_name_tss")
    idx = df.groupby("gene_name_tss")["p_age"].idxmin()
    rnk = df.loc[idx, ["gene_name_tss", "ranking", "p_age", "beta_age", "bin_id"]]
    rnk = rnk.rename(columns={"gene_name_tss": "gene"}).sort_values(
        "ranking", ascending=False
    )
    return rnk.reset_index(drop=True)


def run_region(region, tss_bed):
    reg_path = C.OUT_REGRESSION / f"{region}.tsv"
    if not reg_path.exists():
        logging.warning("skip %s (no regression table)", region)
        return
    reg_df = pd.read_csv(reg_path, sep="\t")
    reg_df = annotate_bins(reg_df, tss_bed)
    rnk = gene_ranking(reg_df)
    out_region = C.OUT_GSEA_PATHWAY / region
    out_region.mkdir(parents=True, exist_ok=True)
    rnk.to_csv(out_region / f"{region}.ranked_genes.tsv", sep="\t", index=False)
    logging.info("%s: %d genes ranked", region, len(rnk))
    if len(rnk) < 100:
        return
    rnk_gp = rnk.set_index("gene")[["ranking"]]
    for lib in C.PATHWAY_LIBRARIES:
        try:
            res = gp.prerank(
                rnk=rnk_gp,
                gene_sets=lib,
                outdir=str(out_region / lib.replace(" ", "_")),
                permutation_num=C.PATHWAY_PERMUTATIONS,
                seed=C.GSEA_SEED,
                verbose=False,
                min_size=C.MIN_SET_SIZE,
                max_size=C.PATHWAY_MAX_SIZE,
            )
            res.res2d.to_csv(
                out_region / f"{lib.replace(' ', '_')}.report.csv", index=False
            )
            n_sig = (
                (res.res2d["NOM p-val"] < C.GSEA_P_CUTOFF)
                & (res.res2d["FDR q-val"] < C.GSEA_FDR_CUTOFF)
            ).sum()
            logging.info("  %s: %d terms, %d significant", lib, len(res.res2d), n_sig)
        except Exception as exc:  # keep going: network / library failures
            logging.error("  %s failed: %s", lib, exc)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(C.LOG_DIR / "02_gsea_pathway.log", mode="w"),
                  logging.StreamHandler()],
    )
    logging.info("gseapy %s", gp.__version__)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as tmp:
        tss_bed = tmp.name
    build_tss_bed(tss_bed)
    for reg_path in sorted(C.OUT_REGRESSION.glob("*.tsv")):
        run_region(reg_path.stem, tss_bed)
    logging.info("DONE")


if __name__ == "__main__":
    main()
