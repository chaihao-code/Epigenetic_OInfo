#!/usr/bin/env python3
"""04_summarize.py

Collect significant GSEA terms (NOM p < 0.05 AND FDR q < 0.25, the dual rule)
across all regions into two summary tables.

Outputs : results/tables/gsea_pathway_significant.tsv
          results/tables/gsea_marginal_significant.tsv
"""

import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

# Lead_genes lists every member gene/bin of a term and can reach hundreds of
# thousands of entries for the flat marginal sets -- drop it from summaries
# (full per-term reports stay in the per-region outputs).
DROP_COLUMNS = ["Lead_genes"]


def _drop(df):
    return df.drop(columns=[c for c in DROP_COLUMNS if c in df.columns])


def summarize_pathway():
    rows = []
    for region_dir in sorted(C.OUT_GSEA_PATHWAY.iterdir()):
        if not region_dir.is_dir():
            continue
        for report in region_dir.glob("*.report.csv"):
            df = pd.read_csv(report)
            df = df[
                (df["NOM p-val"] < C.GSEA_P_CUTOFF)
                & (df["FDR q-val"] < C.GSEA_FDR_CUTOFF)
            ].copy()
            if df.empty:
                continue
            df.insert(0, "library", report.stem.replace(".report", ""))
            df.insert(0, "region", region_dir.name)
            rows.append(df)
    out = (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=["region", "library"])
    )
    out = _drop(out)
    out.to_csv(C.OUT_TABLES / "gsea_pathway_significant.tsv", sep="\t", index=False)
    logging.info("pathway: %d significant terms from %d regions", len(out), len(rows))
    return out


def summarize_marginal():
    rows = []
    for region_dir in sorted(C.OUT_GSEA_MARGINAL.iterdir()):
        report = region_dir / "gsea_report.tsv"
        if not report.exists():
            continue
        df = pd.read_csv(report, sep="\t")
        df = df[
            (df["NOM p-val"] < C.GSEA_P_CUTOFF)
            & (df["FDR q-val"] < C.GSEA_FDR_CUTOFF)
        ].copy()
        if df.empty:
            continue
        df.insert(0, "region", region_dir.name)
        rows.append(df)
    out = (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(columns=["region"])
    )
    out = _drop(out)
    out.to_csv(C.OUT_TABLES / "gsea_marginal_significant.tsv", sep="\t", index=False)
    logging.info("marginal: %d significant terms from %d regions", len(out), len(rows))
    return out


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    summarize_pathway()
    summarize_marginal()
