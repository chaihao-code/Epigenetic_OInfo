#!/usr/bin/env python3
"""03_gsea_marginal.py

Signal-change (bin-marginal) preranked GSEA on the immune-adjusted age effect.

Custom gene sets are built exactly as in the original unadjusted pipeline:
each bin is a "gene", and for every mark the bin is classified by the Pearson
correlation of its discrete signal (0-4) across the 6 samples with age:
    up    : r >  0.3 and p < 0.05
    down  : r < -0.3 and p < 0.05
    flat  : everything else
Each mark x direction combination is one gene set (sets with < 5 bins are
dropped). The ranking of the bins is the immune-adjusted metric
    -log10(p_age) * sign(beta_age)
from 01_regress_oinfo_immune.py.

Inputs  : results/regression/<region>.tsv
          <REGION_BASE>/<region>/extracted_signals.tsv
Outputs : results/gsea_marginal/<region>/gsea_report.tsv (+ gseapy outdir)
          logs/03_gsea_marginal.log
"""

import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

AGES = np.array(C.AGES, dtype=float)


def classify_bins(signals_path):
    """Return DataFrame bin_id x {mark}_class in {up, flat, down}."""
    sig = pd.read_csv(signals_path, sep="\t", usecols=["bin_id"] +
                      [f"{s}_{m}" for s in C.SAMPLES for m in C.MARKS])
    out = pd.DataFrame({"bin_id": sig["bin_id"]})
    n = len(AGES)
    denom = np.sqrt(((AGES - AGES.mean()) ** 2).sum())
    for m in C.MARKS:
        X = sig[[f"{s}_{m}" for s in C.SAMPLES]].values.astype(float)
        Xc = X - X.mean(axis=1, keepdims=True)
        r = (Xc * (AGES - AGES.mean())).sum(axis=1) / (denom * np.sqrt((Xc**2).sum(axis=1)))
        with np.errstate(invalid="ignore", divide="ignore"):
            t = r * np.sqrt((n - 2) / (1 - r**2))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
        cls = np.full(len(sig), "flat", dtype=object)
        cls[(r > C.TREND_R_THRESHOLD) & (p < C.TREND_P_THRESHOLD)] = "up"
        cls[(r < -C.TREND_R_THRESHOLD) & (p < C.TREND_P_THRESHOLD)] = "down"
        out[f"{m}_class"] = cls
    return out


def build_gmt(cls_df, gmt_path):
    """Write mark x direction gene sets (bins as genes); return set names."""
    sets = {}
    for m in C.MARKS:
        for direction in ("up", "flat", "down"):
            name = f"{m}_{direction}"
            members = cls_df.loc[cls_df[f"{m}_class"] == direction, "bin_id"].tolist()
            if len(members) >= C.MIN_SET_SIZE:
                sets[name] = members
    with open(gmt_path, "w") as fh:
        for name, members in sets.items():
            fh.write(f"{name}\tna\t" + "\t".join(members) + "\n")
    return list(sets)


def run_region(region):
    reg_path = C.OUT_REGRESSION / f"{region}.tsv"
    sig_path = C.REGION_BASE / region / "extracted_signals.tsv"
    if not (reg_path.exists() and sig_path.exists()):
        logging.warning("skip %s (missing input)", region)
        return
    reg_df = pd.read_csv(reg_path, sep="\t", usecols=["bin_id", "p_age", "beta_age"])
    reg_df["ranking"] = (
        -np.log10(reg_df["p_age"].clip(lower=1e-300)) * np.sign(reg_df["beta_age"])
    )
    cls_df = classify_bins(sig_path)
    rnk = cls_df.merge(reg_df[["bin_id", "ranking"]], on="bin_id", how="inner")
    rnk = rnk.set_index("bin_id")[["ranking"]].dropna()

    out_region = C.OUT_GSEA_MARGINAL / region
    out_region.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="gmt_"))
    gmt = tmp / "marginal.gmt"
    set_names = build_gmt(cls_df, gmt)
    logging.info("%s: %d bins ranked, %d gene sets", region, len(rnk), len(set_names))
    if not set_names or len(rnk) < C.MIN_SET_SIZE:
        return
    res = gp.prerank(
        rnk=rnk,
        gene_sets=str(gmt),
        outdir=str(out_region / "gsea_prerank"),
        permutation_num=C.MARGINAL_PERMUTATIONS,
        seed=C.GSEA_SEED,
        verbose=False,
        min_size=C.MIN_SET_SIZE,
        max_size=C.MARGINAL_MAX_SIZE,
    )
    res.res2d.to_csv(out_region / "gsea_report.tsv", sep="\t", index=False)
    sig = res.res2d[
        (res.res2d["NOM p-val"] < C.GSEA_P_CUTOFF)
        & (res.res2d["FDR q-val"] < C.GSEA_FDR_CUTOFF)
    ]
    logging.info("  significant terms: %d / %d", len(sig), len(res.res2d))


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(C.LOG_DIR / "03_gsea_marginal.log", mode="w"),
                  logging.StreamHandler()],
    )
    logging.info("gseapy %s", gp.__version__)
    for reg_path in sorted(C.OUT_REGRESSION.glob("*.tsv")):
        run_region(reg_path.stem)
    logging.info("DONE")


if __name__ == "__main__":
    main()
