#!/usr/bin/env python3
"""Regress per-block O-information (and TC/DTC) on donor age for neutrophil.

Reads results_oinfo_n_bins_5/<sample>/<sample>_oinfo.csv (full system,
n_missing == 0) and merges donor age from donor_metadata.tsv.
G-CSF/Dex treated samples are excluded. Age ranges (e.g. "60-65") are
converted to midpoints. Per block: OLS metric ~ age (slope, R2, p) plus
Spearman.

Outputs:
  oinfo_age_regression.tsv         all untreated samples (n=9, incl. children)
  oinfo_age_regression_adults.tsv  adults only (n=7)
  oinfo_permutation_block_summary.tsv  per-block permutation-null significance
                                       across donors (median p_dev etc.)
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
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path(_release_path('@ROOT@/cases/aging/neutrophil'))
RESULTS = PROJECT / "results_oinfo_n_bins_5"

meta = pd.read_csv(PROJECT / "donor_metadata.tsv", sep="\t")
meta = meta[meta["treatment"] == "untreated"].copy()


def midpoint(r):
    lo, hi = str(r).split("-")
    return (float(lo) + float(hi)) / 2


meta["age_mid"] = meta["donor_age"].map(midpoint)

# ---- observed O-info (full system) per donor x block ----
rows = []
for f in sorted(glob.glob(str(RESULTS / "*" / "*_oinfo.csv"))):
    if "permutation" in Path(f).name:
        continue
    acc = f.split("/")[-2]
    df = pd.read_csv(f)
    df = df[df["n_missing"] == 0]
    for _, r in df.iterrows():
        rows.append({"epirr": acc, "block": r["block"],
                     "hoi": r["hoi"], "tc": r["tc"], "dtc": r["dtc"]})

# ---- permutation-null p values per donor x block ----
prows = []
for f in sorted(glob.glob(str(RESULTS / "*" / "*_permutation_summary.csv"))):
    acc = f.split("/")[-2]
    df = pd.read_csv(f)
    for _, r in df.iterrows():
        prows.append({"epirr": acc, "block": r["block"],
                      "p_dev": r["p_dev"], "p_up": r["p_up"], "p_down": r["p_down"],
                      "null_mean": r["null_mean"], "null_sd": r["null_sd"]})

oinfo = pd.DataFrame(rows)
perm = pd.DataFrame(prows)
r = oinfo.merge(meta[["epirr_accession", "age_mid", "life_stage"]],
                left_on="epirr", right_on="epirr_accession")

# ---- per-block permutation significance across donors ----
ps = perm.groupby("block").agg(
    n_donors=("epirr", "nunique"),
    p_dev_median=("p_dev", "median"),
    p_dev_min=("p_dev", "min"),
    frac_p_dev_lt_005=("p_dev", lambda x: float(np.mean(x < 0.05))),
    p_up_median=("p_up", "median"),
    p_down_median=("p_down", "median"),
).reset_index()
out_p = PROJECT / "oinfo_permutation_block_summary.tsv"
ps.to_csv(out_p, sep="\t", index=False)
print(f"saved {out_p} ({len(ps)} blocks)")


def regress(df, out_name):
    out_rows = []
    for block, g in df.groupby("block"):
        age = g["age_mid"].to_numpy()
        row = {"block": block, "n": len(g)}
        for metric in ("hoi", "tc", "dtc"):
            v = g[metric].to_numpy()
            lr = stats.linregress(age, v)
            rho, sp = stats.spearmanr(age, v)
            row[f"{metric}_slope"] = lr.slope
            row[f"{metric}_r2"] = lr.rvalue ** 2
            row[f"{metric}_p"] = lr.pvalue
            row[f"{metric}_spearman_rho"] = rho
            row[f"{metric}_spearman_p"] = sp
        out_rows.append(row)
    res = pd.DataFrame(out_rows)
    out = PROJECT / out_name
    res.to_csv(out, sep="\t", index=False)
    print(f"saved {out}\n")
    pd.set_option("display.width", 200)
    show = res[["block", "n", "hoi_slope", "hoi_r2", "hoi_p",
                "hoi_spearman_rho", "hoi_spearman_p"]].copy()
    for c in show.columns[2:]:
        show[c] = show[c].map(lambda x: f"{x:.4g}")
    print(show.sort_values("hoi_p").to_string(index=False))


print(f"n untreated samples: {len(meta)}")
regress(r, "oinfo_age_regression.tsv")
adults = r[r["life_stage"] == "adult"]
print(f"n untreated adults: {adults['epirr'].nunique()}")
regress(adults, "oinfo_age_regression_adults.tsv")
