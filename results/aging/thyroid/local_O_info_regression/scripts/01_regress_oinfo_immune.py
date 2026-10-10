#!/usr/bin/env python3
"""01_regress_oinfo_immune.py

Per-bin linear regression of local O-information on donor age, controlling
for the estimated total immune-cell fraction:

    O_info(bin, sample) ~ age + immune_fraction(sample)

One regression per bin; samples are the 6 healthy thyroid donors
(30/37/37/46/55/67 y). Both an age+immune model (primary) and an age-only
model (for comparison with the original unadjusted analysis) are fitted with
a shared, fully vectorized least-squares routine.

Inputs  : <REGION_BASE>/<region>/oinfo_per_sample.bed  (from run_analysis_generic_excluded.py)
          PROPORTIONS_TSV (Panel A / Model A immune fractions)
Outputs : results/regression/<region>.tsv  (one row per bin)
          logs/01_regress.log
"""

import logging
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import config as C


def load_immune_fractions():
    """Sample -> total immune fraction, restricted to the chosen panel/model."""
    df = pd.read_csv(C.PROPORTIONS_TSV, sep="\t")
    df = df[(df["panel"] == C.PANEL) & (df["model"] == C.MODEL)]
    available = set(df["cell_type"].unique())
    missing = set(C.IMMUNE_CELL_TYPES) - available
    if missing:
        raise ValueError(f"Immune cell types missing from {C.PROPORTIONS_TSV}: {missing}")
    imm = df[df["cell_type"].isin(C.IMMUNE_CELL_TYPES)]
    frac = imm.groupby("sample_id")["fraction"].sum()
    frac = frac.reindex(C.SAMPLES)
    if frac.isna().any():
        raise ValueError(f"Samples without immune fraction: {frac[frac.isna()].index.tolist()}")
    logging.info(
        "Immune fraction (%s/%s) = sum of %s",
        C.PANEL, C.MODEL, "+".join(C.IMMUNE_CELL_TYPES),
    )
    for s, v in frac.items():
        logging.info("  %s: %.4f", s, v)
    return frac


def ols_fit(Y, X):
    """Vectorized OLS of every row of Y on the shared design matrix X.

    Y : (n_bins, n_samples) response matrix
    X : (n_samples, n_params) design matrix (columns already include intercept)
    Returns beta, se, tstat, pval, r2  -- each (n_bins, n_params) except r2.
    """
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = Y @ X @ XtX_inv
    resid = Y - beta @ X.T
    dof = X.shape[0] - X.shape[1]
    sigma2 = (resid**2).sum(axis=1) / dof                      # (n_bins,)
    se = np.sqrt(np.outer(sigma2, np.diag(XtX_inv)))
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = beta / se
    pval = 2 * stats.t.sf(np.abs(tstat), dof)
    ss_res = (resid**2).sum(axis=1)
    ss_tot = ((Y - Y.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    with np.errstate(invalid="ignore"):
        r2 = np.where(ss_tot > 0, 1 - ss_res / ss_tot, np.nan)
    return beta, se, tstat, pval, r2


def discover_regions():
    regions = []
    for d in sorted(C.REGION_BASE.iterdir()):
        if d.is_dir() and (d / "oinfo_per_sample.bed").exists() and (d / "extracted_signals.tsv").exists():
            regions.append(d.name)
    return regions


def regress_region(region, immune):
    d = C.REGION_BASE / region
    df = pd.read_csv(d / "oinfo_per_sample.bed", sep="\t")
    sample_cols = [c for c in df.columns if c.endswith("y") and c.split("_")[0] in C.SAMPLES]
    sample_cols = sorted(sample_cols, key=lambda c: int(c.split("_")[1].rstrip("y")))
    assert [c.split("_")[0] for c in sample_cols] == C.SAMPLES, f"unexpected samples in {region}"

    ages = np.array([int(c.split("_")[1].rstrip("y")) for c in sample_cols], dtype=float)
    sample_ids = [c.split("_")[0] for c in sample_cols]
    immune_v = np.array([immune[s] for s in sample_ids], dtype=float)
    Y = df[sample_cols].values.astype(float)

    # age + immune model (primary)
    X2 = np.column_stack([np.ones_like(ages), ages, immune_v])
    b2, se2, t2, p2, r2 = ols_fit(Y, X2)
    # age-only model (comparison with original analysis)
    X1 = np.column_stack([np.ones_like(ages), ages])
    b1, se1, t1, p1, _ = ols_fit(Y, X1)

    out = df[["chr", "start", "end", "bin_id"]].copy()
    for c in sample_cols:
        out[c] = df[c]
    out["beta_age"] = b2[:, 1]
    out["se_age"] = se2[:, 1]
    out["t_age"] = t2[:, 1]
    out["p_age"] = p2[:, 1]
    out["beta_immune"] = b2[:, 2]
    out["se_immune"] = se2[:, 2]
    out["t_immune"] = t2[:, 2]
    out["p_immune"] = p2[:, 2]
    out["r2_model"] = r2
    out["slope_age_unadj"] = b1[:, 1]
    out["se_age_unadj"] = se1[:, 1]
    out["p_age_unadj"] = p1[:, 1]
    # attenuation of the age slope after immune adjustment
    with np.errstate(invalid="ignore", divide="ignore"):
        out["age_slope_attenuation"] = 1 - out["beta_age"] / out["slope_age_unadj"]

    out_path = C.OUT_REGRESSION / f"{region}.tsv"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    n_sig = (out["p_age"] < 0.05).sum()
    logging.info(
        "%s: bins=%d | adj model r2 mean=%.3f | p_age<0.05: %d (unadj: %d)",
        region, len(out), np.nanmean(r2), n_sig, (out["p_age_unadj"] < 0.05).sum(),
    )
    return len(out)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(C.LOG_DIR / "01_regress.log", mode="w"),
                  logging.StreamHandler()],
    )
    logging.info("Python %s | numpy %s | pandas %s | scipy %s",
                 sys.version.split()[0], np.__version__, pd.__version__,
                 stats.__version__ if hasattr(stats, "__version__") else "na")
    immune = load_immune_fractions()
    regions = discover_regions()
    logging.info("Discovered %d regions under %s", len(regions), C.REGION_BASE)
    total = 0
    for region in regions:
        total += regress_region(region, immune)
    logging.info("DONE. %d regions, %d bins total.", len(regions), total)


if __name__ == "__main__":
    main()
