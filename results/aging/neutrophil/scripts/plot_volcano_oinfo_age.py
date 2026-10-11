#!/usr/bin/env python3
"""Volcano plot: O-information age slope vs -log10(p) per block (adults, n=7).

Reads oinfo_age_regression_adults.tsv produced by oinfo_age_regression.py.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path("/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil")
PLOTS = PROJECT / "plots"
PLOTS.mkdir(exist_ok=True)

res = pd.read_csv(PROJECT / "oinfo_age_regression_adults.tsv", sep="\t")
res["neglog10_p"] = -np.log10(res["hoi_p"].clip(lower=1e-300))

fig, ax = plt.subplots(figsize=(8, 6))
colors = np.where(res["hoi_slope"] > 0, "tomato", "steelblue")
ax.scatter(res["hoi_slope"], res["neglog10_p"], c=colors, s=28, alpha=0.8)

# label blocks with p < 0.05
for _, r in res[res["hoi_p"] < 0.05].iterrows():
    ax.annotate(r["block"].replace("_purified_bins", ""),
                (r["hoi_slope"], r["neglog10_p"]),
                fontsize=6, alpha=0.85,
                xytext=(3, 3), textcoords="offset points")

ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=1)
ax.axvline(0, ls="--", c="grey", lw=1)
ax.set_xlabel("O-information ~ age slope (bits/year)")
ax.set_ylabel(r"$-\log_{10}(p)$")
ax.set_title("Volcano: per-block O-information vs donor age (neutrophil adults, n=7)")
fig.tight_layout()
out = PLOTS / "volcano_oinfo_slope_vs_p_adults.png"
fig.savefig(out, dpi=200)
print(f"saved {out}")
