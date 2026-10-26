#!/usr/bin/env python3
"""Scatter plots of per-block O-information (full system) vs donor age.

neutrophil, untreated donors only. Reads results_oinfo_n_bins_5/<sample>/<sample>_oinfo.csv
(n_missing == 0) and donor_metadata.tsv.

Outputs (plots/):
  oinfo_vs_age_all_blocks.png        40 blocks, all untreated (n=9)
  oinfo_vs_age_all_blocks_adults.png 40 blocks, adults only (n=7)
  oinfo_vs_age_cpg_blocks.png        CpG-island blocks, all untreated
  oinfo_vs_age_cpg_blocks_adults.png CpG-island blocks, adults only
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

PROJECT = Path(_release_path('@ROOT@/cases/aging/neutrophil'))
RESULTS = PROJECT / "results_oinfo_n_bins_5"
PLOTS = PROJECT / "plots"
PLOTS.mkdir(exist_ok=True)

meta = pd.read_csv(PROJECT / "donor_metadata.tsv", sep="\t")
meta = meta[meta["treatment"] == "untreated"].copy()
meta["age_mid"] = meta["donor_age"].map(
    lambda r: (float(str(r).split("-")[0]) + float(str(r).split("-")[1])) / 2)

rows = []
for f in sorted(glob.glob(str(RESULTS / "*" / "*_oinfo.csv"))):
    if "permutation" in Path(f).name:
        continue
    acc = f.split("/")[-2]
    df = pd.read_csv(f)
    df = df[df["n_missing"] == 0]
    for _, r in df.iterrows():
        rows.append({"epirr": acc, "block": r["block"], "hoi": r["hoi"]})

r = pd.DataFrame(rows).merge(
    meta[["epirr_accession", "age_mid", "life_stage"]],
    left_on="epirr", right_on="epirr_accession")


def plot_grid(df, blocks, out_path, title):
    n = len(blocks)
    ncols = 5
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows),
                             sharex=True)
    axes = axes.flatten()
    for ax, block in zip(axes, blocks):
        g = df[df["block"] == block]
        ax.scatter(g["age_mid"], g["hoi"], s=25, color="steelblue")
        if g["age_mid"].nunique() > 1:
            b, a = __import__("numpy").polyfit(g["age_mid"], g["hoi"], 1)
            xs = __import__("numpy").linspace(g["age_mid"].min(),
                                              g["age_mid"].max(), 10)
            ax.plot(xs, a + b * xs, color="tomato", lw=1)
        short = (block.replace("_purified_bins", "").replace("_bins", "")
                 .replace("_AND_", "\nAND "))
        ax.set_title(short, fontsize=7)
        ax.tick_params(labelsize=7)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=11)
    fig.supxlabel("donor age (years, range midpoint)")
    fig.supylabel("O-information (bits)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"saved {out_path}")


all_blocks = sorted(r["block"].unique())
cpg_blocks = [b for b in all_blocks if b.endswith("_AND_CpG_island")]
adults = r[r["life_stage"] == "adult"]

plot_grid(r, all_blocks, PLOTS / "oinfo_vs_age_all_blocks.png",
          f"O-information vs age, all blocks (untreated, n={r['epirr'].nunique()})")
plot_grid(adults, all_blocks, PLOTS / "oinfo_vs_age_all_blocks_adults.png",
          f"O-information vs age, all blocks (adults, n={adults['epirr'].nunique()})")
plot_grid(r, cpg_blocks, PLOTS / "oinfo_vs_age_cpg_blocks.png",
          f"O-information vs age, CpG-island blocks (untreated, n={r['epirr'].nunique()})")
plot_grid(adults, cpg_blocks, PLOTS / "oinfo_vs_age_cpg_blocks_adults.png",
          f"O-information vs age, CpG-island blocks (adults, n={adults['epirr'].nunique()})")
