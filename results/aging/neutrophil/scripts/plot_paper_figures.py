#!/usr/bin/env python3
"""Publication-style figures for the neutrophil aging analysis.

Outputs (plots/):
  heatmap_oinfo_age.png        blocks x untreated donors, O-information
  scatter_selected_blocks.png  hoi vs age for the top blocks by |slope| (adults)
"""
import glob
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT = Path("/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil")
RESULTS = PROJECT / "results_oinfo_n_bins_5"
PLOTS = PROJECT / "plots"
PLOTS.mkdir(exist_ok=True)

meta = pd.read_csv(PROJECT / "donor_metadata.tsv", sep="\t")
meta = meta[meta["treatment"] == "untreated"].copy()
meta["age_mid"] = meta["donor_age"].map(
    lambda r: (float(str(r).split("-")[0]) + float(str(r).split("-")[1])) / 2)
meta = meta.sort_values("age_mid")

rows = []
for f in sorted(glob.glob(str(RESULTS / "*" / "*_oinfo.csv"))):
    if "permutation" in Path(f).name:
        continue
    acc = f.split("/")[-2]
    df = pd.read_csv(f)
    df = df[df["n_missing"] == 0]
    for _, r in df.iterrows():
        rows.append({"epirr": acc, "block": r["block"], "hoi": r["hoi"]})

r = pd.DataFrame(rows)

# ---- heatmap: blocks x donors ----
pivot = r.pivot_table(index="block", columns="epirr", values="hoi")
pivot = pivot[meta["epirr_accession"]]  # order donors by age
pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

fig, ax = plt.subplots(figsize=(7, 10))
vmax = np.nanpercentile(np.abs(pivot.values), 98)
im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
               vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(pivot.columns)))
ax.set_xticklabels([f"{a}\n({m})" for a, m in
                    zip(pivot.columns,
                        meta.set_index("epirr_accession")["age_mid"])],
                   fontsize=7)
ax.set_yticks(range(len(pivot.index)))
ax.set_yticklabels([b.replace("_purified_bins", "") for b in pivot.index],
                   fontsize=6)
fig.colorbar(im, ax=ax, shrink=0.4, label="O-information (bits)")
ax.set_title("O-information by block and donor (neutrophil, untreated)")
fig.tight_layout()
out = PLOTS / "heatmap_oinfo_age.png"
fig.savefig(out, dpi=200)
plt.close(fig)
print(f"saved {out}")

# ---- scatter: top blocks by |age slope| (adults) ----
adult_epirr = set(meta.loc[meta["life_stage"] == "adult", "epirr_accession"])
ra = r[r["epirr"].isin(adult_epirr)].merge(
    meta[["epirr_accession", "age_mid"]], left_on="epirr",
    right_on="epirr_accession")

slopes = []
for block, g in ra.groupby("block"):
    b, a = np.polyfit(g["age_mid"], g["hoi"], 1)
    slopes.append({"block": block, "slope": b})
slopes = pd.DataFrame(slopes)
top = slopes.reindex(slopes["slope"].abs().sort_values(ascending=False).index)
top = top[top["slope"].abs() > 1e-9].head(6)["block"].tolist()

fig, axes = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
for ax, block in zip(axes.flatten(), top):
    g = ra[ra["block"] == block]
    ax.scatter(g["age_mid"], g["hoi"], s=40, color="steelblue")
    b, a = np.polyfit(g["age_mid"], g["hoi"], 1)
    xs = np.linspace(g["age_mid"].min(), g["age_mid"].max(), 10)
    ax.plot(xs, a + b * xs, color="tomato", lw=1.5)
    ax.set_title(f"{block.replace('_purified_bins','')}\nslope={b:.4g} bits/yr",
                 fontsize=8)
    ax.tick_params(labelsize=8)
fig.supxlabel("donor age (years)")
fig.supylabel("O-information (bits)")
fig.suptitle("Top blocks by |O-info age slope| (adults, n=7)")
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = PLOTS / "scatter_selected_blocks.png"
fig.savefig(out, dpi=200)
print(f"saved {out}")
