#!/usr/bin/env python3
"""Scan mi_trend csvs (n3/n5/n10) for blocks where the joint-total curve differs
a lot across n and stays below 1.0 even at n=10."""
from pathlib import Path
import numpy as np
import pandas as pd

PLOTS = Path("/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/mutual_information/plots")
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM']

rows = []
for f in sorted((PLOTS / "n10").glob("*_mi_trend.csv")):
    name = f.name.replace("_mi_trend.csv", "")
    anchor = name.rsplit("_", 1)[-1]
    block = name.rsplit("_", 1)[0]
    means = {}
    ok = True
    for n in (3, 5, 10):
        fn = PLOTS / f"n{n}" / f"{name}_mi_trend.csv"
        if not fn.exists():
            ok = False
            break
        df = pd.read_csv(fn)
        t = df[df['target'] == 'total'].set_index('stage').reindex(STAGES)['normalized_mi']
        means[n] = t.values
    if not ok:
        continue
    m3, m5, m10 = means[3], means[5], means[10]
    rows.append({
        'block': block,
        'anchor': anchor,
        'mean3': np.nanmean(m3), 'mean5': np.nanmean(m5), 'mean10': np.nanmean(m10),
        'max10': np.nanmax(m10),
        'spread_10_minus_3': np.nanmean(m10) - np.nanmean(m3),
    })

df = pd.DataFrame(rows)
# candidate: n10 stays below 1.0, big spread between n3 and n10
cand = df[(df['max10'] < 0.98) & (df['mean10'] > 0.5)].copy()
cand = cand.sort_values('spread_10_minus_3', ascending=False)
print(f"total (block,anchor) pairs: {len(df)}; candidates (max10<0.98, mean10>0.5): {len(cand)}")
print(cand.head(30).to_string(index=False))
print("\n--- also: largest spread regardless of n10 cap ---")
print(df.sort_values('spread_10_minus_3', ascending=False).head(10).to_string(index=False))
