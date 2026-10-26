#!/usr/bin/env python3
"""Scan mi_trend csvs (n3/n5/n10) for blocks where the joint-total curve differs
a lot across n and stays below 1.0 even at n=10."""

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
from pathlib import Path
import numpy as np
import pandas as pd

PLOTS = Path(_release_path('@ROOT@/cases/mouse_preimp/mutual_information/plots'))
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
