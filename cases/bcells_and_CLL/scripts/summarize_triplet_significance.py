#!/usr/bin/env python3
"""Extract triplet permutation tests; BH families are all triplets/blocks per sample."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "results_oinfo_triplet_n_bins_5"


def bh(p):
    p = np.asarray(p, dtype=float)
    if not np.isfinite(p).all():
        raise ValueError("Nonfinite P values")
    order = np.argsort(p)
    q = np.empty(len(p))
    q[order] = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def summarize(sample):
    folder = ROOT / sample
    df = pd.read_csv(folder / f"{sample}_null_summary.csv")
    df = df.loc[df["order"] == 3].copy()
    if len(df) != 1400 or (df["n_perm"].nunique() != 1 or not (df["n_perm"] > 0).all()):
        raise ValueError(f"{sample}: expected 1400 triplets with one positive permutation count (new runs default to 10000)")
    df.insert(0, "sample", sample)
    for p in ("p_dev", "p_up", "p_down"):
        df[p.replace("p_", "q_") + "_sample_BH"] = bh(df[p])
    df.to_csv(folder / f"{sample}_triplet_significance.tsv", sep="\t", index=False)
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", help="Summarize one completed sample")
    args = parser.parse_args()
    if args.sample:
        summarize(args.sample)
    else:
        samples = sorted(p.name for p in ROOT.iterdir() if p.is_dir() and p.name.startswith("CEEHRC."))
        if len(samples) != 17:
            raise ValueError(f"Expected 17 samples, got {len(samples)}")
        df = pd.concat([summarize(s) for s in samples], ignore_index=True)
        for p in ("p_dev", "p_up", "p_down"):
            df[p.replace("p_", "q_") + "_global_BH"] = bh(df[p])
        path = ROOT / "all_samples_triplet_significance.tsv"
        df.to_csv(path, sep="\t", index=False)
        print(f"Wrote {path}: {len(df)} rows")


if __name__ == "__main__":
    main()
