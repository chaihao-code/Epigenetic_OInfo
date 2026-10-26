#!/usr/bin/env python3
"""
Discretize genomic signals using global quantile-based binning.

Input modes:
  1) bigWig + BED bins: --input signal.bigWig --bins bins.bed
  2) 4-column BED:      --input signal.bed (chr start end signal)

Output: 4-column BED with discretized integer bin values.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyBigWig


def assign_bins(signal_series: pd.Series, n_bins: int) -> pd.Series:
    """Global quantile-based discretization.
    Bin 0: zero signal + lowest non-zero quantile (merged).
    Bin 1..(n_bins-1): remaining quantile bins of non-zero signal.
    """
    bins = pd.Series(0, index=signal_series.index, dtype=int)
    nonzero = signal_series[signal_series > 0]
    if nonzero.empty:
        return bins
    K = min(n_bins, nonzero.nunique())
    edges = np.quantile(nonzero, np.linspace(0, 1, K + 1))
    bin_ids = np.digitize(nonzero, edges[1:-1], right=False)
    bin_ids = np.clip(bin_ids, 0, n_bins - 1)
    bins.loc[nonzero.index] = bin_ids
    return bins


def extract_from_bigwig(bw_path: str, bed_df: pd.DataFrame) -> pd.DataFrame:
    """Extract mean signal per bin from bigWig."""
    bw = pyBigWig.open(bw_path)
    if bw is None:
        raise RuntimeError(f"Cannot open bigWig: {bw_path}")

    vals = []
    for _, row in bed_df.iterrows():
        chrom, start, end = row["chr"], int(row["start"]), int(row["end"])
        try:
            mean_val = bw.stats(chrom, start, end, type="mean", exact=True)
            val = mean_val[0] if mean_val and mean_val[0] is not None else 0.0
        except Exception:
            val = 0.0
        vals.append(val)

    bw.close()
    result = bed_df[["chr", "start", "end"]].copy()
    result["signal"] = vals
    return result


def load_bed_signal(path: str) -> pd.DataFrame:
    """Load 4-column BED signal file."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "start", "end", "signal"])
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Discretize genomic signals (bigWig or 4-col BED) into quantile bins."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input signal file (.bigWig/.bw or 4-column BED)")
    parser.add_argument("-b", "--bins", default=None,
                        help="BED file with bin coordinates (required for bigWig input)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output 4-column BED file with discretized values")
    parser.add_argument("-n", "--n-bins", type=int, default=5,
                        help="Number of discrete bins (default: 5)")
    parser.add_argument("--keep-autosomes", action="store_true",
                        help="Keep only autosomes (chr1-chr22)")
    args = parser.parse_args()

    in_path = Path(args.input)

    # ------------------------------------------------------------------
    # Load raw signals
    # ------------------------------------------------------------------
    if in_path.suffix.lower() in (".bigwig", ".bw"):
        if args.bins is None:
            print("Error: --bins is required when input is bigWig.", file=sys.stderr)
            sys.exit(1)
        bed_df = pd.read_csv(args.bins, sep="\t", header=None,
                             names=["chr", "start", "end"])
        print(f"Loading bigWig: {in_path}")
        print(f"Bins: {len(bed_df)} regions from {args.bins}")
        df = extract_from_bigwig(str(in_path), bed_df)
    else:
        print(f"Loading BED signal: {in_path}")
        df = load_bed_signal(str(in_path))

    # ------------------------------------------------------------------
    # Optional autosome filter
    # ------------------------------------------------------------------
    if args.keep_autosomes:
        n_before = len(df)
        df = df[df["chr"].str.match(r"^chr\d+$")].copy()
        print(f"Autosome filter: {n_before} -> {len(df)} bins")

    # ------------------------------------------------------------------
    # Discretize
    # ------------------------------------------------------------------
    print(f"Discretizing with n_bins={args.n_bins} ...")
    df["discrete"] = assign_bins(df["signal"], args.n_bins)

    # Print distribution
    dist = df["discrete"].value_counts().sort_index()
    print("Bin distribution:")
    for b, c in dist.items():
        print(f"  bin {b}: {c} ({c/len(df)*100:.1f}%)")

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    out_df = df[["chr", "start", "end", "discrete"]].copy()
    out_df.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Saved: {args.output} ({len(out_df)} bins)")


if __name__ == "__main__":
    main()
