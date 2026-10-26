#!/usr/bin/env python3
"""
Compute pairwise Pearson correlation between genomic signals.

Simplified from compute_hoi_entropy.py:
- Loads signal matrix using signal_loader (discretized bins) or raw continuous signals
- Optionally processes blocks via block_processor (Cartesian product)
- Outputs 3-column TSV: feature1 feature2 correlation
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from signal_loader import (
    load_signals_from_list,
    build_signal_matrix,
)
from block_processor import (
    parse_block_groups_from_args,
    load_block_group,
    cartesian_product_blocks,
    extract_block_data,
    sample_block_data,
)


def load_continuous_signal_matrix(
    file_paths: List[str],
    feature_names: List[str],
    verbose: bool = True
) -> pd.DataFrame:
    """
    Load continuous (non-discretized) signals from .bp files.

    Parameters
    ----------
    file_paths : List[str]
        List of .bp file paths
    feature_names : List[str]
        Names for each signal column
    verbose : bool
        Print progress

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with index=bin_id, columns=feature_names
    """
    if verbose:
        print(f"\nLoading continuous signals (no binning)...")

    signal_series_list = []
    for path, name in zip(file_paths, feature_names):
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["chr", "start", "end", "signal"])
        # Keep only autosomes (consistent with signal_loader)
        df = df[df["chr"].str.match(r"^chr\d+$")]
        df["bin_id"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        signal_series = df.set_index("bin_id")["signal"].rename(name)
        signal_series_list.append(signal_series)

        if verbose:
            print(f"  Loaded: {name} ({len(df)} bins)")

    signal_matrix = pd.concat(signal_series_list, axis=1)

    if verbose:
        print(f"Signal matrix shape: {signal_matrix.shape}")
        print(f"  - Rows (bin_ids): {len(signal_matrix)}")
        print(f"  - Columns (signals): {len(signal_matrix.columns)}")

    return signal_matrix

def compute_pairwise_correlation(
    df: pd.DataFrame,
    method: str = "pearson",
    drop_na: bool = True
) -> pd.DataFrame:
    """
    Compute pairwise correlation and return long-format 3-column DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Signal matrix (samples x features)
    method : str
        Correlation method: pearson, spearman, kendall
    drop_na : bool
        Drop rows with NaN correlation (e.g., zero-variance features)

    Returns
    -------
    pd.DataFrame
        Columns: feature1, feature2, correlation
    """
    corr_matrix = df.corr(method=method)

    rows = []
    features = corr_matrix.columns.tolist()
    for i, feat1 in enumerate(features):
        for j, feat2 in enumerate(features):
            if i < j:  # upper triangle only, skip self-correlation
                rows.append({
                    "feature1": feat1,
                    "feature2": feat2,
                    "correlation": corr_matrix.iloc[i, j]
                })

    result = pd.DataFrame(rows)
    if drop_na:
        result = result.dropna(subset=["correlation"])
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compute pairwise signal correlations (3-column TSV output)."
    )
    parser.add_argument(
        "-i", "--input-file", required=True,
        help="Path to input file containing list of .bp signal file paths"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file path (3 columns: feature1 feature2 correlation)"
    )
    parser.add_argument(
        "-n", "--n-bins", type=int, default=5,
        help="Number of bins for signal discretization (default: 5)"
    )
    parser.add_argument(
        "--method", default="pearson",
        choices=["pearson", "spearman", "kendall"],
        help="Correlation method (default: pearson)"
    )
    parser.add_argument(
        "--continuous", action="store_true",
        help="Use continuous raw signals instead of discretized bins (default: discrete)"
    )
    parser.add_argument(
        "--block-files", nargs="+", default=None,
        help="Block annotation files (mutually exclusive)"
    )
    parser.add_argument(
        "--block-groups", type=str, default=None,
        help="Block groups string for Cartesian product (separated by --group-separator)"
    )
    parser.add_argument(
        "--group-separator", type=str, default="::BLOCKGROUP::",
        help="Separator between block groups (default: ::BLOCKGROUP::)"
    )
    parser.add_argument(
        "--per-block", action="store_true",
        help="Compute correlation per block combination (adds block column to output)"
    )
    parser.add_argument(
        "--n-max", type=int, default=None,
        help="Maximum number of samples per block (only used with --per-block)"
    )
    parser.add_argument(
        "--random-seed", type=int, default=42,
        help="Random seed for sampling (only used with --per-block)"
    )
    parser.add_argument(
        "--with-header", action="store_true",
        help="Include header row in output (default: no header)"
    )
    parser.add_argument(
        "-q", "--quiet", action="store_true",
        help="Suppress verbose output"
    )

    args = parser.parse_args()
    verbose = not args.quiet

    # Load signal matrix
    if verbose:
        print(f"Loading signals from: {args.input_file}")
    file_paths, feature_names = load_signals_from_list(args.input_file)
    if args.continuous:
        signal_matrix = load_continuous_signal_matrix(
            file_paths=file_paths,
            feature_names=feature_names,
            verbose=verbose
        )
    else:
        signal_matrix = build_signal_matrix(
            file_paths=file_paths,
            feature_names=feature_names,
            n_bins=args.n_bins,
            verbose=verbose
        )

    if args.per_block:
        # Parse block groups and compute Cartesian product
        if not args.block_groups and not args.block_files:
            print("Error: --block-groups or --block-files required for --per-block", file=sys.stderr)
            sys.exit(1)

        block_files = parse_block_groups_from_args(
            args.block_files, args.block_groups, args.group_separator
        )
        block_groups = [load_block_group(group) for group in block_files]
        block_combinations = cartesian_product_blocks(
            block_groups, all_bin_ids=set(signal_matrix.index)
        )

        if verbose:
            print(f"\nProcessing {len(block_combinations)} block combinations...")

        all_results = []
        for block_label, bin_ids in block_combinations:
            if len(bin_ids) == 0:
                continue

            df_block = extract_block_data(
                signal_matrix, bin_ids, block_label, verbose=verbose
            )
            if df_block is None or len(df_block) < 2:
                continue

            if args.n_max is not None and len(df_block) > args.n_max:
                df_block = sample_block_data(df_block, args.n_max, args.random_seed)
                if verbose:
                    print(f"  Sampled {block_label} to {len(df_block)} rows")

            corr_df = compute_pairwise_correlation(df_block, method=args.method)
            if corr_df.empty:
                continue
            corr_df.insert(0, "block", block_label)
            all_results.append(corr_df)

        if all_results:
            result_df = pd.concat(all_results, ignore_index=True)
        else:
            result_df = pd.DataFrame(
                columns=["block", "feature1", "feature2", "correlation"]
            )
    else:
        # Overall correlation across all bins
        if verbose:
            print(f"\nComputing overall {args.method} correlation across {len(signal_matrix)} bins...")
        result_df = compute_pairwise_correlation(signal_matrix, method=args.method)

    # Save output
    os.makedirs(Path(args.output).parent or Path("."), exist_ok=True)
    result_df.to_csv(
        args.output,
        sep="\t",
        index=False,
        header=args.with_header,
        float_format="%.6f"
    )

    if verbose:
        print(f"\nSaved {len(result_df)} correlations to: {args.output}")


if __name__ == "__main__":
    main()
