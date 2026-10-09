#!/usr/bin/env python3
"""
Compute O-information from combined entropy_all_subsets.csv files.

Each cell type has a single summary file:
    {cell_type}/{cell_type}_entropy_all_subsets.csv

This file contains entropy subsets for multiple blocks, distinguished by the
`block` column. The script groups rows by `block` and computes O-information
for each block separately, producing the same output format as
compute_oinfo_from_entropy.py.

Output layout:
    {output_dir}/{cell_type}/{block}/{cell_type}_{block}_oinfo_overall.csv
    {output_dir}/{cell_type}/{block}/{cell_type}_{block}_oinfo_missing_1.csv
    ...
"""

import argparse
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

import pandas as pd

from compute_oinfo_from_entropy import compute_oinfo_from_entropy_df


def process_combined_entropy_file(
    combined_csv_path: Path,
    output_base_dir: Path,
    verbose: bool = True
) -> bool:
    """
    Process one cell-type combined entropy file and write per-block O-info.

    Parameters
    ----------
    combined_csv_path : Path
        Path to {cell_type}/{cell_type}_entropy_all_subsets.csv.
    output_base_dir : Path
        Root directory for O-info outputs.
    verbose : bool
        Print progress.

    Returns
    -------
    bool
        True if at least one block was processed successfully.
    """
    combined_csv_path = Path(combined_csv_path)
    if not combined_csv_path.exists():
        if verbose:
            print(f"[ERROR] File not found: {combined_csv_path}", file=sys.stderr)
        return False

    # cell_type is the parent directory name (e.g., "4cell")
    cell_type = combined_csv_path.parent.name

    try:
        entropy_df = pd.read_csv(combined_csv_path)
    except Exception as e:
        if verbose:
            print(f"[ERROR] Reading {combined_csv_path}: {e}", file=sys.stderr)
        return False

    if entropy_df.empty:
        if verbose:
            print(f"[ERROR] Empty file: {combined_csv_path}", file=sys.stderr)
        return False

    if "block" not in entropy_df.columns:
        if verbose:
            print(
                f"[ERROR] Missing 'block' column in {combined_csv_path}",
                file=sys.stderr,
            )
        return False

    blocks = entropy_df["block"].unique()
    if verbose:
        print(f"\n[CELL] {cell_type}: {len(blocks)} blocks")

    n_success = 0
    for block_label in blocks:
        block_df = entropy_df[entropy_df["block"] == block_label].copy()

        if block_df.empty:
            continue

        # Drop the block column before passing to the entropy->O-info function
        block_df = block_df.drop(columns=["block"])

        output_dir = output_base_dir / cell_type / block_label
        output_dir.mkdir(parents=True, exist_ok=True)

        prefix = f"{cell_type}_{block_label}"
        overall_file = output_dir / f"{prefix}_oinfo_overall.csv"

        if overall_file.exists():
            if verbose:
                print(f"  [SKIP] {cell_type}/{block_label}: already exists")
            n_success += 1
            continue

        try:
            results = compute_oinfo_from_entropy_df(block_df)
        except Exception as e:
            if verbose:
                print(
                    f"  [ERROR] {cell_type}/{block_label}: {e}",
                    file=sys.stderr,
                )
            continue

        for n_missing, df in results.items():
            if n_missing == 0:
                out_name = f"{prefix}_oinfo_overall.csv"
            else:
                out_name = f"{prefix}_oinfo_missing_{n_missing}.csv"
            df.to_csv(output_dir / out_name, index=False)

        if verbose:
            n_features = block_df[block_df["order"] == 1].shape[0]
            print(f"  [DONE] {cell_type}/{block_label} (n_features={n_features})")
        n_success += 1

    return n_success > 0


def main():
    parser = argparse.ArgumentParser(
        description="Compute O-information from combined entropy_all_subsets.csv files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir", "-i",
        required=True,
        help="Root directory containing cell-type result folders (e.g. results_entropy_n_bins_5).",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Root directory for O-info outputs (e.g. results_oinfo_n_bins_5).",
    )
    parser.add_argument(
        "--n-jobs", "-j",
        type=int,
        default=-1,
        help="Number of parallel workers (-1 = all CPUs).",
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    verbose = not args.quiet

    # Find combined entropy files: {cell_type}/{cell_type}_entropy_all_subsets.csv
    entropy_files = sorted(input_dir.rglob("*_entropy_all_subsets.csv"))
    combined_files = [
        f for f in entropy_files
        if f.parent.name == f.stem.replace("_entropy_all_subsets", "")
    ]

    if not combined_files:
        print(f"No combined entropy files found in {input_dir}")
        sys.exit(1)

    if verbose:
        print(f"Found {len(combined_files)} combined entropy files")

    n_jobs = args.n_jobs if args.n_jobs > 0 else cpu_count()
    process_fn = partial(
        process_combined_entropy_file,
        output_base_dir=output_dir,
        verbose=verbose,
    )

    if n_jobs == 1:
        successes = [process_fn(f) for f in combined_files]
    else:
        with Pool(processes=n_jobs) as pool:
            successes = pool.map(process_fn, [str(f) for f in combined_files])

    n_success = sum(bool(s) for s in successes)
    if verbose:
        print(f"\nCompleted: {n_success}/{len(combined_files)} cell types processed successfully")


if __name__ == "__main__":
    main()
