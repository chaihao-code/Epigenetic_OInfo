#!/usr/bin/env python3
"""
Compute O-information from pre-computed joint entropy subsets.

Reads entropy_all_subsets.csv (output of compute_hoi_entropy.py) and
calculates O-information, TC, DTC for the full system and all missing-variable
subsets, matching the format of compute_oinfo_missing_variables().
"""

import argparse
import os
import sys
import json
from pathlib import Path
from itertools import combinations
from multiprocessing import Pool, cpu_count
from functools import partial

import numpy as np
import pandas as pd


def compute_oinfo_from_entropy_df(entropy_df: pd.DataFrame) -> dict:
    """
    Compute O-information results from an entropy_all_subsets DataFrame.
    
    Returns dict: {n_missing: DataFrame}
    """
    # Build lookup: subset_names -> entropy
    entropy_map = {}
    for _, row in entropy_df.iterrows():
        entropy_map[row['subset_names']] = float(row['entropy'])
    
    # Extract ordered feature list from order-1 rows
    order1 = entropy_df[entropy_df['order'] == 1].copy()
    # Parse subset_indices to get proper ordering
    def parse_indices(s):
        if isinstance(s, str):
            return eval(s)
        return list(s)
    order1['idx'] = order1['subset_indices'].apply(parse_indices)
    order1 = order1.sort_values('idx')
    feature_names = order1['subset_names'].tolist()
    n_features = len(feature_names)
    
    results = {}
    
    # Order 0: full system
    full_names = '/'.join(feature_names)
    h_full = entropy_map.get(full_names, 0.0)
    
    tc_full = sum(entropy_map.get(f, 0.0) for f in feature_names) - h_full
    dtc_full = sum(
        entropy_map.get('/'.join([f for f in feature_names if f != fj]), 0.0)
        for fj in feature_names
    ) - (n_features - 1) * h_full
    oinfo_full = tc_full - dtc_full
    
    # Match existing output format: entropy = marginal entropy of first remaining feature
    entropy_first_full = entropy_map.get(feature_names[0], 0.0)
    
    df_full = pd.DataFrame({
        'order': [n_features],
        'n_missing': [0],
        'missing_indices': [[]],
        'missing_names': [''],
        'remaining_indices': [list(range(n_features))],
        'remaining_names': [full_names],
        'hoi': [oinfo_full],
        'tc': [tc_full],
        'dtc': [dtc_full],
        'entropy': [entropy_first_full],
        'interpretation': ['Redundancy' if oinfo_full > 0 else 'Synergy']
    })
    results[0] = df_full
    
    # Orders 1 to n-2 (need at least 2 remaining variables for non-trivial O-info)
    max_missing = n_features - 2
    for n_missing in range(1, max_missing + 1):
        rows = []
        for excluded_indices in combinations(range(n_features), n_missing):
            remaining_indices = [i for i in range(n_features) if i not in excluded_indices]
            n_remaining = len(remaining_indices)
            if n_remaining < 2:
                continue
            
            remaining_features = [feature_names[i] for i in remaining_indices]
            remaining_names = '/'.join(remaining_features)
            h_remaining = entropy_map.get(remaining_names, 0.0)
            
            tc_val = sum(entropy_map.get(f, 0.0) for f in remaining_features) - h_remaining
            dtc_val = sum(
                entropy_map.get('/'.join([f for f in remaining_features if f != fj]), 0.0)
                for fj in remaining_features
            ) - (n_remaining - 1) * h_remaining
            oinfo_val = tc_val - dtc_val
            
            excluded_names = [feature_names[i] for i in excluded_indices]
            
            # Match existing output format: entropy = marginal entropy of first remaining feature
            entropy_first = entropy_map.get(remaining_features[0], 0.0)
            
            rows.append({
                'order': n_remaining,
                'n_missing': n_missing,
                'missing_indices': list(excluded_indices),
                'missing_names': '/'.join(excluded_names),
                'remaining_indices': remaining_indices,
                'remaining_names': remaining_names,
                'hoi': oinfo_val,
                'tc': tc_val,
                'dtc': dtc_val,
                'entropy': entropy_first,
                'interpretation': 'Redundancy' if oinfo_val > 0 else 'Synergy'
            })
        
        if rows:
            results[n_missing] = pd.DataFrame(rows)
    
    return results


def process_single_entropy_file(entropy_csv_path: str, output_base_dir: str, verbose: bool = False) -> bool:
    """
    Process a single entropy_all_subsets.csv file and write O-info results.
    
    Path structure:
      input:  .../{cell_type}/{block}/{sample}_{block}_entropy_all_subsets.csv
      output: .../{cell_type}/{block}/{sample}_{block}_oinfo_*.csv
    """
    entropy_csv_path = Path(entropy_csv_path)
    if not entropy_csv_path.exists():
        return False
    
    # Derive sample name and block from path
    block_dir = entropy_csv_path.parent
    cell_type_dir = block_dir.parent
    
    sample_name = cell_type_dir.name
    block_label = block_dir.name
    
    stem = entropy_csv_path.stem  # e.g. E13.5_Liver_proximal_ctcf_bound_enhancer_AND_None_entropy_all_subsets
    # Remove '_entropy_all_subsets' suffix to get prefix
    if '_entropy_all_subsets' in stem:
        prefix = stem.replace('_entropy_all_subsets', '')
    else:
        prefix = stem
    
    output_dir = Path(output_base_dir) / cell_type_dir.name / block_label
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if already done (look for oinfo_overall.csv)
    overall_file = output_dir / f"{prefix}_oinfo_overall.csv"
    if overall_file.exists():
        if verbose:
            print(f"  [SKIP] Already exists: {overall_file}")
        return True
    
    try:
        entropy_df = pd.read_csv(entropy_csv_path)
        if entropy_df.empty:
            return False
        
        results = compute_oinfo_from_entropy_df(entropy_df)
        
        for n_missing, df in results.items():
            if n_missing == 0:
                out_name = f"{prefix}_oinfo_overall.csv"
            else:
                out_name = f"{prefix}_oinfo_missing_{n_missing}.csv"
            df.to_csv(output_dir / out_name, index=False)
        
        if verbose:
            print(f"  [DONE] {prefix} (n_features={entropy_df[entropy_df['order']==1].shape[0]})")
        return True
    except Exception as e:
        if verbose:
            print(f"  [ERROR] {entropy_csv_path}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Compute O-information from entropy subset results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir", "-i",
        help="Root directory containing entropy results (e.g. results_entropy_n_bins_5). Mutually exclusive with --input-file."
    )
    parser.add_argument(
        "--input-file", "-f",
        help="Single entropy_all_subsets.csv file to process. Mutually exclusive with --input-dir."
    )
    parser.add_argument(
        "--output-dir", "-o", required=True,
        help="Root directory for O-info outputs (e.g. results_n_bins_5)"
    )
    parser.add_argument(
        "--n-jobs", type=int, default=-1,
        help="Number of parallel workers (-1 = all CPUs)"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress progress output"
    )
    args = parser.parse_args()
    
    if not args.input_dir and not args.input_file:
        parser.error("Either --input-dir or --input-file must be provided.")
    if args.input_dir and args.input_file:
        parser.error("--input-dir and --input-file are mutually exclusive.")
    
    verbose = not args.quiet
    
    if args.input_file:
        # Single file mode
        entropy_files = [Path(args.input_file)]
        if not entropy_files[0].exists():
            print(f"File not found: {args.input_file}")
            sys.exit(1)
        if verbose:
            print(f"Processing single file: {args.input_file}")
        n_jobs = 1
    else:
        # Batch directory mode
        input_dir = Path(args.input_dir)
        entropy_files = sorted(input_dir.rglob("*_entropy_all_subsets.csv"))
        if not entropy_files:
            print(f"No entropy_all_subsets.csv files found in {input_dir}")
            sys.exit(1)
        if verbose:
            print(f"Found {len(entropy_files)} entropy files to process")
        n_jobs = args.n_jobs if args.n_jobs > 0 else cpu_count()
    
    process_fn = partial(process_single_entropy_file, output_base_dir=args.output_dir, verbose=verbose)
    
    if n_jobs == 1:
        successes = [process_fn(f) for f in entropy_files]
    else:
        with Pool(processes=n_jobs) as pool:
            successes = pool.map(process_fn, [str(f) for f in entropy_files])
    
    n_success = sum(successes)
    if verbose:
        print(f"\nCompleted: {n_success}/{len(entropy_files)} files processed successfully")


if __name__ == "__main__":
    main()
