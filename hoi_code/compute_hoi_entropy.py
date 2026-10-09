#!/usr/bin/env python3
"""
HOI (Higher-Order Interactions) Analysis Script V3

Modular version with clear separation of concerns:
- signal_loader.py: Genome signal loading and binning
- block_processor.py: Block definitions and Cartesian product

This script combines both modules for end-to-end HOI analysis.
"""

import argparse
import os
import sys
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from itertools import combinations
from multiprocessing import Pool, cpu_count, shared_memory
from functools import partial

import numpy as np
import pandas as pd

# Import local modules
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

# HOI imports
try:
    from hoi.metrics import Oinfo, TC, DTC
    from hoi.core.entropies import get_entropy
    import jax
    import jax.numpy as jnp
except ImportError as e:
    print(f"Error: hoi package is required. Please activate hoi_env environment.")
    print(f"Import error: {e}")
    sys.exit(1)


# ============================================================================
# Configuration
# ============================================================================

DEFAULT_N_BINS = 3
DEFAULT_MAX_MISSING = 3
DEFAULT_N_MAX = 100000
DEFAULT_RANDOM_SEED = 42
DEFAULT_N_JOBS = -1  # -1 means use all available CPUs


# ============================================================================
# HOI Computation
# ============================================================================

class SharedMemoryManager:
    """Manages shared memory for signal matrix (DataFrame) across processes."""
    
    def __init__(self, df: pd.DataFrame):
        """
        Create shared memory from pandas DataFrame.
        
        Parameters
        ----------
        df : pd.DataFrame
            The DataFrame to share across processes
        """
        # Store DataFrame metadata
        self.index = df.index
        self.columns = df.columns
        
        # Convert to numpy array for sharing
        array = df.values
        self.shape = array.shape
        self.dtype = array.dtype
        self.nbytes = array.nbytes
        
        # Create shared memory
        self.shm = shared_memory.SharedMemory(create=True, size=self.nbytes)
        
        # Copy data to shared memory
        shm_array = np.ndarray(self.shape, dtype=self.dtype, buffer=self.shm.buf)
        shm_array[:] = array[:]
        
        self.name = self.shm.name
        print(f"  Shared memory created: {self.name}, size: {self.nbytes / 1024 / 1024:.1f} MB")
    
    def get_dataframe(self) -> pd.DataFrame:
        """Get DataFrame from shared memory (for child processes)."""
        existing_shm = shared_memory.SharedMemory(name=self.name)
        array = np.ndarray(self.shape, dtype=self.dtype, buffer=existing_shm.buf)
        return pd.DataFrame(array, index=self.index, columns=self.columns)
    
    def cleanup(self):
        """Clean up shared memory (call from parent process)."""
        self.shm.close()
        self.shm.unlink()
        print(f"  Shared memory cleaned up: {self.name}")


def process_single_block(
    block_idx: int,
    block_label: str,
    bin_ids: np.ndarray,
    shm_name: Optional[str],
    shm_shape: tuple,
    shm_dtype: np.dtype,
    shm_index: Optional[list],
    shm_columns: Optional[list],
    feature_names: List[str],
    sample_name: str,
    output_dir: str,
    n_bins: int,
    n_max: Optional[int],
    random_seed: int,
    max_missing: int,
    method: str,
    resume: bool,
    verbose: bool
) -> Optional[Dict]:
    """
    Process a single block - designed to run in a separate process.
    Uses shared memory to access signal matrix.
    
    Parameters
    ----------
    shm_name : str or None
        Shared memory name. If None, shm_shape is the actual DataFrame.
    shm_shape : tuple
        Shape of shared memory array
    shm_dtype : np.dtype
        Dtype of shared memory array
    shm_index : list or None
        DataFrame index (list of bin_ids)
    shm_columns : list or None
        DataFrame columns (feature names)
    
    Returns a dict with results if successful, None otherwise.
    """
    # Attach to shared memory (or use direct DataFrame for single process mode)
    if shm_name is None:
        # Single process mode: shm_shape is actually the DataFrame
        signal_matrix = shm_shape
    else:
        # Multi-process mode: reconstruct DataFrame from shared memory
        try:
            existing_shm = shared_memory.SharedMemory(name=shm_name)
            array = np.ndarray(shm_shape, dtype=shm_dtype, buffer=existing_shm.buf)
            signal_matrix = pd.DataFrame(array, index=shm_index, columns=shm_columns)
        except Exception as e:
            print(f"Error attaching to shared memory for block {block_label}: {e}")
            return None
    
    # Check if block already completed (for resume mode)
    block_output_dir = os.path.join(output_dir, block_label.replace('/', '_'))
    if resume:
        summary_file = os.path.join(
            block_output_dir, 
            f"{sample_name}_{block_label.replace('/', '_').replace('&', '_AND_')}_summary.json"
        )
        if os.path.exists(summary_file):
            if verbose:
                print(f"  [RESUME] Block already completed: {block_label}")
            return None
    
    if len(bin_ids) == 0:
        if verbose:
            print(f"  Skipping empty block: {block_label}")
        return None
    
    # Extract data for this block
    df_block = extract_block_data(
        signal_matrix=signal_matrix,
        bin_ids=bin_ids,
        block_label=block_label,
        verbose=verbose
    )
    
    if df_block is None:
        if verbose:
            print(f"  No valid data for block: {block_label}")
        return None
    
    # Sample if needed
    if n_max is not None and len(df_block) > n_max:
        df_block = sample_block_data(df_block, n_max, random_seed)
        if verbose:
            print(f"  Sampled to {len(df_block)} samples (n_max={n_max})")
    
    # Prepare data matrix for HOI
    data = df_block.values
    if method != 'binning':
        data = data.astype(float)
    
    if verbose:
        print(f"  Final data shape: {data.shape}")
    
    # Run entropy analysis for all subsets
    try:
        entropy_df = compute_all_subset_entropies(
            data=data,
            feature_names=feature_names,
            method=method,
            verbose=verbose
        )
    except Exception as e:
        if verbose:
            print(f"  Error computing entropy for block {block_label}: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Add block label to results
    entropy_df['block'] = block_label
    
    # Save intermediate results
    os.makedirs(block_output_dir, exist_ok=True)
    
    metadata = {
        'sample_name': sample_name,
        'block_label': block_label,
        'n_bins': n_bins,
        'n_max': n_max,
        'random_seed': random_seed,
        'max_missing': max_missing,
        'method': method,
        'n_samples': len(data),
    }
    
    save_results_v3(
        block_output_dir,
        f"{sample_name}_{block_label.replace('/', '_').replace('&', '_AND_')}",
        {0: entropy_df},
        feature_names,
        metadata
    )
    
    return {
        'block_label': block_label,
        'results': {0: entropy_df},
        'feature_names': feature_names
    }



def compute_all_subset_entropies(
    data: np.ndarray,
    feature_names: List[str],
    method: str = 'binning',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Compute joint entropy H(S) for all non-empty subsets S of features.
    
    For n features, computes:
    - Order 1: H(1), H(2), ..., H(n)
    - Order 2: H(1,2), H(1,3), ..., H(n-1,n)
    - ...
    - Order n: H(1,2,...,n)
    """
    n_features = data.shape[1]
    rows = []
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing all subset joint entropies")
        print(f"  n_features={n_features}")
        print(f"{'='*60}")
    
    entropy_fn = None if method == 'binning' else get_entropy(method=method)
    
    for r in range(1, n_features + 1):
        subsets_r = list(combinations(range(n_features), r))
        n_combs = len(subsets_r)
        if verbose:
            print(f"\n[Order {r}] Computing {n_combs} subsets")
        
        for subset_indices in subsets_r:
            subset_indices = list(subset_indices)
            data_subset = data[:, subset_indices]
            
            # Compute joint entropy
            if method == 'binning':
                # For discrete data, use frequency counting (exact)
                if r == 1:
                    values, counts = np.unique(data_subset[:, 0], return_counts=True)
                else:
                    values, counts = np.unique(data_subset, axis=0, return_counts=True)
                probs = counts / counts.sum()
                entropy_val = float(-np.sum(probs * np.log2(probs + 1e-12)))
            else:
                # For continuous methods, use hoi's get_entropy
                data_t = data_subset.transpose(1, 0)  # (r, n_samples)
                entropy_val = float(entropy_fn(data_t))
            
            subset_names = [feature_names[i] for i in subset_indices]
            
            rows.append({
                'order': r,
                'subset_indices': subset_indices,
                'subset_names': '/'.join(subset_names),
                'entropy': entropy_val,
            })
            
            if verbose and len(rows) <= 3 and r <= 2:
                print(f"  H({'/'.join(subset_names)}) = {entropy_val:.6f}")
        
        if verbose:
            print(f"  Completed {n_combs} subsets of order {r}")
    
    df = pd.DataFrame(rows)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Total subsets computed: {len(df)}")
        print(f"{'='*60}")
    
    return df


def compute_oinfo_missing_variables(
    data: np.ndarray,
    feature_names: List[str],
    max_missing: int = DEFAULT_MAX_MISSING,
    method: str = 'gc',
    verbose: bool = True
) -> Dict[int, pd.DataFrame]:
    """
    Compute O-information for subsets with missing variables.
    
    For max_missing=3, computes:
    - order 0: Ω(Xⁿ) - full system
    - order 1: Ω(Xⁿ₋ᵢ) for each excluded variable i
    - order 2: Ω(Xⁿ₋ᵢ₋ⱼ) for each excluded pair (i,j)
    - order 3: Ω(Xⁿ₋ᵢ₋ⱼ₋ₖ) for each excluded triplet (i,j,k)
    """
    n_features = data.shape[1]
    results = {}
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Computing O-information with missing variables")
        print(f"  n_features={n_features}, max_missing={max_missing}")
        print(f"{'='*60}")
    
    # Order 0: Full system
    if verbose:
        print(f"\n[Order 0] Full system Ω(Xⁿ) - {n_features} variables")
    print("data.shape", data.shape)
    # Compute O-information for full system
    model_full = Oinfo(data, verbose=False)
    hoi_full = model_full.fit(minsize=n_features, maxsize=n_features, method=method)
    oinfo_full = float(hoi_full[0, 0]) if hoi_full.ndim > 1 else float(hoi_full[0])
    
    # Compute TC (Total Correlation) for full system
    model_tc = TC(data, verbose=False)
    hoi_tc = model_tc.fit(minsize=n_features, maxsize=n_features, method=method)
    tc_full = float(hoi_tc[0, 0]) if hoi_tc.ndim > 1 else float(hoi_tc[0])
    
    # Compute DTC (Dual Total Correlation) for full system
    model_dtc = DTC(data, verbose=False)
    hoi_dtc = model_dtc.fit(minsize=n_features, maxsize=n_features, method=method)
    dtc_full = float(hoi_dtc[0, 0]) if hoi_dtc.ndim > 1 else float(hoi_dtc[0])
    
    # Compute total entropy H(Xⁿ) for full system
    entropy_fn = jax.vmap(get_entropy(method=method))
    data_t = data.transpose(1, 0)  # (n_features, n_samples)
    if data_t.ndim == 2:
        data_t = data_t[:, jnp.newaxis, :]  # (n_features, 1, n_samples)
    total_entropy_full = float(entropy_fn(data_t)[0])
    
    df_full = pd.DataFrame({
        'order': [n_features],
        'n_missing': [0],
        'missing_indices': [[]],
        'missing_names': [''],
        'remaining_indices': [list(range(n_features))],
        'remaining_names': ['/'.join(feature_names)],
        'hoi': [oinfo_full],
        'tc': [tc_full],
        'dtc': [dtc_full],
        'entropy': [total_entropy_full],
        'interpretation': ['Redundancy' if oinfo_full > 0 else 'Synergy']
    })
    results[0] = df_full
    
    if verbose:
        print(f"  Ω(Xⁿ) = {oinfo_full:.6f} ({'Redundancy' if oinfo_full > 0 else 'Synergy'})")
        print(f"  TC(Xⁿ) = {tc_full:.6f}, DTC(Xⁿ) = {dtc_full:.6f}")
        print(f"  H(Xⁿ) = {total_entropy_full:.6f}")
    
    # Orders 1 to max_missing
    for n_missing in range(1, min(max_missing + 1, n_features)):
        if verbose:
            print(f"\n[Order {n_missing}] Missing {n_missing} variable(s)")
        
        rows = []
        for excluded_indices in combinations(range(n_features), n_missing):
            mask = np.ones(n_features, dtype=bool)
            mask[list(excluded_indices)] = False
            remaining_indices = np.where(mask)[0]
            
            if len(remaining_indices) < 2:
                continue
            
            data_subset = data[:, mask]
            n_remaining = data_subset.shape[1]
            
            print(f"  Excluding indices {excluded_indices} ({[feature_names[i] for i in excluded_indices]}), "
                  f"remaining {n_remaining} variables")
            print(f"  data_subset.shape: {data_subset.shape}")
            
            # Compute O-information for subset
            model = Oinfo(data_subset, verbose=False)
            hoi = model.fit(minsize=n_remaining, maxsize=n_remaining, method=method)
            oinfo_val = float(hoi[0, 0]) if hoi.ndim > 1 else float(hoi[0])
            
            # Compute TC for subset
            model_tc = TC(data_subset, verbose=False)
            hoi_tc = model_tc.fit(minsize=n_remaining, maxsize=n_remaining, method=method)
            tc_val = float(hoi_tc[0, 0]) if hoi_tc.ndim > 1 else float(hoi_tc[0])
            
            # Compute DTC for subset
            model_dtc = DTC(data_subset, verbose=False)
            hoi_dtc = model_dtc.fit(minsize=n_remaining, maxsize=n_remaining, method=method)
            dtc_val = float(hoi_dtc[0, 0]) if hoi_dtc.ndim > 1 else float(hoi_dtc[0])
            
            # Compute total entropy H(X^{n-k}) for subset
            data_subset_t = data_subset.transpose(1, 0)  # (n_remaining, n_samples)
            if data_subset_t.ndim == 2:
                data_subset_t = data_subset_t[:, jnp.newaxis, :]  # (n_remaining, 1, n_samples)
            total_entropy_val = float(entropy_fn(data_subset_t)[0])
            
            excluded_names = [feature_names[i] for i in excluded_indices]
            remaining_names = [feature_names[i] for i in remaining_indices]
            
            rows.append({
                'order': n_remaining,
                'n_missing': n_missing,
                'missing_indices': list(excluded_indices),
                'missing_names': '/'.join(excluded_names),
                'remaining_indices': list(remaining_indices),
                'remaining_names': '/'.join(remaining_names),
                'hoi': oinfo_val,
                'tc': tc_val,
                'dtc': dtc_val,
                'entropy': total_entropy_val,
                'interpretation': 'Redundancy' if oinfo_val > 0 else 'Synergy'
            })
        
        if rows:
            df_missing = pd.DataFrame(rows)
            results[n_missing] = df_missing
            
            if verbose:
                print(f"  Computed {len(rows)} combinations")
                print(f"  Mean Ω(Xⁿ₋{n_missing}) = {df_missing['hoi'].mean():.6f}")
    
    if verbose:
        print(f"\n{'='*60}")
        total_combinations = sum(len(df) for df in results.values())
        print(f"Total combinations computed: {total_combinations}")
        print(f"{'='*60}")
    
    return results


# ============================================================================
# Output Functions
# ============================================================================

def convert_to_serializable(obj):
    """Convert numpy types to Python native types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


def save_results_v3(
    output_dir: str,
    sample_name: str,
    oinfo_results: Dict[int, pd.DataFrame],
    feature_names: List[str],
    metadata: Dict
):
    """Save all results to output directory (Entropy version)."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Save each result set to separate file
    for n_missing, df in oinfo_results.items():
        if 'entropy' in df.columns and 'hoi' not in df.columns:
            filename = f"{sample_name}_entropy_all_subsets.csv"
        elif n_missing == 0:
            filename = f"{sample_name}_oinfo_overall.csv"
        else:
            filename = f"{sample_name}_oinfo_missing_{n_missing}.csv"
        
        filepath = os.path.join(output_dir, filename)
        df.to_csv(filepath, index=False)
        print(f"  - {filename} ({len(df)} records)")
    
    # Save summary JSON
    summary = {
        'feature_names': feature_names,
        'n_features': len(feature_names),
        'total_records': sum(len(df) for df in oinfo_results.values()),
    }
    
    # Entropy statistics by order
    subsets_by_order = {}
    for n_missing, df in oinfo_results.items():
        if 'order' in df.columns and 'entropy' in df.columns:
            for order, group in df.groupby('order'):
                order_key = f'order_{order}'
                subsets_by_order[order_key] = {
                    'n_subsets': len(group),
                    'mean_entropy': float(group['entropy'].mean()),
                    'min_entropy': float(group['entropy'].min()),
                    'max_entropy': float(group['entropy'].max()),
                }
    if subsets_by_order:
        summary['subsets_by_order'] = subsets_by_order
    
    # Overall statistics
    all_dfs = [df for df in oinfo_results.values()]
    if all_dfs:
        all_df = pd.concat(all_dfs, ignore_index=True)
        if 'entropy' in all_df.columns:
            summary['entropy'] = {
                'mean': float(all_df['entropy'].mean()),
                'min': float(all_df['entropy'].min()),
                'max': float(all_df['entropy'].max()),
            }
    
    summary_path = os.path.join(output_dir, f"{sample_name}_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(convert_to_serializable(summary), f, indent=2)
    print(f"  - {sample_name}_summary.json")
    
    # Save metadata
    metadata_path = os.path.join(output_dir, f"{sample_name}_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump({
            'sample_name': sample_name,
            'feature_names': feature_names,
            'n_features': len(feature_names),
            **metadata
        }, f, indent=2)
    print(f"  - {sample_name}_metadata.json")
    
    print(f"\nResults saved to: {output_dir}")


# ============================================================================
# Main Function
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Compute HOI with missing variable subsets (V3 - Modular)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--input-file", "-i", required=True,
        help="Path to input file containing list of signal files"
    )
    parser.add_argument(
        "--output-dir", "-o", required=True,
        help="Output directory for results"
    )
    parser.add_argument(
        "--sample-name", "-s", default=None,
        help="Sample name (default: derived from input filename)"
    )
    
    # Data parameters
    parser.add_argument(
        "--n-bins", type=int, default=DEFAULT_N_BINS,
        help="Number of bins for signal discretization"
    )
    parser.add_argument(
        "--n-max", type=int, default=DEFAULT_N_MAX,
        help="Maximum number of samples per block (None for all)"
    )
    parser.add_argument(
        "--random-seed", type=int, default=DEFAULT_RANDOM_SEED,
        help="Random seed for sampling"
    )
    parser.add_argument(
        "--block-files", nargs="+", default=None,
        help="Block annotation files for partitioning"
    )
    parser.add_argument(
        "--block-groups", type=str, default=None,
        help="Block groups string for Cartesian product"
    )
    parser.add_argument(
        "--group-separator", type=str, default="::BLOCKGROUP::",
        help="Separator between block groups"
    )
    
    # HOI parameters
    parser.add_argument(
        "--max-missing", type=int, default=DEFAULT_MAX_MISSING,
        help="Maximum number of variables to exclude (1 to max_missing)"
    )
    parser.add_argument(
        "--method", default='binning',
        choices=['gc', 'binning', 'knn', 'kernel', 'gauss'],
        help="Method for entropy estimation"
    )
    
    # Other options
    parser.add_argument(
        "--quiet", "-q", action='store_true',
        help="Suppress verbose output"
    )
    parser.add_argument(
        "--resume", action='store_true',
        help="Resume from previous run: skip blocks that already have output directories"
    )
    parser.add_argument(
        "--n-jobs", type=int, default=DEFAULT_N_JOBS,
        help="Number of parallel workers (-1 means use all CPUs, 1 means sequential)"
    )
    
    args = parser.parse_args()
    
    verbose = not args.quiet
    
    if args.sample_name is None:
        args.sample_name = Path(args.input_file).stem
    
    # Parse block groups
    block_files = parse_block_groups_from_args(
        args.block_files, 
        args.block_groups, 
        args.group_separator
    )
    
    if verbose:
        print(f"{'='*60}")
        print(f"HOI Analysis V3 (Modular) for: {args.sample_name}")
        print(f"{'='*60}")
        print(f"Input file: {args.input_file}")
        print(f"Output directory: {args.output_dir}")
        print(f"Block groups ({len(block_files)} groups):")
        for i, group in enumerate(block_files, 1):
            if len(group) == 1:
                print(f"  Group {i}: {os.path.basename(group[0])}")
            else:
                print(f"  Group {i} (mutually exclusive):")
                for f in group:
                    print(f"    - {os.path.basename(f)}")
        print(f"n_bins={args.n_bins}, n_max={args.n_max}")
        print(f"max_missing={args.max_missing}, method={args.method}")
        print(f"{'='*60}\n")
    
    # ================================================================
    # Part 1: Load and bin signals (signal_loader module)
    # ================================================================
    file_paths, feature_names = load_signals_from_list(args.input_file)
    
    if verbose:
        print(f"Features ({len(feature_names)}):")
        for name in feature_names:
            print(f"  - {name}")
    
    signal_matrix = build_signal_matrix(
        file_paths=file_paths,
        feature_names=feature_names,
        n_bins=args.n_bins,
        verbose=verbose
    )
    
    # ================================================================
    # Part 2: Process blocks (block_processor module)
    # ================================================================
    block_groups = [load_block_group(group) for group in block_files]
    block_combinations = cartesian_product_blocks(block_groups, all_bin_ids=set(signal_matrix.index))
    
    if verbose:
        print(f"\nGenerated {len(block_combinations)} block combinations:")
        for label, _ in block_combinations:
            print(f"  - {label}")
        print()
    
    # ================================================================
    # Part 3: Run HOI analysis for each block (Parallel processing with Shared Memory)
    # ================================================================
    
    # Determine number of workers
    n_jobs = args.n_jobs if args.n_jobs > 0 else cpu_count()
    if verbose:
        print(f"\nParallel processing: Using {n_jobs} workers out of {cpu_count()} CPUs")
    
    # Prepare arguments for each block (filter first)
    blocks_to_process = []
    for block_idx, (block_label, bin_ids) in enumerate(block_combinations):
        # Pre-check for resume mode to filter out completed blocks
        if args.resume:
            block_output_dir = os.path.join(args.output_dir, block_label.replace('/', '_'))
            summary_file = os.path.join(
                block_output_dir, 
                f"{args.sample_name}_{block_label.replace('/', '_').replace('&', '_AND_')}_summary.json"
            )
            if os.path.exists(summary_file):
                if verbose:
                    print(f"  [RESUME] Block already completed: {block_label}")
                continue
        
        if len(bin_ids) == 0:
            if verbose:
                print(f"  Skipping empty block: {block_label}")
            continue
        
        blocks_to_process.append((block_idx, block_label, bin_ids))
    
    all_results = []
    
    if len(blocks_to_process) == 0:
        if verbose:
            print("  No blocks to process (all completed or empty)")
    elif n_jobs == 1:
        # Sequential processing for single worker (no shared memory needed)
        if verbose:
            print(f"\nProcessing {len(blocks_to_process)} blocks sequentially...")
        for block_idx, block_label, bin_ids in blocks_to_process:
            # In single process mode, pass the actual DataFrame as shm_shape (5th param)
            task = (
                block_idx, block_label, bin_ids, None, signal_matrix, None, None, None,
                feature_names, args.sample_name, args.output_dir, args.n_bins, args.n_max, args.random_seed,
                args.max_missing, args.method, args.resume, verbose
            )
            result = process_single_block(*task)
            if result is not None:
                all_results.append(result)
    else:
        # Parallel processing with shared memory
        if verbose:
            print(f"\nProcessing {len(blocks_to_process)} blocks in parallel using shared memory...")
            print(f"  Signal matrix size: {signal_matrix.values.nbytes / 1024 / 1024:.1f} MB")
        
        # Create shared memory for signal matrix
        shm_manager = SharedMemoryManager(signal_matrix)
        
        try:
            # Prepare tasks with shared memory info
            block_tasks = []
            for block_idx, block_label, bin_ids in blocks_to_process:
                task = (
                    block_idx, block_label, bin_ids, shm_manager.name, 
                    shm_manager.shape, shm_manager.dtype,
                    list(shm_manager.index), list(shm_manager.columns),
                    feature_names, args.sample_name, args.output_dir, args.n_bins, args.n_max, args.random_seed,
                    args.max_missing, args.method, args.resume, verbose
                )
                block_tasks.append(task)
            
            with Pool(processes=n_jobs) as pool:
                results = pool.starmap(process_single_block, block_tasks)
                all_results = [r for r in results if r is not None]
            
            if verbose:
                print(f"  Completed {len(all_results)}/{len(blocks_to_process)} blocks successfully")
        finally:
            # Clean up shared memory
            shm_manager.cleanup()
            if verbose:
                print(f"  Shared memory cleaned up")
    
    # ================================================================
    # Part 4: Combine all results
    # ================================================================
    if all_results:
        combined_results = {}
        max_missing_global = max(max(r['results'].keys()) for r in all_results)
        
        for n_missing in range(max_missing_global + 1):
            dfs = []
            for r in all_results:
                if n_missing in r['results']:
                    dfs.append(r['results'][n_missing])
            if dfs:
                combined_results[n_missing] = pd.concat(dfs, ignore_index=True)
        
        # Save combined results
        os.makedirs(args.output_dir, exist_ok=True)
        combined_metadata = {
            'sample_name': args.sample_name,
            'n_blocks': len(all_results),
            'block_labels': [r['block_label'] for r in all_results],
            'n_bins': args.n_bins,
            'n_max': args.n_max,
            'max_missing': args.max_missing,
            'method': args.method,
        }
        
        save_results_v3(
            args.output_dir,
            args.sample_name,
            combined_results,
            all_results[0]['feature_names'],
            combined_metadata
        )
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Analysis complete! Processed {len(all_results)} blocks.")
        print(f"{'='*60}")


if __name__ == "__main__":
    main()
