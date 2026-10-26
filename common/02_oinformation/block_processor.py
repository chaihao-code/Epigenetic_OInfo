"""
Block Processor Module

Part 2: Handle block definitions and Cartesian product.
Takes a signal matrix and produces block-specific data subsets.
"""

from itertools import product
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np


def load_block_group(file_list: List[str]) -> Dict[str, set]:
    """
    Load a group of block files, return dict of {block_name: set of bin_ids}
    
    Parameters
    ----------
    file_list : List[str]
        List of block file paths
    
    Returns
    -------
    Dict[str, set]
        Dictionary mapping block names to sets of bin_ids
    """
    block_dict = {}
    for bf in file_list:
        df_block = pd.read_csv(bf, sep="\t")
        block_name = df_block.columns[1]
        bin_ids = set(df_block[df_block[block_name] == 1]["bin_id"])
        block_dict[block_name] = bin_ids
    return block_dict


def cartesian_product_blocks(
    block_groups: List[Dict[str, set]],
    all_bin_ids: Optional[set] = None
) -> List[Tuple[str, set]]:
    """
    Compute Cartesian product of block groups.
    
    Within each group, blocks are mutually exclusive.
    Across groups, compute all combinations.
    
    Parameters
    ----------
    block_groups : List[Dict[str, set]]
        List of block groups (each group is a dict of {block_name: bin_ids})
    all_bin_ids : Optional[set]
        If provided, use this as the universal set of bin_ids (e.g. all bins
        in the signal matrix). Otherwise defaults to the union of all flag=1
        bins across block_groups.
    
    Returns
    -------
    List[Tuple[str, set]]
        List of (label, bin_ids) tuples for each combination
    """
    if not block_groups:
        return []
    
    # Compute union of all bin_ids
    if all_bin_ids is None:
        all_bin_ids = set()
        for group in block_groups:
            for bin_ids in group.values():
                all_bin_ids |= bin_ids
    
    if not all_bin_ids:
        return []
    
    # Add "None" option for each group (bins not covered by any block in the group)
    group_with_none = []
    for group in block_groups:
        group_bins = set()
        for bin_ids in group.values():
            group_bins |= bin_ids
        none_bins = all_bin_ids - group_bins
        
        extended_group = dict(group)
        if none_bins:
            extended_group["None"] = none_bins
        group_with_none.append(extended_group)
    
    # Compute Cartesian product
    result = []
    group_names = [list(g.keys()) for g in group_with_none]
    
    for combo in product(*group_names):
        bin_sets = [group_with_none[i][combo[i]] for i in range(len(combo))]
        intersect_bins = all_bin_ids.copy()
        for bin_set in bin_sets:
            intersect_bins &= bin_set
        
        if intersect_bins:
            label = "_AND_".join(combo)
            result.append((label, intersect_bins))
    
    return result


def extract_block_data(
    signal_matrix: pd.DataFrame,
    bin_ids: set,
    block_label: str,
    verbose: bool = True
) -> Optional[pd.DataFrame]:
    """
    Extract data for a specific block from the signal matrix.
    
    Parameters
    ----------
    signal_matrix : pd.DataFrame
        Wide-format DataFrame with index=bin_id
    bin_ids : set
        Set of bin_ids for this block
    block_label : str
        Label for this block (for logging)
    verbose : bool
        Print progress
        
    Returns
    -------
    pd.DataFrame or None
        DataFrame with block data, or None if no valid data
    """
    # Filter to bins in this block
    df_block = signal_matrix.loc[signal_matrix.index.isin(bin_ids)]
    
    # Check for missing data
    if df_block.isnull().any().any():
        missing_signals = df_block.columns[df_block.isnull().any()].tolist()
        if verbose:
            print(f"  Warning: Missing signals in block '{block_label}': {missing_signals}")
        return None
    
    if verbose:
        print(f"  Block '{block_label}': {len(df_block)} bins")
    
    return df_block


def process_all_blocks(
    signal_matrix: pd.DataFrame,
    block_combinations: List[Tuple[str, set]],
    verbose: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    Process all block combinations and return a dictionary of block data.
    
    Parameters
    ----------
    signal_matrix : pd.DataFrame
        Wide-format DataFrame with index=bin_id
    block_combinations : List[Tuple[str, set]]
        List of (label, bin_ids) from cartesian_product_blocks
    verbose : bool
        Print progress
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary mapping block labels to their DataFrames
    """
    results = {}
    
    if verbose:
        print(f"\nProcessing {len(block_combinations)} block combinations...")
    
    for label, bin_ids in block_combinations:
        if len(bin_ids) == 0:
            if verbose:
                print(f"  Skipping empty block: {label}")
            continue
        
        df_block = extract_block_data(signal_matrix, bin_ids, label, verbose)
        if df_block is not None:
            results[label] = df_block
    
    if verbose:
        print(f"Successfully processed {len(results)} blocks")
    
    return results


def sample_block_data(
    df: pd.DataFrame,
    n_max: Optional[int],
    random_state: Optional[int] = None
) -> pd.DataFrame:
    """
    Sample data if needed.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    n_max : Optional[int]
        Maximum samples to keep
    random_state : Optional[int]
        Random seed
        
    Returns
    -------
    pd.DataFrame
        Sampled DataFrame
    """
    if n_max is None or len(df) <= n_max:
        return df
    
    return df.sample(n=n_max, random_state=random_state)


def parse_block_groups_from_args(
    block_files: Optional[List[str]] = None,
    block_groups_str: Optional[str] = None,
    group_separator: str = "::BLOCKGROUP::"
) -> List[List[str]]:
    """
    Parse block files/groups from command line arguments.
    
    Parameters
    ----------
    block_files : Optional[List[str]]
        List of block files (mutually exclusive)
    block_groups_str : Optional[str]
        Block groups string for Cartesian product
    group_separator : str
        Separator between block groups
        
    Returns
    -------
    List[List[str]]
        List of block file groups
    """
    if block_groups_str is not None:
        group_strings = block_groups_str.split(group_separator)
        block_files = []
        for group_str in group_strings:
            if group_str.strip():
                files = [f.strip() for f in group_str.split(';') if f.strip()]
                block_files.append(files)
        return block_files
    elif block_files is not None:
        return [[f] for f in block_files]
    else:
        raise ValueError("Either block_files or block_groups_str must be specified")


def main():
    """CLI for testing block processing independently."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Process blocks for signal matrix")
    parser.add_argument("-s", "--signal-matrix", required=True, 
                        help="Input signal matrix (parquet)")
    parser.add_argument("--block-files", nargs="+", 
                        help="Block annotation files")
    parser.add_argument("--block-groups", type=str, default=None,
                        help="Block groups string for Cartesian product")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for block data")
    parser.add_argument("-n", "--n-max", type=int, default=None,
                        help="Max samples per block")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Suppress output")
    
    args = parser.parse_args()
    
    # Load signal matrix
    signal_matrix = pd.read_parquet(args.signal_matrix)
    
    # Parse block groups
    block_groups_list = parse_block_groups_from_args(
        args.block_files, args.block_groups
    )
    
    # Load and compute Cartesian product
    block_groups = [load_block_group(group) for group in block_groups_list]
    block_combinations = cartesian_product_blocks(block_groups)
    
    # Process blocks
    results = process_all_blocks(
        signal_matrix, block_combinations, verbose=not args.quiet
    )
    
    # Save results
    import os
    os.makedirs(args.output_dir, exist_ok=True)
    
    for label, df in results.items():
        # Sample if needed
        df_sampled = sample_block_data(df, args.n_max, args.seed)
        
        # Save
        safe_label = label.replace('/', '_').replace('&', '_AND_')
        output_path = os.path.join(args.output_dir, f"{safe_label}.parquet")
        df_sampled.to_parquet(output_path)
        
        if not args.quiet:
            print(f"  Saved: {safe_label} ({len(df_sampled)} rows)")
    
    if not args.quiet:
        print(f"\nAll blocks saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
