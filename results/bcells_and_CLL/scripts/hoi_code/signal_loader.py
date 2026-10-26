"""
Genome Signal Loader Module

Part 1: Load genomic signals and discretize into bins.
Outputs a wide-format DataFrame (rows=bin_ids, cols=signal_names).
"""

import os
from pathlib import Path
from typing import List, Tuple, Optional
import pandas as pd
import numpy as np


def load_signal(path: str) -> pd.DataFrame:
    """Load signal file (bed format) and return DataFrame with bin_id and signal columns."""
    df = pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "signal"])
    # Keep only autosomes
    df = df[df["chr"].str.match(r"^chr\d+$")]
    df["bin_id"] = df["chr"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    return df[["bin_id", "signal"]]


def parse_bin_id(bin_id: str) -> tuple:
    """Parse bin_id string 'chr:start-end' into (chr, start, end)."""
    chrom, coords = bin_id.split(":")
    start, end = map(int, coords.split("-"))
    return chrom, start, end


def compute_bin_center(chrom: str, start: int, end: int) -> tuple:
    """Compute center position of a bin."""
    return (chrom, (start + end) // 2)


def map_coarse_to_fine(
    coarse_df: pd.DataFrame,
    fine_bin_ids: pd.Index,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Map coarse-grained signal to fine-grained bins using nearest center matching.
    
    For each fine-grained bin, find the coarse-grained bin whose center is closest
    to the fine bin's center, and assign the coarse bin's signal value.
    
    Parameters
    ----------
    coarse_df : pd.DataFrame
        DataFrame with 'bin_id' and 'signal' columns from coarse-grained file
    fine_bin_ids : pd.Index
        Index of fine-grained bin_ids to map to
    verbose : bool
        Print progress information
        
    Returns
    -------
    pd.DataFrame
        DataFrame with same shape as fine_bin_ids, containing mapped signal values
    """
    if verbose:
        print(f"    Mapping coarse ({len(coarse_df)} bins) to fine ({len(fine_bin_ids)} bins)...")
    
    # Parse coarse bins
    coarse_parsed = coarse_df["bin_id"].apply(lambda x: pd.Series(parse_bin_id(x), index=["chr", "start", "end"]))
    coarse_df = coarse_df.copy()
    coarse_df["chr"] = coarse_parsed["chr"]
    coarse_df["center"] = (coarse_parsed["start"] + coarse_parsed["end"]) // 2
    
    # Parse fine bins
    fine_parsed = pd.DataFrame(
        [parse_bin_id(bid) for bid in fine_bin_ids],
        columns=["chr", "start", "end"],
        index=fine_bin_ids
    )
    fine_parsed["center"] = (fine_parsed["start"] + fine_parsed["end"]) // 2
    fine_parsed["bin_id"] = fine_parsed.index
    
    # For each chromosome, find nearest coarse bin for each fine bin
    mapped_signals = []
    
    for chrom in fine_parsed["chr"].unique():
        fine_chr = fine_parsed[fine_parsed["chr"] == chrom]
        coarse_chr = coarse_df[coarse_df["chr"] == chrom]
        
        if len(coarse_chr) == 0:
            if verbose:
                print(f"    Warning: No coarse bins found for {chrom}, filling with 0")
            for bin_id in fine_chr["bin_id"]:
                mapped_signals.append({"bin_id": bin_id, "signal": 0})
            continue
        
        # For each fine bin, find nearest coarse bin center
        for _, fine_row in fine_chr.iterrows():
            fine_center = fine_row["center"]
            # Find coarse bin with closest center
            distances = np.abs(coarse_chr["center"] - fine_center)
            nearest_idx = distances.idxmin()
            nearest_signal = coarse_chr.loc[nearest_idx, "signal"]
            mapped_signals.append({"bin_id": fine_row["bin_id"], "signal": nearest_signal})
    
    result_df = pd.DataFrame(mapped_signals)
    
    if verbose:
        print(f"    Mapping complete: {len(result_df)} fine bins mapped")
    
    return result_df


def assign_bins(signal_series: pd.Series, n_bins: int) -> pd.Series:
    """
    Global quantile-based discretization.
    Bin 0: zero signal + lowest non-zero bin (merged)
    Bin 1..(n_bins-1): remaining quantile bins of non-zero signal

    Manually corrected by H. C. on 2026.3.28 to ensure zero signal is always in bin 0,
        and non-zero signals are binned separately.
    """
    bins = pd.Series(0, index=signal_series.index)
    nonzero = signal_series[signal_series > 0]

    if nonzero.empty:
        return bins

    K = min(n_bins, nonzero.nunique())
    edges = np.quantile(nonzero, np.linspace(0, 1, K + 1))
    bin_ids = np.digitize(nonzero, edges[1:-1], right=False)

    bin_ids_adjusted = bin_ids.copy()
    bin_ids_adjusted = np.clip(bin_ids_adjusted, 0, n_bins - 1)

    bins.loc[nonzero.index] = bin_ids_adjusted
    return bins


def load_signals_from_list(file_list_path: str) -> Tuple[List[str], List[str]]:
    """Load signal file paths from a list file."""
    with open(file_list_path, 'r') as f:
        file_paths = [line.strip() for line in f if line.strip()]
    
    feature_names = []
    for path in file_paths:
        basename = os.path.basename(path)
        name = basename.replace('.bin1000.bp', '').replace('.bp', '')
        feature_names.append(name)
    
    return file_paths, feature_names


def build_signal_matrix(
    file_paths: List[str],
    feature_names: List[str],
    n_bins: int,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Load all signals and build a wide-format DataFrame.
    
    Parameters
    ----------
    file_paths : List[str]
        List of signal file paths
    feature_names : List[str]
        Names for each signal (column names in output)
    n_bins : int
        Number of bins for discretization
    verbose : bool
        Print progress
        
    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with index=bin_id, columns=feature_names
    """
    if verbose:
        print(f"\nLoading and binning {len(file_paths)} signals...")
    
    signal_series_list = []
    
    for path, name in zip(file_paths, feature_names):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Signal file not found: {path}")
        
        df = load_signal(path)
        df['signal_bin'] = assign_bins(df['signal'], n_bins)
        
        # Convert to Series indexed by bin_id
        signal_series = df.set_index('bin_id')['signal_bin'].rename(name)
        signal_series_list.append(signal_series)
        
        if verbose:
            print(f"  Loaded & binned: {name} ({len(df)} bins)")
    
    # Merge into wide-format DataFrame
    signal_matrix = pd.concat(signal_series_list, axis=1)
    
    if verbose:
        print(f"Signal matrix shape: {signal_matrix.shape}")
        print(f"  - Rows (bin_ids): {len(signal_matrix)}")
        print(f"  - Columns (signals): {len(signal_matrix.columns)}")
    
    return signal_matrix


def main():
    """CLI for testing signal loading independently."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Load and bin genomic signals")
    parser.add_argument("-i", "--input", required=True, help="File list of signal paths")
    parser.add_argument("-o", "--output", required=True, help="Output parquet file path")
    parser.add_argument("-n", "--n-bins", type=int, default=3, help="Number of bins")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress output")
    
    args = parser.parse_args()
    
    file_paths, feature_names = load_signals_from_list(args.input)
    signal_matrix = build_signal_matrix(
        file_paths, feature_names, 
        n_bins=args.n_bins, 
        verbose=not args.quiet
    )
    
    signal_matrix.to_parquet(args.output)
    if not args.quiet:
        print(f"\nSaved signal matrix to: {args.output}")


if __name__ == "__main__":
    main()
