#!/usr/bin/env python3
"""
Consolidated HOI pipeline: one pass from signal list to two big CSVs.

Pipeline (per sample):
  1. Load + discretize signals and partition into blocks, reusing
     signal_loader.py / block_processor.py unchanged, so preprocessing is
     bit-identical to the legacy compute_hoi_entropy.py runs.
  2. Observed: all-subset joint entropies per block -> a single
     <sample>_entropy_all_subsets.csv (legacy schema, no per-block files);
     O-information (legacy arithmetic via compute_oinfo_from_entropy.py) for
     the full system plus an optional reduced subsystem (e.g. 7-signal)
     -> a single <sample>_oinfo.csv.
  3. Permutation null (--n-permutations, default 1000): for each replicate,
     every signal track is independently permuted across its own genomic
     positions (marginal distribution, sparsity and discretized bin
     occupancy preserved; locus-wise cross-signal correspondence destroyed),
     and O-information is recomputed on the same per-block sampled loci.
     -> <sample>_permutation_oinfo.csv + <sample>_permutation_summary.csv.
     Optional --save-permutation-entropy also dumps the full all-subset
     entropy table of every replicate (very large; off by default).
     With --permutation-level all, the permutation instead recomputes
     O-information for EVERY subsystem (all subsets of >= 2 signals) per
     block; the null is stored as a compact float64 .npy of shape
     (n_permutations, n_blocks, n_subsystems) plus mapping tables. For wide
     panels (many signals) --permutation-max-order caps the null to small
     subsystems, since the full-system entropy table grows as n_bins^n_signals.

Entropy estimator: exact frequency counting via base-n_bins integer row
encoding + np.bincount, mathematically identical to the legacy
np.unique-based counting.
"""

import argparse
import ast
import json
import os
import sys
import time
from itertools import combinations
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd

HOI_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
if HOI_CODE_DIR not in sys.path:
    sys.path.insert(0, HOI_CODE_DIR)

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
from compute_oinfo_from_entropy import compute_oinfo_from_entropy_df

NAN_CODE = -1  # int8 sentinel for missing values (real codes are 0..n_bins-1)


# ---------------------------------------------------------------------------
# Entropy helpers (exact frequency counting)
# ---------------------------------------------------------------------------

def encode_subset(cols):
    """Encode a list of 1-D integer code arrays as a single base-n_bins int."""
    enc = cols[0].astype(np.int64)
    for c in cols[1:]:
        enc = enc * _N_BINS[0] + c.astype(np.int64)
    return enc


_N_BINS = [5]  # mutable cell so worker processes see the value via fork


def entropy_from_codes(enc, n_samples, table_size):
    counts = np.bincount(enc, minlength=table_size)
    probs = counts[counts > 0] / n_samples
    return float(-np.sum(probs * np.log2(probs + 1e-12)))


def all_subset_entropies(data, feature_names, n_bins):
    """
    Joint entropy H(S) for all non-empty subsets S (legacy loop order).
    data: (n_samples, n_features) integer codes in [0, n_bins).
    Returns DataFrame with columns order/subset_indices/subset_names/entropy.
    """
    n_features = data.shape[1]
    cols = [data[:, i] for i in range(n_features)]
    n_samples = data.shape[0]
    rows = []
    for r in range(1, n_features + 1):
        for subset in combinations(range(n_features), r):
            enc = encode_subset([cols[i] for i in subset])
            h = entropy_from_codes(enc, n_samples, n_bins ** r)
            rows.append({
                'order': r,
                'subset_indices': list(subset),
                'subset_names': '/'.join(feature_names[i] for i in subset),
                'entropy': h,
            })
    return pd.DataFrame(rows)


def all_subsystems(n_features):
    """All subsets of size >= 2, ordered by (size, lexicographic indices)."""
    return [sub for r in range(2, n_features + 1)
            for sub in combinations(range(n_features), r)]


def entropy_map_dict(x, n_features, n_bins, max_order=None):
    """{(i,): H_i, (i,j): H_ij, ...} for non-empty subsets of x's columns.

    max_order: if given, only compute subsets of size <= max_order (used to
    bound the permutation-null cost for wide panels, where the full table
    n_bins**n_features would be infeasible).
    """
    cols = [x[:, i] for i in range(n_features)]
    m = x.shape[0]
    h = {}
    r_max = n_features if max_order is None else min(max_order, n_features)
    for r in range(1, r_max + 1):
        size = n_bins ** r
        for sub in combinations(range(n_features), r):
            h[sub] = entropy_from_codes(
                encode_subset([cols[i] for i in sub]), m, size)
    return h


def oinfo_from_entropy_map(h, sub):
    """O-information of subsystem sub from an entropy_map_dict."""
    k = len(sub)
    h_n = h[sub]
    tc = sum(h[(j,)] for j in sub) - h_n
    dtc = sum(h[sub[:i] + sub[i + 1:]] for i in range(k)) - (k - 1) * h_n
    return tc - dtc


def validate_against_legacy(oinfo_df, ref_root, sample_name,
                            tolerance=1e-9):
    """Compare observed O-info against legacy per-block oinfo csv files.

    Legacy layout: <ref_root>/<sample>/<block>/<sample>_<block>_oinfo_*.csv
    Returns (n_checked, max_diff, missing_keys, extra_keys).
    """
    ref_dir = Path(ref_root) / sample_name
    system = oinfo_df['system'].iloc[0]
    obs = {}
    for _, row in oinfo_df[oinfo_df['system'] == system].iterrows():
        key = (row['block'], tuple(_as_index_list(row['remaining_indices'])))
        obs[key] = (float(row['hoi']), float(row['tc']), float(row['dtc']))
    n_checked = 0
    max_diff = 0.0
    missing = []
    for block_dir in sorted(ref_dir.iterdir()):
        if not block_dir.is_dir():
            continue
        for csv_path in sorted(block_dir.glob("*_oinfo*.csv")):
            df = pd.read_csv(csv_path)
            for _, row in df.iterrows():
                key = (block_dir.name,
                       tuple(_as_index_list(row['remaining_indices'])))
                if key not in obs:
                    missing.append(key)
                    continue
                o = obs.pop(key)
                d = max(abs(float(row['hoi']) - o[0]),
                        abs(float(row['tc']) - o[1]),
                        abs(float(row['dtc']) - o[2]))
                max_diff = max(max_diff, d)
                n_checked += 1
    extra = sorted(obs.keys())
    return n_checked, max_diff, missing, extra


# ---------------------------------------------------------------------------
# Subsystem (e.g. 7-signal) filtering
# ---------------------------------------------------------------------------

def _as_index_list(idx):
    if isinstance(idx, str):
        return ast.literal_eval(idx)
    return list(idx)


def find_drop_indices(entropy_df, tokens):
    """Original feature indices whose order-1 name contains any token."""
    drop = set()
    order1 = entropy_df[entropy_df['order'] == 1]
    for _, row in order1.iterrows():
        name = str(row['subset_names'])
        if any(tok and tok in name for tok in tokens):
            drop.update(_as_index_list(row['subset_indices']))
    return drop


def filter_subsystem(entropy_df, drop):
    mask = entropy_df['subset_indices'].apply(
        lambda idx: not any(i in drop for i in _as_index_list(idx)))
    return entropy_df[mask].copy()


# ---------------------------------------------------------------------------
# Observed phase workers (fork-shared globals)
# ---------------------------------------------------------------------------

_G_CODES = None          # (N, n_features) int8 discretized signal matrix
_G_FEATURE_NAMES = None  # list of str
_G_BLOCK_TASKS = None    # list of (block_label, positions) for permutation phase
_G_BLOCK_MATS = None     # per-block observed sample matrices (block scope)
_G_SUBSET_STRUCT = None  # dict with precomputed index lists per system
_G_SUBSYSTEMS = None     # list of subsystem index-tuples (permutation-level all)
_G_MAX_ORDER = [None]    # optional cap on permutation subsystem order
_PERM_SCOPE = ['block']  # 'genome' or 'block'


def _observed_block_worker(task):
    block_label, positions = task
    data = _G_CODES[positions]
    df = all_subset_entropies(data, _G_FEATURE_NAMES, _N_BINS[0])
    df['block'] = block_label
    return df


def _permutation_worker(args):
    perm_id, save_entropy_path = args
    rng = np.random.default_rng(_PERM_SEED[0] * 1000003 + perm_id)

    if _PERM_SCOPE[0] == 'genome':
        n_rows, n_tracks = _G_CODES.shape
        pm = np.full((n_rows, n_tracks), NAN_CODE, dtype=np.int8)
        for j in range(n_tracks):
            pos = _G_COL_POS[j]
            val = _G_COL_VAL[j]
            pm[pos, j] = val[rng.permutation(val.shape[0])]
        block_sources = [pm[positions] for _, positions in _G_BLOCK_TASKS]
    else:
        block_sources = _G_BLOCK_MATS

    rows = []
    entropy_rows = []
    for block_idx, (block_label, _) in enumerate(_G_BLOCK_TASKS):
        x0 = block_sources[block_idx]
        m = x0.shape[0]
        if _PERM_SCOPE[0] == 'block':
            # independent within-block permutation of each track's values
            # across the block's sampled loci: the per-track value multiset
            # within the block is preserved exactly
            x = x0.copy()
            for j in range(x0.shape[1]):
                x[:, j] = x0[rng.permutation(m), j]
        else:
            x = x0
        for system, struct in _G_SUBSET_STRUCT.items():
            all_idx = struct['all']
            h1 = {}
            for j in all_idx:
                enc = encode_subset([x[:, j]])
                h1[j] = entropy_from_codes(enc, m, _N_BINS[0])
            h_n = entropy_from_codes(
                encode_subset([x[:, j] for j in all_idx]),
                m, _N_BINS[0] ** len(all_idx))
            h_nm1 = []
            for js in struct['minus1']:
                h_nm1.append(entropy_from_codes(
                    encode_subset([x[:, j] for j in js]),
                    m, _N_BINS[0] ** len(js)))
            n = len(all_idx)
            tc = sum(h1.values()) - h_n
            dtc = sum(h_nm1) - (n - 1) * h_n
            rows.append({
                'perm': perm_id, 'block': block_label, 'system': system,
                'hoi': tc - dtc, 'tc': tc, 'dtc': dtc,
                'mean_h1': float(np.mean(list(h1.values()))),
            })
        if save_entropy_path is not None:
            df = all_subset_entropies(x, _G_FEATURE_NAMES, _N_BINS[0])
            df.insert(0, 'block', block_label)
            df.insert(0, 'perm', perm_id)
            entropy_rows.append(df)

    if save_entropy_path is not None and entropy_rows:
        pd.concat(entropy_rows, ignore_index=True).to_csv(
            save_entropy_path, index=False)
    return rows


def _permutation_worker_all(perm_id):
    """All-subsystem permutation null: returns (n_blocks, n_subsystems) array
    of O-information values for this replicate."""
    rng = np.random.default_rng(_PERM_SEED[0] * 1000003 + perm_id)
    n_blocks = len(_G_BLOCK_TASKS)
    n_sub = len(_G_SUBSYSTEMS)
    out = np.empty((n_blocks, n_sub), dtype=np.float64)
    for b_idx, (block_label, positions) in enumerate(_G_BLOCK_TASKS):
        if _PERM_SCOPE[0] == 'genome':
            # genome-scope is only implemented for permutation-level system
            raise NotImplementedError(
                "permutation-level all requires --permutation-scope block")
        x0 = _G_BLOCK_MATS[b_idx]
        m = x0.shape[0]
        x = x0.copy()
        for j in range(x0.shape[1]):
            x[:, j] = x0[rng.permutation(m), j]
        h = entropy_map_dict(x, x0.shape[1], _N_BINS[0], _G_MAX_ORDER[0])
        for s_idx, sub in enumerate(_G_SUBSYSTEMS):
            out[b_idx, s_idx] = oinfo_from_entropy_map(h, sub)
    return perm_id, out


_PERM_SEED = [1234]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Consolidated entropy + O-information pipeline with permutation null",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input-file", "-i", required=True,
                        help="File listing signal track paths (one per line)")
    parser.add_argument("--output-dir", "-o", required=True)
    parser.add_argument("--sample-name", "-s", default=None)
    parser.add_argument("--n-bins", type=int, default=5)
    parser.add_argument("--n-max", type=int, default=100000)
    parser.add_argument("--random-seed", type=int, default=42,
                        help="Seed for per-block n_max sampling")
    parser.add_argument("--block-files", nargs="+", default=None)
    parser.add_argument("--block-groups", type=str, default=None)
    parser.add_argument("--group-separator", type=str, default="::BLOCKGROUP::")
    parser.add_argument("--drop-signals", type=str, default="H2A.Z,H3K36me3",
                        help="Comma-separated name tokens defining the reduced "
                             "subsystem (e.g. the two signals absent from the "
                             "gastrulation data). Empty string disables it.")
    parser.add_argument("--n-permutations", type=int, default=1000,
                        help="Number of locus-permutation replicates (0 = off)")
    parser.add_argument("--permutation-seed", type=int, default=1234)
    parser.add_argument("--permutation-scope", choices=['block', 'genome'],
                        default='block',
                        help="'block': independently permute each track within "
                             "each block's sampled loci (block-level marginals "
                             "preserved exactly); 'genome': permute each track "
                             "across its own genomic positions (genome-wide "
                             "marginals preserved)")
    parser.add_argument("--permutation-level", choices=['system', 'all'],
                        default='system',
                        help="'system': null only for the full system (and the "
                             "drop-signals subsystem); 'all': null for EVERY "
                             "subsystem (all subsets of >= 2 signals) per "
                             "block, stored as a compact .npy + summary")
    parser.add_argument("--permutation-min-order", type=int, default=2,
                        help="Minimum subsystem order included in the permutation null")
    parser.add_argument("--permutation-max-order", type=int, default=None,
                        help="With --permutation-level all, restrict the null "
                             "to subsystems of at most this order (the entropy "
                             "table of the full system grows as n_bins^n_signals, "
                             "so wide panels should cap this, e.g. 3 for "
                             "triplet analysis). Observed O-info is still "
                             "computed for every subsystem.")
    parser.add_argument("--validate-against", type=str, default=None,
                        help="Optional legacy results root "
                             "(<root>/<sample>/<block>/<sample>_<block>_oinfo_*.csv) "
                             "to validate observed O-info against; exits "
                             "non-zero on mismatch > 1e-9")
    parser.add_argument("--save-permutation-entropy", action="store_true",
                        help="Also dump the full all-subset entropy of every "
                             "permutation replicate (very large)")
    parser.add_argument("--n-jobs", type=int, default=-1)
    parser.add_argument("--quiet", "-q", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet
    if args.sample_name is None:
        args.sample_name = Path(args.input_file).stem
    _N_BINS[0] = args.n_bins
    _PERM_SEED[0] = args.permutation_seed

    os.makedirs(args.output_dir, exist_ok=True)

    t0 = time.time()
    if verbose:
        print(f"=== compute_hoi_oinfo: {args.sample_name} ===")

    # ----- Part 0: load + discretize (identical to legacy pipeline) -----
    file_paths, feature_names = load_signals_from_list(args.input_file)
    signal_matrix = build_signal_matrix(
        file_paths, feature_names, n_bins=args.n_bins, verbose=verbose)

    block_files = parse_block_groups_from_args(
        args.block_files, args.block_groups, args.group_separator)
    block_groups = [load_block_group(g) for g in block_files]
    block_combinations = cartesian_product_blocks(
        block_groups, all_bin_ids=set(signal_matrix.index))

    # Integer-code matrix (NaN -> sentinel) + per-track non-NaN positions
    mat = signal_matrix.values.astype(np.float64)
    codes = np.full(mat.shape, NAN_CODE, dtype=np.int8)
    valid = ~np.isnan(mat)
    codes[valid] = mat[valid].astype(np.int8)
    n_nan = int((codes == NAN_CODE).sum())
    if n_nan:
        raise ValueError(
            f"{n_nan} NaN-coded cells in the discretized matrix; the "
            f"entropy counters cannot handle them. Resolve missing data "
            f"before running (historical datasets had complete coverage).")

    col_pos = [np.flatnonzero(valid[:, j]) for j in range(mat.shape[1])]
    col_val = [codes[col_pos[j], j] for j in range(mat.shape[1])]

    # Per-block extraction + sampling, exactly as the legacy code path
    block_tasks = []
    for block_label, bin_ids in block_combinations:
        df_block = extract_block_data(signal_matrix, bin_ids, block_label,
                                      verbose=False)
        if df_block is None:
            if verbose:
                print(f"  [SKIP] block with missing signals: {block_label}")
            continue
        df_sampled = sample_block_data(df_block, args.n_max, args.random_seed)
        positions = signal_matrix.index.get_indexer(df_sampled.index)
        block_tasks.append((block_label, positions))

    if verbose:
        print(f"Blocks to process: {len(block_tasks)}")

    # ----- Part 1: observed all-subset entropies (parallel over blocks) -----
    global _G_CODES, _G_FEATURE_NAMES
    _G_CODES = codes
    _G_FEATURE_NAMES = feature_names

    n_jobs = args.n_jobs if args.n_jobs > 0 else cpu_count()
    if len(block_tasks) > 1 and n_jobs > 1:
        with Pool(processes=min(n_jobs, len(block_tasks))) as pool:
            block_dfs = pool.map(_observed_block_worker, block_tasks)
    else:
        block_dfs = [_observed_block_worker(t) for t in block_tasks]

    entropy_df = pd.concat(block_dfs, ignore_index=True)
    entropy_path = os.path.join(
        args.output_dir, f"{args.sample_name}_entropy_all_subsets.csv")
    entropy_df.to_csv(entropy_path, index=False)
    if verbose:
        print(f"Wrote {entropy_path} ({len(entropy_df)} rows) "
              f"[{time.time() - t0:.1f}s]")

    # ----- Part 2: observed O-information (9-sig + optional subsystem) -----
    tokens = [t.strip() for t in args.drop_signals.split(',') if t.strip()]
    drop = find_drop_indices(block_dfs[0], tokens) if tokens else set()
    if tokens:
        if verbose:
            print(f"Subsystem drop indices: {sorted(drop)} ({len(drop)} signals)")
        if len(drop) not in (0, 2):
            raise ValueError(
                f"drop-signals matched {len(drop)} signals, expected 0 or 2")
    compute_7sig = len(drop) == 2

    oinfo_frames = []
    for (block_label, _), block_ent in zip(block_tasks, block_dfs):
        res = compute_oinfo_from_entropy_df(block_ent)
        for _, df in res.items():
            df.insert(0, 'block', block_label)
            df.insert(0, 'system', '9sig')
            oinfo_frames.append(df)
        if compute_7sig:
            sub = filter_subsystem(block_ent, drop)
            res = compute_oinfo_from_entropy_df(sub)
            for _, df in res.items():
                df.insert(0, 'block', block_label)
                df.insert(0, 'system', '7sig')
                oinfo_frames.append(df)

    oinfo_df = pd.concat(oinfo_frames, ignore_index=True)
    oinfo_path = os.path.join(args.output_dir, f"{args.sample_name}_oinfo.csv")
    oinfo_df.to_csv(oinfo_path, index=False)
    if verbose:
        print(f"Wrote {oinfo_path} ({len(oinfo_df)} rows) "
              f"[{time.time() - t0:.1f}s]")

    # ----- Part 2b: optional validation against legacy results -----
    if args.validate_against:
        n_checked, max_diff, missing, extra = validate_against_legacy(
            oinfo_df, args.validate_against, args.sample_name)
        print(f"[validate] {args.sample_name}: {n_checked} subsystem rows "
              f"checked vs {args.validate_against}, max|diff| = {max_diff:.3e} "
              f"({'OK' if max_diff < 1e-9 else 'FAIL'}), "
              f"{len(missing)} legacy rows unmatched, {len(extra)} computed "
              f"rows unmatched")
        if max_diff >= 1e-9 or missing or extra:
            if missing[:5]:
                print(f"  missing e.g. {missing[:5]}")
            if extra[:5]:
                print(f"  extra e.g. {extra[:5]}")
            sys.exit(1)

    # ----- Part 3: permutation null -----
    summary_path = None
    if args.n_permutations > 0:
      if args.permutation_level == 'all':
        # Null for EVERY subsystem (all subsets of >= 2 signals) per block.
        if args.permutation_scope != 'block':
            raise ValueError(
                "--permutation-level all requires --permutation-scope block")
        subsystems = all_subsystems(len(feature_names))
        if args.permutation_max_order is not None:
            subsystems = [s for s in subsystems
                          if len(s) <= args.permutation_max_order]
        subsystems = [sub for sub in subsystems if len(sub) >= args.permutation_min_order]
        if not subsystems:
            raise ValueError("No subsystems in the requested permutation order range")
        n_sub = len(subsystems)
        if verbose:
            print(f"Permutation level: all ({n_sub} subsystems x "
                  f"{len(block_tasks)} blocks"
                  f"{', order <= ' + str(args.permutation_max_order) if args.permutation_max_order else ''})")

        # mapping tables
        blocks_df = pd.DataFrame({
            'block_idx': range(len(block_tasks)),
            'block': [bl for bl, _ in block_tasks],
            'n_samples': [int(len(pos)) for _, pos in block_tasks],
        })
        blocks_df.to_csv(os.path.join(
            args.output_dir, f"{args.sample_name}_blocks.csv"), index=False)
        subsets_df = pd.DataFrame({
            'subset_idx': range(n_sub),
            'order': [len(s) for s in subsystems],
            'subset_indices': [list(s) for s in subsystems],
            'subset_names': ['/'.join(feature_names[i] for i in s)
                             for s in subsystems],
        })
        subsets_df.to_csv(os.path.join(
            args.output_dir, f"{args.sample_name}_subsets.csv"), index=False)

        global _G_SUBSYSTEMS, _G_BLOCK_TASKS, _G_BLOCK_MATS, _PERM_SCOPE
        _G_SUBSYSTEMS = subsystems
        _G_BLOCK_TASKS = block_tasks
        _G_BLOCK_MATS = [codes[pos] for _, pos in block_tasks]
        _PERM_SCOPE[0] = 'block'
        _G_MAX_ORDER[0] = args.permutation_max_order

        if verbose:
            print(f"Running {args.n_permutations} permutations on "
                  f"{min(n_jobs, args.n_permutations)} workers...")
        t_perm = time.time()
        if n_jobs > 1 and args.n_permutations > 1:
            with Pool(processes=min(n_jobs, args.n_permutations)) as pool:
                results = pool.map(_permutation_worker_all,
                                   range(args.n_permutations))
        else:
            results = [_permutation_worker_all(p)
                       for p in range(args.n_permutations)]
        results.sort(key=lambda r: r[0])
        null_arr = np.stack([r[1] for r in results])  # (perm, block, subset)
        null_path = os.path.join(
            args.output_dir, f"{args.sample_name}_permutation_null.npy")
        np.save(null_path, null_arr)
        if verbose:
            print(f"Wrote {null_path} {null_arr.shape} "
                  f"[perm phase {time.time() - t_perm:.1f}s]")

        # summary vs observed (full system rows of oinfo_df)
        system = oinfo_df['system'].iloc[0]
        obs_lookup = {}
        for _, row in oinfo_df[oinfo_df['system'] == system].iterrows():
            key = (row['block'], tuple(_as_index_list(row['remaining_indices'])))
            obs_lookup[key] = float(row['hoi'])
        block_labels = [bl for bl, _ in block_tasks]
        n_perm = null_arr.shape[0]
        rows = []
        for b_idx, bl in enumerate(block_labels):
            for s_idx, sub in enumerate(subsystems):
                null = null_arr[:, b_idx, s_idx]
                obs = obs_lookup[(bl, sub)]
                p_down = (1 + int(np.sum(null <= obs))) / (n_perm + 1)
                p_up = (1 + int(np.sum(null >= obs))) / (n_perm + 1)
                rows.append({
                    'block': bl, 'subset_idx': s_idx, 'order': len(sub),
                    'subset_names': '/'.join(feature_names[i] for i in sub),
                    'n_perm': n_perm, 'observed_hoi': obs,
                    'null_mean': float(np.mean(null)),
                    'null_sd': float(np.std(null)),
                    'null_median': float(np.median(null)),
                    'p_down': p_down, 'p_up': p_up,
                    'p_dev': min(1.0, 2 * min(p_up, p_down)),
                })
        summary_df = pd.DataFrame(rows)
        summary_path = os.path.join(
            args.output_dir, f"{args.sample_name}_null_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        if verbose:
            print(f"Wrote {summary_path} ({len(summary_df)} rows)")
      else:
        if len(drop) == 2:
            sub_idx = [i for i in range(len(feature_names)) if i not in drop]
        else:
            sub_idx = None
        struct9 = {'all': list(range(len(feature_names))),
                   'minus1': [[j for j in range(len(feature_names)) if j != i]
                              for i in range(len(feature_names))]}
        subset_struct = {'9sig': struct9}
        if sub_idx is not None:
            subset_struct['7sig'] = {
                'all': sub_idx,
                'minus1': [[j for j in sub_idx if j != i] for i in sub_idx]}
        if verbose:
            print(f"Permutation systems: {list(subset_struct.keys())}")

        global _G_SUBSET_STRUCT, _G_COL_POS, _G_COL_VAL
        _G_BLOCK_TASKS = block_tasks
        _G_SUBSET_STRUCT = subset_struct
        _G_COL_POS = col_pos
        _G_COL_VAL = col_val
        _PERM_SCOPE[0] = args.permutation_scope
        if args.permutation_scope == 'block':
            _G_BLOCK_MATS = [codes[pos] for _, pos in block_tasks]
        if verbose:
            print(f"Permutation scope: {args.permutation_scope}")

        perm_args = []
        for p in range(args.n_permutations):
            chunk_path = None
            if args.save_permutation_entropy:
                chunk_path = os.path.join(
                    args.output_dir, f"_perm_entropy_{args.sample_name}_{p:05d}.csv")
            perm_args.append((p, chunk_path))

        if verbose:
            print(f"Running {args.n_permutations} permutations on "
                  f"{min(n_jobs, args.n_permutations)} workers...")
        t_perm = time.time()
        if n_jobs > 1 and args.n_permutations > 1:
            chunksize = max(1, args.n_permutations // (n_jobs * 4))
            with Pool(processes=min(n_jobs, args.n_permutations)) as pool:
                perm_rows = pool.map(_permutation_worker, perm_args,
                                     chunksize=chunksize)
        else:
            perm_rows = [_permutation_worker(a) for a in perm_args]
        perm_df = pd.DataFrame([r for batch in perm_rows for r in batch])
        perm_path = os.path.join(
            args.output_dir, f"{args.sample_name}_permutation_oinfo.csv")
        perm_df.to_csv(perm_path, index=False)
        if verbose:
            print(f"Wrote {perm_path} ({len(perm_df)} rows) "
                  f"[perm phase {time.time() - t_perm:.1f}s]")

        if args.save_permutation_entropy:
            import glob as globmod
            chunks = sorted(globmod.glob(os.path.join(
                args.output_dir, f"_perm_entropy_{args.sample_name}_*.csv")))
            with open(os.path.join(
                    args.output_dir,
                    f"{args.sample_name}_permutation_entropy.csv"), 'w') as out:
                for i, c in enumerate(chunks):
                    with open(c) as f:
                        for line_no, line in enumerate(f):
                            if i > 0 and line_no == 0:
                                continue
                            out.write(line)
            for c in chunks:
                os.remove(c)
            if verbose:
                print("Concatenated permutation entropy chunks")

        # Summary: observed vs the permutation null DISTRIBUTION.
        # p_up   = (1 + #{null >= obs})/(N+1): observed shifted UP relative to
        #          the null (redundancy-type; under dependence the estimation
        #          bias can only shrink, so this implies true redundancy).
        # p_down = (1 + #{null <= obs})/(N+1): synergy-type deviation.
        # p_dev  = min(1, 2*min(p_up, p_down)): two-sided deviation.
        obs_lookup = oinfo_df[oinfo_df['n_missing'] == 0].set_index(
            ['system', 'block'])['hoi'].to_dict()
        summary_rows = []
        for (system, block), group in perm_df.groupby(['system', 'block']):
            null = group['hoi'].to_numpy()
            obs = float(obs_lookup[(system, block)])
            n = len(null)
            mean_null = float(np.mean(null))
            p_down = (1 + int(np.sum(null <= obs))) / (n + 1)
            p_up = (1 + int(np.sum(null >= obs))) / (n + 1)
            p_dev = min(1.0, 2 * min(p_up, p_down))
            summary_rows.append({
                'system': system, 'block': block, 'n_perm': n,
                'observed_hoi': obs,
                'null_mean': mean_null,
                'null_sd': float(np.std(null)),
                'p_dev': p_dev, 'p_up': p_up, 'p_down': p_down,
            })
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(
            args.output_dir, f"{args.sample_name}_permutation_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        if verbose:
            print(f"Wrote {summary_path} ({len(summary_df)} rows)")

    # ----- Metadata -----
    metadata = {
        'sample_name': args.sample_name,
        'input_file': os.path.abspath(args.input_file),
        'feature_names': feature_names,
        'n_features': len(feature_names),
        'n_bins': args.n_bins,
        'n_max': args.n_max,
        'random_seed': args.random_seed,
        'blocks': {bl: int(len(pos)) for bl, pos in block_tasks},
        'drop_signals': args.drop_signals,
        'n_permutations': args.n_permutations,
        'permutation_seed': args.permutation_seed,
        'permutation_scope': args.permutation_scope,
        'permutation_level': args.permutation_level,
        'permutation_min_order': args.permutation_min_order,
        'permutation_max_order': args.permutation_max_order,
        'validate_against': args.validate_against,
        'save_permutation_entropy': bool(args.save_permutation_entropy),
        'elapsed_seconds': round(time.time() - t0, 1),
    }
    with open(os.path.join(args.output_dir,
                           f"{args.sample_name}_metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)

    if verbose:
        print(f"=== done in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
