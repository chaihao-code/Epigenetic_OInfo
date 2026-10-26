"""Single supported exact-count O-information engine.

Loads binned tracks, discretizes globally per track, partitions annotation classes,
samples at most n_max loci/class, computes all-subset entropies and Omega=TC-DTC,
and runs independent-track permutations. Triplet inference defaults to 10,000
within-class permutations and lower-tail BH correction within each sample/class.
Optional local values and sampled-locus IDs use precisely the same state matrix.
Historical observations are not recalculated merely by installing this release.
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
def measures_from_subset_entropies(entropy_df: pd.DataFrame) -> dict:
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
            return ast.literal_eval(s)
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

NAN_CODE = -1  # int8 sentinel for missing values (real codes are 0..n_bins-1)


# ---------------------------------------------------------------------------
# Entropy helpers (exact frequency counting)
# ---------------------------------------------------------------------------

def encode_subset(cols, n_bins=None):
    """Encode a list of 1-D integer code arrays as a single base-n_bins int."""
    radix = _N_BINS[0] if n_bins is None else n_bins
    enc = cols[0].astype(np.int64)
    for c in cols[1:]:
        enc = enc * radix + c.astype(np.int64)
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
            enc = encode_subset([cols[i] for i in subset], n_bins)
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
                encode_subset([cols[i] for i in sub], n_bins), m, size)
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

def _validated_states(x, n_states):
    x = np.asarray(x)
    if x.ndim != 2 or x.shape[0] == 0 or x.shape[1] < 3:
        raise ValueError('Expected a nonempty loci-by-signals matrix with at least three signals')
    if not np.isfinite(x).all() or not np.equal(x, np.floor(x)).all():
        raise ValueError('Discrete states must be finite integers')
    if n_states < 2 or np.any(x < 0) or np.any(x >= n_states):
        raise ValueError('State codes must be between zero and n_states - 1')
    return np.ascontiguousarray(x, dtype=np.int64)


def entropy_of_states(x, n_states=5):
    x = np.asarray(x, dtype=np.int64)
    if x.ndim != 2 or x.shape[0] == 0:
        raise ValueError('Expected a nonempty two-dimensional state matrix')
    return entropy_from_codes(encode_subset([x[:, i] for i in range(x.shape[1])], n_states),
                              len(x), n_states ** x.shape[1])


def information_from_codes(x, n_states=5):
    """Exact global Omega, TC and DTC, sharing the CLI entropy functions."""
    x = _validated_states(x, n_states)
    n_features = x.shape[1]
    h = entropy_map_dict(x, n_features, n_states)
    full = tuple(range(n_features))
    tc = sum(h[(j,)] for j in full) - h[full]
    dtc = sum(h[full[:i] + full[i + 1:]] for i in full) - (n_features - 1) * h[full]
    return dict(hoi=float(tc-dtc), tc=float(tc), dtc=float(dtc), n=len(x))


def local_oinfo_from_codes(x, n_states=5):
    """Local Omega per locus. Mean equals global Omega on the identical matrix."""
    x = _validated_states(x, n_states)
    n_loci, n_features = x.shape
    cols = [x[:, j] for j in range(n_features)]
    def surprisal(columns):
        enc = encode_subset(columns, n_states)
        counts = np.bincount(enc, minlength=n_states ** len(columns))
        return -np.log2(counts[enc] / n_loci + 1e-12)
    values = (n_features - 2) * surprisal(cols)
    for j in range(n_features):
        values += surprisal([cols[j]]) - surprisal([cols[k] for k in range(n_features) if k != j])
    return values


def permutation_null_from_codes(x, n_perm=10000, n_states=5, seed=1234):
    """Single-region global-O null for the existing local-regression adapter."""
    x = _validated_states(x, n_states)
    if n_perm < 1:
        raise ValueError('n_perm must be positive')
    obs = information_from_codes(x, n_states)['hoi']
    rng = np.random.default_rng(seed)
    values = np.empty(n_perm)
    for i in range(n_perm):
        xp = x.copy()
        for j in range(x.shape[1]):
            xp[:, j] = x[rng.permutation(len(x)), j]
        values[i] = information_from_codes(xp, n_states)['hoi']
    down = (1 + int(np.sum(values <= obs))) / (n_perm + 1)
    up = (1 + int(np.sum(values >= obs))) / (n_perm + 1)
    return dict(observed_hoi=obs, null_mean=float(values.mean()), null_sd=float(values.std()),
                null_median=float(np.median(values)), p_down=down, p_up=up,
                p_dev=min(1., 2 * min(down, up)), n_perm=n_perm)


def benjamini_hochberg(p_values):
    """BH across the entire supplied family, including nonsignificant/positive-O rows."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('P-values must be a finite one-dimensional array in [0, 1]')
    if p.size == 0:
        return p.copy()
    order = np.argsort(p, kind='stable')
    q = np.empty(len(p))
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0., 1.)
    return q


def correct_triplet_fdr(frame, family='sample-block', alpha=0.05):
    """Lower-tail BH, then require negative observed Omega; never prefilter by sign/P."""
    required = {'sample', 'block', 'order', 'observed_hoi', 'p_down'}
    if not required.issubset(frame.columns):
        raise ValueError(f'Missing columns: {sorted(required - set(frame.columns))}')
    if family not in ['sample-block', 'sample', 'global'] or not 0 < alpha < 1:
        raise ValueError('Invalid FDR family or alpha')
    out = frame.loc[frame['order'].eq(3)].copy().reset_index(drop=True)
    keys = {'sample-block': ['sample', 'block'], 'sample': ['sample'], 'global': []}[family]
    groups = out.groupby(keys, sort=False).indices.values() if keys else [np.arange(len(out))]
    out['q_down_BH'] = np.nan
    out['fdr_family_size'] = 0
    for idx in groups:
        out.loc[idx, 'q_down_BH'] = benjamini_hochberg(out.loc[idx, 'p_down'])
        out.loc[idx, 'fdr_family_size'] = len(idx)
    out['fdr_family'] = family
    out['fdr_alpha'] = alpha
    out['significant_synergy'] = out['observed_hoi'].lt(0) & out['q_down_BH'].le(alpha)
    return out


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
    parser.add_argument("--drop-signals", type=str, default="",
                        help="Comma-separated name tokens defining the reduced "
                             "subsystem (e.g. the two signals absent from the "
                             "gastrulation data). Empty string disables it.")
    parser.add_argument("--n-permutations", type=int, default=10000,
                        help="Positive number of locus-permutation replicates")
    parser.add_argument("--permutation-seed", type=int, default=1234)
    parser.add_argument("--permutation-scope", choices=['block', 'genome'],
                        default='block',
                        help="'block': independently permute each track within "
                             "each block's sampled loci (block-level marginals "
                             "preserved exactly); 'genome': permute each track "
                             "across its own genomic positions (genome-wide "
                             "marginals preserved)")
    parser.add_argument("--permutation-level", choices=['system', 'all'],
                        default='all',
                        help="'system': null only for the full system (and the "
                             "drop-signals subsystem); 'all': null for EVERY "
                             "subsystem (all subsets of >= 2 signals) per "
                             "block, stored as a compact .npy + summary")
    parser.add_argument("--permutation-min-order", type=int, default=3,
                        help="Minimum subsystem order included in the permutation null")
    parser.add_argument("--permutation-max-order", type=int, default=3,
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
    parser.add_argument('--fdr-family', choices=['sample-block', 'sample'], default='sample-block',
                        help='Lower-tail BH family for triplets; global correction is a separate multi-sample step')
    parser.add_argument('--fdr-alpha', type=float, default=0.05)
    parser.add_argument('--save-local-oinfo', action='store_true',
                        help='Save full-panel local Omega on the same sampled loci used for global Omega')
    parser.add_argument('--save-sampled-bins', action='store_true')
    args = parser.parse_args()
    if args.n_permutations < 1:
        parser.error('--n-permutations must be positive; the obsolete no-permutation CLI has been removed')
    if args.n_bins < 2 or args.n_bins > 127:
        parser.error('--n-bins must be between 2 and 127 for int8 state storage')
    if not 0 < args.fdr_alpha < 1:
        parser.error('--fdr-alpha must lie between zero and one')
    return args


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
        res = measures_from_subset_entropies(block_ent)
        for _, df in res.items():
            df.insert(0, 'block', block_label)
            df.insert(0, 'system', f'{len(feature_names)}sig')
            oinfo_frames.append(df)
        if compute_7sig:
            sub = filter_subsystem(block_ent, drop)
            res = measures_from_subset_entropies(sub)
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

    if args.save_local_oinfo or args.save_sampled_bins:
        sampled_frames, local_frames, state_frames = [], [], []
        checks = []
        for block_label, positions in block_tasks:
            ids = signal_matrix.index[positions]
            sampled_frames.append(pd.DataFrame({'sample': args.sample_name, 'block': block_label, 'bin_id': ids}))
            if args.save_local_oinfo:
                x = codes[positions]
                values = local_oinfo_from_codes(x, args.n_bins)
                local_frames.append(pd.DataFrame({'sample': args.sample_name, 'block': block_label,
                                                 'bin_id': ids, 'local_oinfo': values}))
                states = pd.DataFrame(x, columns=feature_names)
                states.insert(0, 'bin_id', ids)
                states.insert(0, 'block', block_label)
                state_frames.append(states)
                global_value = float(oinfo_df.loc[(oinfo_df['block'].eq(block_label)) &
                                                  oinfo_df['n_missing'].eq(0)].iloc[0]['hoi'])
                error = abs(float(values.mean()) - global_value)
                if error > 1e-8:
                    raise ValueError(f'Local/global mismatch for {block_label}: {error}')
                checks.append({'block': block_label, 'n_loci': len(ids), 'global_oinfo': global_value,
                               'mean_local_oinfo': float(values.mean()), 'absolute_difference': error})
        pd.concat(sampled_frames).to_csv(os.path.join(args.output_dir,
                        f'{args.sample_name}_sampled_bins.tsv.gz'), sep='\t', index=False)
        if args.save_local_oinfo:
            pd.concat(local_frames).to_csv(os.path.join(args.output_dir,
                        f'{args.sample_name}_local_oinfo.tsv.gz'), sep='\t', index=False)
            pd.concat(state_frames).to_csv(os.path.join(args.output_dir,
                        f'{args.sample_name}_discrete_states.tsv.gz'), sep='\t', index=False)
            pd.DataFrame(checks).to_csv(os.path.join(args.output_dir,
                        f'{args.sample_name}_local_global_check.tsv'), sep='\t', index=False)

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
        triplets = summary_df.loc[summary_df['order'].eq(3)].copy()
        if not triplets.empty:
            triplets.insert(0, 'sample', args.sample_name)
            triplets = correct_triplet_fdr(triplets, family=args.fdr_family, alpha=args.fdr_alpha)
            triplet_path = os.path.join(args.output_dir, f'{args.sample_name}_triplet_significance.tsv')
            triplets.to_csv(triplet_path, sep='\t', index=False)
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
        subset_struct = {f'{len(feature_names)}sig': struct9}
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
        'fdr_family': args.fdr_family,
        'fdr_alpha': args.fdr_alpha,
        'save_local_oinfo': args.save_local_oinfo,
        'elapsed_seconds': round(time.time() - t0, 1),
    }
    with open(os.path.join(args.output_dir,
                           f"{args.sample_name}_metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)

    if verbose:
        print(f"=== done in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
