#!/usr/bin/env python3
"""Existing thyroid API, backed exclusively by the shared compute_hoi_oinfo engine."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.epigenetic_oinfo_root').is_file())
sys.path.insert(0, str(ROOT / 'common/02_oinformation'))
from compute_hoi_oinfo import local_oinfo_from_codes, information_from_codes, permutation_null_from_codes

def oi_local_hoi(X, n_states=5):
    return local_oinfo_from_codes(X, n_states)

def global_oinfo(X, n_states=5):
    return information_from_codes(X, n_states)['hoi']

def permutation_null(X, n_perm=10000, n_states=5, seed=1234):
    return permutation_null_from_codes(X, n_perm=n_perm, n_states=n_states, seed=seed)
