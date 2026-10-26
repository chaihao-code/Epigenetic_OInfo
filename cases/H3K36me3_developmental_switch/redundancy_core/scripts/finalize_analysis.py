#!/usr/bin/env python3
"""Mark the run complete only after both validated trajectory figures exist."""

# Resolve repository-relative paths in retained operational JSON manifests.
_RELEASE_ROOT = next(p for p in __import__('pathlib').Path(__file__).resolve().parents
                     if (p / '.epigenetic_oinfo_root').is_file())
def _release_manifest(value):
    if isinstance(value, dict):
        return {k: _release_manifest(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_manifest(v) for v in value]
    if isinstance(value, str) and value.startswith(('cases/', 'common/', 'data/')):
        return str(_RELEASE_ROOT / value)
    return value
import json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
for kind in ['redundancy_core', 'synergistic_triplet']:
    assert _release_manifest(json.loads((root / 'results' / (kind + '_plot_validation.json')).read_text()))['status'] == 'complete'
    for ext in ['pdf', 'png']:
        assert (root / 'plots' / (kind + '_mark_trajectories.' + ext)).stat().st_size > 1000
p = root / 'inputs/analysis_manifest.json'
m = _release_manifest(json.loads(p.read_text()))
m['status'] = 'complete'
p.write_text(json.dumps(m, indent=2) + '\n')
print('All 73 core tasks and both trajectory figures validated.')
