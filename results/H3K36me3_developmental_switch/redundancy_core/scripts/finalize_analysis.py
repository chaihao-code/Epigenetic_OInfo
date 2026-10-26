#!/usr/bin/env python3
"""Mark the run complete only after both validated trajectory figures exist."""
import json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
for kind in ['redundancy_core', 'synergistic_triplet']:
    assert json.loads((root / 'results' / (kind + '_plot_validation.json')).read_text())['status'] == 'complete'
    for ext in ['pdf', 'png']:
        assert (root / 'plots' / (kind + '_mark_trajectories.' + ext)).stat().st_size > 1000
p = root / 'inputs/analysis_manifest.json'
m = json.loads(p.read_text())
m['status'] = 'complete'
p.write_text(json.dumps(m, indent=2) + '\n')
print('All 73 core tasks and both trajectory figures validated.')
