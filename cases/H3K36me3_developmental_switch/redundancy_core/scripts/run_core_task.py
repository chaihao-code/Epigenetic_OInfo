#!/usr/bin/env python3
"""Run one manifest sample through the shared redundancy-core CLI."""

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
import os
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest = _release_manifest(json.loads((root / 'inputs/analysis_manifest.json').read_text()))
task = manifest['tasks'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
output = root / 'results' / task['sample'] / 'redundancy_cores.csv'
command = [sys.executable, manifest['core_program'], '--input', task['input'],
           '--output', str(output), '--sample', task['sample'], '--signals', *manifest['signals'],
           '--min-order', str(manifest['min_order']), '--z-threshold', str(manifest['z_threshold']),
           '--eta-threshold', str(manifest['eta_threshold']), '--expected-blocks', '40']
subprocess.run(command, check=True)
