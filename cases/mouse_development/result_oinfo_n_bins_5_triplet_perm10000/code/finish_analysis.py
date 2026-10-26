#!/usr/bin/env python3
"""Require all 82 results, run per-block BH35, and regenerate the plots."""

# Portable repository paths; no analysis is executed by this helper.
from pathlib import Path as _ReleasePath
_RELEASE_ROOT = next(p for p in _ReleasePath(__file__).resolve().parents
                     if (p / ".epigenetic_oinfo_root").is_file())
def _release_path(value):
    return value.replace("@ROOT@", str(_RELEASE_ROOT))
def _release_config(value):
    if isinstance(value, dict):
        return {k: _release_config(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_config(v) for v in value]
    return _release_path(value) if isinstance(value, str) else value
import csv
import os
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parent.parent
    entries = list(csv.DictReader((root / 'input_manifest.tsv').open(), delimiter='\t'))
    assert len(entries) == 82
    for entry in entries:
        sample = entry['sample']
        for suffix in ('null_summary.csv', 'metadata.json'):
            if not (root / sample / f'{sample}_{suffix}').exists():
                raise RuntimeError(f'Missing completed result: {sample}/{suffix}')
    subprocess.run([sys.executable, '-u', str(root / 'negative_fdr' / 'analyze.py')], check=True)
    env = os.environ.copy()
    env.update(HOI_NULL_DIR=str(root), HOI_OUT_DIR=str(root / 'negative_fdr' / 'plots'),
               MPLCONFIGDIR='/tmp/mouse-development-mpl', XDG_CACHE_HOME='/tmp/mouse-development-cache')
    plot = root.parent / 'script' / 'plot_triplet_counts_corrected_by_stage.py'
    subprocess.run([__import__('sys').executable, '-u', str(plot)], env=env, check=True)


if __name__ == '__main__':
    main()
