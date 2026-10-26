#!/usr/bin/env python3
"""Reproduce the complete figure audit from the restored release intermediates."""
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
env = os.environ.copy()
env.update(EPINFO_ROOT=str(ROOT), MPLBACKEND='Agg', OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
# The mouse MI runner itself runs its three bin settings in parallel. Keeping
# the groups sequential avoids multiplying their own concurrency limits.
commands = [
    HERE / 'bcells/run.py',
    HERE / 'bcells/make_report.py',
    HERE / 'aging/run_aging_plots.py',
    HERE / 'aging/make_report.py',
    HERE / 'gsea/run.py',
    HERE / 'switch/run_all.py',
    HERE / 'mouse/run_development.py',
    HERE / 'mouse/run_mi.py',
    HERE / 'mouse/make_report.py',
    HERE / 'make_report.py',
    ROOT / 'tools/index_case_figures.py',
]
failed = []
for path in commands:
    print(f'Running {path.relative_to(ROOT)}', flush=True)
    result = subprocess.run([sys.executable, str(path)], cwd=ROOT, env=env)
    if result.returncode:
        failed.append(str(path.relative_to(ROOT)))
if failed:
    print('Runner failures: ' + ', '.join(failed), file=sys.stderr)
    sys.exit(1)
