#!/usr/bin/env python3
"""Compatibility path forwarding all computation to common/02_oinformation."""
import importlib.util
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.epigenetic_oinfo_root').is_file())
CORE_DIR = ROOT / 'common/02_oinformation'
sys.path.insert(0, str(CORE_DIR))
NAME = '_epio_canonical_compute_hoi_oinfo'
if NAME not in sys.modules:
    spec = importlib.util.spec_from_file_location(NAME, CORE_DIR / 'compute_hoi_oinfo.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[NAME] = module
    spec.loader.exec_module(module)
core = sys.modules[NAME]
for name in dir(core):
    if not name.startswith('__'):
        globals()[name] = getattr(core, name)
if __name__ == '__main__':
    core.main()
