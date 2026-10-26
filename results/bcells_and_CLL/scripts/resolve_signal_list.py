#!/usr/bin/env python3
"""Resolve the retained signal manifests against a configurable data directory."""
import argparse
from pathlib import Path


def resolve_manifest(manifest, data_root):
    paths = []
    for line in Path(manifest).read_text().splitlines():
        if not line.strip():
            continue
        # Retain filenames and their data-relative hierarchy from historical manifests.
        relative = line.strip().split('/data/', 1)[-1]
        path = Path(data_root) / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        paths.append(path.resolve())
    if len(paths) != 7:
        raise ValueError(f'Expected seven signal files, found {len(paths)}')
    return paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('data_root', type=Path)
    args = parser.parse_args()
    print('\n'.join(map(str, resolve_manifest(args.manifest, args.data_root))))
