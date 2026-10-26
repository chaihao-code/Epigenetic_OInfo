#!/usr/bin/env python3
"""Export a source-only directory using explicit source/input allowlists."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

PROJECT = Path(__file__).resolve().parents[1]


def source_files():
    paths = {PROJECT / name for name in ['README.md', 'CODE_MAP.md', 'requirements.txt', '.gitignore']}
    for folder in ['scripts', 'dmr_repro_cemt/scripts']:
        paths.update(p for p in (PROJECT / folder).rglob('*')
                     if p.is_file() and p.suffix in {'.py', '.R', '.sh', '.slurm', '.json', '.md'})
    paths.update((PROJECT / 'global_oinfo_and_triplets').glob('*.py'))
    paths.add(PROJECT / 'global_oinfo_and_triplets/reference_config.json')
    paths.update((PROJECT / 'data_info').glob('CEEHRC.CEMT*.txt'))
    paths.add(PROJECT / 'data_info/cll_metadata.tsv')
    paths.add(PROJECT / 'data_info/external_input_provenance.json')
    lookup = PROJECT / 'data_info/ihec_rna_lookup'
    for name in ['rna_sample_metadata.tsv', 'RNA_download_guide_zh.md', 'download_rna_cll13.sh',
                 'rna_genes_cll13.urls.txt', 'rna_isoforms_cll13.urls.txt', 'rna_qc_cll13.urls.txt']:
        paths.add(lookup / name)
    for name in ['gene_annotation.tsv', 'sample_metadata.tsv', 'sample_qc.tsv']:
        paths.add(PROJECT / 'results_rna_de_slurm' / name)
    for name in ['ABC_triplet_classes/README_zh.md', 'global_oinfo_and_triplets/README_zh.md',
                 'dmr_repro_cemt/README.md', 'results_rna_de_slurm/README.md',
                 'results_oinfo_triplet_n_bins_5/README.md', 'results_rna_de_female_validation/README_zh.md']:
        paths.add(PROJECT / name)
    return sorted(paths)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, nargs='?', default=PROJECT.with_name('bcells_github'))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f'Output already exists; choose a new directory: {output}')
    if output.is_relative_to(PROJECT):
        parser.error('Export outside the working analysis directory')
    paths = source_files()
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    manifest = {}
    output.mkdir(parents=True)
    for path in paths:
        relative = path.relative_to(PROJECT)
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        manifest[str(relative)] = {'bytes': path.stat().st_size,
                                   'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    (output / 'SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Exported {len(paths)} source/documentation/input files to {output}')


if __name__ == '__main__':
    main()
