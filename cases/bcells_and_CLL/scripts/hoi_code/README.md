# cases/bcells_and_CLL/scripts/hoi_code

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Use the retained sample metadata, imported counts/TPM in the B-cell case and the frozen GO/GOA annotation snapshot under `data/annotations/go_human/`. O-information programs instead use that case's seven-track signal manifests and original hg38 masks.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `block_processor.py`: `python cases/bcells_and_CLL/scripts/hoi_code/block_processor.py --help`; options include `--signal-matrix`, `--block-files`, `--block-groups`, `--output-dir`, `--n-max`, `--seed`, `--quiet`.
- `compute_hoi_oinfo.py` forwards imports and command-line execution to the single shared engine.
- `signal_loader.py`: `python cases/bcells_and_CLL/scripts/hoi_code/signal_loader.py --help`; options include `--input`, `--output`, `--n-bins`, `--quiet`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `block_processor.py`
- `compute_hoi_oinfo.py`
- `signal_loader.py`
