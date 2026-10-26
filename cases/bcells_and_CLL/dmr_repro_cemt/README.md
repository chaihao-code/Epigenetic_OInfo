# cases/bcells_and_CLL/dmr_repro_cemt

Directory index and supporting files for the parent workflow.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Prepared four-coordinate/sample methylation matrix at `input/all_bins_1kb_cemt_g2_g3.csv`, the retained 13-donor age/sex metadata and optional GENCODE annotation. Beta values are on the original 0–1 scale.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

The canonical case wrapper writes limma results, BH q-values, delta-beta classifications and run diagnostics to its `output/` directory. Adjacent 1-kb windows are not merged.

## Subdirectories

- [input](input/README.md)
- [scripts](scripts/README.md)

