# common

Directory index and supporting files for the parent workflow.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [01_preprocessing](01_preprocessing/README.md)
- [02_oinformation](02_oinformation/README.md)
- [03_triplets_and_redundancy](03_triplets_and_redundancy/README.md)
- [04_local_oinformation](04_local_oinformation/README.md)
- [05_dmr](05_dmr/README.md)
- [06_de_and_gsea](06_de_and_gsea/README.md)

