# cases

Directory index and supporting files for the parent workflow.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Generated result images: [case-by-case locations](FIGURES.md).

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [H3K36me3_developmental_switch](H3K36me3_developmental_switch/README.md)
- [aging](aging/README.md)
- [bcells_and_CLL](bcells_and_CLL/README.md)
- [mouse_after_E11.5](mouse_after_E11.5/README.md)
- [mouse_development](mouse_development/README.md)
- [mouse_preimp](mouse_preimp/README.md)

