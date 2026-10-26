# cases/H3K36me3_developmental_switch/redundancy_core/full_signal_occurrence

Directory index and supporting files for the parent workflow.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [data](data/README.md)
- [inputs](inputs/README.md)

