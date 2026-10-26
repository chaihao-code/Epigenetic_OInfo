# cases/bcells_and_CLL/results_oinfo_triplet_n_bins_5

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [CEEHRC.CEMT0004](CEEHRC.CEMT0004/README.md)
- [CEEHRC.CEMT0005](CEEHRC.CEMT0005/README.md)
- [CEEHRC.CEMT0006](CEEHRC.CEMT0006/README.md)
- [CEEHRC.CEMT0025](CEEHRC.CEMT0025/README.md)
- [CEEHRC.CEMT0026](CEEHRC.CEMT0026/README.md)
- [CEEHRC.CEMT0027](CEEHRC.CEMT0027/README.md)
- [CEEHRC.CEMT0028](CEEHRC.CEMT0028/README.md)
- [CEEHRC.CEMT0029](CEEHRC.CEMT0029/README.md)
- [CEEHRC.CEMT0030](CEEHRC.CEMT0030/README.md)
- [CEEHRC.CEMT0094](CEEHRC.CEMT0094/README.md)
- [CEEHRC.CEMT0095](CEEHRC.CEMT0095/README.md)
- [CEEHRC.CEMT0096](CEEHRC.CEMT0096/README.md)
- [CEEHRC.CEMT0097](CEEHRC.CEMT0097/README.md)
- [CEEHRC.CEMT0129](CEEHRC.CEMT0129/README.md)
- [CEEHRC.CEMT0134](CEEHRC.CEMT0134/README.md)
- [CEEHRC.CEMT0135](CEEHRC.CEMT0135/README.md)
- [CEEHRC.CEMT0138](CEEHRC.CEMT0138/README.md)

