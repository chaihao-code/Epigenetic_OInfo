# Code, inputs and outputs

| Stage | Supported implementation | Input | Output |
|---|---|---|---|
| Fixed-grid preprocessing | `demo/scripts/prepare_from_bigwig.py`; retained common preprocessing utilities | Full bigWigs, unchanged original BED windows and matching Input | Original four-column `.bp`; six histone `.divInput.bp`; methylation `.bp` |
| Global O-information | `common/02_oinformation/compute_hoi_oinfo.py` | Complete binned tracks, unchanged masks, explicit output directory | Subset entropies, Omega/TC/DTC, independent-track permutation summaries |
| Triplet FDR | Same engine; `common/03_triplets_and_redundancy/correct_triplet_fdr.py` for pooling samples | Complete lower-tail triplet test families | BH q-values, family sizes, negative-O plus q-threshold calls |
| Redundancy cores | `common/03_triplets_and_redundancy/compute_redundancy_core.py` | Observed subsystem CSV, named signal universe | Replacement-background core calls and diagnostics |
| Local O-information | Same shared engine; `common/04_local_oinformation/oi_local_hoi.py` is an adapter | Identical state matrix and locus set | Per-locus values and mean-local/global identity checks |
| Thyroid age/immune regression | `cases/aging/thyroid/local_O_info_regression/scripts/` | Retained six-donor local tables, ages and proportions | Age/immune coefficients, P values, existing pathway/marginal GSEA outputs |
| DMR | `cases/bcells_and_CLL/dmr_repro_cemt/scripts/run_dmr.R` | Prepared 13-donor beta matrix and age/sex design | limma M-value test, BH, delta-beta-selected 1-kb windows |
| DE/GSEA | B-cell `scripts/analyze_rna_de_male.py`, `gsea_rna_male.py` | Imported RNA matrices, method/IGHV covariates, frozen GO/GOA | Male R/S DE and ranked GSEA |

The only global/local O-information arithmetic is in `compute_hoi_oinfo.py`; former case engines forward to it and thyroid recomputation helpers call its API. Core identification and group/age regression are separate statistical analyses of its results.

New triplet defaults are 10,000 permutations. Preserve explicit family boundaries: mouse seven-signal sample/class BH35; shared-eight-signal BH56; preimplantation BH84; B-cell sample BH1400 and global BH23800. Older reference runs are not relabelled as new 10,000-permutation results. Full-signal sign-only occurrence branches retain their historical definition and are labelled accordingly.

The retained bigWig averaging implementation uses the sixth UCSC output field, covered-base `mean`. The demo follows that actual implementation; it does not silently substitute `mean0`. The legacy observed `entropy` field is a first-feature marginal entropy; all joint entropies are supplied separately.
