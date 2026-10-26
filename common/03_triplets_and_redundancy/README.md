# Synergistic triplets, FDR and redundancy cores

## Triplet inputs and computation

Use the shared `common/02_oinformation/compute_hoi_oinfo.py`, defaulting to 10,000 within-class track permutations. For each triplet, lower-tail P = `(1 + number(null <= observed))/(N + 1)`.

BH correction uses every lower-tail triplet P in the chosen family, then a synergistic call requires `observed_hoi < 0` and `q_down_BH <= 0.05`. Never filter negative/significant-looking rows before BH. Default family: sample × class (35 tests for seven signals, 56 for eight, 84 for nine). The B-cell per-sample family is 1,400 tests and the 17-sample global family is 23,800 tests; choose those explicitly.

```bash
python common/03_triplets_and_redundancy/correct_triplet_fdr.py \
  --input path/to/sampleA_null_summary.csv path/to/sampleB_null_summary.csv \
  --output results/triplets_global_BH.tsv --family global
```

Inputs must contain the complete selected testing families once, including positive-O rows; duplicate sample/class/triplet rows are rejected. Outputs include original observations/P values, q-values, family labels/sizes and the final sign-plus-FDR call. Metadata records the input files and definition.

## Redundancy-core inputs and computation

`compute_redundancy_core.py` reads an observed all-subsystem O-information CSV and an explicit signal universe. This is the retained replacement-background procedure (positive Omega, modified-z, GMM/BIC and incoming-signal modality criteria), not the triplet permutation test.

```bash
python common/03_triplets_and_redundancy/compute_redundancy_core.py --help
```

Outputs are core calls, candidate diagnostics and metadata. Existing case core settings are preserved; no full cohort is rerun by opening this folder.
