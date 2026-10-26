# Human thyroid aging

Six healthy donors: CEMT0098, CEMT0130, CEMT0132, CEMT0042, CEMT0044 and CEMT0040. Seven signals, hg38, the original fixed 1-kb grid and 40 annotation classes. The complete CEMT0098 donor is the release demo.

## Inputs

Use the retained case manifests, original genome grid and unchanged annotation masks. Raw/processed local paths are in `config/local_input_links.tsv`. Four-column `.bp` tracks have no header. Reference tables under this case retain their original scientific settings.

## Usage

Start from the release root. For a configured O-information command, use `python tools/plan_oinfo.py --case aging/thyroid --sample <sample>`; the H3K36me3 comparison uses its input manifests rather than a new independent cohort. Read `demo/README.md` for a complete thyroid run.

All new O-information values come from `common/02_oinformation/compute_hoi_oinfo.py`; same-named case files are forwarding adapters. Old entropy-only programs and launchers were removed. New triplet workflows default to 10,000 permutations; historical 1,000-permutation tables remain historical. Use the FDR family described above rather than pooling cases.

## Outputs

Generated thyroid images: [figure locations and index](FIGURES.md).

Existing observed tables are grouped in `tables/per_sample/` and `tables/per_class/` where this is an independent dataset. The H3K36me3 comparison reuses the early/late datasets. New previewed commands write to `outputs/oinfo/<sample>/`. Shared-engine local exports write selected locus IDs, local values, states and local/global checks; original class membership files remain separate and unchanged.

Read `docs/CODE_MAP.md`, `docs/PATHS.md` and subdirectory READMEs for DMR, RNA/GSEA, regression, core and plotting entrypoints. These analyses keep their original models; only the explicit demo is executed during this update.

## Scope of reproduction: archived thyroid values

This demo validates the currently available complete inputs and their end-to-end calculation, not a bytewise reproduction of the archived paper-era thyroid values. The upstream QC reports already document that the original histone divInput files used for those stored values were on an unavailable volume. See `cases/aging/thyroid/deconvolution/results/qc/hoi_recompute_validation.md` and `step9_recomputed_slope_divergence.md`.

The demo records its measured difference from the archived table in `validation.json`; it does not overwrite that table or use the new donor example to claim replication of the archived six-donor age effect. Agreement between the bigWig and supplied-bin routes is checked separately.
