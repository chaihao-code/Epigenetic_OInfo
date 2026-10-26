# B cells and CLL

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Use the retained case manifests, original genome grid and unchanged annotation masks. Raw/processed local paths are in `config/local_input_links.tsv`. Four-column `.bp` tracks have no header. Reference tables under this case retain their original scientific settings.

## Usage

Start from the release root. For a configured O-information command, use `python tools/plan_oinfo.py --case bcells_and_CLL --sample <sample>`; the H3K36me3 comparison uses its input manifests rather than a new independent cohort. Read `demo/README.md` for a complete thyroid run.

All new O-information values come from `common/02_oinformation/compute_hoi_oinfo.py`; same-named case files are forwarding adapters. Old entropy-only programs and launchers were removed. New triplet workflows default to 10,000 permutations; historical 1,000-permutation tables remain historical. Use the FDR family described above rather than pooling cases.

## Outputs

Generated result images: [figure locations and file index](FIGURES.md).

Existing observed tables are grouped in `tables/per_sample/` and `tables/per_class/` where this is an independent dataset. The H3K36me3 comparison reuses the early/late datasets. New previewed commands write to `outputs/oinfo/<sample>/`. Shared-engine local exports write selected locus IDs, local values, states and local/global checks; original class membership files remain separate and unchanged.

Read `docs/CODE_MAP.md`, `docs/PATHS.md` and subdirectory READMEs for DMR, RNA/GSEA, regression, core and plotting entrypoints. These analyses keep their original models; only the explicit demo is executed during this update.
