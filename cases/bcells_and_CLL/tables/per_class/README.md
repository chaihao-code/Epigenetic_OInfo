# cases/bcells_and_CLL/tables/per_class

Observed O-information tables grouped by the original genomic annotation class. Rows retain individual sample/subsystem values; no new across-sample mean is computed.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/bcells_and_CLL/tables/per_class/None_func_AND_CpG_island.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `None_func_AND_CpG_island.tsv.gz`
- `None_func_AND_LINE.tsv.gz`
- `None_func_AND_LTR.tsv.gz`
- `None_func_AND_None_rep.tsv.gz`
- `None_func_AND_SINE.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_None_rep.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_None_rep.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `exon_purified_bins_AND_CpG_island.tsv.gz`
- `exon_purified_bins_AND_LINE.tsv.gz`
- `exon_purified_bins_AND_LTR.tsv.gz`
- `exon_purified_bins_AND_None_rep.tsv.gz`
- `exon_purified_bins_AND_SINE.tsv.gz`
- `intron_purified_bins_AND_CpG_island.tsv.gz`
- `intron_purified_bins_AND_LINE.tsv.gz`
- `intron_purified_bins_AND_LTR.tsv.gz`
- `intron_purified_bins_AND_None_rep.tsv.gz`
- `intron_purified_bins_AND_SINE.tsv.gz`
- `promoter_purified_bins_AND_CpG_island.tsv.gz`
- `promoter_purified_bins_AND_LINE.tsv.gz`
- `promoter_purified_bins_AND_LTR.tsv.gz`
- `promoter_purified_bins_AND_None_rep.tsv.gz`
- `promoter_purified_bins_AND_SINE.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_None_rep.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_None_rep.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
