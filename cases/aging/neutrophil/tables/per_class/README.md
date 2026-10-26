# cases/aging/neutrophil/tables/per_class

Observed O-information tables grouped by the original genomic annotation class. Rows retain individual sample/subsystem values; no new across-sample mean is computed.

12 retained samples; original O-information and age-regression reference tables. This cleanup does not invent a new neutrophil triplet study.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/aging/neutrophil/tables/per_class/None_AND_CpG_island.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `None_AND_CpG_island.tsv.gz`
- `None_AND_LINE.tsv.gz`
- `None_AND_LTR.tsv.gz`
- `None_AND_None.tsv.gz`
- `None_AND_SINE.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_None.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_None.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `exon_purified_bins_AND_CpG_island.tsv.gz`
- `exon_purified_bins_AND_LINE.tsv.gz`
- `exon_purified_bins_AND_LTR.tsv.gz`
- `exon_purified_bins_AND_None.tsv.gz`
- `exon_purified_bins_AND_SINE.tsv.gz`
- `intron_purified_bins_AND_CpG_island.tsv.gz`
- `intron_purified_bins_AND_LINE.tsv.gz`
- `intron_purified_bins_AND_LTR.tsv.gz`
- `intron_purified_bins_AND_None.tsv.gz`
- `intron_purified_bins_AND_SINE.tsv.gz`
- `promoter_purified_bins_AND_CpG_island.tsv.gz`
- `promoter_purified_bins_AND_LINE.tsv.gz`
- `promoter_purified_bins_AND_LTR.tsv.gz`
- `promoter_purified_bins_AND_None.tsv.gz`
- `promoter_purified_bins_AND_SINE.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_None.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_LINE.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_LTR.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_None.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins_AND_SINE.tsv.gz`
