# cases/aging/thyroid/tables/per_class

Observed O-information tables grouped by the original genomic annotation class. Rows retain individual sample/subsystem values; no new across-sample mean is computed.

Six healthy donors: CEMT0098, CEMT0130, CEMT0132, CEMT0042, CEMT0044 and CEMT0040. Seven signals, hg38, the original fixed 1-kb grid and 40 annotation classes. The complete CEMT0098 donor is the release demo.

## Inputs

The original hg38 1-kb grid, seven binned tracks per donor, unchanged genomic masks and retained donor ages/cell fractions. The complete input-to-output example is `demo/README.md`; raw/processed local input locations are listed in `config/local_input_links.tsv`.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/aging/thyroid/tables/per_class/None_AND_CpG_island.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

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
- `promoter_AND_CpG_island.tsv.gz`
- `promoter_AND_LINE.tsv.gz`
- `promoter_AND_LTR.tsv.gz`
- `promoter_AND_None.tsv.gz`
- `promoter_AND_SINE.tsv.gz`
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
