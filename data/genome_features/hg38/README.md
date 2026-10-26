# data/genome_features/hg38

Directory index and supporting files for the parent workflow.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `data/genome_features/hg38/CpG_island_bins.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `CpG_island_bins.tsv.gz`
- `LINE_bins.tsv.gz`
- `LTR_bins.tsv.gz`
- `None_func_bins.tsv.gz`
- `None_rep_bins.tsv.gz`
- `SINE_bins.tsv.gz`
- `distal_ctcf_bound_enhancer_purified_bins.tsv.gz`
- `distal_not_ctcf_bound_enhancer_purified_bins.tsv.gz`
- `exon_purified_bins.tsv.gz`
- `intron_purified_bins.tsv.gz`
- `promoter_bins.tsv.gz`
- `promoter_purified_bins.tsv.gz`
- `proximal_ctcf_bound_enhancer_purified_bins.tsv.gz`
- `proximal_not_ctcf_bound_enhancer_purified_bins.tsv.gz`
