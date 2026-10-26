# demo/input/binned

Directory index and supporting files for the parent workflow.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Follow the parent workflow; these files are inputs or reference outputs rather than executable entrypoints.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `98308.CEEHRC.CEMT0098.H3K4me1.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98311.CEEHRC.CEMT0098.H3K4me3.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98314.CEEHRC.CEMT0098.H3K9me3.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98317.CEEHRC.CEMT0098.H3K27me3.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98320.CEEHRC.CEMT0098.H3K36me3.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98323.CEEHRC.CEMT0098.H3K27ac.signal_unstranded.bin1000.step1000.divInput.bp.gz`
- `98642.CEEHRC.CEMT0098.WGB-Seq.methylation_profile.bin1000.step1000.bp.gz`
