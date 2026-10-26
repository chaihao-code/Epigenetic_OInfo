# cases/aging/thyroid/local_oinfo_age_correlation

Directory index and supporting files for the parent workflow.

Six healthy donors: CEMT0098, CEMT0130, CEMT0132, CEMT0042, CEMT0044 and CEMT0040. Seven signals, hg38, the original fixed 1-kb grid and 40 annotation classes. The complete CEMT0098 donor is the release demo.

## Inputs

The original hg38 1-kb grid, seven binned tracks per donor, unchanged genomic masks and retained donor ages/cell fractions. The complete input-to-output example is `demo/README.md`; raw/processed local input locations are listed in `config/local_input_links.tsv`.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [None_AND_CpG_island](None_AND_CpG_island/README.md)
- [None_AND_LINE](None_AND_LINE/README.md)
- [None_AND_LTR](None_AND_LTR/README.md)
- [None_AND_None](None_AND_None/README.md)
- [None_AND_SINE](None_AND_SINE/README.md)
- [distal_ctcf_bound_enhancer_purified_AND_CpG_island](distal_ctcf_bound_enhancer_purified_AND_CpG_island/README.md)
- [distal_ctcf_bound_enhancer_purified_AND_LINE](distal_ctcf_bound_enhancer_purified_AND_LINE/README.md)
- [distal_ctcf_bound_enhancer_purified_AND_LTR](distal_ctcf_bound_enhancer_purified_AND_LTR/README.md)
- [distal_ctcf_bound_enhancer_purified_AND_None](distal_ctcf_bound_enhancer_purified_AND_None/README.md)
- [distal_ctcf_bound_enhancer_purified_AND_SINE](distal_ctcf_bound_enhancer_purified_AND_SINE/README.md)
- [distal_not_ctcf_bound_enhancer_purified_AND_CpG_island](distal_not_ctcf_bound_enhancer_purified_AND_CpG_island/README.md)
- [distal_not_ctcf_bound_enhancer_purified_AND_LINE](distal_not_ctcf_bound_enhancer_purified_AND_LINE/README.md)
- [distal_not_ctcf_bound_enhancer_purified_AND_LTR](distal_not_ctcf_bound_enhancer_purified_AND_LTR/README.md)
- [distal_not_ctcf_bound_enhancer_purified_AND_None](distal_not_ctcf_bound_enhancer_purified_AND_None/README.md)
- [distal_not_ctcf_bound_enhancer_purified_AND_SINE](distal_not_ctcf_bound_enhancer_purified_AND_SINE/README.md)
- [exon_purified_AND_CpG_island](exon_purified_AND_CpG_island/README.md)
- [exon_purified_AND_LINE](exon_purified_AND_LINE/README.md)
- [exon_purified_AND_LTR](exon_purified_AND_LTR/README.md)
- [exon_purified_AND_None](exon_purified_AND_None/README.md)
- [exon_purified_AND_SINE](exon_purified_AND_SINE/README.md)
- [intron_purified_AND_CpG_island](intron_purified_AND_CpG_island/README.md)
- [intron_purified_AND_LINE](intron_purified_AND_LINE/README.md)
- [intron_purified_AND_LTR](intron_purified_AND_LTR/README.md)
- [intron_purified_AND_None](intron_purified_AND_None/README.md)
- [intron_purified_AND_SINE](intron_purified_AND_SINE/README.md)
- [promoter_purified_AND_CpG_island](promoter_purified_AND_CpG_island/README.md)
- [promoter_purified_AND_LINE](promoter_purified_AND_LINE/README.md)
- [promoter_purified_AND_LTR](promoter_purified_AND_LTR/README.md)
- [promoter_purified_AND_None](promoter_purified_AND_None/README.md)
- [promoter_purified_AND_SINE](promoter_purified_AND_SINE/README.md)
- [proximal_ctcf_bound_enhancer_purified_AND_CpG_island](proximal_ctcf_bound_enhancer_purified_AND_CpG_island/README.md)
- [proximal_ctcf_bound_enhancer_purified_AND_LINE](proximal_ctcf_bound_enhancer_purified_AND_LINE/README.md)
- [proximal_ctcf_bound_enhancer_purified_AND_LTR](proximal_ctcf_bound_enhancer_purified_AND_LTR/README.md)
- [proximal_ctcf_bound_enhancer_purified_AND_None](proximal_ctcf_bound_enhancer_purified_AND_None/README.md)
- [proximal_ctcf_bound_enhancer_purified_AND_SINE](proximal_ctcf_bound_enhancer_purified_AND_SINE/README.md)
- [proximal_not_ctcf_bound_enhancer_purified_AND_CpG_island](proximal_not_ctcf_bound_enhancer_purified_AND_CpG_island/README.md)
- [proximal_not_ctcf_bound_enhancer_purified_AND_LINE](proximal_not_ctcf_bound_enhancer_purified_AND_LINE/README.md)
- [proximal_not_ctcf_bound_enhancer_purified_AND_LTR](proximal_not_ctcf_bound_enhancer_purified_AND_LTR/README.md)
- [proximal_not_ctcf_bound_enhancer_purified_AND_None](proximal_not_ctcf_bound_enhancer_purified_AND_None/README.md)
- [proximal_not_ctcf_bound_enhancer_purified_AND_SINE](proximal_not_ctcf_bound_enhancer_purified_AND_SINE/README.md)

