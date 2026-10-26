# cases/aging/thyroid/local_O_info_regression/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

Six healthy donors: CEMT0098, CEMT0130, CEMT0132, CEMT0042, CEMT0044 and CEMT0040. Seven signals, hg38, the original fixed 1-kb grid and 40 annotation classes. The complete CEMT0098 donor is the release demo.

## Inputs

The original hg38 1-kb grid, seven binned tracks per donor, unchanged genomic masks and retained donor ages/cell fractions. The complete input-to-output example is `demo/README.md`; raw/processed local input locations are listed in `config/local_input_links.tsv`.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `01_regress_oinfo_immune.py`
- `02_gsea_pathway.py`
- `03_gsea_marginal.py`
- `04_summarize.py`
- `config.py`
- `run_all.sh`
