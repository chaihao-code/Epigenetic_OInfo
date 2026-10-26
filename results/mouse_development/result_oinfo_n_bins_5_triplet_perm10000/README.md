# Mouse development: 10,000-permutation triplets

Replacement for the cancelled 100,000-permutation run. All 82 samples run in a new directory, preserving existing 100,000-permutation outputs.

- Five bins, at most 100,000 loci per annotation block; sampling seed 42.
- 10,000 independent within-block permutations; permutation seed 1234.
- 40 annotation blocks × 35 triplets per sample; triplet-only permutation tests.
- Same validated scientific computing snapshot as the previous run.
- 32 CPUs per task, partition C064M0256G, 24-hour time limit, no array concurrency throttle.
- Null tensor per sample: (10000, 40, 35), approximately 112 MB.
- Minimum attainable one-sided p-value: 1/10001, approximately 0.0001.

Run `python submit.py` to submit the array and an afterok-dependent analysis job. Accepted job IDs are recorded in `submission.json`; submission is idempotent. The dependent job requires all 82 results, applies BH separately to all 35 `p_down` values in each sample × block, then selects negative O-information with q <= 0.05. It regenerates the A/B/C plots in `negative_fdr/plots/`, including `triplet_ABC_stacked.png` and `.pdf`.

`code_provenance.json` records the new snapshot hashes. `validation.json` is inherited evidence for the unchanged computing code; `run_settings.json` records this lineage and the updated permutation count. Input panels and their hashes are listed in `input_manifest.tsv`.
