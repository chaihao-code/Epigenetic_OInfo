# Three-signal O-information permutation significance

17 samples, 40 genomic annotation blocks, 35 triplets per block (23,800 tests).
Same preprocessing as results_oinfo_n_bins_5: 5 bins, at most 100,000 loci
per block, sampling seed 42. Each track is independently shuffled within
each sampled block, preserving its marginal distribution (1,000 permutations,
seed 1234). The generic pipeline also outputs two-signal nulls; triplet tables
contain only order 3.

Per-sample outputs: `<sample>/<sample>_triplet_significance.tsv`.
Combined output after all samples finish: `all_samples_triplet_significance.tsv`.

- `observed_hoi`: observed triplet O-information in bits.
- `p_up`, `p_down`: upper/lower-tail empirical P values, with +1 correction.
- `p_dev`: min(1, 2 * min(p_up, p_down)).
- `q_*_sample_BH`: BH correction separately for each P-value type across
  all 1,400 triplet/block tests in the same sample.
- `q_*_global_BH` (combined table only): BH correction separately for each
  P-value type across all 23,800 sample/block/triplet tests.

These P values assess departure from independently shuffled tracks;
they do not test the difference between clinical groups.
The sign of observed O-information and direction relative to the permutation
null are distinct, since finite-sample null O-information need not be zero.

Run `scripts/submit_oinfo_triplets.slurm` via sbatch. Submit
`scripts/aggregate_oinfo_triplets.slurm` with an afterok dependency on the
array job to produce the combined table automatically.
