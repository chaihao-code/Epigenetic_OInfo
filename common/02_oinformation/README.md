# Single O-information engine

`compute_hoi_oinfo.py` is the only numerical implementation of global/local Omega, exact-count entropy and permutation nulls in this release. Case copies with the same filename are forwarding adapters. The former no-permutation calculation programs and their launch scripts have been deleted.

## Inputs

`--input-file`: one headerless four-column `.bp` signal path per line. `--block-groups`: unchanged two-column `bin_id`/0–1 annotation masks, separated by `;` within groups and `::BLOCKGROUP::` between groups. Paths may be absolute or repository-root-relative. The original `signal_loader.py` five-state rule and `block_processor.py` membership construction are retained. Compressed TSV masks are accepted by pandas.

## Usage

```bash
python common/02_oinformation/compute_hoi_oinfo.py --help
python tools/plan_oinfo.py --case mouse_development --sample 2cell
python demo/run_demo.py --jobs 32
```

Triplet defaults: 10,000 permutations, orders 3–3, within-class independent-track shuffling, seeds 42/1234, lower-tail BH within each sample/class, alpha 0.05. `--permutation-level system` instead tests the full signal system; the thyroid demo uses this explicitly. The CLI rejects zero permutations. `--n-max 100000` is the original cap on genomic loci and is a separate parameter from the permutation count.

`--fdr-family sample` pools triplets/classes within a sample. Global correction across samples is provided by `common/03_triplets_and_redundancy/correct_triplet_fdr.py`; it cannot be honestly calculated by a single-sample run. Correction includes all triplets in the chosen family before selecting negative Omega.

## Outputs

Observed `*_entropy_all_subsets.csv`, `*_oinfo.csv`, null summaries/tensors and `*_metadata.json`. Triplet nulls additionally produce `*_triplet_significance.tsv`. Add `--save-local-oinfo --save-sampled-bins` to save local values, state codes, selected locus IDs and local/global mean checks. The local export uses the same sampled matrix as the global calculation; fixed class memberships remain separate.

## Python API

`information_from_codes(X, n_states=5)` returns global `hoi/tc/dtc/n`; `local_oinfo_from_codes(X, n_states=5)` returns local values; `permutation_null_from_codes` supports the existing thyroid adapter. Rows of X are loci and columns are signals. Integer codes must lie in `[0, n_states)`. Pure numerical APIs support invariant checks and adapters; all standalone CLI workflows include permutation testing.

The legacy `entropy` column in observed O-information rows retains the first remaining signal's marginal entropy for compatibility. Joint entropies are in `*_entropy_all_subsets.csv`; do not relabel that compatibility column as joint system entropy.

The local calculation uses the same `log2(p + 1e-12)` numerical stabilization as the retained global entropy estimator, so the local/global mean identity also holds for sparse seven-signal state tables. This preserves the original global O-information arithmetic.
