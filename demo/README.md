# Complete human thyroid O-information demo

## Dataset and fixed inputs

This example uses the complete healthy thyroid donor **CEMT0098**, a 30-year-old female, from the existing thyroid-aging case. There are six histone tracks (H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K27ac, H3K36me3), one DNA-methylation track and a matching Input control for histone normalization.

The original hg38 grid is supplied unchanged as `input/hg38_bins_1kb.bed.gz`: **3,209,513 windows** across the supplied chromosome/contig list. The original loader keeps the **2,875,012 autosomal windows** for information calculations. These are the original 1-kb genomic windows, including original chromosome-end partial windows. They are distinct from the five discrete signal states and the 10,000 permutation repetitions.

All original annotation masks are unchanged. The same 8 functional categories × 5 CpG/repeat categories yield 40 classes; the `None` groups retain their original complement definition. See `config/demo.json` for exact mask paths and original hashes. No coordinates are cropped and no new annotation category is fitted to this demo.

## Two entry routes

Run from the repository root:

```bash
export EPINFO_ROOT="$PWD"
python -m pip install -r demo/requirements.txt
# Route A: supplied original full prepared tracks; local links or compressed copies.
python demo/run_demo.py --jobs 32
# Route B: full original/derived bigWigs -> the same formats -> global/local Omega.
python demo/run_demo.py --from-bigwig --jobs 32
```

The second route needs the UCSC `bigWigAverageOverBed` utility on PATH and the full eight bigWigs listed in `input/bigwig_sources.tsv`. Histone/Input URLs point to bigWigs. The methylation URL points to its original CpG BED.gz; the `.bigWig` used here is the retained conversion, which must be restored as a derived input. Do not rename that CpG download to `.bigWig`. Local workspace links are already configured.

For Slurm:

```bash
export EPINFO_ROOT="$PWD"
export EPINFO_PYTHON="$(command -v python3)"
sbatch demo/submit_demo.slurm
```

The supplied cluster disables `--mem`, so the example omits that directive. Partition and time limits must be adapted on another cluster. The job runs the scientific unit checks and the full bigWig route. Only this donor/demo is executed; the six-donor age-regression study is not automatically rerun.

## What each step computes

1. **Fixed-window binning:** decompress the original BED grid, add a temporary window ID, run `bigWigAverageOverBed` on each whole bigWig and retain its sixth output field (`mean`, over covered bases). Output: headerless `chr start end signal` `.bp` files, with the original filenames, coordinate order and boundaries. The grid is not regenerated or resized.
2. **Input normalization:** divide each histone's binned signal by the matching Input bin plus 0.01. The exact retained thyroid formula is `signal/(Input+0.01)`, with awk output precision; a zero signal therefore stays zero even when Input is zero. DNA methylation is not divided by Input. Its chromosome/start row order is retained from the original bedGraph-derived prepared file; sorting changes no window or membership. The six histone filenames end in `.divInput.bp`; the methylation filename remains `.bp`.
3. **Global signal discretization:** the shared `signal_loader` processes complete autosomal tracks using the original five-state quantile rule. This is done before genomic class partitioning, not independently within a cropped region.
4. **Class membership and sampling:** use the unchanged masks and original 40-class construction. Preserve the original `n_max=100000` and sampling seed 42. All classes remain defined in full; the exact selected loci are written separately, so sampled loci cannot be mistaken for class definitions.
5. **Global O-information:** `compute_hoi_oinfo.py` computes all-subsystem entropies and Omega = TC − DTC. The seven-signal full-system table has one row per original class; all subsystem observations are also retained.
6. **Permutation null:** by default the thyroid demonstration tests its seven-signal full system in each class using 10,000 independent-track permutations, seed 1234, plus-one directional P values. It does not interpret a sign alone as significance.
7. **Local O-information:** the same shared engine calculates a local value at each selected locus on the same state matrix. The mean local value must equal that class's global value within numerical tolerance. This demo's local file uses the global-run sampled loci; the full annotation memberships are supplied separately and unchanged. The original case's full-class local regression inputs remain in its own case directory.
8. **Checks:** verify coordinate/NA agreement of regenerated versus original prepared tracks, record measured signal differences, validate all 40 classes and the local/global identity, and compare observed class values with the original thyroid reference table. The validation JSON reports actual differences rather than claiming equality in advance.

## Outputs

Outputs live under `results/from_bins/` or `results/from_bigwig/`:

| File | Meaning |
|---|---|
| `preprocessing/raw_binned/*.bp` | Eight full fixed-grid averages, before Input normalization (bigWig route) |
| `preprocessing/input_normalized/*.divInput.bp` | Six Input-normalized histone tracks |
| `oinfo/CEEHRC.CEMT0098_entropy_all_subsets.csv` | Joint entropies for every feature subset in each class |
| `oinfo/CEEHRC.CEMT0098_oinfo.csv` | Observed subsystem Omega, TC and DTC; use `order=7`, `n_missing=0` for the full system |
| `global_oinfo_by_class.tsv` | The 40 full-system observations |
| `oinfo/CEEHRC.CEMT0098_permutation_summary.csv` | Observed seven-signal Omega versus its 10,000-permutation null |
| `oinfo/CEEHRC.CEMT0098_sampled_bins.tsv.gz` | Actual selected loci, separate from fixed annotation-class memberships |
| `oinfo/CEEHRC.CEMT0098_discrete_states.tsv.gz` | The seven states at each selected locus |
| `oinfo/CEEHRC.CEMT0098_local_oinfo.tsv.gz` | One full-panel local value per selected locus |
| `oinfo/CEEHRC.CEMT0098_local_global_check.tsv` | Forty checks of mean(local Omega) = global Omega |
| `validation.json`, `run.log` | Actual inputs, commands, software versions, job ID and check results |

## Optional triplet/FDR route

```bash
python demo/run_demo.py --triplets --output demo/results/triplets --jobs 32
```

This uses the same fixed inputs and the generic shared engine: 35 triplets × 40 classes, 10,000 permutations, and BH across all 35 lower-tail tests within each sample/class. It writes `CEEHRC.CEMT0098_triplet_significance.tsv` with `p_down`, `q_down_BH`, `fdr_family_size=35` and `significant_synergy`. Correction precedes selection by Omega < 0 and q <= 0.05. This optional methodological demonstration is not presented as a newly established thyroid biological result.

The local calculation uses the same `log2(p + 1e-12)` numerical stabilization as the retained global entropy estimator, so the local/global mean identity also holds for sparse seven-signal state tables. This preserves the original global O-information arithmetic.

To verify the supplied prepared-track route on Slurm, use `sbatch demo/submit_from_bins.slurm` with the same exported root/interpreter. It uses the complete tracks and all 40 classes, not a reduced window set.

## Scope of reproduction: archived thyroid values

This demo validates the currently available complete inputs and their end-to-end calculation, not a bytewise reproduction of the archived paper-era thyroid values. The upstream QC reports already document that the original histone divInput files used for those stored values were on an unavailable volume. See `cases/aging/thyroid/deconvolution/results/qc/hoi_recompute_validation.md` and `step9_recomputed_slope_divergence.md`.

The demo records its measured difference from the archived table in `validation.json`; it does not overwrite that table or use the new donor example to claim replication of the archived six-donor age effect. Agreement between the bigWig and supplied-bin routes is checked separately.

After both routes finish, run `python demo/compare_runs.py` to compare decoded observed, null, state, locus and local tables and write `demo/results/verification.json`.

## Verified execution

Both complete routes were executed successfully on Slurm: bigWig job **16241245** and prepared-input job **16241246**. Each used all 40 original classes, seven signals and 10,000 permutations. All observed, null, selected-locus, discrete-state and local result files agree exactly after gzip decoding. The maximum mean-local/global difference is 3.141931159689193e-14 bits. The exported selected-locus count is 1955604. Read `results/verification.json` for the full comparison and the separate archived-result limitation.
