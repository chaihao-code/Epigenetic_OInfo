#!/usr/bin/env bash
# Submit consolidated entropy, O-information and all-subsystem permutation jobs.
# The historical filename is retained for compatibility.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
REPO_DIR="$(cd "$PROJECT_DIR/../.." && pwd)"
INPUT_DIR="$PROJECT_DIR/data/input_files"
RESULTS_DIR="$PROJECT_DIR/result_oinfo_n_bins_5"
PYTHON_SCRIPT="$REPO_DIR/script/hoi_code/compute_hoi_oinfo.py"
PYTHON_BIN="${HOI_PYTHON:-/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3}"
FEATURE_DIR="$REPO_DIR/data/genome_features/mm10"

PARTITION="C064M0256G"
TIME_LIMIT="04:00:00"
CPUS=32
N_MAX=100000
N_PERMUTATIONS=1000
PREPARE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --partition) PARTITION="${2:?Missing partition}"; shift 2 ;;
        --time) TIME_LIMIT="${2:?Missing time limit}"; shift 2 ;;
        --cpus) CPUS="${2:?Missing CPU count}"; shift 2 ;;
        --n-max) N_MAX="${2:?Missing maximum bin count}"; shift 2 ;;
        --n-permutations) N_PERMUTATIONS="${2:?Missing permutation count}"; shift 2 ;;
        --prepare-only) PREPARE_ONLY=1; shift ;;
        --help|-h)
            cat <<'HELP'
Usage: bash submit_entropy_computation.sh [OPTIONS]

Run compute_hoi_oinfo.py for every seven-signal data/input_files/*.txt panel.
Outputs: result_oinfo_n_bins_5/<sample>/.

  --prepare-only       Validate inputs and prepare jobs without submission
  --partition NAME     Slurm partition (default: C064M0256G)
  --time HH:MM:SS       Time limit (default: 04:00:00)
  --cpus N             CPUs and workers per job (default: 32)
  --n-max N            Maximum bins per block (default: 100000)
  --n-permutations N   Permutation replicates (default: 1000)

Fixed: five states, sampling seed 42, permutation seed 1234,
block-scope permutations for all subsystems, no dropped signals.
Set HOI_PYTHON to override Python. No --mem directive is used.
Recorded submissions are skipped to avoid duplicate jobs on rerun.
HELP
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
for value in "$CPUS" "$N_MAX"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Expected positive integer: $value" >&2; exit 2; }
done
[[ "$N_PERMUTATIONS" =~ ^[0-9]+$ ]] || { echo "Invalid permutation count" >&2; exit 2; }
[[ -x "$PYTHON_BIN" && -f "$PYTHON_SCRIPT" ]] || { echo "Missing Python executable or pipeline" >&2; exit 1; }
if (( ! PREPARE_ONLY )); then
    command -v sbatch >/dev/null || { echo "sbatch is not available" >&2; exit 1; }
fi

functional=(promoter_bins.tsv proximal_ctcf_bound_enhancer_purified_bins.tsv
    proximal_not_ctcf_bound_enhancer_purified_bins.tsv
    distal_ctcf_bound_enhancer_purified_bins.tsv
    distal_not_ctcf_bound_enhancer_purified_bins.tsv exon_purified_bins.tsv intron_purified_bins.tsv)
sequence=(CpG_island_bins.tsv LINE_bins.tsv SINE_bins.tsv LTR_bins.tsv)
groups=()
for group in functional sequence; do
    declare -n names="$group"
    joined=""
    for name in "${names[@]}"; do
        path="$FEATURE_DIR/$name"
        [[ -f "$path" ]] || { echo "Missing annotation: $path" >&2; exit 1; }
        joined+="${joined:+;}$path"
    done
    groups+=("$joined")
done
BLOCK_GROUPS="${groups[0]}::BLOCKGROUP::${groups[1]}"
shopt -s nullglob
inputs=("$INPUT_DIR"/*.txt)
(( ${#inputs[@]} > 0 )) || { echo "No input panels found" >&2; exit 1; }
# Check all inputs before preparing or submitting any job.
for input in "${inputs[@]}"; do
    count=0
    while IFS= read -r signal || [[ -n "$signal" ]]; do
        signal="${signal%$'\r'}"
        [[ -n "$signal" ]] || continue
        [[ -f "$signal" ]] || { echo "Missing signal in $input: $signal" >&2; exit 1; }
        count=$((count + 1))
    done < "$input"
    (( count == 7 )) || { echo "Expected seven signals in $input; got $count" >&2; exit 1; }
done

mkdir -p "$RESULTS_DIR"
manifest="$RESULTS_DIR/input_manifest.tsv"
printf 'sample\tinput_file\toutput_dir\n' > "$manifest"
submitted=0
for input in "${inputs[@]}"; do
    sample="$(basename "$input" .txt)"
    output="$RESULTS_DIR/$sample"
    mkdir -p "$output"
    printf '%s\t%s\t%s\n' "$sample" "$input" "$output" >> "$manifest"
    cmd=("$PYTHON_BIN" "$PYTHON_SCRIPT" --input-file "$input" --output-dir "$output"
        --sample-name "$sample" --n-bins 5 --n-max "$N_MAX" --random-seed 42
        --n-jobs "$CPUS" --block-groups "$BLOCK_GROUPS" --drop-signals ""
        --n-permutations "$N_PERMUTATIONS" --permutation-seed 1234
        --permutation-scope block --permutation-level all)
    {
        printf '#!/usr/bin/env bash\n'
        printf '#SBATCH --job-name=oinfo_%s\n' "$sample"
        printf '#SBATCH --partition=%s\n#SBATCH --time=%s\n' "$PARTITION" "$TIME_LIMIT"
        printf '#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=%s\n' "$CPUS"
        printf '#SBATCH --output=%s/slurm_%%j.out\n#SBATCH --error=%s/slurm_%%j.err\n' "$output" "$output"
        printf 'set -euo pipefail\nexport OMP_NUM_THREADS=1\nexport MKL_NUM_THREADS=1\nexport OPENBLAS_NUM_THREADS=1\n'
        printf 'echo "Started: $(date --iso-8601=seconds), host: $(hostname), job: ${SLURM_JOB_ID:-local}"\n'
        printf '%q ' "${cmd[@]}"
        printf '\necho "Completed: $(date --iso-8601=seconds)"\n'
    } > "$output/job.sh"
    bash -n "$output/job.sh"
    if (( PREPARE_ONLY )); then
        continue
    fi
    record="$RESULTS_DIR/submissions.tsv"
    if [[ ! -e "$record" ]]; then
        printf 'sample\tjob_id\tsubmitted_at\tjob_script\n' > "$record"
    fi
    if awk -F '\t' -v name="$sample" 'NR>1 && $1==name{found=1}END{exit !found}' "$record"; then
        echo "Already submitted: $sample"
        continue
    fi
    job_id="$(sbatch --parsable "$output/job.sh")"
    job_id="${job_id%%;*}"
    [[ "$job_id" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $job_id" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\n' "$sample" "$job_id" "$(date --iso-8601=seconds)" "$output/job.sh" >> "$record"
    submitted=$((submitted + 1))
    echo "Submitted $sample: $job_id"
done
if (( PREPARE_ONLY )); then
    echo "Prepared ${#inputs[@]} jobs; none submitted. Output: $RESULTS_DIR"
else
    echo "Submitted $submitted jobs. Output: $RESULTS_DIR"
fi
