#!/bin/bash
#SBATCH --job-name=triplet_summary
#SBATCH --partition=C064M0256G
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
set -euo pipefail
PROJECT_DIR="${BCELLS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
export BCELLS_DATA_ROOT="${BCELLS_DATA_ROOT:-$(cd "$PROJECT_DIR/../.." && pwd)/data}"
if [[ -z "${BCELLS_PYTHON:-}" && -x "$PROJECT_DIR/.rna_de_env/bin/python" ]]; then
    BCELLS_PYTHON="$PROJECT_DIR/.rna_de_env/bin/python"
fi
PY="${BCELLS_PYTHON:-python3}"
"$PY" "$PROJECT_DIR/scripts/summarize_triplet_significance.py"
