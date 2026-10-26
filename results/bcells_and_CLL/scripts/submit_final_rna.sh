#!/usr/bin/env bash
set -euo pipefail
PROJECT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export BCELLS_PROJECT_DIR="$PROJECT"
stage=${1:-all}
case "$stage" in de|gsea|figures|all) ;; *) echo "Usage: $0 {de|gsea|figures|all}" >&2; exit 2;; esac
mkdir -p "$PROJECT/logs"
sbatch --output="$PROJECT/logs/final_rna_%j.out" --error="$PROJECT/logs/final_rna_%j.err" "$PROJECT/scripts/run_final_rna.slurm" "$stage"
