#!/usr/bin/env bash

set -euo pipefail

IN_DIR="$1"
OUT_DIR="$2"
THRESHOLD="${3:-8}"

mkdir -p logs

for f in "$IN_DIR"/*.bed.gz
do
    echo "Submitting $f"
    sbatch slurm_convert_cpg_to_bedgraph.sh "$f" "$OUT_DIR" "$THRESHOLD"
done
