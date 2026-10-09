#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR="/lustre/home/2101110354/sequence_entropy/data/mouse_preimplementation/original/methylation_rep_averaged"
OUT_DIR="/lustre/home/2101110354/sequence_entropy/data/mouse_preimplementation/original/methylation_stage_averaged"
LOG_DIR="logs"


mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "[INFO] Detecting stage × experiment schema..."

find "${INPUT_DIR}" -type f -name "*_embryo*.bedGraph" \
| sed -E 's#.*/##' \
| awk -F'.' '
{
    # Example:
    # 2-Cell_embryo3.ACG.TCG.mm10.bedGraph
    # 2-Cell_embryo3.GCA.GCC.GCT.mm10.bedGraph

    split($1, a, "_embryo")
    stage = a[1]

    # experiment schema = everything between embryoX. and .mm10
    schema = ""
    for (i = 2; i <= NF - 2; i++) {
        if (schema == "")
            schema = $i
        else
            schema = schema "." $i
    }

    if (stage != "" && schema != "")
        print stage "\t" schema
}' \
| sort -u \
| while read -r STAGE SCHEMA; do

    echo "[INFO] Submitting ${STAGE} / ${SCHEMA}"

    sbatch \
      --job-name=merge_${STAGE}_${SCHEMA} \
      --output=${LOG_DIR}/merge_${STAGE}_${SCHEMA}.%j.out \
      --error=${LOG_DIR}/merge_${STAGE}_${SCHEMA}.%j.err \
      merge_stage_schema.slurm \
      "${STAGE}" "${SCHEMA}" "${INPUT_DIR}" "${OUT_DIR}"

done
