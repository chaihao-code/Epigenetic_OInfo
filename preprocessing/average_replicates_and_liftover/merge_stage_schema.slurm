#!/bin/bash
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00

set -euo pipefail

STAGE="$1"
SCHEMA="$2"
INPUT_DIR="$3"
OUT_DIR="$4"

module load bedtools 2>/dev/null || true

TMPDIR=$(mktemp -d)
trap "rm -rf ${TMPDIR}" EXIT

echo "[INFO] Stage: ${STAGE}"
echo "[INFO] Experiment schema: ${SCHEMA}"

ls "${INPUT_DIR}/${STAGE}_embryo"*.${SCHEMA}.mm10.bedGraph \
> "${TMPDIR}/files.txt"

NFILES=$(wc -l < "${TMPDIR}/files.txt")

if [[ "${NFILES}" -lt 2 ]]; then
    echo "[WARN] ${STAGE} ${SCHEMA}: only ${NFILES} file(s), skipping"
    exit 0
fi

echo "[INFO] Merging ${NFILES} files"

bedtools unionbedg \
  -i $(cat "${TMPDIR}/files.txt") \
| awk -v OFS="\t" '
{
    sum = 0;
    n = 0;
    for (i = 4; i <= NF; i++) {
        if ($i != "." && $i != "NA") {
            sum += $i;
            n += 1;
        }
    }
    if (n > 0) {
        print $1, $2, $3, sum / n;
    }
}' \
> "${OUT_DIR}/${STAGE}.${SCHEMA}.stage_mean.bedGraph"

echo "[DONE] ${STAGE} ${SCHEMA}"
