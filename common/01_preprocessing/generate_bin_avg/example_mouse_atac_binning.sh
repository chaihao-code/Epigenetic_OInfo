#!/usr/bin/env bash
set -euo pipefail
EPINFO_ROOT="${EPINFO_ROOT:?Set EPINFO_ROOT to the release root}"
bash "$EPINFO_ROOT/common/01_preprocessing/generate_bin_avg/batch_compute_bin_average.sh" \
  "$EPINFO_ROOT/data/mouse_preimplementation/original/atac" \
  "$EPINFO_ROOT/cases/mouse_preimp/outputs/binned_atac" \
  "$EPINFO_ROOT/data/genome_features/mm10/mm10.chrom.sizes" 1000 1000
