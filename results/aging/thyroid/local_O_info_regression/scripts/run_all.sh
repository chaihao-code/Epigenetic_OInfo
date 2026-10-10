#!/bin/bash
# Run the full immune-adjusted O-information pipeline.
#   01  per-bin regression  O-info ~ age + immune_fraction
#   02  pathway GSEA (GO BP 2021 + MSigDB Hallmark 2020, TSS gene mapping)
#   03  signal-change GSEA (mark x trend custom gene sets)
#   04  collect significant terms across regions
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

python3 01_regress_oinfo_immune.py
python3 02_gsea_pathway.py
python3 03_gsea_marginal.py
python3 04_summarize.py

echo "ALL DONE"
