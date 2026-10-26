#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
python validation/plot_rerun_20261026/mouse/run_development.py
python validation/plot_rerun_20261026/mouse/run_mi.py
python validation/plot_rerun_20261026/mouse/make_report.py
