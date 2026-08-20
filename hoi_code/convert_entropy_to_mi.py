#!/usr/bin/env python3
"""
Convert entropy_all_subsets.csv (from HOI analysis) into MI-format CSVs
required by plot_by_block_all_anchors.py and plot_joint_vs_pairwise_ratio_heatmap_v3.py.

Input: directory containing stage subdirectories, each with {stage}_entropy_all_subsets.csv
Output: directory with {stage}/{stage}_mi_self.csv, {stage}_mi_pairwise.csv, {stage}_mi_joint.csv

Usage:
    python convert_entropy_to_mi.py -i /path/to/results_entropy_n_bins_5 -o /path/to/pairwise_mi_n5
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert entropy subsets to MI-format CSVs"
    )
    parser.add_argument(
        "-i", "--input-dir", required=True,
        help="Input directory containing stage subdirs (e.g., results_entropy_n_bins_5)"
    )
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory for MI-format CSVs (e.g., pairwise_mi_n5)"
    )
    return parser.parse_args()


def process_stage(input_dir: Path, output_dir: Path, stage: str):
    """Process a single stage's entropy_all_subsets.csv and write MI CSVs."""
    entropy_path = input_dir / stage / f"{stage}_entropy_all_subsets.csv"
    metadata_path = input_dir / stage / f"{stage}_metadata.json"

    if not entropy_path.exists():
        print(f"  [SKIP] {entropy_path} not found")
        return

    # Read metadata for n_bins
    n_bins = 5
    if metadata_path.exists():
        with open(metadata_path) as f:
            meta = json.load(f)
        n_bins = meta.get("n_bins", 5)
    else:
        pass

    # Read entropy subsets
    df = pd.read_csv(entropy_path)

    # Get list of blocks
    blocks = df['block'].unique()
    print(f"  Stage {stage}: {len(blocks)} blocks, n_bins={n_bins}")

    self_rows = []
    pairwise_rows = []
    joint_rows = []

    for block in blocks:
        block_df = df[df['block'] == block].copy()

        # Build lookup: frozenset of feature names -> entropy
        ent_map = {}
        for _, row in block_df.iterrows():
            names = frozenset(row['subset_names'].split('/'))
            ent_map[names] = row['entropy']

        # Self entropy (order 1)
        self_map = {}
        for names, H in ent_map.items():
            if len(names) == 1:
                feat = list(names)[0]
                self_map[feat] = H
                self_rows.append({
                    'feature': feat,
                    'H_X': H,
                    'n_bins': n_bins,
                    'block': block,
                })

        # Pairwise MI (order 2)
        for names, H_xy in ent_map.items():
            if len(names) == 2:
                f1, f2 = sorted(list(names))
                H_x = self_map.get(f1, np.nan)
                H_y = self_map.get(f2, np.nan)
                I_xy = H_x + H_y - H_xy if not (np.isnan(H_x) or np.isnan(H_y)) else np.nan
                pairwise_rows.append({
                    'feature_1': f1,
                    'feature_2': f2,
                    'H_X': H_x,
                    'H_Y': H_y,
                    'H_XY': H_xy,
                    'I_XY': I_xy,
                    'n_bins': n_bins,
                    'block': block,
                })

        # Joint MI: I(X_i; X_rest)
        all_features = sorted(self_map.keys())
        n_features = len(all_features)
        if n_features < 2:
            continue

        # Full joint entropy (order n)
        full_set = frozenset(all_features)
        H_full = ent_map.get(full_set, np.nan)

        for feat in all_features:
            H_x = self_map.get(feat, np.nan)
            # Rest set
            rest = frozenset([f for f in all_features if f != feat])
            H_rest = ent_map.get(rest, np.nan)
            I_joint = H_x + H_rest - H_full if not (np.isnan(H_x) or np.isnan(H_rest) or np.isnan(H_full)) else np.nan

            joint_rows.append({
                'feature': feat,
                'H_X': H_x,
                'H_X_minus_i': H_rest,
                'H_X_joint': H_full,
                'I_X_joint': I_joint,
                'n_other_features': n_features - 1,
                'n_joint_states_observed': np.nan,  # not available from entropy CSV
                'n_bins': n_bins,
                'block': block,
            })

    # Write output
    stage_out = output_dir / stage
    stage_out.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(self_rows).to_csv(stage_out / f"{stage}_mi_self.csv", index=False)
    pd.DataFrame(pairwise_rows).to_csv(stage_out / f"{stage}_mi_pairwise.csv", index=False)
    pd.DataFrame(joint_rows).to_csv(stage_out / f"{stage}_mi_joint.csv", index=False)

    print(f"    Written: {len(self_rows)} self, {len(pairwise_rows)} pairwise, {len(joint_rows)} joint rows")


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect stages: subdirectories that contain {subdir}_entropy_all_subsets.csv
    stages = []
    for subdir in sorted(input_dir.iterdir()):
        if subdir.is_dir() and (subdir / f"{subdir.name}_entropy_all_subsets.csv").exists():
            stages.append(subdir.name)

    if not stages:
        print(f"No stage directories found in {input_dir}")
        return

    print(f"Found stages: {stages}")

    for stage in stages:
        process_stage(input_dir, output_dir, stage)

    print(f"\nAll done! Output: {output_dir}")


if __name__ == "__main__":
    main()
