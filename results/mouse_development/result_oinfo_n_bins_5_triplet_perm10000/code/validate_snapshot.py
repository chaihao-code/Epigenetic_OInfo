#!/usr/bin/env python3
"""Check that selecting triplet nulls leaves their numerical values unchanged."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="triplet-null-validation-") as tmp:
        root = Path(tmp)
        rng = np.random.default_rng(781)
        tracks = []
        for mark in ["ATAC", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "CpG_methyl"]:
            path = root / f"fixture.{mark}.bin1000.bp"
            values = rng.integers(0, 10, size=24)
            path.write_text("".join(f"chr1\t{i*1000}\t{(i+1)*1000}\t{v}\n" for i, v in enumerate(values)))
            tracks.append(str(path))
        inputs = root / "inputs.txt"
        inputs.write_text("\n".join(tracks) + "\n")
        mask = root / "mask.tsv"
        mask.write_text("bin_id\tfixture_region\n" + "".join(
            f"chr1:{i*1000}-{(i+1)*1000}\t{int(i<12)}\n" for i in range(24)))
        for name, script in [("original", args.original), ("snapshot", args.snapshot)]:
            cmd = [sys.executable, str(script), "--input-file", str(inputs),
                   "--output-dir", str(root / name), "--sample-name", "fixture",
                   "--n-bins", "5", "--n-max", "20", "--random-seed", "42",
                   "--block-files", str(mask), "--drop-signals", "",
                   "--n-permutations", "3", "--permutation-seed", "1234",
                   "--permutation-level", "all", "--permutation-scope", "block",
                   "--permutation-max-order", "3", "--n-jobs", "1"]
            if name == "snapshot":
                cmd += ["--permutation-min-order", "3"]
            subprocess.run(cmd, check=True)
        old, new = root / "original", root / "snapshot"
        old_subsets = pd.read_csv(old / "fixture_subsets.csv")
        new_subsets = pd.read_csv(new / "fixture_subsets.csv")
        assert len(old_subsets) == 56 and len(new_subsets) == 35
        assert new_subsets["order"].eq(3).all()
        expected = np.load(old / "fixture_permutation_null.npy")[:, :, old_subsets["order"].to_numpy() == 3]
        actual = np.load(new / "fixture_permutation_null.npy")
        np.testing.assert_array_equal(actual, expected)
        assert (old / "fixture_oinfo.csv").read_bytes() == (new / "fixture_oinfo.csv").read_bytes()
        assert (old / "fixture_entropy_all_subsets.csv").read_bytes() == (new / "fixture_entropy_all_subsets.csv").read_bytes()
        summary = pd.read_csv(new / "fixture_null_summary.csv")
        assert len(summary) == 70 and summary["order"].eq(3).all()
        report = dict(status="passed", n_signals=7, n_blocks=2, n_permutations=3,
                      original_null_subsystems=56, snapshot_null_subsystems=35,
                      null_values_bitwise_equal=True, observed_outputs_bytewise_equal=True,
                      triplet_only_summary=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
