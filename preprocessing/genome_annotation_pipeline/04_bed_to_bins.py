#!/usr/bin/env python3
"""把 BED 区间转成 1-kb bin 的 0/1 mask 表。

用法:
    python3 04_bed_to_bins.py <bed_file> <bins_bed> <feature_name> <output_dir>

输出: <output_dir>/<feature_name>_bins.tsv
      列: bin_id, <feature_name>（1 = 该 bin 与 BED 任一区间有重叠）
"""
import os
import sys

import pandas as pd


def mark_bins(bed_file, bins_file, name, out_dir):
    bins = pd.read_csv(bins_file, sep="\t", header=None, names=["chrom", "start", "end"])
    bins["start"] = bins["start"].astype(int)
    bins["end"] = bins["end"].astype(int)
    bins["bin_id"] = (bins["chrom"] + ":" + bins["start"].astype(str)
                      + "-" + bins["end"].astype(str))
    bins[name] = 0

    bed = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chrom", "start", "end"])
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
    bed = bed.dropna(subset=["start", "end"])

    for chrom in bins["chrom"].unique():
        bins_chrom = bins[bins["chrom"] == chrom]
        bed_chrom = bed[bed["chrom"] == chrom]
        for _, row in bed_chrom.iterrows():
            idx = bins_chrom[(bins_chrom["end"] > row["start"])
                             & (bins_chrom["start"] < row["end"])].index
            bins.loc[idx, name] = 1

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{name}_bins.tsv")
    bins[["bin_id", name]].to_csv(out, sep="\t", index=False)
    print(f"[done] {out}: {int(bins[name].sum())} positive bins")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 04_bed_to_bins.py <bed> <bins_bed> <name> <out_dir>")
        sys.exit(1)
    mark_bins(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
