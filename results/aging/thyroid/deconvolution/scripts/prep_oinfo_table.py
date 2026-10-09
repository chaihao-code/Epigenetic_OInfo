#!/usr/bin/env python3
"""prep_oinfo_table.py — 汇总原始 O-information 结果为 sample × region 长表.

输出 data/processed/oinfo_long.tsv:
  sample_id, region, hoi, tc, dtc, n_samples(entropy 实际使用), entropy_used
hoi 取自 results_oinfo_n_bins_5 每区域 order=7, n_missing=0 行。
"""
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CFG = yaml.safe_load(open(ROOT / "config.yaml"))


def main():
    odir = Path(CFG["data"]["oinfo_results_dir"])
    rows = []
    for sample in CFG["samples"]["included"]:
        sdir = odir / f"CEEHRC.{sample}"
        for rdir in sorted(sdir.iterdir()):
            f = rdir / f"CEEHRC.{sample}_{rdir.name}_oinfo_overall.csv"
            if not f.exists():
                continue
            df = pd.read_csv(f)
            full = df[(df["order"] == 7) & (df["n_missing"] == 0)]
            if full.empty:
                continue
            r = full.iloc[0]
            meta_f = Path(CFG["data"]["oinfo_results_dir"]).parent / \
                f"results_entropy_n_bins_5/CEEHRC.{sample}/" / \
                f"CEEHRC.{sample}_{rdir.name}_metadata.json"
            n_used = None
            if meta_f.exists():
                import json
                n_used = json.load(open(meta_f)).get("n_samples")
            rows.append(dict(sample_id=sample, region=rdir.name,
                             hoi=float(r["hoi"]), tc=float(r["tc"]), dtc=float(r["dtc"]),
                             n_entropy_samples=n_used))
    out = ROOT / "data/processed/oinfo_long.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    print(f"{len(df)} rows -> {out}")
    print(df.pivot(index="sample_id", columns="region", values="hoi").shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
