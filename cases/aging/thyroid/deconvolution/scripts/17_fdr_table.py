#!/usr/bin/env python3
"""17_fdr_table.py — 给 Q4 回归结果加多重检验校正列, 生成最终 FDR 表.

输入: results/tables/oinfo_age_adjustment.tsv (07 的输出, 主分析 Panel A/Model A)
输出: results/tables/oinfo_age_fdr_table.tsv
  按 p_adj 排序, 含 nominal p 与两种 FDR (BH / two-stage BKY),
  供 methods_zh 3.2 节表格与 fig3/fig4 使用。
"""
import sys
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"

COLS = ["region", "beta_age_unadj", "p_unadj", "q_bh_unadj", "q_bky_unadj",
        "beta_age_adj", "ci_lo_adj", "ci_hi_adj", "p_adj", "q_bh_adj",
        "q_bky_adj", "attenuation"]


def main():
    oa = pd.read_csv(TBL / "oinfo_age_adjustment.tsv", sep="\t")
    m = oa[(oa.panel == "A") & (oa.model == "A")].copy()
    _, m["q_bh_unadj"], _, _ = multipletests(m.p_unadj.values, 0.05, "fdr_bh")
    _, m["q_bky_unadj"], _, _ = multipletests(m.p_unadj.values, 0.05, "fdr_tsbky")
    _, m["q_bh_adj"], _, _ = multipletests(m.p_adj.values, 0.05, "fdr_bh")
    _, m["q_bky_adj"], _, _ = multipletests(m.p_adj.values, 0.05, "fdr_tsbky")
    m = m.sort_values("p_adj")
    m[COLS].to_csv(TBL / "oinfo_age_fdr_table.tsv", sep="\t", index=False,
                   float_format="%.6g")
    print("wrote oinfo_age_fdr_table.tsv | rows:", len(m),
          "| FDR<0.05 (adj, BH):", int((m.q_bh_adj < 0.05).sum()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
