#!/usr/bin/env python3
"""01_inventory.py — Step 1/2: 本地输入数据盘点.

扫描甲状腺 6 样本的全部输入: 七轨道 bigWig、1kb bins、注释 mask、
离散化矩阵、O-info 结果, 用 pyBigWig 实测验证基因组版本/数值范围/零值语义,
输出 input_inventory.tsv 与 inventory_summary.md。
"""

# Portable repository paths; no analysis is executed by this helper.
from pathlib import Path as _ReleasePath
_RELEASE_ROOT = next(p for p in _ReleasePath(__file__).resolve().parents
                     if (p / ".epigenetic_oinfo_root").is_file())
def _release_path(value):
    return value.replace("@ROOT@", str(_RELEASE_ROOT))
def _release_config(value):
    if isinstance(value, dict):
        return {k: _release_config(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_config(v) for v in value]
    return _release_path(value) if isinstance(value, str) else value
import csv
import json
import logging
import sys
from pathlib import Path

import pyBigWig
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/01_inventory.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("inventory")

ROOT = Path(__file__).resolve().parents[1]
CFG = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
BW_DIR = Path(CFG["data"]["bigwig_dir"])
OUT = ROOT / "results/tables/input_inventory.tsv"

MARKS = {
    "DNA_methylation": "WGB-Seq.methylation_profile.bigWig",
    "H3K4me1": "H3K4me1.signal_unstranded.bigWig",
    "H3K4me3": "H3K4me3.signal_unstranded.bigWig",
    "H3K9me3": "H3K9me3.signal_unstranded.bigWig",
    "H3K27me3": "H3K27me3.signal_unstranded.bigWig",
    "H3K27ac": "H3K27ac.signal_unstranded.bigWig",
    "H3K36me3": "H3K36me3.signal_unstranded.bigWig",
}

# GRCh38 chr21 长度用于判定基因组版本 (hg19 chr21 = 48129895)
GRCH38_CHR21 = 46709983
PROBE_REGIONS = [("chr1", 1_000_000, 1_100_000), ("chr21", 10_000_000, 10_100_000)]


def probe_bigwig(path: Path, is_methylation: bool):
    """返回 dict: 染色体命名/版本判定/数值范围/零与缺失统计。
    pyBigWig 打开失败(截断文件)时尝试 .1 后缀的完整副本。"""
    info = {"exists": path.exists()}
    if not info["exists"]:
        return info
    bw = None
    try:
        bw = pyBigWig.open(str(path))
        if not bw.isBigWig():
            raise RuntimeError("not bigWig")
    except Exception as e:
        info["open_error"] = f"{type(e).__name__}: {e}"
        alt = path.parent / (path.name + ".1")
        if alt.exists():
            info["fallback_file"] = str(alt)
            try:
                bw = pyBigWig.open(str(alt))
                info["used_fallback"] = 1
            except Exception:
                bw = None
        if bw is None:
            log.error("cannot open %s (%s)", path, info["open_error"])
            return info
    chroms = bw.chroms()
    info["n_chroms"] = len(chroms)
    info["chr_prefix"] = all(c.startswith("chr") for c in list(chroms)[:25])
    info["chr21_len"] = chroms.get("chr21")
    info["build_inferred"] = "hg38" if chroms.get("chr21") == GRCH38_CHR21 else "unknown/other"
    # 抽两个 100kb 窗口统计值分布、零值与 nan
    n_vals = n_zero = n_nan = 0
    vmin, vmax = float("inf"), float("-inf")
    for chrom, s, e in PROBE_REGIONS:
        if chrom not in chroms or chroms[chrom] <= s:
            continue
        vals = bw.values(chrom, s, min(e, chroms[chrom]))
        for v in vals:
            if v is None or v != v:
                n_nan += 1
                continue
            n_vals += 1
            vmin = min(vmin, v)
            vmax = max(vmax, v)
            if v == 0:
                n_zero += 1
    bw.close()
    info.update(
        n_probe_vals=n_vals,
        n_zero=n_zero,
        n_nan=n_nan,
        zero_frac=round(n_zero / n_vals, 4) if n_vals else None,
        probe_min=round(vmin, 3) if n_vals else None,
        probe_max=round(vmax, 3) if n_vals else None,
    )
    if is_methylation and n_vals:
        # 0-100 百分制判定
        info["unit_inferred"] = "percent_0_100" if vmax > 1.5 else "fraction_0_1"
    return info


def main():
    rows = []
    meta_path = Path(CFG["samples"]["metadata"])
    meta = {}
    with open(meta_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            meta.setdefault(r["sample_id"], {"age": r["age"], "sex": r["sex"]})

    for sample in CFG["samples"]["included"]:
        for mark, suffix in MARKS.items():
            hits = sorted(BW_DIR.glob(f"*CEEHRC.{sample}.{suffix}"))
            if not hits:
                hits = sorted(BW_DIR.glob(f"*{sample}.*{suffix}"))
            if not hits:
                log.warning("MISSING %s %s", sample, mark)
                rows.append(dict(sample_id=sample, signal=mark, path="", exists=0))
                continue
            p = hits[0]
            inf = probe_bigwig(p, mark == "DNA_methylation")
            rows.append(
                dict(
                    sample_id=sample,
                    age=meta.get(sample, {}).get("age"),
                    sex=meta.get(sample, {}).get("sex"),
                    signal=mark,
                    path=str(p),
                    **inf,
                )
            )
            log.info("%s %s build=%s unit=%s zero_frac=%s",
                     sample, mark, inf.get("build_inferred"), inf.get("unit_inferred"), inf.get("zero_frac"))

    # 结构文件
    structural = [
        ("1kb_bins", CFG["data"]["bins_1kb"]),
        ("gencode_gtf_hg38", CFG["data"]["gencode_gtf"]),
        ("oinfo_results_dir", CFG["data"]["oinfo_results_dir"]),
        ("local_oinfo_dir", CFG["data"]["local_oinfo_dir"]),
        ("go_bp_summary", CFG["data"]["go_summary"]),
        ("discretized_dir", CFG["data"]["discretized_dir"]),
        ("genome_features_dir", CFG["data"]["genome_features_dir"]),
    ]
    for name, p in structural:
        p = Path(p)
        rows.append(dict(sample_id="(structural)", signal=name, path=str(p),
                         exists=1 if p.exists() else 0))

    # 离散化文件清点
    ddir = Path(CFG["data"]["discretized_dir"])
    if ddir.exists():
        for sample in CFG["samples"]["included"]:
            n = len(list(ddir.glob(f"*CEEHRC.{sample}.*discrete.bed*")))
            rows.append(dict(sample_id=sample, signal="discrete_bed_files", path=str(ddir),
                             exists=1, n_files=n))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sorted({k for r in rows for k in r}),
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log.info("wrote %s (%d rows)", OUT, len(rows))

    # 汇总 md
    builds = {r.get("build_inferred") for r in rows if r.get("build_inferred")}
    units = {r.get("unit_inferred") for r in rows if r.get("unit_inferred")}
    missing = [(r["sample_id"], r["signal"]) for r in rows if r.get("exists") == 0 and r["sample_id"] != "(structural)"]
    md = ROOT / "results/qc/inventory_summary.md"
    md.write_text(
        "# Step 1/2 数据盘点摘要\n\n"
        f"- 样本: {CFG['samples']['included']}(原分析排除 {CFG['samples']['excluded']})\n"
        f"- bigWig 基因组版本(实测 chr21 长度): {builds}\n"
        f"- 甲基化单位(实测): {units}\n"
        f"- 缺失文件: {missing if missing else '无'}\n"
        "- 关键限制: methylation_profile bigWig 为 0-100 百分制且无 coverage,"
        " 无法运行 UXM_deconv; 采用 constrained NNLS beta 混合近似并单独报告。\n"
        "- 零值语义: zero_frac 列见 input_inventory.tsv; bigWig 无值区间返回 NaN(缺失),"
        " 0 为实测信号零值; WGBS 未覆盖 CpG 处以 NaN/无记录表示。\n"
    )
    print(md.read_text())


if __name__ == "__main__":
    sys.exit(main())
