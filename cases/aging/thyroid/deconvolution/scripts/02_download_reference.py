#!/usr/bin/env python3
"""02_download_reference.py — Step 3: 下载 GSE186458 (Loyfer 2023) 参考数据.

下载内容(均为 hg38, 与目标样本同构建, 无需 liftOver):
  1. 每细胞类型 .hg38.beta (wgbstools 二进制, #meth/#coverage uint8 对, 保留覆盖度)
  2. 系列级 atlas blocks GSE186458_blocks.s207.hg38.bed.gz (免疫特异性区域分析用)
  3. UXM 现成 marker 列表 Markers.U250.hg38.tsv (UXM_deconv GitHub)
  4. UCSC hg38.2bit (用于把 beta 数组索引映射为基因组坐标)

输出: data/reference/, reference_download_manifest.tsv, logs/02_download_reference.log
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
import logging
import subprocess
import sys
import time
from pathlib import Path

import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/02_download_reference.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("download")

ROOT = Path(__file__).resolve().parents[1]
CFG = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
REF = ROOT / CFG["reference"]["dir"]
REF.mkdir(parents=True, exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/{gsm}/suppl/{fname}"

# 细粒度 Model A 所需细胞类型 -> (GEO cell-type label, [(gsm, title)])
CELL_TYPES = {
    "Thyroid_epithelial": [
        ("GSM5652264", "Thyroid-Epithelial-Z0000042S"),
        ("GSM5652265", "Thyroid-Epithelial-Z0000042T"),
        ("GSM5652266", "Thyroid-Epithelial-Z0000042U"),
    ],
    "T_CD4": [("GSM5652279", "Blood-T-CD4-Z000000P"), ("GSM5652280", "Blood-T-CD4-Z000000Q"),
              ("GSM5652281", "Blood-T-CD4-Z000000R")],
    "T_CD8": [("GSM5652282", "Blood-T-CD8-Z000000S"), ("GSM5652283", "Blood-T-CD8-Z000000T"),
              ("GSM5652284", "Blood-T-CD8-Z000000U")],
    "B_cell": [("GSM5652316", "Blood-B-Z000000TX"), ("GSM5652317", "Blood-B-Z000000TY"),
               ("GSM5652318", "Blood-B-Z000000TZ")],
    "Monocyte": [("GSM5652302", "Blood-Monocytes-Z000000U8"), ("GSM5652303", "Blood-Monocytes-Z000000U9"),
                 ("GSM5652304", "Blood-Monocytes-Z000000UA")],
    "NK": [("GSM5652299", "Blood-NK-Z000000U5"), ("GSM5652300", "Blood-NK-Z000000U6"),
           ("GSM5652301", "Blood-NK-Z000000U7")],
    "Granulocyte": [("GSM5652313", "Blood-Granulocytes-Z000000TU"), ("GSM5652314", "Blood-Granulocytes-Z000000TV"),
                    ("GSM5652315", "Blood-Granulocytes-Z000000TW")],
    "Endothelial": [("GSM5652179", "Aorta-Endothel-Z0000043A"), ("GSM5652180", "Aorta-Endothel-Z0000043B"),
                    ("GSM5652181", "Aorta-Endothel-Z0000043C")],
    "Fibroblast": [("GSM5652198", "Colon-Fibroblast-Z0000043P"), ("GSM5652199", "Colon-Fibroblast-Z0000043Q"),
                   ("GSM5652200", "Colon-Fibroblast-Z0000043R")],
}
# 注意: 部分 title 为调查时记录的推测命名, 下载时若 404 则列目录修正文件名(见 download 函数)。

SERIES_FILES = {
    "GSE186458_blocks.s207.hg38.bed.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/suppl/GSE186458_blocks.s207.hg38.bed.gz",
}

GITHUB_FILES = {
    "Markers.U250.hg38.tsv":
        "https://raw.githubusercontent.com/nloyfer/UXM_deconv/main/supplemental/Markers.U250.hg38.tsv",
    "Atlas.U250.l4.hg38.full.tsv":
        "https://raw.githubusercontent.com/nloyfer/UXM_deconv/main/supplemental/Atlas.U250.l4.hg38.full.tsv",
}

HG38_2BIT = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"


def download(url: str, out: Path, retries: int = 3) -> dict:
    rec = {"url": url, "path": str(out), "status": "failed"}
    if out.exists() and out.stat().st_size > 0:
        rec["status"] = "exists"
        rec["size_bytes"] = out.stat().st_size
        return rec
    for i in range(retries):
        try:
            t0 = time.time()
            subprocess.run(["curl", "-fSL", "--retry", "3", "-o", str(out), url],
                           check=True, capture_output=True, timeout=3600)
            rec.update(status="downloaded", size_bytes=out.stat().st_size,
                       seconds=round(time.time() - t0, 1))
            log.info("OK %s (%s bytes)", out.name, rec["size_bytes"])
            return rec
        except subprocess.CalledProcessError as e:
            log.warning("attempt %d failed for %s: %s", i + 1, url,
                        e.stderr.decode()[-200:] if e.stderr else e)
            time.sleep(10 * (i + 1))
    return rec


def main():
    manifest = []
    beta_dir = REF / "beta_hg38"
    beta_dir.mkdir(exist_ok=True)

    for ct, samples in CELL_TYPES.items():
        for gsm, title in samples:
            fname = f"{gsm}_{title}.hg38.beta"
            url = FTP.format(gsm=gsm, fname=fname)
            out = beta_dir / f"{ct}__{gsm}.hg38.beta"
            rec = download(url, out)
            if rec["status"] == "failed":
                # 列目录找真实文件名
                try:
                    listing = subprocess.run(
                        ["curl", "-fsSL", f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5652nnn/{gsm}/suppl/"],
                        capture_output=True, timeout=120).stdout.decode()
                    import re
                    cands = re.findall(r'href="([^"]+\.hg38\.beta)"', listing)
                    if cands:
                        real = cands[0]
                        log.info("retry %s with real filename %s", gsm, real)
                        rec = download(FTP.format(gsm=gsm, fname=real), out)
                        rec["note"] = f"filename corrected: {real}"
                except Exception as e:
                    log.error("listing failed for %s: %s", gsm, e)
            rec["cell_type"] = ct
            rec["gsm"] = gsm
            rec["geo_label"] = title
            manifest.append(rec)

    for fname, url in SERIES_FILES.items():
        rec = download(url, REF / fname)
        rec["cell_type"] = "(series)"; rec["gsm"] = ""
        manifest.append(rec)

    for fname, url in GITHUB_FILES.items():
        rec = download(url, REF / fname)
        rec["cell_type"] = "(github)"; rec["gsm"] = ""
        manifest.append(rec)

    rec = download(HG38_2BIT, REF / "hg38.2bit")
    rec["cell_type"] = "(ucsc)"; rec["gsm"] = ""
    manifest.append(rec)

    import csv
    mpath = ROOT / "results/tables/reference_download_manifest.tsv"
    with open(mpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_type", "gsm", "geo_label", "url", "path",
                                          "status", "size_bytes", "seconds", "note"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
    ok = sum(1 for r in manifest if r["status"] in ("downloaded", "exists"))
    log.info("manifest: %d/%d ok -> %s", ok, len(manifest), mpath)
    return 0 if ok == len(manifest) else 1


if __name__ == "__main__":
    sys.exit(main())
