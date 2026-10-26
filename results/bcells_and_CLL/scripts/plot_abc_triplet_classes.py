#!/usr/bin/env python3
"""Reference-style ABC heatmaps, retaining every R/S/pB donor as a row."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/bcells-triplet-mpl")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/bcells-plot-cache")
from itertools import combinations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "ABC_triplet_classes"
PLOTS = OUT / "plots"
SIGNALS = ["H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac", "WGB-Seq"]
ANCHORS = {"H3K9me3", "H3K27me3", "H3K36me3"}
REGIONS = ["promoter_purified_bins", "proximal_ctcf_bound_enhancer_purified_bins",
    "proximal_not_ctcf_bound_enhancer_purified_bins", "distal_ctcf_bound_enhancer_purified_bins",
    "distal_not_ctcf_bound_enhancer_purified_bins", "exon_purified_bins", "intron_purified_bins", "None_func"]
FEATURES = ["CpG_island", "LINE", "SINE", "LTR", "None_rep"]
BLOCKS = [r + "_AND_" + f for r in REGIONS for f in FEATURES]
REGION_LABELS = ["Promoter", "Proximal\nenh. CTCF+", "Proximal\nenh. CTCF−", "Distal\nenh. CTCF+", "Distal\nenh. CTCF−", "Exon", "Intron", "None"]


def classify(signals):
    signals = set(signals)
    if signals & ANCHORS:
        return "A" if "WGB-Seq" in signals else "B"
    return "C"


def canonical(names):
    found = [next(s for s in SIGNALS if s in name) for name in names.split("/")]
    assert len(set(found)) == 3, names
    return tuple(s for s in SIGNALS if s in found)


def render(matrix, row_labels, class_name, capacity, stem, *, pooled=False, criterion="nominal"):
    n = len(row_labels)
    height_cm = 14 if not pooled else 7.5
    fig = plt.figure(figsize=(18/2.54, height_cm/2.54))
    bottom, top = (.30, .90) if not pooled else (.51, .78)
    ax = fig.add_axes([.17, bottom, .74, top-bottom])
    im = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=capacity, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(n), row_labels)
    ax.set_xticks([])
    ax.tick_params(axis="y", length=0)
    if not pooled:
        ax.axhline(12.5, color="black", lw=.8)
    for k in range(1, 8):
        ax.axvline(k*5-.5, color="gray", lw=.5)
    for i in range(n):
        for j in range(40):
            value = matrix[i, j]
            label = f"{value:.1f}" if pooled else str(int(value))
            ax.text(j, i, label, ha="center", va="center", fontsize=3.5,
                color="white" if value >= capacity*.62 else "black")
    arrow_y, region_y, feature_y = (-.025, -.06, -.17) if not pooled else (-.06, -.16, -.48)
    for k, label in enumerate(REGION_LABELS):
        ax.annotate("", xy=(k*5+4.4, arrow_y), xytext=(k*5-.4, arrow_y),
            xycoords=ax.get_xaxis_transform(), arrowprops=dict(arrowstyle="<->", color="gray", lw=.5))
        ax.text(k*5+2, region_y, label, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5)
    for j in range(40):
        ax.text(j, feature_y, ["CpG island", "LINE", "SINE", "LTR", "None"][j % 5],
            transform=ax.get_xaxis_transform(), ha="center", va="top", rotation=90, fontsize=5)
    cax = fig.add_axes([.935, bottom, .013, top-bottom])
    ticks = sorted(set([0, round(capacity/3), round(2*capacity/3), capacity]))
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    cb.set_label(("Mean number" if pooled else "Number") + " of synergistic triplets", fontsize=5)
    cb.ax.tick_params(labelsize=5)
    title = f"Class {class_name} synergistic triplets across B-cell samples"
    if pooled:
        title = f"Class {class_name} synergistic triplets: R, S and pB group means"
    fig.text(.17, .95 if not pooled else .91, title, fontsize=7)
    footer = ("R, S and pB pooled separately; within-group unweighted donor means."
        if pooled else "R and S treated as parallel CLL samples; each row is one donor.")
    fig.text(.17, .055, footer, fontsize=5)
    threshold = "nominal permutation P < 0.05" if criterion == "nominal" else "global BH FDR < 0.05"
    fig.text(.17, .030, f"O-information < 0; lower-tail {threshold}; {capacity} possible class-{class_name} triplets per region.", fontsize=4.7)
    for ext in ["pdf", "png"]:
        fig.savefig(PLOTS / f"{stem}.{ext}", dpi=300)
    plt.close(fig)


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(PROJECT / "results_oinfo_triplet_n_bins_5/all_samples_triplet_significance.tsv", sep="\t")
    clinical = pd.read_csv(PROJECT / "data_info/cll_metadata.tsv", sep="\t", comment="#")
    info = pd.read_csv(PROJECT / "data_info/ihec_rna_lookup/rna_sample_metadata.tsv", sep="\t", keep_default_na=False)
    group = dict(zip(clinical["sample"], clinical["group"]))
    samples = sorted(df["sample"].unique())
    pb = ["CEMT0129", "CEMT0134", "CEMT0135", "CEMT0138"]
    mapping = []
    for sample in samples:
        short = sample.removeprefix("CEEHRC.")
        if short in group:
            label_group = group[short].replace("CLL-", "")
            pool = label_group
        else:
            assert short in pb
            metadata = info[info.sample_id.eq(short)]
            assert len(metadata) == 1 and metadata.iloc[0]["cell_type"] == "B cells"
            assert metadata.iloc[0]["disease"] == "None"
            label_group, pool = "pB", "pB"
        mapping.append(dict(sample=sample, sample_id=short, label_group=label_group,
            pooled_group=pool, display_label=f"{label_group}  {short}"))
    mapping = pd.DataFrame(mapping)
    mapping["display_order"] = mapping.label_group.map({"R": 0, "S": 1, "pB": 2})
    mapping = mapping.sort_values(["display_order", "sample_id"]).reset_index(drop=True)
    assert mapping.label_group.value_counts().to_dict() == {"R": 8, "S": 5, "pB": 4}
    mapping.to_csv(OUT / "sample_labels.tsv", sep="\t", index=False)
    capacities = pd.Series([classify(t) for t in combinations(SIGNALS, 3)]).value_counts().to_dict()
    assert capacities == {"B": 19, "A": 12, "C": 4}
    df["triplet_signals"] = df.subset_names.map(canonical)
    df["triplet"] = df.triplet_signals.map(lambda x: "+".join(x))
    df["class"] = df.triplet_signals.map(classify)
    assert len(df) == 23800 and df.order.eq(3).all() and df.n_perm.eq(1000).all()
    assert not df.duplicated(["sample", "block", "triplet"]).any()
    assert set(df.block) == set(BLOCKS)
    df["synergistic_nominal"] = df.observed_hoi.lt(0) & df.p_down.lt(.05)
    df["synergistic_FDR005"] = df.observed_hoi.lt(0) & df.q_down_global_BH.lt(.05)
    counts = df.groupby(["sample", "block", "class"]).agg(
        n_tested=("triplet", "size"), count_nominal=("synergistic_nominal", "sum"),
        count_FDR005=("synergistic_FDR005", "sum")).reset_index()
    assert counts.n_tested.eq(counts["class"].map(capacities)).all()
    assert len(counts) == 17*40*3
    counts = counts.merge(mapping, on="sample", validate="many_to_one")
    counts["fraction_nominal"] = counts.count_nominal / counts.n_tested
    counts["fraction_FDR005"] = counts.count_FDR005 / counts.n_tested
    counts.to_csv(OUT / "counts_by_sample_region_class.tsv", sep="\t", index=False)
    df.drop(columns="triplet_signals").to_csv(OUT / "all_triplet_classes_and_tests.tsv", sep="\t", index=False)
    pd.DataFrame([dict(triplet="+".join(t), class_name=classify(t))
        for t in combinations(SIGNALS, 3)]).to_csv(OUT / "triplet_class_definitions.tsv", sep="\t", index=False)

    font = Path(os.environ.get("BCELLS_ARIAL_FONT", str(Path.home() / ".local/share/fonts/arial.ttf")))
    if font.exists():
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({"font.family": "Arial" if font.exists() else "DejaVu Sans",
        "font.size": 5, "pdf.fonttype": 42, "ps.fonttype": 42})
    for cl in "ABC":
        sub = counts[counts["class"].eq(cl)]
        for metric, suffix, criterion in [("count_nominal", "", "nominal"), ("count_FDR005", "_FDR005", "FDR")]:
            matrix = sub.pivot(index="sample", columns="block", values=metric).reindex(index=mapping["sample"], columns=BLOCKS)
            assert matrix.notna().all().all()
            export = matrix.copy()
            export.index = mapping.display_label
            export.index.name = "sample_label"
            export.to_csv(PLOTS / f"{cl}_sample_counts{suffix}.tsv", sep="\t")
            render(matrix.to_numpy(), mapping.display_label.tolist(), cl, capacities[cl], f"{cl}_sample_heatmap{suffix}", criterion=criterion)
        mean_counts = sub.groupby(["pooled_group", "block"]).agg(
            n_samples=("sample", "nunique"), mean_count=("count_nominal", "mean"),
            std_count=("count_nominal", "std"), mean_count_FDR005=("count_FDR005", "mean")).reset_index()
        mean_counts["group_order"] = mean_counts.pooled_group.map({"R": 0, "S": 1, "pB": 2})
        mean_counts["block_order"] = mean_counts.block.map({b: i for i, b in enumerate(BLOCKS)})
        mean_counts = mean_counts.sort_values(["group_order", "block_order"]).drop(columns=["group_order", "block_order"])
        assert mean_counts.n_samples.eq(mean_counts.pooled_group.map({"R": 8, "S": 5, "pB": 4})).all()
        mean_counts.to_csv(PLOTS / f"{cl}_pooled_average_counts.tsv", sep="\t", index=False)
        mean_matrix = mean_counts.pivot(index="pooled_group", columns="block", values="mean_count").reindex(index=["R", "S", "pB"], columns=BLOCKS)
        render(mean_matrix.to_numpy(), ["R (n=8)", "S (n=5)", "pB (n=4)"], cl, capacities[cl], f"{cl}_pooled_average_heatmap", pooled=True)

    config = dict(source="results_oinfo_triplet_n_bins_5/all_samples_triplet_significance.tsv",
        reference="project/H3K36me3_triplet_robustness/ABC_triplet_classes/plots",
        signals=SIGNALS, anchor_signals=sorted(ANCHORS), class_capacities=capacities,
        classes={"A": "DNA methylation present AND at least one anchor signal",
            "B": "DNA methylation absent AND at least one anchor signal", "C": "No anchor signals"},
        samples={"R": 8, "S": 5, "pB": 4}, R_S_handling="Parallel donors; no R/S group significance comparison. R, S and pB pooled separately with equal donor weights within each group.",
        pB_mapping=pb, pB_annotation_source="data_info/ihec_rna_lookup/rna_sample_metadata.tsv: B cells, disease None, origin Tonsil",
        primary_criterion="observed_hoi < 0 AND p_down < 0.05, same nominal lower-tail criterion as reference",
        supplementary_criterion="observed_hoi < 0 AND q_down_global_BH < 0.05; inherited global BH family of 23800 tests",
        pooled_combined_criterion="observed_hoi < 0 AND q_down_global_BH < 0.05; global BH family of 23800 tests",
        pooled_combined_script="../scripts/plot_abc_pooled_combined.py",
        capacity_note="Seven-signal B-cell universe lacks ATAC present in eight-signal reference; capacities are 12/19/4 instead of 15/31/10.")
    (OUT / "plot_config.json").write_text(json.dumps(config, indent=2))
    (OUT / "README_zh.md").write_text("""# B 细胞 ABC 协同 triplet 热图

主图在 plots/A_sample_heatmap、B_sample_heatmap、C_sample_heatmap（PNG/PDF）。
纵轴依次为 8 个 R、5 个 S、4 个 pB 样本，每行一个样本；R/S 仅为标签，作为平行 CLL 样本处理，未做 R/S 差异检验。
横轴为 8 种功能区域 × 5 种重复序列/CpG 注释，沿用参考图的顺序、配色及单元格计数。
黑色横线仅分隔 13 个 CLL 样本与 4 个 pB 样本。

主图与参考图一致，统计 Ω<0 且下尾置换 p_down<0.05 的 triplet 数量；每项有 1000 次置换。
后缀 _FDR005 的图改为 Ω<0 且 q_down_global_BH<0.05；其 FDR 沿用全部 23800 项检验，主图没有把名义 P 值写成 FDR。
_pooled_average_heatmap 分别展示 R（n=8）、S（n=5）、pB（n=4）三组均值，各组内部按样本等权平均。
R 和 S 分开汇总，不合并为一个 CLL 均值；这里仅做描述性绘图，不进行 R/S 组间显著性比较。
合并图 plots/ABC_pooled_average_heatmap 已更新为 Ω<0 且全局 BH FDR<0.05，
覆盖全部 23800 项检验，在按符号、类别、组别筛选前完成校正。
其绘图数值由原始检验表重新汇总，ABC_pooled_average_counts.tsv 的 mean_count 为 FDR 后计数均值，
mean_count_nominal 保留名义 P 值均值供对照。该图单独运行 ../scripts/plot_abc_pooled_combined.py 重画；
该脚本独立复核 BH 值及样本级计数，保留原有 18 cm 宽度、正方形格子和配色。
合并图行标签仅显示 R、S、pB，不显示 n；格内数字保持 3.5 pt，其余文字统一为 8 pt。
三个色条共用一个文字标签，画布为 18 × 8.035 cm。

ABC 分类沿用参考定义，令 T={H3K9me3,H3K27me3,H3K36me3}：
- A：含 DNA 甲基化，且至少包含 T 中一个标记。
- B：不含 DNA 甲基化，且至少包含 T 中一个标记。
- C：不含 T 中任何标记。

本数据只有 7 个信号，ABC 的组合上限分别为 12、19、4，总计 35；
参考图使用 8 个信号，上限为 15、31、10，不能直接套用。
各类原始计数受组合上限影响，counts_by_sample_region_class.tsv 同时提供各类比例。

pB 对应 CEMT0129、CEMT0134、CEMT0135、CEMT0138；本地原始 metadata 注释为 B cells、disease=None，来源为 Tonsil。
sample_labels.tsv 保留完整样本映射。所有单项检验、分类和计数均已导出，可从现有数据重画，无需重新提交计算任务。
复现脚本：../scripts/plot_abc_triplet_classes.py。
""")
    print("Class capacities:", capacities)
    print("Sample rows:", mapping[["display_label", "pooled_group"]].to_string(index=False))
    print("Saved plots to", PLOTS)


if __name__ == "__main__":
    main()
