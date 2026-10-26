# H3K36me3作图复跑记录

入口16个，成功14个，缺输入2个。新画24张图（每张PNG和PDF）。

全部入口均读取中间结果与元数据；未读取bigWig、FASTQ或原始信号，也未重跑置换或冗余核心分析。部分原入口会从已有统计全表重汇总并核对BH值。

复跑：`python validation/plot_rerun_20261026/switch/run_all.py`。日志和逐文件路径见`results.json`。

|入口|画了什么|状态|张数|输出目录|
|---|---|---|---:|---|
|plot_batch_effect.py|16种信号来源替换组合的同阶段标记参与比例|success|1|cases/H3K36me3_developmental_switch/batch_effect_figures/plots|
|plot_support_without_h3k27me3_replacement.py|去除H3K27me3替换后的8组合支持度分布（两个FDR阈值）|success|2|cases/H3K36me3_developmental_switch/batch_effect_figures/without_h3k27me3_replacement|
|plot_triplet_support_distribution.py|16组合的triplet支持度分布（两个FDR阈值）|success|2|cases/H3K36me3_developmental_switch/batch_effect_figures/plots|
|plot_combined_trajectories.py|triplet与冗余核心标记轨迹并排合图|success|1|cases/H3K36me3_developmental_switch/redundancy_core/plots|
|plot_core_occurrence_frequency_e11.5.py|历史发育后期冗余核心Top25频率、累计线和点阵|blocked|0|FileNotFoundError: [Errno 2] No such file or directory: '/lustre/home/2101110354/Epigenetic_OInfo_release/cases/mouse_after_E11.5/supreme_cores_batch_v4_replacement_bg_mad_z2.5_bic_k1k2_eta05/all_supreme_cores_filtered.csv'|
|plot_core_occurrence_frequency_preimp.py|历史植入前冗余核心Top25频率、累计线和点阵|blocked|0|FileNotFoundError: [Errno 2] No such file or directory: '/lustre/home/2101110354/Epigenetic_OInfo_release/cases/mouse_preimp/supreme_cores_batch_v4_replacement_bg_mad_z2.5_bic_k1k2_eta05/all_supreme_cores_filtered.csv'|
|plot_full_signal_combination_occurrence.py|全信号triplet/冗余核心Top25频率（4单图和2×2合图）|success|5|cases/H3K36me3_developmental_switch/redundancy_core/full_signal_occurrence/plots|
|plot_redundancy_core_trajectories.py|8标记在冗余核心中的发育轨迹|success|1|cases/H3K36me3_developmental_switch/redundancy_core/plots|
|plot_synergistic_triplet_trajectories.py|8标记在协同triplet中的发育轨迹|success|1|cases/H3K36me3_developmental_switch/redundancy_core/plots|
|plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2.py|含H3K9ac/H3K4me2全信号A/B/C组合频次点阵（仅符号筛选）|success|3|cases/H3K36me3_developmental_switch/redundancy_core/plots/legacy_triplet_combinations|
|plot_mark_occurrence_per_tissue_by_stage.py|标记参与次数：植入前、发育后期及合图|success|3|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|
|plot_original_development_heatmap.py|原始8信号14阶段标记参与比例热图|success|1|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|
|plot_triplet_ABC_selected_regions_post_E11p5.py|8信号选定区域A/B/C相对出现次数轨迹|success|1|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|
|plot_triplet_ABC_stacked_post_E11p5.py|8信号发育后期A/B/C三层热图|success|1|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|
|plot_triplet_ABC_stacked_preimplantation.py|8信号植入前A/B/C三层热图|success|1|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|
|plot_triplet_A_B_C_post_E11p5.py|8信号4区域A/B/C相对出现次数轨迹|success|1|cases/H3K36me3_developmental_switch/synergistic_triplets/plots|

修复记录：

- 从provenance指定的本机upstream备份和原工作目录恢复6个缺失的统计中间文件。SHA256、源路径和大小见`restored_inputs.json`。这些文件现已复制到release中，复跑不依赖备份目录。
- 为原脚本输出创建plots目录；trajectory和full-signal入口增加`--plots-only`，从已有汇总表重画。
- 旧triplet组合脚本默认D没有分类定义，原本退出0但不出图。本次保留ABC分类规则、显式跑A/B/C并读取现有符号筛选计数；默认修为合法A，非法类会报错。
- 植入前堆叠热图的Zygote标签超出画布；调整同一绘图函数左右布局后通过原有文字边界和单元格尺寸断言。后期共享同一布局函数。
- 两个historical core Top25入口缺原始的`all_supreme_cores_filtered.csv`；未用其他核心结果替代。
- `trajectory_utils.py`是无main的共享绘图模块，由两个trajectory入口调用，不单独计入口。

输入表可能记录历史统计设置；本次仅重画，不将旧结果标记为重新完成10000次置换。
