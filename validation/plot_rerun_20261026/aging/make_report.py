#!/usr/bin/env python3
"""Verify generated aging figures and produce the review inventory."""
import hashlib
import json
from pathlib import Path
from PIL import Image, ImageDraw, ImageOps

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
DESCRIPTIONS = {
 '04_select_markers': ('参考细胞甲基化marker热图：Panel A/B', '原入口从参考beta/npz重选marker；本次从reference_matrix_A/B中间表调用heatmap()'),
 '05_validate_deconvolution': ('模拟混合比例与估计免疫比例验证图（4面板）', '原入口重新模拟/NNLS；本次读取simulation_results.tsv执行原绘图段'),
 '06_deconvolve_samples': ('6供体细胞组成与模型敏感性图（4面板）', '原入口读取甲基化bigWig并重新去卷积/bootstrap；本次读取细胞比例表调用figures()'),
 '06_immune_age_association': ('免疫细胞比例随年龄变化及回归（4面板）', '中间细胞比例表+样本年龄；本次执行轻量回归及绘图'),
 '07_oinfo_age_adjustment': ('免疫调整前后年龄斜率、置信区间及稳健性（6面板）', '原入口从O-info及比例表重新回归/bootstrap；本次读取oinfo_age_adjustment.tsv及散点源表执行绘图段'),
 '08_region_filtering': ('排除免疫特异区域与matched downsampling敏感性（3面板）', '原入口读离散轨迹重新计算O-information；本次读取filtering_values/region_filtering中间表执行绘图段'),
 '09_marker_go_analysis': ('免疫/甲状腺marker表达、年龄关联和免疫区间重叠（3面板）', '原入口读RNA bigWig、GTF、离散轨迹/注释并做GO分析；本次只读immune_marker_analysis/immune_region_enrichment表执行绘图段'),
 '10_histone_association': ('6种组蛋白信号的年龄斜率；信号与免疫比例关系（各6面板）', '原入口优先bin_means_cache缓存，缺缓存读取组蛋白bigWig；本次只读histone_immune_association及细胞比例表执行绘图段'),
 '12_fig_q1_composition': ('细胞组成及Panel/Model稳健性（3面板）', '中间细胞比例表+样本元数据；纯绘图'),
 '15_fig_volcano': ('免疫调整后的年龄斜率/P值火山图', '中间oinfo_age_fdr_table.tsv；纯绘图'),
 '16_fig_lollipop': ('年龄斜率调整前后相对差的棒棒糖图', '中间oinfo_age_fdr_table.tsv；纯绘图'),
 'plot_oinfo_vs_age': ('40区域及CpG区域O-info与年龄散点：全部/仅成人，共4图', '中间O-information CSV+供体元数据；拟合图中趋势线'),
 'plot_paper_figures': ('O-info供体×区域热图；年龄斜率最大的6区域散点', '中间O-information CSV+供体元数据；重新选择斜率最大的区域'),
 'plot_volcano_oinfo_age': ('成人中性粒细胞O-info年龄斜率/P值火山图', '中间oinfo_age_regression_adults.tsv；纯绘图'),
}

rows = []
previews = []
for stem, (description, nature) in DESCRIPTIONS.items():
    p = BASE / f'{stem}.json'
    if not p.exists():
        raise RuntimeError(f'Missing execution record: {p}')
    row = json.loads(p.read_text())
    row.update(description_zh=description, input_nature_zh=nature)
    row['source_sha256'] = hashlib.sha256((ROOT / row['script']).read_bytes()).hexdigest()
    row['figure_count'] = len({str(Path(o['path']).with_suffix('')) for o in row['outputs']})
    row['file_count'] = len(row['outputs'])
    for o in row['outputs']:
        path = ROOT / o['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == o['sha256'], path
        if path.suffix == '.png':
            with Image.open(path) as im:
                im.verify()
            previews.append((path, stem))
        elif path.suffix == '.pdf':
            data = path.read_bytes()
            assert data.startswith(b'%PDF-') and b'%%EOF' in data[-2048:], path
    rows.append(row)

report = dict(entries=rows, successful_entries=sum(r['status']=='success' for r in rows),
              removed_entries=sum(r['status']=='removed_by_request' for r in rows),
              figure_count=sum(r['figure_count'] for r in rows), file_count=sum(r['file_count'] for r in rows))
(BASE / 'report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False))
(BASE / 'results.json').write_text(json.dumps(rows, indent=2, ensure_ascii=False))

lines = ['# 衰老案例作图重跑', '',
         f"14个历史入口，当前保留成功{report['successful_entries']}个，按用户要求移除{report['removed_entries']}个；当前{report['figure_count']}幅、{report['file_count']}个图文件。", '',
         '甲状腺去卷积仅保留figQ1_composition.png；其他去卷积图片已按要求删除，默认重跑不会再生成。甲状腺GSEA和其他项目图片保留。', '',
         '本次所有数值输入均为release内已有中间表和元数据，未从原工作目录恢复数据；没有重算O-information、模拟混合或重新去卷积。混合入口通过run_aging_plots.py调用原绘图函数或执行原main中的绘图AST段，详细命令与实际读取文件哈希见report.json。', '',
         '| 程序 | 内容 | 幅/文件 | 输入核实与本次方式 | 输出 |', '|---|---|---:|---|---|']
for r in rows:
    folder = (ROOT / r['outputs'][0]['path']).parent if r['outputs'] else BASE / 'logs'
    lines.append(f"| [{Path(r['script']).name}]({ROOT / r['script']}) | {r['description_zh']} | {r['figure_count']}/{r['file_count']} | {r['input_nature_zh']} | [{r['status']}]({folder}) |")
lines += ['', '## 运行与核验', '',
          '```bash', 'python validation/plot_rerun_20261026/aging/run_aging_plots.py', 'python validation/plot_rerun_20261026/aging/make_report.py', '```', '',
          'PNG逐张解码，PDF检查文件头和EOF，全部图记录SHA256。原程序05的模拟验证图原来固定标注sigma=0.03，但实际选sigma=0.05数据；已改为按实际数据标注。Panel A/B heatmap保留中间表marker原始顺序。Arial字体已作为真实文件保存到data/fonts/arial.ttf。', '',
          '原入口中06_deconvolve_samples、09_marker_go_analysis会读取原始/参考bigWig等文件，10_histone_association缺缓存时也会读取bigWig；08_region_filtering从离散轨迹重算O-info。不能把这些原入口当成单纯作图命令。', '',
          '保留的图片位于各案例结果目录，路径见上表；日志在logs/。[机器清单](report.json)，[缩略图](contact_sheet.jpg)。']
(BASE / 'REPORT.md').write_text('\n'.join(lines) + '\n')

if previews:
    cols, width, height = 4, 420, 330
    sheet = Image.new('RGB', (cols*width, ((len(previews)+cols-1)//cols)*height), 'white')
    draw = ImageDraw.Draw(sheet)
    for i,(path,label) in enumerate(previews):
        im = Image.open(path).convert('RGB')
        im.thumbnail((width-12,height-45))
        x, y = (i%cols)*width, (i//cols)*height
        sheet.paste(im,(x+(width-im.width)//2,y+28))
        draw.text((x+5,y+5),label,fill='black')
    sheet.save(BASE/'contact_sheet.jpg',quality=85)
print(json.dumps({k:v for k,v in report.items() if k!='entries'},ensure_ascii=False))
