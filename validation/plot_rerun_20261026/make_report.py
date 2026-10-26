#!/usr/bin/env python3
"""Merge execution records and verify every newly generated figure file."""
import csv
import hashlib
import json
from pathlib import Path
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = [('bcells', 'bcells/report.json'), ('aging', 'aging/report.json'),
           ('switch', 'switch/results.json'), ('mouse', 'mouse/report.json'),
           ('gsea', 'gsea/report.json')]
GROUPS = {
    'bcells': 'B细胞 / CLL', 'thyroid': '甲状腺去卷积与年龄',
    'neutrophil': '中性粒细胞衰老', 'switch': 'H3K36me3发育转换',
    'mouse_preimp': '小鼠植入前MI', 'mouse_development': '小鼠整体发育',
    'gsea': '甲状腺局部O-info GSEA',
}
rows, artifact_map = [], {}
for group, name in SOURCES:
    path = HERE / name
    data = json.loads(path.read_text())
    entries = data['entries'] if isinstance(data, dict) else data
    for original in entries:
        r = dict(original)
        script = r['script']
        if group == 'aging':
            category = 'neutrophil' if '/neutrophil/' in script else 'thyroid'
        elif group == 'mouse':
            category = 'mouse_preimp' if '/mouse_preimp/' in script else 'mouse_development'
        else:
            category = group
        outputs = r.get('outputs', [])
        figures = []
        for output in outputs:
            relative = output['path'] if isinstance(output, dict) else output
            p = Path(relative)
            if p.is_absolute():
                p = p.relative_to(ROOT)
            if p.suffix.lower() not in {'.png', '.pdf', '.svg'}:
                continue
            target = ROOT / p
            assert target.resolve().is_relative_to(ROOT.resolve()), target
            assert target.is_file(), target
            relative = str(p)
            figures.append(relative)
            if relative in artifact_map:
                continue
            binary = target.read_bytes()
            if p.suffix == '.png':
                with Image.open(target) as image:
                    dimensions = list(image.size)
                    image.verify()
            else:
                dimensions = None
                if p.suffix == '.pdf':
                    assert binary.startswith(b'%PDF-') and b'%%EOF' in binary[-2048:], target
            artifact_map[relative] = dict(path=relative, bytes=len(binary),
                sha256=hashlib.sha256(binary).hexdigest(), dimensions=dimensions,
                group=category, script=script)
        figures = sorted(set(figures))
        status = r.get('status', 'success' if r.get('returncode') == 0 and figures else 'failed')
        if status in {'ok', 'passed', 'successful'}:
            status = 'success'
        error = r.get('error', r.get('reason', r.get('note', '') if status!='success' else ''))
        if status!='success' and not error:
            log_path = r.get('log')
            if log_path and (ROOT/log_path).exists():
                messages = (ROOT/log_path).read_text(errors='replace').strip().splitlines()
                error = messages[-1] if messages else ''
        rows.append(dict(group=category, group_zh=GROUPS[category], script=script,
            description_zh=r.get('description_zh', r.get('title', r.get('description', r.get('id', Path(script).stem)))),
            status=status, figure_count=len({str(Path(x).with_suffix('')) for x in figures}),
            file_count=len(figures), outputs=figures,
            input_nature_zh=r.get('input_nature_zh',r.get('input_type',r.get('input_nature',''))),
            command=r.get('command',r.get('commands',[])),
            error=error,
            detail_report=str(Path(name).parent / 'REPORT.md')))

# Three discretization settings are executions of the same MI entrypoint.
# Count independent source entrypoints and combine all their variant outputs.
merged = {}
for r in rows:
    if r['script'] not in merged:
        r['commands'] = [r['command']]
        merged[r['script']] = r
    else:
        target = merged[r['script']]
        target['outputs'] = sorted(set(target['outputs'] + r['outputs']))
        target['commands'].append(r['command'])
        if r['status']!='success':
            target['status']=r['status']; target['error']=r['error']
for r in merged.values():
    r['figure_count']=len({str(Path(x).with_suffix('')) for x in r['outputs']})
    r['file_count']=len(r['outputs'])
    if len(r['commands'])>1:
        r['command']='\n'.join(' '.join(map(str,x)) if isinstance(x,list) else x for x in r['commands'])
    if r['script'].endswith('/plot_by_block_all_anchors.py'):
        r['description_zh']='早胚MI轨迹：3/5/10状态，每种40区域×9个锚点，共1080幅'
    elif r['script'].endswith('/plot_mi_ratio_hist.py'):
        r['description_zh']='3/5/10状态的成对/联合归一化MI密度直方图，共3幅'
rows = list(merged.values())
inventory = json.loads((HERE/'source_inventory.json').read_text())
expected = set(inventory['independent_entrypoints'])
actual = {r['script'] for r in rows}
assert expected == actual, dict(missing=sorted(expected-actual), unexpected=sorted(actual-expected))

groups = []
for key, label in GROUPS.items():
    historical_entries = [r for r in rows if r['group']==key]
    entries = [r for r in historical_entries if r['status']!='removed_by_request']
    outputs = {x for r in entries for x in r['outputs']}
    groups.append(dict(group=key, title=label, entries=len(entries),
        removed_entries=len(historical_entries)-len(entries),
        success=sum(r['status']=='success' for r in entries),
        figure_count=len({str(Path(x).with_suffix('')) for x in outputs}), file_count=len(outputs)))
summary = dict(entries=len(rows), successful_entries=sum(r['status']=='success' for r in rows),
    failed_entries=sum(r['status'] not in ('success','removed_by_request') for r in rows),
    removed_entries=sum(r['status']=='removed_by_request' for r in rows),
    figure_count=len({str(Path(x).with_suffix('')) for x in artifact_map}), file_count=len(artifact_map))
(HERE/'report.json').write_text(json.dumps(dict(summary=summary,groups=groups,entries=rows),ensure_ascii=False,indent=2))
(HERE/'artifacts.json').write_text(json.dumps(list(artifact_map.values()),ensure_ascii=False,indent=2))
with (HERE/'plot_inventory.tsv').open('w') as f:
    fields=['group_zh','script','description_zh','status','figure_count','file_count','input_nature_zh','command','error']
    writer=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');writer.writeheader()
    for r in rows:
        row=dict(r)
        if isinstance(row['command'],list):row['command']=' '.join(map(str,row['command']))
        writer.writerow(row)
with (HERE/'figure_files.tsv').open('w') as f:
    writer=csv.DictWriter(f,['path','bytes','sha256','group','script'],delimiter='\t',extrasaction='ignore')
    writer.writeheader();writer.writerows(artifact_map.values())

lines=['# 全仓作图重跑检查表','',
    f"历史检查覆盖 **{summary['entries']}个独立绘图入口**；当前保留成功 **{summary['successful_entries']}个**，未成功 **{summary['failed_entries']}个**，按用户要求移除 **{summary['removed_entries']}个**。当前保留 **{summary['figure_count']}幅图、{summary['file_count']}个图文件**（同名PNG/PDF按一幅计）。",'',
    '甲状腺去卷积仅保留figQ1_composition.png。其他去卷积图已按用户要求删除，默认重跑会跳过这些入口；甲状腺GSEA和其他项目保留。','',
    '所有新输出都在当前release目录下。下面按案例汇总，随后逐项列出实际内容、数量和输出链接。没有把复制来的历史图片计入本次成果。','',
    '| 案例 | 成功 / 尝试入口 | 新画幅数 | 图文件数 | 详细记录 |',
    '|---|---:|---:|---:|---|']
for g in groups:
    subgroup = 'aging' if g['group'] in {'thyroid','neutrophil'} else ('mouse' if g['group'].startswith('mouse_') else g['group'])
    lines.append(f"| {g['title']} | {g['success']}/{g['entries']} | {g['figure_count']} | {g['file_count']} | [{subgroup}]({subgroup}/REPORT.md) |")
lines += ['','## 每个程序画了什么','',
          '| 案例 / 程序 | 内容 | 状态 | 幅数 / 文件数 | 本次图文件 |', '|---|---|---|---:|---|']
for r in rows:
    folders = sorted({str(Path(p).parent) for p in r['outputs']})
    links = '<br>'.join(f'[{Path(d).name}]({ROOT / d})' for d in folders)
    lines.append(f"| {r['group_zh']}<br>[{Path(r['script']).name}]({ROOT / r['script']}) | {r['description_zh']} | {r['status']} | {r['figure_count']} / {r['file_count']} | {links or '—'} |")
lines += ['','## 未成功项目','']
for r in rows:
    if r['status'] not in ('success','removed_by_request'):
        lines.append(f"- `{r['script']}`：{r['error'] or '缺少历史输入，见分组日志。'} [{r['group_zh']}详情]({r['detail_report']})。")
lines += ['','## 是否都只读中间文件','',
    '**不能把所有原入口都视为纯作图程序。** 本次尽量使用已有中间结果调用原绘图段，以下例外已核实：','',
    '| 程序 | 原入口还会做什么 / 读取什么 | 本次实际方式 |', '|---|---|---|',
    '| 甲状腺04_select_markers / 05_validate_deconvolution | 读取参考beta/npz筛marker；重新模拟混合与NNLS | 从reference_matrix和simulation_results中间表调用原绘图代码 |',
    '| 甲状腺06_deconvolve_samples | 读取甲基化bigWig，重新去卷积/bootstrap | 从thyroid_cell_proportions表绘图 |',
    '| 甲状腺07_oinfo_age_adjustment / 08_region_filtering | 重新回归/bootstrap；从离散轨迹重新计算O-information | 从保留的调整/过滤结果表绘图，未重算O-info |',
    '| 甲状腺09_marker_go_analysis / 10_histone_association | 读取RNA或组蛋白bigWig、GTF；10缺缓存时读取bigWig | 从marker/enrichment/histone结果表绘图 |',
    '| 甲状腺local回归02_gsea_pathway / 03_gsea_marginal | 从回归、GTF或提取的离散信号构建排名与基因集，运行GSEA并由库内部出图 | 从已有排名/GMT/report重建富集曲线，不重跑置换，见gsea报告 |',
    '| B细胞plot_igh_tracks_corrected | 读取prepared beta矩阵和GENCODE GTF参考注释，计算局部组均值 | 实际执行原入口；输入仍有现有外部软链接 |',
    '| B细胞validate_male_de_in_females | 从导入的counts/TPM重做归一化及男性模型拟合 | 实际执行原入口，不是单纯读取DE汇总 |',
    '| 小鼠resize_redundancy_core | 读取既有PDF矢量路径重排版面 | 实际重排源PDF，不读数值中间表 |','',
    '其余成功入口主要读取结果CSV/TSV及元数据。有些会重新聚合、验证BH、拟合图中趋势线或做Welch检验，这些在各组报告标明；均未重跑O-information置换分析。历史1000次置换的中间表仍按历史结果标注。','',
    '## 输入恢复、核验与重跑','',
    '原release并非所有图都能直接一次跑通：有遗漏的中间表、旧路径和绘图错误。本次恢复的中间输入及SHA256记录在各组restored_inputs.json；没有从其他结果替代缺失的历史分析。具体修复见各组报告。', '',
    '所有PNG逐文件解码验证，PDF检查文件头和EOF；完整路径、大小和SHA256见[figure_files.tsv](figure_files.tsv)。各组另保留日志、版面断言、数值核对和抽样缩略图检查。', '',
    '[逐程序TSV](plot_inventory.tsv) · [逐文件TSV](figure_files.tsv) · [机器清单](report.json) · [入口覆盖审计](source_inventory.json)', '',
    '使用本机已有依赖环境，从仓库根目录复跑：', '',
    '```bash', 'python validation/plot_rerun_20261026/run_all.py', '```','',
    '缺少的历史输入仍会在对应子日志中明确失败；其他任务继续运行。甲状腺混合入口通过本次runner使用绘图函数/AST绘图段，直接运行其原main仍可能重做分析。']
(HERE/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(summary,ensure_ascii=False),flush=True)
