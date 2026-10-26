#!/usr/bin/env python3
"""Index retained result images by biological case without duplicating plots."""
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/'validation/plot_rerun_20261026'
CASES={
 'bcells_and_CLL':('B cells and CLL','bcells/report.json'),
 'mouse_preimp':('Mouse preimplantation','mouse/report.json'),
 'mouse_development':('Mouse development','mouse/report.json'),
 'mouse_after_E11.5':('Mouse development after E11.5',None),
 'H3K36me3_developmental_switch':('H3K36me3 developmental switch','switch/results.json'),
 'aging/neutrophil':('Neutrophil aging','aging/report.json'),
 'aging/thyroid':('Thyroid aging','aging/report.json'),
}
FIELDS=['category','figure','region','library','pdf','png','svg']


def read_entries(name):
    data=json.loads((AUDIT/name).read_text())
    return data['entries'] if isinstance(data,dict) else data


def write_table(path,rows,fields):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def main():
    overview=[];all_rows=[]
    for case,(title,manifest) in CASES.items():
        base=ROOT/'cases'/case
        entries=read_entries(manifest) if manifest else []
        if case=='aging/thyroid':entries+=read_entries('gsea/report.json')
        pairs={}
        for entry in entries:
            if not entry['script'].startswith(f'cases/{case}/') or entry.get('status')!='success':continue
            for item in entry.get('outputs',[]):
                value=item['path'] if isinstance(item,dict) else item
                path=Path(value);path=path if path.is_absolute() else ROOT/path
                if path.suffix not in ('.png','.pdf','.svg'):continue
                assert path.is_file(),path
                key=str(path.with_suffix(''))
                row=pairs.setdefault(key,dict(category=Path(entry['script']).stem,figure=path.stem,
                                              region='',library='',pdf='',png='',svg=''))
                row[path.suffix[1:]]=str(path.relative_to(ROOT))
                if 'gsea_pathway' in path.parts or 'gsea_marginal' in path.parts:
                    marker='gsea_pathway' if 'gsea_pathway' in path.parts else 'gsea_marginal'
                    index=path.parts.index(marker)
                    row['region']=path.parts[index+1]
                    row['library']=path.parts[index+2] if marker=='gsea_pathway' else 'signal_change_gene_sets'
        rows=sorted(pairs.values(),key=lambda r:(r['category'],r['figure'],r['pdf'],r['png']))
        write_table(base/'figure_index.tsv',rows,FIELDS)
        all_rows.extend(dict(case=case,**r) for r in rows)
        grouped=defaultdict(list)
        for row in rows:grouped[row['category']].append(row)
        lines=[f'# {title}: figures','',f'{len(rows)} retained result figures. PNG and PDF files with the same stem count as one figure.','',
               'The [file index](figure_index.tsv) lists the exact paths of every retained PDF, PNG or SVG.','']
        if rows:
            lines+=['| Figure group | Figures | Location |','|---|---:|---|']
            for category,members in grouped.items():
                dirs=set()
                for row in members:
                    for ext in ('pdf','png','svg'):
                        if not row[ext]:continue
                        path=ROOT/row[ext]
                        folder=path.parent
                        for marker in ('gsea_pathway','gsea_marginal'):
                            if marker in path.parts:
                                folder=Path(*path.parts[:path.parts.index(marker)+1])
                        dirs.add(folder)
                links='<br>'.join(f'[{str(folder.relative_to(ROOT))}]({os.path.relpath(folder,base)}/)' for folder in sorted(dirs))
                lines.append(f'| `{category}` | {len(members)} | {links} |')
        if case=='aging/thyroid':
            lines+=['','Only [figQ1_composition.png](deconvolution/results/figures/figQ1_composition.png) is retained for deconvolution, by user request. Other deconvolution images, including its PDF version, were removed. Pathway and marginal GSEA curves are retained.']
        if case=='mouse_after_E11.5':
            lines+=['','There are no independent result images in this case. Shared early/late and H3K36me3 comparisons are indexed under [H3K36me3 developmental switch](../H3K36me3_developmental_switch/FIGURES.md). Numerical results in this case are retained.']
        if case=='mouse_development':
            lines+=['','Original-size TC/DTC figures remain in the execution-audit directory and are linked above. Their 18 cm versions are in this case. The source PDF used for vector reconstruction is an input, so it is not counted as a newly rendered result.']
        (base/'FIGURES.md').write_text('\n'.join(lines)+'\n')
        readme=base/'README.md';text=readme.read_text()
        if 'FIGURES.md' not in text:
            anchor='## Outputs\n';addition='\nGenerated result images: [figure locations and file index](FIGURES.md).\n'
            assert anchor in text,readme
            readme.write_text(text.replace(anchor,anchor+addition,1))
        count_files=sum(bool(r[k]) for r in rows for k in ('pdf','png','svg'))
        overview.append(dict(case=case,title=title,figures=len(rows),files=count_files))
        print(case,len(rows),'figures',count_files,'files',flush=True)
    write_table(ROOT/'cases/figure_index.tsv',all_rows,['case',*FIELDS])
    lines=['# Figure locations by case','',
           'Each case has a FIGURES.md navigation page and a figure_index.tsv with exact image paths. Counts reflect currently retained result figures.','',
           '| Case | Figures | Image files | Locations and index |','|---|---:|---:|---|']
    for r in overview:
        lines.append(f"| {r['title']} | {r['figures']} | {r['files']} | [{r['case']}]({r['case']}/FIGURES.md) |")
    lines+=['','Only the thyroid deconvolution images were reduced to figQ1_composition.png. All thyroid GSEA curves and other cases were retained. No narrative reports/ directories remain.','',
            'Full image paths: [figure_index.tsv](figure_index.tsv). Rebuild navigation with `python tools/index_case_figures.py`.']
    (ROOT/'cases/FIGURES.md').write_text('\n'.join(lines)+'\n')
    readme=ROOT/'cases/README.md';text=readme.read_text()
    if 'FIGURES.md' not in text:
        readme.write_text(text.replace('## Outputs\n','## Outputs\n\nGenerated result images: [case-by-case locations](FIGURES.md).\n',1))


if __name__=='__main__':main()
