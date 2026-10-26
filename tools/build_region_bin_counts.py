#!/usr/bin/env python3
"""Export English mm10 and hg38 region-bin counts, checked against fixed BEDs."""
import csv
import gzip
import json
from pathlib import Path

import pandas as pd
from openpyxl import Workbook
from build_oinfo_summary import excel_sheet

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'supplementary_data/final_tables'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    index=ROOT/'supplementary_data/bins_per_class_index.tsv'
    with index.open() as handle:
        entries=list(csv.DictReader(handle,delimiter='\t'))
    checks=[]
    for genome in ('mm10','hg38'):
        selected=[r for r in entries if r['genome']==genome and
                  f'/{genome}_global_annotation/' in r['file']]
        assert len(selected)==40 and len({r['class'] for r in selected})==40
        rows=[]
        for row in selected:
            path=ROOT/row['file']
            count=0;coordinates=set()
            with gzip.open(path,'rt') as handle:
                for line in handle:
                    if not line.strip() or line.startswith(('#','track','browser')):
                        continue
                    chrom,start,end,*_=line.rstrip().split('\t')
                    start,end=int(start),int(end)
                    assert 0<=start<end,(path,line)
                    coordinates.add((chrom,start,end));count+=1
            assert count==int(row['bins']),(path,count,row['bins'])
            assert len(coordinates)==count,('Duplicate coordinates in annotation class',path,count,len(coordinates))
            functional,repeat=row['class'].split('_AND_',1)
            rows.append(dict(Genome=genome,Genome_region=row['class'],N_bins=count,
                             Functional_region=functional,Repeat_context=repeat,Source_BED=row['file']))
            checks.append(dict(genome=genome,region=row['class'],index_count=int(row['bins']),
                               bed_rows=count,unique_coordinates=len(coordinates),source_bed=row['file']))
        table=pd.DataFrame(rows)
        table.to_csv(OUT/f'{genome}_region_bin_counts.tsv',sep='\t',index=False)
        guide=[
            ('Genome',genome),
            ('Scope','The 40 fixed global annotation classes; each combines a functional annotation and a CpG / repeat context.'),
            ('N_bins','Number of unique chromosome/start/end intervals in the corresponding supplied BED.gz file.'),
            ('Verification','Every count was checked against bins_per_class_index.tsv and the actual BED rows; no duplicate intervals occur within any class.'),
            ('Coordinates','BED coordinates are zero-based and half-open. The original genomic grid and class memberships are unchanged.'),
            ('Sampling','These are annotation-class membership counts, not the per-sample valid loci or loci selected for entropy / O-information estimation.'),
            ('Overlapping classes','Different annotation classes can overlap. Summing class counts does not give unique genomic coverage.'),
            ('Local thyroid analysis','Separate target-bin definitions under hg38_thyroid_local are outside this global annotation table.'),
            ('Region names','Original class identifiers are retained for joining to source data. Source-specific promoter / None aliases in analysis outputs may differ.'),
            ('Source index','supplementary_data/bins_per_class_index.tsv'),
            ('Reproduction','python tools/build_region_bin_counts.py'),
        ]
        book=Workbook(write_only=True)
        excel_sheet(book,'Region_bin_counts',table,{'Genome_region':66,'Functional_region':55,'Repeat_context':20,'Source_BED':110})
        excel_sheet(book,'Readme',pd.DataFrame(guide,columns=['Item','Description']),{'Item':27,'Description':110})
        book.save(OUT/f'{genome}_region_bin_counts.xlsx')
        print(genome,'40 regions; all BED counts verified',flush=True)
    (OUT/'bin_count_checks.json').write_text(json.dumps(dict(regions_checked=len(checks),entries=checks),indent=2)+'\n')
    (OUT/'README.md').write_text('''# Final tables

All workbook sheet names, headers, text cells and filenames are in English.

| File | Contents |
|---|---|
| [O_information_summary.xlsx](O_information_summary.xlsx) | Full-system O-information, permutation P values, BH q values, triplet significance, coverage and method notes. |
| [mm10_region_bin_counts.xlsx](mm10_region_bin_counts.xlsx) | Bin counts for the 40 mm10 global annotation classes. |
| [hg38_region_bin_counts.xlsx](hg38_region_bin_counts.xlsx) | Bin counts for the 40 hg38 global annotation classes. |

The two bin-count tables also have matching TSV exports. Bin counts were checked
against the source BED files and index, with duplicate intervals checked within
each class. These counts describe fixed annotation memberships, not sample-specific
valid / sampled loci. Class counts may overlap across regions.

Detailed O-information TSV exports and statistical definitions remain under
[oinfo_summary](../oinfo_summary/README.md). Blank P/q cells indicate missing or
unperformed tests. Full-system q values use BH across 40 regions per sample and
system; triplet q values retain their explicitly recorded original testing families.

Recreate the workbooks:

```bash
python tools/build_oinfo_summary.py --reuse-tables
python tools/build_region_bin_counts.py
```
''')


if __name__=='__main__':
    main()
