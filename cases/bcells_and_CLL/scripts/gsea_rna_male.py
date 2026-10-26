#!/usr/bin/env python3
"""Whole-gene preranked GO GSEA using male RNA-seq Wald statistics."""

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
from pathlib import Path
from collections import defaultdict
import gzip
import json
import numpy as np
import pandas as pd
import gseapy as gp
from functools import lru_cache

import os
PROJECT=Path(__file__).resolve().parents[1]
DATA_ROOT=Path(os.environ.get('BCELLS_DATA_ROOT', _RELEASE_ROOT / 'data'))
OUT=PROJECT/'results_gsea_rna_male'
ANN=DATA_ROOT/'annotations/go_human'


def ontology():
    terms={};alt={};record=None;headers=[]
    def finish():
        if record and 'id' in record and not record.get('obsolete',False):
            terms[record['id']]=record.copy()
            for x in record['alt']:alt[x]=record['id']
    with (ANN/'go-basic.obo').open() as h:
        for line in h:
            line=line.strip()
            if line.startswith('['):
                finish();record={'parents':set(),'alt':[]} if line=='[Term]' else None
            elif record is None:
                if line.startswith(('data-version:','date:','format-version:')):headers.append(line)
            elif line.startswith('id: '):record['id']=line[4:]
            elif line.startswith('name: '):record['name']=line[6:]
            elif line.startswith('namespace: '):record['namespace']=line[11:]
            elif line.startswith('alt_id: '):record['alt'].append(line[8:])
            elif line.startswith('is_a: '):record['parents'].add(line.split()[1])
            elif line.startswith('relationship: part_of '):record['parents'].add(line.split()[2])
            elif line=='is_obsolete: true':record['obsolete']=True
    finish()
    @lru_cache(None)
    def ancestors(go):
        go=alt.get(go,go)
        if go not in terms:return frozenset()
        out={go}
        for parent in terms[go]['parents']:out.update(ancestors(parent))
        return frozenset(out)
    return terms,alt,ancestors,headers


def gene_sets(background):
    terms,alt,ancestors,obo_headers=ontology()
    by_symbol=defaultdict(set)
    for row in background.itertuples():
        if isinstance(row.gene_name,str):by_symbol[row.gene_name].add(row.gene_id)
    term_genes=defaultdict(set)
    with gzip.open(ANN/'goa_human.gaf.gz','rt') as h:
        for line in h:
            if line.startswith('!'):continue
            x=line.rstrip('\n').split('\t')
            if len(x)<15 or 'taxon:9606' not in x[12].split('|') or 'NOT' in x[3].split('|') or x[6]=='ND':continue
            go=alt.get(x[4],x[4])
            if go not in terms:continue
            ids=by_symbol.get(x[2],set())
            if len(ids)!=1:ids=set().union(*(by_symbol.get(s,set()) for s in x[10].split('|')))
            if len(ids)!=1:continue
            gid=next(iter(ids))
            for term in ancestors(go):term_genes[term].add(gid)
    return terms,term_genes,obo_headers


def main():
    OUT.mkdir(exist_ok=True)
    inputs={'male_RNA_method_IGHV_adjusted':PROJECT/'results_rna_de_male_ighv/rna_method_ighv_adjusted_all_genes.tsv'}
    background=pd.read_csv(inputs['male_RNA_method_IGHV_adjusted'],sep='\t')
    terms,term_genes,obo_headers=gene_sets(background)
    names=dict(zip(background.gene_id,background.gene_name))
    summaries=[];model_config={}
    for label,path in inputs.items():
        table=pd.read_csv(path,sep='\t')
        # Keep every eligible finite Wald statistic, irrespective of significance.
        # Exclude NA p-values, including those flagged by Cook's filtering.
        rank=table[np.isfinite(table.stat)&np.isfinite(table.pvalue)][['gene_id','stat']].copy()
        assert not rank.gene_id.duplicated().any()
        rank=rank.sort_values(['stat','gene_id'],ascending=[False,True],kind='stable')
        ranked_ids=set(rank.gene_id)
        model_out=OUT/label;model_out.mkdir(exist_ok=True)
        rank.to_csv(model_out/'all_genes_Wald.rnk',sep='\t',index=False,header=False)
        model_config[label]=dict(input=str(path),n_ranked_genes=len(rank),n_exact_ties=int(rank.stat.duplicated().sum()),n_positive=int(rank.stat.gt(0).sum()),n_negative=int(rank.stat.lt(0).sum()))
        for namespace,short in [('biological_process','BP'),('molecular_function','MF'),('cellular_component','CC')]:
            collections={go+'|'+terms[go]['name']:sorted(gs&ranked_ids) for go,gs in term_genes.items() if terms[go]['namespace']==namespace and 15<=len(gs&ranked_ids)<=500}
            gmt=model_out/(short+'.gmt')
            if gmt.exists():
                saved={}
                with gmt.open() as h:
                    for line in h:
                        fields=line.rstrip('\n').split('\t')
                        genes=sorted(set(fields[2:]) & ranked_ids)
                        if 15 <= len(genes) <= 500:
                            saved[fields[0]]=genes
                collections=saved
            assert collections
            with gmt.open('w') as h:
                for name,genes in sorted(collections.items()):h.write(name+'\tGOA_human\t'+'\t'.join(genes)+'\n')
            print(label,short,'genes',len(rank),'gene sets',len(collections),flush=True)
            gsea=gp.prerank(rnk=rank,gene_sets=collections,outdir=str(model_out/short),min_size=15,max_size=500,permutation_num=1000,weight=1.,ascending=False,threads=4,seed=42,graph_num=0,no_plot=True,verbose=True)
            result=gsea.res2d.copy()
            for col in ['ES','NES','NOM p-val','FDR q-val','FWER p-val']:result[col]=pd.to_numeric(result[col],errors='coerce')
            result['GO_ID']=result.Term.str.split('|',n=1).str[0]
            result['GO_term']=result.Term.str.split('|',n=1).str[1]
            result['enriched_group']=np.where(result.NES.gt(0),'R','S')
            result['leading_edge_gene_names']=result.Lead_genes.fillna('').map(lambda s:';'.join(names.get(g,g) for g in s.split(';') if g))
            result.sort_values('FDR q-val').to_csv(model_out/(short+'_all_terms.tsv'),sep='\t',index=False)
            for alpha in [.05,.1,.25]:
                result[result['FDR q-val'].lt(alpha)].sort_values('FDR q-val').to_csv(model_out/(short+f'_FDR_{alpha:.2f}.tsv'),sep='\t',index=False)
            summaries.append(dict(model=label,ontology=short,n_ranked_genes=len(rank),n_sets_tested=len(result),n_FDR_005=int(result['FDR q-val'].lt(.05).sum()),n_FDR_010=int(result['FDR q-val'].lt(.1).sum()),n_FDR_025=int(result['FDR q-val'].lt(.25).sum()),n_R_FDR010=int((result['FDR q-val'].lt(.1)&result.NES.gt(0)).sum()),n_S_FDR010=int((result['FDR q-val'].lt(.1)&result.NES.lt(0)).sum())))
            pd.DataFrame(summaries).to_csv(OUT/'gsea_summary.tsv',sep='\t',index=False)
    config=dict(method='GSEApy prerank',gseapy_version=gp.__version__,ranking='All finite RNA DE Wald statistics with non-NA pvalue; no FDR/significance cutoff; no restriction to 21 or 120 candidates',positive_NES='R-enriched',negative_NES='S-enriched',models=model_config,gene_sets='GO BP/MF/CC; same official human GOA and go-basic annotation as previous ORA',annotation_propagation='is_a and part_of; excludes NOT and ND',identifier='stable Ensembl IDs; GOA symbols and unambiguous synonyms mapped via local GENCODE v47',min_size=15,max_size=500,permutation_num=1000,seed=42,weight=1,threads=4,FDR='GSEA empirical FDR q-value calculated within each model/ontology; report 0.05, 0.1 and exploratory 0.25',obo_headers=obo_headers,sources=json.loads((ANN/'sources.json').read_text()),limitations=['Preranked gene-set permutations, not donor/phenotype permutations','Only 4 R and 5 S male donors; original covariate imbalance persists','Enrichment is an association, not proof of pathway activation or causal regulation','Gene ranking can still be affected by individual donors and uncertain technical effects','GO terms overlap; parent/child terms are not independent biological findings'])
    (OUT/'analysis_config.json').write_text(json.dumps(config,indent=2))
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)


if __name__=='__main__':main()
