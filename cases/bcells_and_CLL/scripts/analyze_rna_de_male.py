#!/usr/bin/env python3
"""Final male-only R/S RNA comparison adjusting RNA method and IGHV."""

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
import json
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests

import os
PROJECT=Path(__file__).resolve().parents[1]
DATA_ROOT=Path(os.environ.get('BCELLS_DATA_ROOT', _RELEASE_ROOT / 'data'))
SOURCE=PROJECT/'results_rna_de_slurm'
OUT=PROJECT/'results_rna_de_male_ighv'
GENES='EZH2 EZH1 EED SUZ12 KDM6A KDM6B HDAC1 HDAC2 EP300 CREBBP KMT2A SETD1A KDM5A KMT2C KMT2D KDM1A SUV39H1 SUV39H2 SETDB1 SETD2 PAX5 EBF1 TCF4 NFATC1 BRD4'.split()


def main():
    OUT.mkdir(exist_ok=True)
    meta=pd.read_csv(SOURCE/'sample_metadata.tsv',sep='\t',index_col=0)
    meta=meta[meta.sex.eq('Male')].copy()
    assert meta.group.value_counts().to_dict()=={'CLL-S':5,'CLL-R':4}
    meta.to_csv(OUT/'sample_metadata.tsv',sep='\t')
    counts=pd.read_csv(SOURCE/'length_scaled_counts.tsv',sep='\t',index_col=0).T.loc[meta.index]
    counts=np.rint(counts).astype(np.int64)
    # Refilter within the male subset, using its smaller group size (4).
    counts=counts.loc[:,counts.ge(10).sum(axis=0).ge(4)]
    tpm=pd.read_csv(SOURCE/'gene_tpm_tximport.tsv',sep='\t',index_col=0).loc[:,meta.index]
    ann=pd.read_csv(SOURCE/'gene_annotation.tsv',sep='\t',index_col='gene_id')
    is_r=meta.group.eq('CLL-R')
    base=pd.DataFrame({'intercept':1.,'unstranded':(~meta.analyzed_as_stranded).astype(float),'IGHV_M':meta.ighv_class.eq('IGHV-M').astype(float),'R_vs_S':is_r.astype(float)},index=meta.index)
    models={'rna_method_ighv_adjusted':['intercept','unstranded','IGHV_M','R_vs_S']}
    summaries=[]
    for name,cols in models.items():
        design=base[cols]
        assert np.linalg.matrix_rank(design)==len(cols)
        design.to_csv(OUT/(name+'_design.tsv'),sep='\t')
        print('Fitting',name,'samples',len(meta),'genes',counts.shape[1],flush=True)
        dds=DeseqDataSet(counts=counts,metadata=meta,design=design,n_cpus=4,refit_cooks=True)
        dds.deseq2()
        contrast=np.zeros(len(cols));contrast[-1]=1
        stats=DeseqStats(dds,contrast=contrast,alpha=.1,n_cpus=4)
        stats.summary()
        res=stats.results_df.copy();res.index.name='gene_id_version'
        res['gene_id']=res.index.str.split('.').str[0]
        res=res.join(ann,on='gene_id')
        res['mean_TPM_R']=tpm.loc[res.index,meta.index[is_r]].mean(axis=1)
        res['mean_TPM_S']=tpm.loc[res.index,meta.index[~is_r]].mean(axis=1)
        res.sort_values('padj').to_csv(OUT/(name+'_all_genes.tsv'),sep='\t')
        for alpha in [.05,.1]:
            hit=res[res.padj.lt(alpha)]
            hit.sort_values('padj').to_csv(OUT/(name+f'_FDR_{alpha:.2f}.tsv'),sep='\t')
            summaries.append(dict(model=name,FDR=alpha,n_R=4,n_S=5,n_prefiltered_genes=counts.shape[1],n_tested_p=int(res.pvalue.notna().sum()),n_tested_padj=int(res.padj.notna().sum()),n_significant=len(hit),n_absLFC_ge1=int(hit.log2FoldChange.abs().ge(1).sum()),n_higher_R=int(hit.log2FoldChange.gt(0).sum()),n_higher_S=int(hit.log2FoldChange.lt(0).sum())))
        candidate=pd.DataFrame({'gene_name':GENES}).merge(res[res.gene_name.isin(GENES)][['gene_name','gene_id','baseMean','log2FoldChange','pvalue','padj','mean_TPM_R','mean_TPM_S']],on='gene_name',how='left',validate='one_to_one')
        candidate['q_candidate25']=multipletests(candidate.pvalue.fillna(1),method='fdr_bh')[1]
        candidate['tested']=candidate.pvalue.notna()
        candidate.to_csv(OUT/(name+'_candidate25.tsv'),sep='\t',index=False)
        pd.DataFrame(summaries).to_csv(OUT/'model_summary.tsv',sep='\t',index=False)
    config=dict(n_R=4,n_S=5,samples_R=meta.index[is_r].tolist(),samples_S=meta.index[~is_r].tolist(),models=models,IGHV_correction=True,age_correction=False,counts='Reuse verified tximport lengthScaledTPM counts with common gene-length reference from all 13 donors; subset to 9 males and refit all parameters',prefilter='>=10 counts in >=4 male donors',contrast='R/S; positive log2FoldChange=higher in R',FDR='genome-wide BH with PyDESeq2 independent filtering optimized at alpha=0.1; report both 0.05 and 0.1',limits=['4 R vs 5 S donors; small sample size','All 4 R are stranded; only 2 S are stranded','No verified wet-lab batch correction', 'IGHV imbalance: R has 1 M / 3 UM; S has 4 M / 1 UM', 'Few replicates within covariate combinations; Cook filtering does not cover every donor'])
    (OUT/'analysis_config.json').write_text(json.dumps(config,indent=2))
    print(pd.DataFrame(summaries).to_string(index=False),flush=True)


if __name__=='__main__': main()
