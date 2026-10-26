#!/usr/bin/env python3
"""Project four held-out female CLL-R donors onto the fixed male DE signatures.

No female CLL-S donors exist. This is expression-pattern concordance, not
an independent female R/S differential-expression or diagnostic validation.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/bcells-rna-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bcells-rna-cache')
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / 'results_rna_de_slurm'
MALE = PROJECT / 'results_rna_de_male_ighv'
OUTPUT = PROJECT / 'results_rna_de_female_validation'
PRIMARY = 'logcounts_method_IGHV_adjusted'


def correlation(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def projection(expression, metadata, train, covariates):
    """Fit male-only log-expression model, then freeze nuisance effects/scaling."""
    is_r = metadata['group'].eq('CLL-R').astype(float)
    design = pd.DataFrame({'intercept': 1.}, index=metadata.index)
    for covariate in covariates:
        design[covariate] = metadata[covariate].astype(float)
    design['R_vs_S'] = is_r
    x = design.loc[train].to_numpy(float)
    assert np.linalg.matrix_rank(x) == x.shape[1]
    beta = np.linalg.lstsq(x, expression.loc[train].to_numpy(float), rcond=None)[0]
    nuisance = np.zeros_like(expression.to_numpy(float))
    for j, covariate in enumerate(covariates, 1):
        nuisance += design[covariate].to_numpy(float)[:, None] * beta[j]
    adjusted = pd.DataFrame(expression.to_numpy(float)-nuisance, index=expression.index, columns=expression.columns)
    center = adjusted.loc[train].mean()
    scale = adjusted.loc[train].std(ddof=1).clip(lower=1e-6)
    z = (adjusted-center)/scale
    fitted_beta = pd.DataFrame(beta.T, index=expression.columns, columns=design.columns)
    return adjusted, z, fitted_beta, dict(n_train=len(train), n_parameters=x.shape[1],
                                        residual_degrees_of_freedom=len(train)-x.shape[1],
                                        condition_number=float(np.linalg.cond(x)), covariates=covariates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    OUTPUT.mkdir(exist_ok=True)
    meta = pd.read_csv(SOURCE/'sample_metadata.tsv', sep='\t', index_col=0)
    male = meta.index[meta.sex.eq('Male')].tolist()
    female = meta.index[meta.sex.eq('Female')].tolist()
    assert len(male) == 9 and len(female) == 4
    assert meta.loc[female, 'group'].eq('CLL-R').all()
    meta['unstranded'] = ~meta.analyzed_as_stranded.astype(str).str.lower().eq('true')
    meta['IGHV_M'] = meta.ighv_class.eq('IGHV-M')
    raw = pd.read_csv(SOURCE/'length_scaled_counts.tsv', sep='\t', index_col=0).loc[:, meta.index]
    counts = np.rint(raw).astype(np.int64)
    # Freeze DESeq-style ratio normalization on the nine male training donors.
    eligible = counts.loc[:, male].ge(10).sum(axis=1).ge(4)
    reference_genes = eligible & counts.loc[:, male].gt(0).all(axis=1)
    log_reference = np.log(counts.loc[reference_genes, male]).mean(axis=1)
    factors = {}
    for sample in meta.index:
        sample_counts = counts.loc[reference_genes, sample]
        positive = sample_counts.gt(0)
        factors[sample] = float(np.exp(np.median(np.log(sample_counts.loc[positive])-log_reference.loc[positive])))
    factors = pd.Series(factors)
    logcounts = np.log2(counts.div(factors, axis=1)+1).T
    tpm = pd.read_csv(SOURCE/'gene_tpm_tximport.tsv', sep='\t', index_col=0).loc[:, meta.index]
    logtpm = np.log2(tpm+1).T
    lists = {label: pd.read_csv(MALE/f'rna_method_ighv_adjusted_FDR_{alpha}.tsv', sep='\t', index_col=0)
             for label, alpha in [('FDR05', '0.05'), ('FDR10', '0.10')]}
    assert len(lists['FDR05']) == 38 and len(lists['FDR10']) == 60
    lists['FDR05_autosomal'] = lists['FDR05'].loc[~lists['FDR05'].chrom.isin(['chrX', 'chrY', 'X', 'Y'])]
    union = lists['FDR10'].index
    assert union.isin(logcounts.columns).all()
    stranded_male = [s for s in male if not meta.loc[s, 'unstranded']]
    specifications = {
        PRIMARY: (logcounts.loc[:, union], male, ['unstranded', 'IGHV_M']),
        'logcounts_unadjusted': (logcounts.loc[:, union], male, []),
        'logTPM_method_IGHV_adjusted': (logtpm.loc[:, union], male, ['unstranded', 'IGHV_M']),
        'stranded_logcounts_IGHV_adjusted': (logcounts.loc[:, union], stranded_male, ['IGHV_M']),
    }
    sample_results, gene_results, all_scores, model_info = [], [], [], {}
    primary_z = None
    for method, (expression, train, covariates) in specifications.items():
        adjusted, z, beta, info = projection(expression, meta, train, covariates)
        model_info[method] = info
        beta.to_csv(OUTPUT/(method+'_projection_coefficients.tsv'), sep='\t', index_label='gene_id_version')
        train_r = [s for s in train if meta.loc[s, 'group'] == 'CLL-R']
        train_s = [s for s in train if meta.loc[s, 'group'] == 'CLL-S']
        for label, genes in lists.items():
            ids = genes.index
            values = z.loc[:, ids]
            centroid_r = values.loc[train_r].mean().to_numpy()
            centroid_s = values.loc[train_s].mean().to_numpy()
            delta = centroid_r-centroid_s
            direction = np.sign(genes.log2FoldChange.to_numpy())
            for sample in meta.index:
                v = values.loc[sample].to_numpy()
                dist_r = float(np.sqrt(np.mean((v-centroid_r)**2)))
                dist_s = float(np.sqrt(np.mean((v-centroid_s)**2)))
                position = float(np.dot(v-centroid_s, delta)/np.dot(delta, delta))
                result = dict(sample=sample, sex=meta.loc[sample, 'sex'], group=meta.loc[sample, 'group'],
                              ighv_class=meta.loc[sample, 'ighv_class'], method=method, gene_set=label,
                              n_genes=len(ids), distance_R=dist_r, distance_S=dist_s,
                              nearer_male_group='CLL-R' if dist_r < dist_s else 'CLL-S',
                              projection_S0_R1=position, correlation_R=correlation(v, centroid_r),
                              correlation_S=correlation(v, centroid_s),
                              signed_signature_score=float(np.mean(v*direction)),
                              n_genes_in_male_R_direction_vs_S=int(np.sum((v-centroid_s)*direction > 0)),
                              n_genes_closer_to_R=int(np.sum(np.abs(v-centroid_r) < np.abs(v-centroid_s))))
                all_scores.append(result)
                if sample in female:
                    sample_results.append(result)
            if method == PRIMARY and label == 'FDR05':
                primary_z = values.copy()
        if method == PRIMARY:
            genes = lists['FDR10']
            female_average = adjusted.loc[female].mean()
            male_r = adjusted.loc[train_r].mean()
            male_s = adjusted.loc[train_s].mean()
            for gene, entry in genes.iterrows():
                direction = np.sign(entry.log2FoldChange)
                record = dict(gene_id_version=gene, gene_name=entry.gene_name, chrom=entry.chrom,
                              male_DESeq2_log2FoldChange=entry.log2FoldChange, male_DESeq2_padj=entry.padj,
                              in_male_FDR05=gene in lists['FDR05'].index,
                              male_log_expression_R_minus_S=male_r[gene]-male_s[gene],
                              female_mean_minus_male_S=female_average[gene]-male_s[gene],
                              female_mean_minus_male_R=female_average[gene]-male_r[gene],
                              female_mean_same_direction=bool((female_average[gene]-male_s[gene])*direction > 0),
                              male_log_model_same_direction_as_DESeq2=bool(beta.loc[gene, 'R_vs_S']*direction > 0))
                n_direction = 0
                for sample in female:
                    difference = adjusted.loc[sample, gene]-male_s[gene]
                    record[sample+'_minus_male_S'] = difference
                    record[sample+'_same_direction'] = bool(difference*direction > 0)
                    record[sample+'_TPM'] = tpm.loc[gene, sample]
                    n_direction += difference*direction > 0
                record['n_females_same_direction'] = int(n_direction)
                gene_results.append(record)
    scores = pd.DataFrame(sample_results)
    scores.to_csv(OUTPUT/'female_signature_concordance.tsv', sep='\t', index=False)
    pd.DataFrame(all_scores).to_csv(OUTPUT/'all_sample_signature_scores.tsv', sep='\t', index=False)
    gene_table = pd.DataFrame(gene_results)
    gene_table.to_csv(OUTPUT/'gene_direction_validation.tsv', sep='\t', index=False)
    primary_z.to_csv(OUTPUT/'primary_adjusted_male_scaled_expression.tsv', sep='\t', index_label='sample')
    metadata_columns = ['group', 'sex', 'age', 'ighv_class', 'analyzed_as_stranded']
    selected_meta = meta[metadata_columns].copy()
    selected_meta['male_frozen_size_factor'] = factors
    selected_meta.to_csv(OUTPUT/'sample_metadata.tsv', sep='\t', index_label='sample')
    summary = []
    for (method, label), group in scores.groupby(['method', 'gene_set'], sort=False):
        summary.append(dict(method=method, gene_set=label, n_genes=int(group.n_genes.iloc[0]),
                            n_female_R=len(group), n_nearer_male_R=int(group.nearer_male_group.eq('CLL-R').sum()),
                            mean_projection_S0_R1=float(group.projection_S0_R1.mean()),
                            min_projection_S0_R1=float(group.projection_S0_R1.min()),
                            max_projection_S0_R1=float(group.projection_S0_R1.max())))
    pd.DataFrame(summary).to_csv(OUTPUT/'validation_summary.tsv', sep='\t', index=False)
    config = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
                  females_R=female, females_S=[], males_R=meta.index[meta.sex.eq('Male') & meta.group.eq('CLL-R')].tolist(),
                  males_S=meta.index[meta.sex.eq('Male') & meta.group.eq('CLL-S')].tolist(),
                  primary_genes='Fixed male FDR<0.05 set, 38 genes; never selected using female data',
                  normalization='RSEM lengthScaledTPM counts rounded; DESeq-style ratio factors against a frozen male-only geometric-mean reference; log2(normalized counts+1)',
                  n_normalization_reference_genes=int(reference_genes.sum()),
                  projection='OLS log-expression fit using males only, removing RNA-method and IGHV terms; z scaling learned in males only. This projection is distinct from the original DESeq2 test.',
                  models=model_info, interpretation='Distance/projection are descriptive concordance, not validated diagnostic classification; female donors have no S controls.',
                  no_female_DE_pvalues=True, input_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'sample_metadata.tsv', SOURCE/'length_scaled_counts.tsv', SOURCE/'gene_tpm_tximport.tsv', MALE/'rna_method_ighv_adjusted_FDR_0.05.tsv', MALE/'rna_method_ighv_adjusted_FDR_0.10.tsv']})
    (OUTPUT/'validation_config.json').write_text(json.dumps(config, indent=2)+'\n')
    print(scores.loc[(scores.method == PRIMARY) & (scores.gene_set == 'FDR05')].to_string(index=False), flush=True)
    print(pd.DataFrame(summary).to_string(index=False), flush=True)
    for label, ids in [('FDR05', lists['FDR05'].index), ('FDR10', lists['FDR10'].index)]:
        selected = gene_table.set_index('gene_id_version').loc[ids]
        print(label, 'female mean same direction', int(selected.female_mean_same_direction.sum()), '/', len(selected),
              '; all four females same direction', int(selected.n_females_same_direction.eq(4).sum()), flush=True)
    if args.no_plots:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    male_s = meta.index[meta.sex.eq('Male') & meta.group.eq('CLL-S')].tolist()
    male_r = meta.index[meta.sex.eq('Male') & meta.group.eq('CLL-R')].tolist()
    order = male_s+male_r+female
    genes = lists['FDR05'].sort_values('log2FoldChange', ascending=False)
    matrix = primary_z.loc[order, genes.index].T.to_numpy()
    labels = [str(row.gene_name) if pd.notna(row.gene_name) and str(row.gene_name) else str(row.gene_id)
              for _, row in genes.iterrows()]
    fig, axes = plt.subplots(1, 2, figsize=(12, 11), gridspec_kw={'width_ratios': [2.2, 1]})
    ax = axes[0]
    im = ax.imshow(matrix, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(labels)), labels=labels, fontsize=7)
    ax.set_xticks(range(len(order)), labels=[s+'\n'+('Female R' if s in female else 'Male '+meta.loc[s, 'group'][-1]) for s in order], rotation=60, ha='right', fontsize=7)
    ax.axvline(len(male_s)-.5, color='black', lw=1.5)
    ax.axvline(len(male_s)+len(male_r)-.5, color='black', lw=1.5)
    ax.set_title('38 male DE genes: male-trained method/IGHV-adjusted expression', fontsize=10)
    fig.colorbar(im, ax=ax, fraction=.025, pad=.02, label='Gene z-score (male reference)')
    ax = axes[1]
    group = scores[(scores.method == PRIMARY) & (scores.gene_set == 'FDR05')].set_index('sample').loc[female]
    ax.axvspan(-.5, .5, color='#4575b4', alpha=.08)
    ax.axvline(0, color='#4575b4', linestyle='--', label='Male S centroid')
    ax.axvline(1, color='#d73027', linestyle='--', label='Male R centroid')
    ax.axvline(.5, color='gray', linestyle=':', linewidth=1)
    ax.scatter(group.projection_S0_R1, range(len(female)), color='#8e44ad', s=55, zorder=3)
    for i, (sample, row) in enumerate(group.iterrows()):
        ax.annotate(f'{row.projection_S0_R1:.2f}', (row.projection_S0_R1, i), xytext=(6, 6), textcoords='offset points', fontsize=9)
    ax.set_yticks(range(len(female)), labels=female, fontsize=9)
    ax.set_ylim(-.6, len(female)-.3)
    ax.set_xlim(min(-.3, group.projection_S0_R1.min()-.3), max(1.4, group.projection_S0_R1.max()+.4))
    ax.set_xlabel('Projection onto male S → R expression contrast')
    ax.set_title('Held-out female donors (all CLL-R)', fontsize=10)
    ax.legend(loc='upper left', fontsize=8, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.suptitle('Female concordance with the fixed male RNA signature\nNo female CLL-S controls: descriptive expression-pattern validation', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, .95])
    for ext in ('pdf', 'png'):
        path = OUTPUT/('female_male_signature_concordance.'+ext)
        fig.savefig(path, dpi=250, bbox_inches='tight')
        print('Saved', path, flush=True)
    plt.close(fig)


if __name__ == '__main__':
    main()
