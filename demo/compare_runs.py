#!/usr/bin/env python3
"""Compare both complete thyroid demo routes, without recomputing any information measure."""
from pathlib import Path
import argparse,gzip,hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def digest(path, compressed=False):
    value=hashlib.sha256();opener=gzip.open if compressed else open
    with opener(path,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()
def compare():
    results=ROOT/'demo/results';a=results/'from_bins';b=results/'from_bigwig'
    ra=json.loads((a/'validation.json').read_text());rb=json.loads((b/'validation.json').read_text())
    assert ra['status']==rb['status']=='passed'
    assert ra['n_permutations']==rb['n_permutations']==10000
    assert ra['checks']['n_classes']==rb['checks']['n_classes']==40
    stem='CEEHRC.CEMT0098'
    checks={}
    for suffix in ['_oinfo.csv','_entropy_all_subsets.csv','_permutation_oinfo.csv','_permutation_summary.csv',
                   '_sampled_bins.tsv.gz','_discrete_states.tsv.gz','_local_oinfo.tsv.gz','_local_global_check.tsv']:
        pa=a/'oinfo'/(stem+suffix);pb=b/'oinfo'/(stem+suffix)
        checks[suffix]=digest(pa,suffix.endswith('.gz'))==digest(pb,suffix.endswith('.gz'))
    if not all(checks.values()):raise ValueError(f'Thyroid route comparison differs: {checks}')
    maximum=max(ra['checks']['maximum_local_global_difference'],rb['checks']['maximum_local_global_difference'])
    report=dict(status='passed',sample=stem,full_genome_windows=3209513,autosomal_windows=2875012,
        classes=40,signals=7,n_permutations=10000,slurm_jobs={'from_bins':ra['slurm_job_id'],'from_bigwig':rb['slurm_job_id']},
        n_sampled_loci=ra['checks']['n_sampled_loci'],maximum_local_global_difference=maximum,
        identical_decoded_result_files=checks,
        signal_input_comparisons=rb['input_comparisons'],
        historical_values_reproduced=False,
        historical_max_difference=ra['checks'].get('maximum_difference_from_original_oinfo'),
        historical_note='The upstream QC reports already document unavailable paper-era histone inputs. This validates available complete inputs and does not claim reproduction of archived paper-era values.')
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'demo/results/verification.json');a=p.parse_args()
    report=compare();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
