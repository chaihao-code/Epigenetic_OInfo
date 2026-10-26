import shutil,hashlib,json
from pathlib import Path
root=Path.cwd(); case=root/'cases/H3K36me3_developmental_switch'; upstream=Path('/tmp/Epigenetic_OInfo_upstream_20261026/results/H3K36me3_developmental_switch'); original=Path('/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch')
paths=['batch_effect_figures/data/all_16variant_triplets_bh56.tsv.gz','batch_effect_figures/data/late_reference_triplets_bh56.tsv.gz','synergistic_triplets/inputs/source_manifest.json','synergistic_triplets/inputs/preimplantation_additional_stages_bh56.tsv.gz','batch_effect_figures/fdr_lt_0p1/triplet_support_across_variants.tsv','batch_effect_figures/fdr_lt_0p1/threshold_comparison.tsv']
rows=[]
for rel in paths:
 src=upstream/rel
 if not src.exists():src=original/rel
 dst=case/rel;dst.parent.mkdir(exist_ok=True,parents=True)
 if not dst.exists():shutil.copy2(src,dst)
 rows.append(dict(source=str(src),destination=str(dst.relative_to(root)),sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),bytes=dst.stat().st_size))
for d in ['batch_effect_figures/plots','redundancy_core/plots','redundancy_core/full_signal_occurrence/plots','synergistic_triplets/plots']:(case/d).mkdir(exist_ok=True,parents=True)
(root/'validation/plot_rerun_20261026/switch/restored_inputs.json').write_text(json.dumps(rows,indent=2))
print('Restored',len(rows),'intermediates')
