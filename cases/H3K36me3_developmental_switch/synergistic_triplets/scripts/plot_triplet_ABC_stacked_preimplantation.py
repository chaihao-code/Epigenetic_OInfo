#!/usr/bin/env python3
"""Eight-signal A/B/C count heatmaps for seven preimplantation samples."""

# Resolve repository-relative paths in retained operational JSON manifests.
_RELEASE_ROOT = next(p for p in __import__('pathlib').Path(__file__).resolve().parents
                     if (p / '.epigenetic_oinfo_root').is_file())
def _release_manifest(value):
    if isinstance(value, dict):
        return {k: _release_manifest(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_manifest(v) for v in value]
    if isinstance(value, str) and value.startswith(('cases/', 'common/', 'data/')):
        return str(_RELEASE_ROOT / value)
    return value
import csv
import gzip
import hashlib
import itertools
import json
from collections import Counter, defaultdict

from plot_triplet_ABC_stacked_post_E11p5 import (
    ROOT, MARKS, REPRESSIVE, BLOCKS, classify, write_tsv, read_tsv, draw_heatmaps, np,
)

STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']
LABELS = ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE']
PREFIX = 'triplet_ABC_stacked_preimplantation'


def read_gzip(path):
    with gzip.open(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    manifest = _release_manifest(json.loads((ROOT / 'inputs/source_manifest.json').read_text()))
    original_path = ROOT.parent / 'batch_effect_figures/data/all_16variant_triplets_bh56.tsv.gz'
    snapshot_path = ROOT / 'inputs/preimplantation_additional_stages_bh56.tsv.gz'
    assert hashlib.sha256(snapshot_path.read_bytes()).hexdigest() == manifest['snapshot_sha256']
    rows = [r for r in read_gzip(original_path) if int(r['mask']) == 0]
    assert {r['stage'] for r in rows} == {'2cell', '8cell', 'ICM'}
    for r in rows:
        r['q'] = float(r['q_down_block_bh56'])
    original_stages = {r['stage'] for r in rows}
    for r in read_gzip(snapshot_path):
        if r['stage'] in STAGES and r['stage'] not in original_stages:
            r['q'] = float(r['q_down_sample_block_bh56'])
            rows.append(r)
    assert len(rows) == 7 * 40 * 56
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    maxima = dict(Counter(classify(t) for t in universe))
    assert maxima == {'A': 15, 'B': 31, 'C': 10}
    assert classify(['DNA_methylation', 'H3K36me3', 'ATAC']) == 'A'
    assert classify(['H3K36me3', 'H3K4me1', 'ATAC']) == 'B'
    families, sample_stages = defaultdict(list), defaultdict(set)
    for r in rows:
        assert r['stage'] in STAGES and r['block'] in BLOCKS
        assert int(r['n_perm']) == 10000 and frozenset(r['triplet'].split('/')) in universe
        if 'triplet_class' in r:
            assert r['triplet_class'] == classify(r['triplet'].split('/'))
        families[r['sample'], r['block']].append(r)
        sample_stages[r['stage']].add(r['sample'])
    assert len(families) == 7 * 40 and set(sample_stages) == set(STAGES)
    assert all(len(v) == 1 for v in sample_stages.values())
    sample_by_stage = {s: next(iter(sample_stages[s])) for s in STAGES}
    per_sample, summary, lookup, participation = [], [], {}, Counter()
    for stage in STAGES:
        sample = sample_by_stage[stage]
        for block in BLOCKS:
            family = families[sample, block]
            assert len(family) == 56 and {frozenset(r['triplet'].split('/')) for r in family} == universe
            ordered = sorted(family, key=lambda r: float(r['p_down']))
            q = np.minimum(1, np.minimum.accumulate((np.array([float(r['p_down']) for r in ordered])
                                                    * 56 / np.arange(1, 57))[::-1])[::-1])
            assert np.allclose(q, [r['q'] for r in ordered], rtol=0, atol=1e-12)
            significant = [r for r in family if float(r['observed_hoi']) < 0 and r['q'] <= .05]
            counts = Counter(classify(r['triplet'].split('/')) for r in significant)
            assert sum(counts.values()) == len(significant)
            for r in significant:
                for mark in r['triplet'].split('/'):
                    participation[sample, mark] += 1
            for cls in 'ABC':
                count = counts[cls]
                assert 0 <= count <= maxima[cls]
                per_sample.append(dict(sample=sample, stage=stage, block=block, triplet_class=cls,
                                       count=count, theoretical_max=maxima[cls]))
                summary.append(dict(triplet_class=cls, stage=stage, block=block, mean_count=count,
                                    theoretical_max=maxima[cls], n=1))
                lookup[cls, stage, block] = count
    # Check all eight individual mark counts against the existing early analysis.
    occurrences = read_tsv(ROOT / 'data/original_development_mark_participation_by_sample.tsv')
    checked = 0
    for r in occurrences:
        if r['stage'] in STAGES:
            assert r['sample'] == sample_by_stage[r['stage']]
            assert participation[r['sample'], r['mark']] == int(r['n_significant'])
            checked += 1
    assert checked == 7 * 8
    write_tsv(ROOT / 'data' / (PREFIX + '_by_sample.tsv'), per_sample)
    write_tsv(ROOT / 'data' / (PREFIX + '_summary.tsv'), summary)
    titles, figure_metadata = draw_heatmaps(
        lookup, stages=STAGES, stage_labels=[s + ' (n=1)' for s in LABELS], prefix=PREFIX,
        period_title='preimplantation', aggregation_title='one original sample per stage', colorbar_label='Count')
    sources = [dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
               for p in [original_path, snapshot_path]]
    validation = dict(status='complete', sources=sources, signals=MARKS,
                      repressive_marks=sorted(REPRESSIVE), theoretical_max=maxima,
                      class_definitions=titles, stages=STAGES, sample_by_stage=sample_by_stage,
                      n_samples=7, samples_per_stage=1, ICM_TE_separate=True,
                      n_blocks=40, n_permutations=10000, BH_family_size=56,
                      criterion='observed_hoi < 0 and BH56 q <= 0.05',
                      aggregation='One original sample per stage; raw counts, no tissue averaging',
                      original_sources='Mask 0 for 2cell/8cell/ICM; preserved BH56 snapshot for zygote/4cell/morula/TE',
                      color_scale='Independent Reds scale for each class and developmental period',
                      n_summary_rows=len(summary), n_sample_rows=len(per_sample),
                      eight_signal_families_complete=True, BH56_recomputed_and_matched=True,
                      ABC_partition_checked=True, H3K36me3_classification_checked=True,
                      snapshot_sha256_matches_manifest=True, all_56_mark_occurrence_counts_match=True,
                      **figure_metadata)
    (ROOT / 'data' / (PREFIX + '_validation.json')).write_text(json.dumps(validation, indent=2) + '\n')
    print('Saved', ROOT / 'plots' / (PREFIX + '.png'))
    print('Saved', ROOT / 'plots' / (PREFIX + '.pdf'))
    print('Validated 7 stages × 40 blocks × 56 triplets; all 56 mark-occurrence counts match.')


if __name__ == '__main__':
    main()
