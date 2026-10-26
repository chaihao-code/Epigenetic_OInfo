#!/usr/bin/env python3
"""Validate beta scale, coordinates and coverage without changing the input."""
import csv
import hashlib
import json
import math
from pathlib import Path

_ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.epigenetic_oinfo_root').is_file())
ROOT = _ROOT / 'cases/bcells_and_CLL/dmr_repro_cemt'
source = ROOT / 'input/all_bins_1kb_cemt_g2_g3.csv'
canonical = {'chr'+str(i) for i in range(1, 23)}
total = autosomes = kept = invalid = nonfinite = partial = 0
minimum, maximum = math.inf, -math.inf
coverage, group_coverage, chromosomes = {}, {}, {}
with source.open() as handle:
    reader = csv.DictReader(handle)
    samples = reader.fieldnames[3:]
    assert len(samples) == 13 and len(set(samples)) == 13
    for row in reader:
        total += 1
        if row['chrom'] not in canonical:
            continue
        autosomes += 1
        chromosomes[row['chrom']] = chromosomes.get(row['chrom'], 0)+1
        width = int(row['end'])-int(row['start'])
        assert 0 < width <= 1000
        partial += width != 1000
        values = [float(row[s]) if row[s].strip() not in ('', 'NA', 'NaN', 'nan') else None for s in samples]
        for value in values:
            if value is not None:
                if not math.isfinite(value):
                    nonfinite += 1
                elif not 0 <= value <= 1:
                    invalid += 1
                else:
                    minimum, maximum = min(minimum, value), max(maximum, value)
        n = sum(v is not None for v in values)
        coverage[n] = coverage.get(n, 0)+1
        if n >= 8:
            kept += 1
            nr, ns = sum(v is not None for v in values[:8]), sum(v is not None for v in values[8:])
            group = str(nr)+'R_'+str(ns)+'S'
            group_coverage[group] = group_coverage.get(group, 0)+1
assert invalid == 0 and nonfinite == 0
report = dict(n_rows_input=total, n_autosomal_bins=autosomes, n_bins_at_least_8_samples=kept,
              n_invalid_beta=invalid, n_nonfinite_beta=nonfinite, beta_min=minimum, beta_max=maximum,
              n_short_terminal_autosomal_bins=partial, samples=samples, coverage_counts=coverage,
              retained_group_coverage=group_coverage, bins_by_chromosome=chromosomes,
              input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              script_sha256=hashlib.sha256((ROOT/'scripts/dmr_limma_1kb_cemt_g2_g3.R').read_bytes()).hexdigest())
(ROOT/'input_validation.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2), flush=True)
