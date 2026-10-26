#!/usr/bin/env python3
"""Re-render aging figures from retained intermediates using original plot code.

Mixed analysis/plotting entrypoints are imported without running main. Their
original plotting functions or the AST statements in their figure section are
executed with the retained tables as inputs. No raw-track analysis is invoked.
"""
import argparse
import ast
import contextlib
import hashlib
import json
import os
from pathlib import Path
import runpy
import sys
import time
import traceback

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / '.epigenetic_oinfo_root').exists())
OUT = Path(__file__).resolve().parent
(OUT / 'cache').mkdir(exist_ok=True)
(OUT / 'cache/fonts.conf').write_text(
    '<?xml version="1.0"?><!DOCTYPE fontconfig SYSTEM "fonts.dtd">'
    f'<fontconfig><dir>{ROOT / "data/fonts"}</dir><dir>/usr/share/fonts/dejavu</dir>'
    f'<cachedir>{OUT / "cache/fontconfig"}</cachedir></fontconfig>')
os.environ.setdefault('MPLCONFIGDIR', str(OUT / 'cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME', str(OUT / 'cache'))
os.environ.setdefault('FONTCONFIG_FILE', str(OUT / 'cache/fonts.conf'))
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')

THYROID = ROOT / 'cases/aging/thyroid/deconvolution'
SCRIPTS = THYROID / 'scripts'
RETAINED_THYROID_SCRIPT = '12_fig_q1_composition.py'


def figure_tail(path, namespace, marker, end_marker=None):
    """Execute unchanged statements from the selected figure section of main."""
    source = path.read_text()
    marker_line = source[:source.index(marker)].count('\n') + 1
    tree = ast.parse(source, filename=str(path))
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    end_line = source[:source.index(end_marker)].count('\n') + 1 if end_marker else len(source.splitlines()) + 1
    statements = [n for n in main.body if marker_line < n.lineno < end_line and not isinstance(n, ast.Return)]
    if not statements:
        raise RuntimeError(f'No plotting statements found after {marker!r}')
    exec(compile(ast.Module(body=statements, type_ignores=[]), str(path), 'exec'), namespace)


def render(script):
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    plt.rcdefaults()  # Each original script starts with a fresh plotting style.
    font_manager.fontManager.addfont(str(ROOT / 'data/fonts/arial.ttf'))
    src = ROOT / script
    if '/thyroid/' in script:
        if src.name != RETAINED_THYROID_SCRIPT:
            raise ValueError('Only figQ1_composition.png is retained for thyroid deconvolution.')
        folder = THYROID / 'results/figures'
    else:
        folder = ROOT / 'cases/aging/neutrophil/plots'
    folder.mkdir(parents=True, exist_ok=True)
    inputs = []
    original_read_csv = pd.read_csv

    def read_csv(path, *args, **kwargs):
        if isinstance(path, (str, Path)):
            p = Path(path).absolute()
            if p.exists():
                record = dict(path=str(p), resolved=str(p.resolve()), bytes=p.stat().st_size,
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                if record not in inputs:
                    inputs.append(record)
        return original_read_csv(path, *args, **kwargs)

    pd.read_csv = read_csv
    try:
        if 'neutrophil' in script:
            source = src.read_text().replace('PLOTS = PROJECT / "plots"', f'PLOTS = Path({str(folder)!r})')
            ns = dict(__file__=str(src), __name__='__main__')
            exec(compile(source, str(src), 'exec'), ns)
        else:
            (THYROID / 'logs').mkdir(exist_ok=True)
            sys.path.insert(0, str(SCRIPTS))
            ns = runpy.run_path(str(src), run_name='_plot_rerun_')
            ns['FIG'] = folder
            tbl = ns['TBL']
            name = src.name
            if name.startswith('04_'):
                # reference_matrix_* are exactly the selected marker means.
                panels = []
                for panel in ('A', 'B'):
                    df = pd.read_csv(THYROID / f'data/processed/reference_matrix_{panel}.tsv', sep='\t')
                    wide = df.pivot(index='marker_id', columns='cell_type', values='methylation')
                    wide = wide.reindex(df.marker_id.drop_duplicates())
                    wide = wide.rename(columns=lambda c: f'mean_{c}').reset_index()
                    panels.append(wide)
                ns['heatmap'].__globals__.update(FIG=folder)
                ns['heatmap'](*panels)
            elif name.startswith('05_'):
                ns['sim'] = pd.read_csv(tbl / 'simulation_results.tsv', sep='\t')
                figure_tail(src, ns, '# ---- 图 ----', '# ---- identifiability report ----')
            elif name.startswith('06_deconvolve'):
                ns['figures'].__globals__.update(FIG=folder)
                ns['figures'](pd.read_csv(tbl / 'thyroid_cell_proportions.tsv', sep='\t'))
            elif name.startswith('07_'):
                ns.update(out=pd.read_csv(tbl / 'oinfo_age_adjustment.tsv', sep='\t'),
                          oinfo=pd.read_csv(THYROID / 'data/processed/oinfo_long.tsv', sep='\t'),
                          prop=pd.read_csv(tbl / 'thyroid_cell_proportions.tsv', sep='\t'),
                          meta=pd.read_csv(THYROID / 'data/sample_metadata/analysis_samples.tsv', sep='\t'))
                ns['main_df'] = ns['out'][(ns['out'].panel == 'A') & (ns['out'].model == 'A')]
                figure_tail(src, ns, '# ---------------- Figure 2 ----------------')
            elif name.startswith('08_'):
                ns.update(vals=pd.read_csv(tbl / 'oinfo_filtering_values.tsv', sep='\t'),
                          reg=pd.read_csv(tbl / 'oinfo_region_filtering.tsv', sep='\t'))
                figure_tail(src, ns, '# ---------------- Figure 3 ----------------')
            elif name.startswith('09_'):
                ns.update(df=pd.read_csv(tbl / 'immune_marker_analysis.tsv', sep='\t'),
                          ov=pd.read_csv(tbl / 'immune_region_enrichment.tsv', sep='\t'))
                figure_tail(src, ns, '# ---------------- Figure 4 ----------------')
            elif name.startswith('10_'):
                out = pd.read_csv(tbl / 'histone_immune_association.tsv', sep='\t')
                imm = pd.read_csv(tbl / 'thyroid_cell_proportions.tsv', sep='\t')
                imm_frac = (imm[(imm.panel == 'A') & (imm.model == 'A') &
                                imm.cell_type.isin(['T_CD4', 'T_CD8', 'B_cell', 'Monocyte', 'NK', 'Granulocyte'])]
                            .groupby('sample_id').fraction.sum())
                ns.update(main=out[out.metric == 'locus_relative'], imm_frac=imm_frac)
                figure_tail(src, ns, '# ---------------- Figure 5 ----------------')
            else:
                # These four entrypoints already read only retained tables.
                ns['main'].__globals__.update(FIG=folder)
                ns['main']()
    finally:
        pd.read_csv = original_read_csv
        plt.close('all')
    return inputs, folder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', help='Run one script basename')
    args = parser.parse_args()
    scripts = sorted(str(p.relative_to(ROOT)) for p in SCRIPTS.glob('*.py') if 'savefig' in p.read_text())
    scripts += sorted(str(p.relative_to(ROOT)) for p in (ROOT / 'cases/aging/neutrophil/scripts').glob('plot*.py'))
    rows = []
    for script in scripts:
        if args.only and Path(script).name != args.only:
            continue
        stem = Path(script).stem
        if '/thyroid/' in script and Path(script).name != RETAINED_THYROID_SCRIPT:
            record_path = OUT / f'{stem}.json'
            result = json.loads(record_path.read_text()) if record_path.exists() else dict(script=script)
            result.update(status='removed_by_request', outputs=[],
                          note='Only figQ1_composition.png is retained for thyroid deconvolution.')
            record_path.write_text(json.dumps(result, ensure_ascii=False, indent=2))
            rows.append(result)
            print(f'removed_by_request: {script}', flush=True)
            continue
        log = OUT / 'logs' / f'{stem}.log'
        started = time.time()
        result = dict(script=script, command=f'{sys.executable} {Path(__file__).relative_to(ROOT)} --only {Path(script).name}',
                      log=str(log.relative_to(ROOT)), started_unix=started)
        with log.open('w') as fh, contextlib.redirect_stdout(fh), contextlib.redirect_stderr(fh):
            try:
                inputs, folder = render(script)
                figures = [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
                                sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                           for p in sorted(folder.iterdir()) if p.suffix in ('.pdf', '.png', '.svg')]
                result.update(status='success' if figures else 'no_figures', inputs=inputs, outputs=figures)
            except Exception as exc:
                traceback.print_exc()
                result.update(status='failed', error=f'{type(exc).__name__}: {exc}', outputs=[])
        result['elapsed_seconds'] = round(time.time() - started, 2)
        rows.append(result)
        (OUT / f'{stem}.json').write_text(json.dumps(result, ensure_ascii=False, indent=2))
        print(f'{result["status"]}: {script} ({result["elapsed_seconds"]}s)', flush=True)
    (OUT / 'results.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
