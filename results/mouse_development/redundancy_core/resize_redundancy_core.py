"""Reconstruct the plot from its PDF vector paths; no source data required."""
from pathlib import Path
import csv
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'redundancy_core_top25_cdf.pdf'
OUTPUT = HERE / 'redundancy_core_top25_cdf_7.5x4cm'
SVG = '{http://www.w3.org/2000/svg}'

def numbers(value):
    return list(map(float, re.findall(r'-?\d*\.?\d+(?:e[-+]?\d+)?', value)))

with tempfile.TemporaryDirectory() as temp:
    svg_path = Path(temp) / 'source.svg'
    subprocess.run(['pdftocairo', '-svg', str(SOURCE), str(svg_path)], check=True)
    root = ET.parse(svg_path).getroot()

surface = root.find(SVG + 'g')
bars = [p for p in surface if 'stroke-width:0.4;' in p.get('style', '')
        and 'Z' in p.get('d', '')]
assert len(bars) == 25
axis_bottom, axis_top = 103.867188, 18.351562
axis_height = axis_bottom - axis_top
counts = np.array([round((numbers(p.get('d'))[5] - numbers(p.get('d'))[1])
                         * 0.995957 / axis_height * 800) for p in bars])
x_centers = np.array([(numbers(p.get('d'))[0] + numbers(p.get('d'))[2])
                      / 2 * 0.995957 for p in bars])

# Black filled cubic paths below the histogram are the active matrix dots.
dots = []
for p in surface:
    if 'fill:rgb(0%,0%,0%)' in p.get('style', '') and 'C' in p.get('d', ''):
        xy = np.array(numbers(p.get('d'))).reshape(-1, 2)
        center = (xy.min(axis=0) + xy.max(axis=0)) / 2
        if center[1] > axis_bottom:
            dots.append(center)
row_centers = sorted(set(round(y, 2) for x, y in dots))
assert len(row_centers) == 7
matrix = np.zeros((7, 25), dtype=bool)
for x, y in dots:
    col = np.argmin(abs(x_centers - x))
    row = np.argmin(abs(np.array(row_centers) - y))
    assert abs(x_centers[col] - x) < 0.02
    matrix[row, col] = True

# Preserve the displayed CDF, including its denominator across all combinations.
line = next(p for p in root.iter(SVG + 'path')
            if 'stroke:rgb(89.399' in p.get('style', ''))
cdf_y = np.array(numbers(line.get('d'))[1::2])
extracted_cdf = (axis_bottom - (18 + 61.512362 - 0.995957 * cdf_y)) / axis_height * 100
cdf = extracted_cdf
assert len(cdf) == 25 and np.all(np.diff(cdf) > 0)

labels = ['ATAC', 'H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K9me3',
          'H3K27me3', 'DNA methylation']
with OUTPUT.with_suffix('.csv').open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['rank', 'occurrence_count', 'cumulative_percent', *labels])
    for i in range(25):
        writer.writerow([i + 1, counts[i], cdf[i], *matrix[:, i].astype(int)])

arial = Path('/lustre/home/2101110354/.local/share/fonts/arial.ttf')
font_manager.fontManager.addfont(str(arial))
plt.rcParams.update({'font.family': 'Arial', 'font.size': 5,
                     'axes.labelsize': 5, 'xtick.labelsize': 5,
                     'ytick.labelsize': 5, 'pdf.fonttype': 42,
                     'ps.fonttype': 42, 'axes.linewidth': 0.5})
width_pt, height_pt = 7.5 / 2.54 * 72, 4 / 2.54 * 72
fig = plt.figure(figsize=(7.5 / 2.54, 4 / 2.54))
left, right = 45, width_pt - 25
matrix_bottom, matrix_height = 4.5, 42  # Seven rows, exactly 6 pt pitch.
hist_bottom, hist_top = 51, height_pt - 5
ax = fig.add_axes([left / width_pt, hist_bottom / height_pt,
                   (right - left) / width_pt, (hist_top - hist_bottom) / height_pt])
mx = fig.add_axes([left / width_pt, matrix_bottom / height_pt,
                   (right - left) / width_pt, matrix_height / height_pt])
x = np.arange(25)
ax.bar(x, counts, width=0.8, color='#87ceeb', edgecolor='black', linewidth=0.3)
ax.set(xlim=(-0.5, 24.5), ylim=(0, 800), ylabel='Occurrence count')
ax.set_yticks([0, 200, 400, 600, 800])
ax.set_xticks([])
ax.tick_params(axis='y', length=1.5, width=0.5, pad=1.5)
ax.yaxis.labelpad = 2
red = '#e41a1c'
ay = ax.twinx()
ay.plot(x, cdf, '-o', color=red, linewidth=0.65, markersize=2,
        markeredgecolor='black', markeredgewidth=0.3)
ay.set(ylim=(0, 100), ylabel='Cumulative %')
ay.set_yticks([0, 20, 40, 60, 80, 100])
ay.tick_params(axis='y', colors=red, length=1.5, width=0.5, pad=1.5)
ay.yaxis.label.set_color(red)
ay.yaxis.labelpad = 2

for row in range(7):
    mx.axhline(row, color='lightgray', linewidth=0.3, zorder=0)
for col in range(25):
    rows = np.flatnonzero(matrix[:, col])
    mx.plot([col, col], [rows.min(), rows.max()], color='black', linewidth=0.5)
    mx.scatter(np.full(len(rows), col), rows, s=4, color='black', zorder=3)
mx.set(xlim=(-0.5, 24.5), ylim=(6.5, -0.5))
mx.set_yticks(range(7), labels=labels)
mx.set_xticks([])
mx.tick_params(axis='y', length=0, pad=2)
for spine in mx.spines.values():
    spine.set_visible(False)

# Do not use bbox_inches='tight': the PDF page must remain exactly 7.5 x 4 cm.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
for a in (ax, ay, mx):
    for label in [*a.get_yticklabels(), a.yaxis.label]:
        if not label.get_text():
            continue
        box = label.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0
        assert box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
fig.savefig(str(OUTPUT) + '.pdf')
fig.savefig(str(OUTPUT) + '.png', dpi=600)
print('Recovered counts:', counts.tolist())
print('CDF percentages preserved from vector coordinates.')
print('Saved:', str(OUTPUT) + '.pdf')
