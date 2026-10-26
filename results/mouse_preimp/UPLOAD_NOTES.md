# GitHub export scope and reproduction

This directory is part of the four-project mouse source export dated 2026-10-26.

Included: retained analysis/plotting source, per-sample job scripts where needed, signal path lists,
operational manifests, documentation, and selected small numerical tables/test snapshots used by the final figures.
Excluded: all `legacy/` and historical `provenance/` code, TF binding-site/motif analyses, raw genomic tracks,
large entropy/O-information matrices and permutation tensors, installed environments, figures and execution logs.
No TF binding-site/motif analysis was identified among the selected non-legacy source files.
CTCF-bound/unbound annotation strata are retained as part of the existing genomic block definitions.

The export preserves original scientific source bytes, parameters, validation evidence and HPC path strings.
**It is a source archive, not a self-contained dataset or a newly rerun analysis.** Before running on a new
machine, update original workspace/interpreter/annotation paths and Slurm partition/resource directives,
and restore the required upstream data. Commands in the retained guides refer to the original workspace layout.
In this repository these four projects live under `results/`, whereas the original workspace used `project/`.
Historical provenance paths document the original calculations rather than a new calculation in this repository.
The original local workspace and its data/results are unchanged by this upload.

## Shared dependencies and project relationships

Eight required shared calculation modules are archived at
`../mouse_development/shared/hoi_code/`, with a source SHA256 manifest.
Use those modules when adapting the original `script/hoi_code/` references; other versions in the repository
need not match the archived computing snapshot. Python package requirements are listed in
`../mouse_development/requirements.txt`; the archive does not install packages or promise an exact environment lock.
The seven-signal 10,000-permutation development run also retains its own validated `code/` snapshot.

| Project | Analysis panel |
|---|---|
| `mouse_preimp` | Seven early stages; nine-signal full systems, BH84 triplets, 3/5/10-bin entropy and MI analyses |
| `mouse_after_E11.5` | 66 late tissue-stage samples; ten-signal systems and eight-signal triplet permutation inputs |
| `mouse_development` | 82 samples / 17 stages; seven-signal triplets with per-sample × block BH35 and developmental plots |
| `H3K36me3_developmental_switch` | Shared-eight-signal BH56 developmental triplets, redundancy cores and source-replacement support |

The switch project consumes original preimplantation/late results and selected BH56 snapshots.
Its full-signal occurrence branch uses a sign-only triplet definition, unlike the BH56 branch.
The seven-signal development project has a different signal universe and sampling panel;
its plots cannot isolate the effect of adding H3K36me3 by direct comparison with eight-signal plots.

For the TC/DTC developmental figures, the original raw inputs are outside these four projects:
`project/permutation_test/2026-09-28/results_mouse_development/` and
`project/mouse_development_fine_te/`. The exact small differentiation tree is included in
`mouse_development/plot_tc_dtc/inputs/`; the cached 492-row plotting table and validation metadata are included.
The retained plotting code still references its original tree/archive locations, which must be adapted or restored.

`SOURCE_MANIFEST.json` lists the exact files and symlinks shipped in this project.
The source export itself does not submit any jobs, refit models or regenerate scientific outputs.
