# H3K36me3 developmental switch: synergy, redundancy and source robustness

This project examines developmental changes in the information-theoretic participation of H3K36me3,
using original developmental synergistic triplets and redundancy cores as the main analyses.
Source-replacement robustness is a supporting analysis.

Here **developmental switch** names the observed developmental change in information organization;
it does not establish a causal molecular switch or a uniform transition at every stage/region.
The early and late panels use different tissues, and the interval E3.5–E11.5 is not sampled here.
Triplet occurrences and core occurrences are different counts, with different selection rules;
their raw magnitudes are not a common quantitative scale.

## Current analyses

| Directory | Purpose |
|---|---|
| `synergistic_triplets/` | Original eight-signal developmental occurrence trajectories, ABC heatmaps and selected-region trajectories; 10,000 permutations, BH56 per sample × block, negative O-information and q≤0.05 |
| `redundancy_core/` | Eight-signal redundancy cores and paired synergy/core trajectories; cores use the replacement-background algorithm, with no permutation P values or FDR assigned to cores |
| `redundancy_core/full_signal_occurrence/` | Nine-signal early and ten-signal late combination occurrence plots; triplets use sign-only selection, without FDR |
| `batch_effect_figures/` | Four-mark, 16-source-combination replacements at 2cell/8cell/ICM, BH56 participation and stability summaries |

The eight-signal core computation covers 7 early + 66 late samples. Main matched developmental
figures use the seven early samples and five common late tissues (Forebrain, Midbrain, Hindbrain,
Heart and Liver) at each late stage; ICM and TE are combined under E3.5 in the paired core/triplet trajectories.

## Evidence motivating the name

In `redundancy_core/results/synergistic_triplet_mark_trajectories.tsv`, H3K36me3 mean participation
is 34 at E3.5 and 178.6 at E11.5, remaining 195.6–246.8 across E12.5–E16.5 (224 at P0).
In `redundancy_core/results/redundancy_core_mark_trajectories.tsv`, its core participation is 27.5
at E3.5 and 15.2 at E11.5, falling to 6.4 at E16.5 (9.4 at P0).
These descriptive trajectories motivate a focus on a developmental shift between redundant and
synergistic organization; they do not locate a switch time or establish causality.

Main paired figure: `redundancy_core/plots/combined_triplet_core_trajectories.pdf`.
Original occurrence figure: `synergistic_triplets/plots/mark_occurrence_per_tissue_by_stage_combined.pdf`.
Full-signal summary: `redundancy_core/full_signal_occurrence/plots/full_signal_triplets_and_cores_2x2.pdf`.
Robustness support: `batch_effect_figures/plots/same_stage_bulk_replacement_all_marks.pdf`.
See each subdirectory README for inputs, methods and commands.

## Related mouse projects and rename

- `../mouse_preimp/`: original nine-signal preimplantation calculations and MI analyses.
- `../mouse_after_E11.5/`: original ten-signal late systems and shared-eight-signal triplet permutations.
- `../mouse_development/`: separate seven-signal integration across 82 samples, including gastrulation.

Renamed on 2026-10-26 from `H3K36me3_triplet_robustness` to `H3K36me3_developmental_switch`.
Live scripts, input lists and operational manifests were updated, together with external live-code references.
The local workspace retains historical provenance and rename records. These historical folders and logs
are excluded from this GitHub export. See UPLOAD_NOTES.md and SOURCE_MANIFEST.json for the uploaded scope.
No analyses were rerun or scientific results regenerated for the rename.

## GitHub source export

See [UPLOAD_NOTES.md](UPLOAD_NOTES.md) for the uploaded scope, excluded legacy/TF analyses,
shared calculation snapshot, original HPC paths and external data requirements.
