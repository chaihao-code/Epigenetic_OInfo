# Paths and data restoration

Run commands from the release root and export `EPINFO_ROOT="$PWD"`. Python adapters locate `.epigenetic_oinfo_root`; shell/Slurm programs inherit `EPINFO_ROOT`. Set `EPINFO_PYTHON` to the desired interpreter before submitting `demo/submit_demo.slurm`. On the supplied cluster do not add `--mem`.

Local input symlinks are listed in `config/local_input_links.tsv`. They point to existing original inputs and are omitted from the GitHub checkout and portable archives. No result output directory is symlinked to the original workspace. The demo includes compressed byte-for-byte copies of the complete seven prepared thyroid signal tracks and original genomic grid. If local `.bp` links are unavailable, `demo/run_demo.py` restores those supplied copies without changing their contents.

Full original/derived thyroid bigWigs are listed in `demo/input/bigwig_sources.tsv`. Restore those full files at the named paths to execute `--from-bigwig`. Original histone/Input source URLs are direct bigWigs; methylation's URL is a CpG BED.gz and its locally converted bigWig is a separately documented derived input.

The existing annotation masks under `data/genome_features/` are fixed inputs. Compressed counterparts are supplied; `python tools/restore_annotation_masks.py` only decompresses files and does not rebuild annotations. The thyroid demo reads compressed masks directly. Existing class-coordinate exports under `supplementary_data/bins_per_class/` remain unchanged locally, but the supplementary_data directory is excluded from GitHub. Fixed annotation masks remain included under data/genome_features.

New calculation commands use `common/02_oinformation/compute_hoi_oinfo.py`. The case files bearing that same name only forward to the shared implementation. Old entropy-only/entropy-conversion programs and their launchers were removed. `tools/plan_oinfo.py` previews a command and checks paths; it does not execute analysis.

New triplet runs default to 10,000 permutations. Existing reference results retain their actual historical metadata/counts. Local values exported by the shared engine use the sampled matrix also used by global Omega; complete annotation classes and exact selected loci are stored separately.
