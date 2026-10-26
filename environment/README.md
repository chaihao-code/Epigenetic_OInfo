# Dependencies

The single exact-count global/local engine needs Python, NumPy and pandas; it does not require hoi/JAX. The complete thyroid bigWig demo additionally requires the UCSC bigWigAverageOverBed utility on PATH. The included generic bigWig helper uses pyBigWig. Install the Python demo dependencies from `demo/requirements.txt`.

Redundancy cores need scikit-learn/statsmodels; regressions use SciPy/statsmodels; plotting uses matplotlib; case GSEA uses gseapy; B-cell RNA uses its retained PyDESeq2/anndata/formulaic dependencies. DMR/import need R, data.table, limma, tximport and statmod. Case-specific requirements remain original dependency records. Actual demo software versions are saved in its validation JSON. No full-cohort environment upgrade or reanalysis is implied.
