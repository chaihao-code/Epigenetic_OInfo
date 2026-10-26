#!/usr/bin/env Rscript
lib <- Sys.getenv("R_LIBS_USER")
stopifnot(nzchar(lib))
dir.create(lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(lib, .libPaths()))
options(timeout=300)
for (p in c("data.table", "BiocManager")) {
  if (!requireNamespace(p, quietly=TRUE)) {
    install.packages(p, lib=lib, repos="https://cloud.r-project.org", Ncpus=4)
  }
  stopifnot(requireNamespace(p, quietly=TRUE))
}
if (!requireNamespace("limma", quietly=TRUE)) {
  BiocManager::install("limma", lib=lib, ask=FALSE, update=FALSE, Ncpus=4)
}
stopifnot(requireNamespace("limma", quietly=TRUE))
sessionInfo()
