#!/usr/bin/env Rscript
args <- grep("^--file=", commandArgs(), value=TRUE)
if (length(args) != 1) stop("Run with Rscript")
root <- dirname(normalizePath(sub("^--file=", "", args)))
while (!file.exists(file.path(root, ".epigenetic_oinfo_root"))) {
  parent <- dirname(root)
  if (parent == root) stop("Release root not found")
  root <- parent
}
status <- system2(file.path(R.home("bin"), "Rscript"),
                  shQuote(file.path(root, "cases/bcells_and_CLL/dmr_repro_cemt/scripts/run_dmr.R")))
quit(status=status)
