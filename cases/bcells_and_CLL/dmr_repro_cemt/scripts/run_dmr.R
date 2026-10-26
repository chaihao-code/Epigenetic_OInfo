#!/usr/bin/env Rscript
# Execution wrapper; the supplied analysis program is sourced unchanged.
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
if (length(script_arg) != 1) stop("Run this wrapper with Rscript")
setwd(dirname(dirname(normalizePath(sub("^--file=", "", script_arg)))))
dir.create("output", showWarnings=FALSE, recursive=TRUE)
suppressPackageStartupMessages(library(data.table))
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
setDTthreads(threads)
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
writeLines(c(paste("started_at", format(Sys.time(), tz="UTC", usetz=TRUE)),
             paste("R", R.version.string),
             paste("data.table", packageVersion("data.table")),
             paste("limma", packageVersion("limma")),
             paste("threads", threads)), "output/run_environment.txt")
source("scripts/dmr_limma_1kb_cemt_g2_g3.R", echo=FALSE)
stopifnot(qr(design)$rank == ncol(design), nrow(res) == nrow(beta_mat))
stopifnot(all(is.na(res$qvalue) | (res$qvalue >= 0 & res$qvalue <= 1)))
expected_q <- p.adjust(res$pvalue, method="BH")
stopifnot(isTRUE(all.equal(expected_q, res$qvalue, tolerance=1e-12)))
stopifnot(nrow(dmrs) == sum(res$is_dmr, na.rm=TRUE))
fwrite(sample_meta, "output/sample_metadata.tsv", sep="\t")
fwrite(data.table(sample=sample_ids, as.data.table(design)), "output/design_matrix.tsv", sep="\t")
summary <- data.table(
  n_bins_tested=nrow(res), n_finite_p=sum(is.finite(res$pvalue)),
  n_FDR_lt005=sum(res$qvalue < FDR_CUT, na.rm=TRUE),
  n_abs_delta_gt02=sum(abs(res$delta_beta) > DELTA, na.rm=TRUE),
  n_DMR=nrow(dmrs),
  n_Group2_high=sum(dmrs$direction == "group2_high"),
  n_Group3_high=sum(dmrs$direction == "group3_high"),
  minimum_q=if (any(is.finite(res$qvalue))) min(res$qvalue, na.rm=TRUE) else NA_real_,
  design_rank=qr(design)$rank, design_columns=ncol(design),
  BH_recalculation_verified=TRUE,
  delta_threshold_applied_after_original_rounding=TRUE
)
fwrite(summary, "output/run_summary.tsv", sep="\t")
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
message("Wrapper verification passed.")
