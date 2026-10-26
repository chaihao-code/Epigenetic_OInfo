#!/usr/bin/env Rscript
# DMR analysis v2: CEMT Group2 vs Group3 (1 kb bins), limma with covariates
#
# Model per described protocol:
#   M_ij = alpha_j + beta_j^group * G_i + gamma1_j * Age_i + gamma2_j * Sex_i + r_ij
#   M = log2((beta + eps)/(1 - beta + eps)), eps = 1e-4
# Empirical-Bayes moderated t test (limma), BH correction.
# DMR: |delta_beta| > 0.2 AND FDR < 0.05.
#
# Groups (from plot_oinfo_vs_age_grouped.py):
#   Group2 (n=8): CEMT0004/0005/0025/0026/0027/0028/0029/0030  (47-79y)
#   Group3 (n=5): CEMT0006/0094/0095/0096/0097                 (57-90y)
# delta_beta = Group2 - Group3.
#
# NOTE: all Group3 donors are male -> the sex covariate is partially
# confounded with group (sex effects are identified from Group2 alone).

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

EPS     <- 1e-4
DELTA   <- 0.2
FDR_CUT <- 0.05
MIN_N   <- 8    # minimum samples with non-missing beta per bin

IN_CSV  <- "input/all_bins_1kb_cemt_g2_g3.csv"
OUT_DIR <- "output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

sample_meta <- data.table(
  sample = c('CEMT0004','CEMT0005','CEMT0025','CEMT0026','CEMT0027',
             'CEMT0028','CEMT0029','CEMT0030',
             'CEMT0006','CEMT0094','CEMT0095','CEMT0096','CEMT0097'),
  group = factor(c(rep('Group2', 8), rep('Group3', 5)), levels = c('Group3', 'Group2')),
  age   = c(74, 68, 62, 79, 71, 79, 47, 56, 60, 69, 71, 57, 90),
  sex   = factor(c('male','female','female','male','male','female','male','female',
                   'male','male','male','male','male'), levels = c('female','male'))
)
sample_ids <- sample_meta$sample

message("Reading beta matrix ...")
beta <- fread(IN_CSV, select = c("chrom", "start", "end", sample_ids))
message(sprintf("  bins: %s", format(nrow(beta), big.mark = ",")))

canonical <- paste0("chr", c(1:22))
beta <- beta[chrom %in% canonical]
message(sprintf("  canonical autosomes (chr1-22): %s", format(nrow(beta), big.mark = ",")))

beta_mat <- as.matrix(beta[, .SD, .SDcols = sample_ids])
storage.mode(beta_mat) <- "double"

n_ok <- rowSums(!is.na(beta_mat))
keep <- n_ok >= MIN_N
message(sprintf("  bins with >= %d covered samples: %s (%.1f%%)",
                MIN_N, format(sum(keep), big.mark = ","), 100 * mean(keep)))
beta_mat <- beta_mat[keep, ]
beta <- beta[keep, ]

M_mat <- log2((beta_mat + EPS) / (1 - beta_mat + EPS))

design <- model.matrix(~ group + age + sex, data = sample_meta)
colnames(design) <- c("Intercept", "groupGroup2", "age", "sexmale")
message("Design matrix:")
print(design)

message("Fitting linear models (limma) ...")
fit <- lmFit(M_mat, design)
fit <- eBayes(fit, trend = TRUE)

message("Extracting group coefficient results ...")
tt <- topTable(fit, coef = "groupGroup2", number = Inf,
               adjust.method = "BH", sort.by = "none")

g2_cols <- sample_ids[sample_meta$group == "Group2"]
g3_cols <- sample_ids[sample_meta$group == "Group3"]
g2_mean <- rowMeans(beta[, .SD, .SDcols = g2_cols], na.rm = TRUE)
g3_mean <- rowMeans(beta[, .SD, .SDcols = g3_cols], na.rm = TRUE)
g2_mean[is.nan(g2_mean)] <- NA
g3_mean[is.nan(g3_mean)] <- NA
delta_beta <- g2_mean - g3_mean

res <- data.table(
  chrom        = beta$chrom,
  start        = beta$start,
  end          = beta$end,
  group2_mean  = round(g2_mean, 4),
  group3_mean  = round(g3_mean, 4),
  delta_beta   = round(delta_beta, 4),
  M_logFC      = round(tt$logFC, 4),
  t            = round(tt$t, 3),
  pvalue       = tt$P.Value,
  qvalue       = tt$adj.P.Val,
  B            = round(tt$B, 3)
)
res[, direction := fifelse(delta_beta > 0, "group2_high", "group3_high")]
res[, is_dmr := !is.na(delta_beta) & abs(delta_beta) > DELTA & qvalue < FDR_CUT]

message("")
message("=== Summary ===")
dmrs <- res[is_dmr == TRUE]
message(sprintf("bins tested:                  %s", format(nrow(res), big.mark = ",")))
message(sprintf("bins with |dBeta| > %.2f:      %s", DELTA,
                format(sum(!is.na(res$delta_beta) & abs(res$delta_beta) > DELTA), big.mark = ",")))
message(sprintf("bins with FDR < %.2f:          %s", FDR_CUT,
                format(sum(res$qvalue < FDR_CUT, na.rm = TRUE), big.mark = ",")))
message(sprintf("min qvalue observed:          %.4g", min(res$qvalue, na.rm = TRUE)))
message(sprintf("DMRs (|dBeta|>%.2f & FDR<%.2f): %s", DELTA, FDR_CUT, format(nrow(dmrs), big.mark = ",")))
if (nrow(dmrs) > 0) print(dmrs[, .N, by = direction][order(-N)])

out_all <- file.path(OUT_DIR, "all_bins_1kb_cemt_g2_g3_limma_results.csv")
fwrite(res, out_all)
message(sprintf("Saved: %s", out_all))

out_dmr <- file.path(OUT_DIR, "cemt_g2_vs_g3_dmrs_1kb_limma.csv")
dmrs[, abs_delta := abs(delta_beta)]
setorder(dmrs, -abs_delta)
fwrite(dmrs, out_dmr)
message(sprintf("Saved: %s", out_dmr))

message("")
message("Top 10 DMRs by |delta_beta|:")
print(head(dmrs[, .(chrom, start, end, group2_mean, group3_mean,
                    delta_beta, pvalue, qvalue, direction)], 10))
message("Done.")
