#!/usr/bin/env Rscript

# Reproduces the cell-level UU-FF endpoint comparison reported in the article.
# Run from any directory; paths are resolved from this script's location.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
input_path <- file.path(
  repo_root, "data", "derived", "endpoint_thermal_retained_cores_83.csv"
)
output_dir <- file.path(repo_root, "results", "reproduced_endpoint_statistics")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
dat <- dat[dat$core_default %in% c("UU", "FF"), ]

uu <- dat$dlst[dat$core_default == "UU"]
ff <- dat$dlst[dat$core_default == "FF"]

stopifnot(length(uu) == 47L, length(ff) == 17L)

welch <- t.test(uu, ff, var.equal = FALSE, conf.level = 0.95)
wilcoxon <- wilcox.test(uu, ff, exact = FALSE)

n_uu <- length(uu)
n_ff <- length(ff)
pooled_sd <- sqrt(
  ((n_uu - 1) * var(uu) + (n_ff - 1) * var(ff)) /
    (n_uu + n_ff - 2)
)
cohens_d <- (mean(uu) - mean(ff)) / pooled_sd
hedges_correction <- 1 - 3 / (4 * (n_uu + n_ff) - 9)
hedges_g <- hedges_correction * cohens_d

result <- data.frame(
  n_UU = n_uu,
  n_FF = n_ff,
  mean_UU_C = mean(uu),
  mean_FF_C = mean(ff),
  difference_UU_minus_FF_C = mean(uu) - mean(ff),
  welch_ci95_low_C = unname(welch$conf.int[1]),
  welch_ci95_high_C = unname(welch$conf.int[2]),
  welch_p = welch$p.value,
  wilcoxon_p = wilcoxon$p.value,
  hedges_g = hedges_g
)

write.csv(
  result,
  file.path(output_dir, "endpoint_UU_FF_statistics.csv"),
  row.names = FALSE
)

print(result, digits = 8)

