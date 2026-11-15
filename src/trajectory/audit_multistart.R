#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
})

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
tables <- file.path(root, "results", "trajectory", "tables")
models <- file.path(root, "results", "trajectory", "models")
reports <- file.path(root, "results", "trajectory", "reports")

source_file <- file.path(tables, "all_icu_gbmtm_prepared_training.csv")
selected_file <- file.path(models, "gbmtm_all_icu_mimiciv_d1d7_easix_siri_selected.rds")
if (!file.exists(source_file)) stop("Prepared GBMTM training data were not exported.")

d <- fread(source_file)
model_data <- as.data.frame(d[, .(id, time, log2_easix_interp_z, log2_siri_interp_z)])
selected <- readRDS(selected_file)
reference <- selected$assign
selected_k <- length(selected$prior)

adjusted_rand <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  index <- sum(choose2(tab))
  row_sum <- sum(choose2(rowSums(tab)))
  col_sum <- sum(choose2(colSums(tab)))
  expected <- row_sum * col_sum / choose2(n)
  maximum <- (row_sum + col_sum) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

seeds <- c(20260727L, 20260731L, 20260803L, 20260807L, 20260811L)
rows <- vector("list", length(seeds))
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  fit <- gbmt(
    x.names = c("log2_easix_interp_z", "log2_siri_interp_z"),
    unit = "id", time = "time", d = 2L, ng = selected_k,
    data = model_data, scaling = 0, quiet = TRUE
  )
  props <- as.numeric(prop.table(table(factor(fit$assign, levels = seq_len(selected_k)))))
  rows[[i]] <- data.table(
    seed = seeds[i], AIC = unname(fit$ic["aic"]), BIC = unname(fit$ic["bic"]),
    min_class_pct = 100 * min(props), ari_vs_selected = adjusted_rand(reference, fit$assign),
    identical_assignment = identical(as.integer(reference), as.integer(fit$assign))
  )
}
out <- rbindlist(rows)
output_stem <- paste0("all_icu_gbmtm_k", selected_k, "_multistart_stability")
fwrite(out, file.path(tables, paste0(output_stem, ".csv")))
writeLines(c(
  paste0("# GBMTM K=", selected_k, " multistart stability audit"), "",
  paste0("- Number of initializations: ", nrow(out), "."),
  paste0("- BIC range: ", sprintf("%.3f", min(out$BIC)), " to ", sprintf("%.3f", max(out$BIC)), "."),
  paste0("- ARI range relative to the selected solution: ", sprintf("%.3f", min(out$ari_vs_selected)), " to ", sprintf("%.3f", max(out$ari_vs_selected)), "."),
  paste0("- Reproducibility of selected K=", selected_k, " across random initializations.")
), file.path(reports, paste0(output_stem, "_report.md")))
