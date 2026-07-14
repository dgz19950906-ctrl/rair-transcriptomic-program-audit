#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
sample_path <- file.path(root, "phase1_cross_definition", "processed",
                         "GSE151179_primary_preRAI_samples.tsv")
out_dir <- file.path(root, "phase1_cross_definition", "results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.delim(sample_path, check.names = FALSE, stringsAsFactors = FALSE)
group_order <- c("RAI_avid_remission", "RAI_avid_persistent",
                 "RAI_nonavid_persistent")
d$analysis_group <- factor(d$analysis_group, levels = group_order)

summaries <- list()
for (variable in c("tumor_purity_class", "lesion_class", "histological_variant")) {
  tab <- table(d$analysis_group, d[[variable]], useNA = "ifany")
  long <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(long) <- c("analysis_group", "level", "n")
  long$variable <- variable
  if (nrow(tab) > 1L && ncol(tab) > 1L) {
    long$fisher_exact_p <- fisher.test(tab)$p.value
  } else {
    long$fisher_exact_p <- NA_real_
  }
  summaries[[variable]] <- long[, c("variable", "analysis_group", "level", "n",
                                    "fisher_exact_p")]
}
out <- do.call(rbind, summaries)
write.table(out, file.path(out_dir, "GSE151179_primary_confounder_balance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)

