#!/usr/bin/env Rscript

# Locked, label-independent thyroid differentiation score (TDS) analysis.
# No feature selection or outcome-trained weights are used.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))

expr_path <- file.path(root, "phase1_cross_definition", "processed",
                       "GSE151179_primary_preRAI_gene_expression.tsv.gz")
sample_path <- file.path(root, "phase1_cross_definition", "processed",
                         "GSE151179_primary_preRAI_samples.tsv")
table_dir <- file.path(root, "phase1_cross_definition", "results", "tables")
figure_dir <- file.path(root, "phase1_cross_definition", "results", "figures")
qc_dir <- file.path(root, "phase1_cross_definition", "qc")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

tds_genes <- c("DIO1", "DIO2", "DUOX1", "DUOX2", "FOXE1", "GLIS3",
               "NKX2-1", "PAX8", "SLC26A4", "SLC5A5", "SLC5A8", "TG",
               "THRA", "THRB", "TPO", "TSHR")

expr <- read.delim(expr_path, row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)
samples <- read.delim(sample_path, check.names = FALSE,
                      stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), samples$geo_accession))

present <- intersect(tds_genes, rownames(expr))
missing <- setdiff(tds_genes, rownames(expr))
if (length(present) < 12L) {
  stop(sprintf("Only %d/16 locked TDS genes are present", length(present)))
}

# Gene-wise z standardisation over all eligible samples, followed by an
# unweighted arithmetic mean. This is fixed before looking at outcome groups.
z <- t(scale(t(expr[present, , drop = FALSE])))
tds <- colMeans(z, na.rm = TRUE)
samples$TDS <- unname(tds[samples$geo_accession])

group_levels <- c("RAI_avid_remission", "RAI_avid_persistent",
                  "RAI_nonavid_persistent")
samples$analysis_group <- factor(samples$analysis_group, levels = group_levels)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  df <- nx + ny - 2
  pooled <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / df)
  d <- (mean(x) - mean(y)) / pooled
  correction <- 1 - 3 / (4 * df - 1)
  correction * d
}

exact_permutation_p <- function(x, y) {
  pooled <- c(x, y)
  nx <- length(x)
  observed <- mean(x) - mean(y)
  idx <- combn(seq_along(pooled), nx)
  perm_diff <- apply(idx, 2, function(i) mean(pooled[i]) - mean(pooled[-i]))
  # Plus-one convention retains valid behaviour even if later changed to sampled permutations.
  (sum(abs(perm_diff) >= abs(observed) - 1e-12) + 1) / (length(perm_diff) + 1)
}

lopo_stability <- function(x, y) {
  full <- sign(mean(x) - mean(y))
  leave_x <- vapply(seq_along(x), function(i) sign(mean(x[-i]) - mean(y)), numeric(1))
  leave_y <- vapply(seq_along(y), function(i) sign(mean(x) - mean(y[-i])), numeric(1))
  mean(c(leave_x, leave_y) == full)
}

contrasts <- list(
  uptake_failure = c("RAI_nonavid_persistent", "RAI_avid_persistent"),
  response_failure_with_uptake = c("RAI_avid_persistent", "RAI_avid_remission"),
  broad_persistence = c("PERSISTENT_COMBINED", "RAI_avid_remission")
)

get_values <- function(label) {
  if (label == "PERSISTENT_COMBINED") {
    samples$TDS[samples$analysis_group %in%
                  c("RAI_avid_persistent", "RAI_nonavid_persistent")]
  } else {
    samples$TDS[samples$analysis_group == label]
  }
}

contrast_rows <- lapply(names(contrasts), function(name) {
  labels <- contrasts[[name]]
  adverse <- get_values(labels[1])
  reference <- get_values(labels[2])
  data.frame(
    contrast = name,
    adverse_group = labels[1],
    reference_group = labels[2],
    n_adverse = length(adverse),
    n_reference = length(reference),
    mean_adverse = mean(adverse),
    mean_reference = mean(reference),
    mean_difference_adverse_minus_reference = mean(adverse) - mean(reference),
    hedges_g = hedges_g(adverse, reference),
    exact_two_sided_permutation_p = exact_permutation_p(adverse, reference),
    leave_one_patient_out_direction_stability = lopo_stability(adverse, reference),
    stringsAsFactors = FALSE
  )
})
contrast_table <- do.call(rbind, contrast_rows)

gene_scores <- do.call(rbind, lapply(present, function(gene) {
  values <- z[gene, ]
  rows <- lapply(names(contrasts)[1:2], function(name) {
    labels <- contrasts[[name]]
    adverse_ids <- if (labels[1] == "PERSISTENT_COMBINED") {
      samples$geo_accession[samples$analysis_group %in%
                              c("RAI_avid_persistent", "RAI_nonavid_persistent")]
    } else samples$geo_accession[samples$analysis_group == labels[1]]
    reference_ids <- samples$geo_accession[samples$analysis_group == labels[2]]
    adverse <- values[adverse_ids]
    reference <- values[reference_ids]
    data.frame(gene = gene, contrast = name,
               mean_difference_adverse_minus_reference = mean(adverse) - mean(reference),
               hedges_g = hedges_g(adverse, reference),
               exact_two_sided_permutation_p = exact_permutation_p(adverse, reference),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}))
gene_scores$BH_FDR_within_contrast <- ave(
  gene_scores$exact_two_sided_permutation_p, gene_scores$contrast,
  FUN = function(p) p.adjust(p, method = "BH")
)

write.table(samples, file.path(table_dir, "GSE151179_TDS_sample_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contrast_table, file.path(table_dir, "GSE151179_TDS_contrasts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_scores, file.path(table_dir, "GSE151179_TDS_gene_components.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

group_labels <- c("RAI-avid\nremission", "RAI-avid\npersistent",
                  "RAI-nonavid\npersistent")
group_cols <- c("#3B82F6", "#F59E0B", "#DC2626")
pdf(file.path(figure_dir, "GSE151179_locked_TDS_three_groups.pdf"),
    width = 7.2, height = 5.2, useDingbats = FALSE)
par(mar = c(5.4, 4.6, 1.2, 1.0), las = 1)
boxplot(TDS ~ analysis_group, data = samples, names = group_labels,
        col = paste0(group_cols, "35"), border = group_cols,
        ylab = "Locked thyroid differentiation score (mean gene-wise z)",
        xlab = "", outline = FALSE, frame.plot = FALSE)
set.seed(151179)
xpos <- jitter(as.numeric(samples$analysis_group), amount = 0.08)
points(xpos, samples$TDS, pch = 21, bg = group_cols[as.numeric(samples$analysis_group)],
       col = "white", cex = 1.35, lwd = 0.7)
abline(h = 0, lty = 3, col = "grey60")
dev.off()

qc_lines <- c(
  sprintf("locked_tds_genes_requested\t%d", length(tds_genes)),
  sprintf("locked_tds_genes_present\t%d", length(present)),
  sprintf("present_genes\t%s", paste(present, collapse = ",")),
  sprintf("missing_genes\t%s", ifelse(length(missing), paste(missing, collapse = ","), "none")),
  "score_definition\tmean of gene-wise z scores across all 17 eligible samples",
  "feature_selection\tnone",
  "outcome_trained_weights\tnone",
  "primary_contrasts\tuptake failure; response failure despite uptake"
)
writeLines(qc_lines, file.path(qc_dir, "locked_TDS_qc.tsv"))

print(contrast_table, row.names = FALSE)
cat(sprintf("\nLocked TDS genes present: %d/16; missing: %s\n",
            length(present), ifelse(length(missing), paste(missing, collapse = ", "), "none")))
