#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

project <- "/home/dony/ThyroidCancer_Project"
audit <- file.path(project, "rair_audit")
input_path <- file.path(audit, "checkpoints", "03_mnn_clustered_counts.rds")
definition_path <- file.path(audit, "manifests", "frozen_signatures.tsv")
gmt_path <- file.path(audit, "manifests", "h.all.v2025.1.Hs.symbols.gmt")

minimum_cells <- 20L
minimum_donors <- 3L
minimum_estimable_iterations <- 9L
minimum_same_direction <- 9L

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), ..., "\n", sep = "")
  flush.console()
}

if (!file.exists(input_path)) stop("Stage 03 checkpoint missing: ", input_path)
if (!file.exists(definition_path)) stop("Frozen signature registry missing: ", definition_path)
log_line("Loading MNN-clustered counts object")
obj <- readRDS(input_path)
DefaultAssay(obj) <- "RNA"
counts <- LayerData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data[colnames(counts), , drop = FALSE]

eligible_cells <- rownames(meta)[
  !is.na(meta$patient) & !is.na(meta$compartment) & meta$compartment != "Uncertain"
]
counts <- counts[, eligible_cells, drop = FALSE]
meta <- meta[eligible_cells, , drop = FALSE]
meta$group <- paste(meta$patient, meta$compartment, sep = "||")

cell_count_table <- as.data.frame(table(meta$patient, meta$compartment), stringsAsFactors = FALSE)
names(cell_count_table) <- c("donor", "compartment", "cells")
cell_count_table$qc_status <- ifelse(cell_count_table$cells >= minimum_cells, "qualified", "insufficient_cells")
write.table(
  cell_count_table,
  file.path(audit, "results", "04_donor_compartment_cell_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

qualified_groups <- with(
  cell_count_table[cell_count_table$cells >= minimum_cells, ],
  paste(donor, compartment, sep = "||")
)
qualified_cells <- rownames(meta)[meta$group %in% qualified_groups]
counts <- counts[, qualified_cells, drop = FALSE]
meta <- meta[qualified_cells, , drop = FALSE]
group_levels <- sort(unique(meta$group))

indicator <- sparseMatrix(
  i = seq_len(nrow(meta)),
  j = match(meta$group, group_levels),
  x = 1,
  dims = c(nrow(meta), length(group_levels)),
  dimnames = list(rownames(meta), group_levels)
)
ordered_counts <- counts[, rownames(indicator), drop = FALSE]
pseudobulk_counts <- ordered_counts %*% indicator
pseudobulk_library <- Matrix::colSums(pseudobulk_counts)
logcpm <- log1p(t(t(pseudobulk_counts) / pseudobulk_library * 1e6))

pb_meta <- do.call(rbind, strsplit(colnames(logcpm), "\\|\\|"))
pb_meta <- data.frame(
  pseudobulk = colnames(logcpm), donor = pb_meta[, 1], compartment = pb_meta[, 2],
  library_size = pseudobulk_library, stringsAsFactors = FALSE
)
rownames(pb_meta) <- pb_meta$pseudobulk

defs <- read.delim(definition_path, check.names = FALSE)
if (file.exists(gmt_path)) {
  control_ids <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
    "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
  )
  gmt <- strsplit(readLines(gmt_path), "\t", fixed = TRUE)
  names(gmt) <- vapply(gmt, `[`, character(1), 1)
  for (signature_id in intersect(control_ids, names(gmt))) {
    genes <- gmt[[signature_id]][-(1:2)]
    defs <- rbind(
      defs,
      data.frame(
        signature_id = signature_id, gene = genes, direction = 1,
        orientation = "higher_more_aggressive_program",
        source = "MSigDB Hallmark v2025.1.Hs",
        role = "aggressiveness_negative_control", lock_status = "locked"
      )
    )
  }
}

score_programs <- function(expression, definition) {
  z <- t(scale(t(as.matrix(expression))))
  z[!is.finite(z)] <- 0
  score_list <- list()
  coverage_list <- list()
  for (signature_id in unique(definition$signature_id)) {
    current <- definition[definition$signature_id == signature_id, , drop = FALSE]
    present <- current$gene[current$gene %in% rownames(z)]
    minimum_present <- max(3L, ceiling(0.7 * nrow(current)))
    coverage_list[[signature_id]] <- data.frame(
      signature_id = signature_id,
      genes_requested = nrow(current),
      genes_present = length(present),
      genes_required = minimum_present,
      passes_coverage = length(present) >= minimum_present,
      missing_genes = paste(setdiff(current$gene, present), collapse = ";")
    )
    if (length(present) < minimum_present) next
    weights <- current$direction[match(present, current$gene)]
    score <- colSums(z[present, , drop = FALSE] * weights) / sum(abs(weights))
    score_list[[signature_id]] <- data.frame(
      pseudobulk = names(score), signature_id = signature_id,
      score = unname(score), orientation = current$orientation[1],
      role = current$role[1]
    )
  }
  list(scores = do.call(rbind, score_list), coverage = do.call(rbind, coverage_list))
}

paired_preference <- function(score_table, metadata, target) {
  merged <- merge(score_table, metadata, by = "pseudobulk", all.x = TRUE)
  donors <- unique(merged$donor)
  differences <- vapply(
    donors,
    function(donor) {
      current <- merged[merged$donor == donor, , drop = FALSE]
      target_score <- current$score[current$compartment == target]
      other_score <- current$score[current$compartment != target]
      if (length(target_score) != 1L || !length(other_score)) return(NA_real_)
      target_score - mean(other_score)
    },
    numeric(1)
  )
  differences[is.finite(differences)]
}

full <- score_programs(logcpm, defs)
write.table(full$coverage, file.path(audit, "results", "04_program_gene_coverage.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
full_scores <- merge(full$scores, pb_meta, by = "pseudobulk", all.x = TRUE)
write.table(full_scores, file.path(audit, "results", "04_pseudobulk_program_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

preference_rows <- list()
selected_rows <- list()
for (signature_id in unique(full_scores$signature_id)) {
  current <- full_scores[full_scores$signature_id == signature_id, , drop = FALSE]
  compartments <- sort(unique(current$compartment))
  rows <- lapply(compartments, function(compartment) {
    differences <- paired_preference(current[, c("pseudobulk", "score")], pb_meta, compartment)
    data.frame(
      signature_id = signature_id, compartment = compartment,
      donors_with_paired_difference = length(differences),
      mean_paired_preference = if (length(differences)) mean(differences) else NA_real_,
      median_paired_preference = if (length(differences)) median(differences) else NA_real_,
      orientation = current$orientation[1], role = current$role[1]
    )
  })
  rows <- do.call(rbind, rows)
  preference_rows[[signature_id]] <- rows
  eligible <- rows[rows$donors_with_paired_difference >= minimum_donors & is.finite(rows$mean_paired_preference), , drop = FALSE]
  if (nrow(eligible)) selected_rows[[signature_id]] <- eligible[which.max(eligible$mean_paired_preference), , drop = FALSE]
}
full_preferences <- do.call(rbind, preference_rows)
selected <- do.call(rbind, selected_rows)
write.table(full_preferences, file.path(audit, "results", "04_all_compartment_preferences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

donors <- sort(unique(pb_meta$donor))
lodo_rows <- list()
for (omitted_donor in donors) {
  retained_columns <- rownames(pb_meta)[pb_meta$donor != omitted_donor]
  iteration <- score_programs(logcpm[, retained_columns, drop = FALSE], defs)
  iteration_scores <- merge(iteration$scores, pb_meta[retained_columns, , drop = FALSE], by = "pseudobulk", all.x = TRUE)
  for (row_index in seq_len(nrow(selected))) {
    signature_id <- selected$signature_id[row_index]
    target <- selected$compartment[row_index]
    current <- iteration_scores[iteration_scores$signature_id == signature_id, , drop = FALSE]
    target_donors <- unique(current$donor[current$compartment == target])
    differences <- paired_preference(current[, c("pseudobulk", "score")], pb_meta[retained_columns, , drop = FALSE], target)
    estimable <- length(target_donors) >= minimum_donors && length(differences) >= minimum_donors
    effect <- if (estimable) mean(differences) else NA_real_
    lodo_rows[[paste(signature_id, omitted_donor, sep = "||")]] <- data.frame(
      signature_id = signature_id, target_compartment = target,
      omitted_donor = omitted_donor, target_contributing_donors = length(target_donors),
      paired_difference_donors = length(differences),
      estimability = ifelse(estimable, "estimable", "not_estimable"),
      lodo_mean_paired_preference = effect,
      full_mean_paired_preference = selected$mean_paired_preference[row_index],
      same_direction = ifelse(estimable, sign(effect) == sign(selected$mean_paired_preference[row_index]), NA)
    )
  }
}
lodo <- do.call(rbind, lodo_rows)
write.table(lodo, file.path(audit, "results", "04_lodo_compartment_preference.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_rows <- lapply(split(lodo, lodo$signature_id), function(current) {
  estimable <- current[current$estimability == "estimable", , drop = FALSE]
  data.frame(
    signature_id = current$signature_id[1],
    target_compartment = current$target_compartment[1],
    full_mean_paired_preference = current$full_mean_paired_preference[1],
    estimable_iterations = nrow(estimable),
    total_iterations = nrow(current),
    same_direction_estimable_iterations = sum(estimable$same_direction, na.rm = TRUE),
    direction_consistency_denominator = nrow(estimable),
    donor_gate_pass = nrow(estimable) >= minimum_estimable_iterations && sum(estimable$same_direction, na.rm = TRUE) >= minimum_same_direction,
    cross_atlas_gate = "pending_GSE191288_and_GSE281736",
    main_text_status = "pending_cross_atlas_AND_rule"
  )
})
summary_table <- do.call(rbind, summary_rows)
write.table(summary_table, file.path(audit, "results", "04_lodo_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

heatmap_data <- aggregate(score ~ signature_id + compartment, full_scores, mean)
pdf(file.path(audit, "figures", "04_program_compartment_heatmap.pdf"), width = 9, height = 7, useDingbats = FALSE)
print(
  ggplot(heatmap_data, aes(compartment, signature_id, fill = score)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "Mean score")
)
dev.off()

pdf(file.path(audit, "figures", "04_lodo_compartment_preference.pdf"), width = 9, height = 7, useDingbats = FALSE)
print(
  ggplot(lodo, aes(lodo_mean_paired_preference, signature_id, color = estimability)) +
    geom_vline(xintercept = 0, linetype = 3, color = "grey50") +
    geom_point(position = position_jitter(height = 0.12), alpha = 0.75) +
    geom_point(
      data = summary_table,
      aes(full_mean_paired_preference, signature_id),
      inherit.aes = FALSE, shape = 23, size = 3, fill = "black"
    ) +
    theme_classic() +
    labs(x = "Target compartment minus within-donor other compartments", y = NULL, color = NULL)
)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(audit, "manifests", "04_sessionInfo.txt"))
writeLines("Stage 04 complete", file.path(audit, "checkpoints", "04_COMPLETE"))
log_line("Stage 04 complete; cross-atlas gate remains pending")
