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
definition_path <- file.path(audit, "manifests", "frozen_signatures.tsv")
gmt_path <- file.path(audit, "manifests", "h.all.v2025.1.Hs.symbols.gmt")
primary_summary_path <- file.path(audit, "results", "04_lodo_summary.tsv")
minimum_cells <- 20L
minimum_donors <- 3L

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), ..., "\n", sep = "")
  flush.console()
}

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
    defs <- rbind(defs, data.frame(
      signature_id = signature_id, gene = genes, direction = 1,
      orientation = "higher_more_aggressive_program",
      source = "MSigDB Hallmark v2025.1.Hs",
      role = "aggressiveness_negative_control", lock_status = "locked"
    ))
  }
}

score_programs <- function(expression, definition) {
  z <- t(scale(t(as.matrix(expression))))
  z[!is.finite(z)] <- 0
  output <- list()
  for (signature_id in unique(definition$signature_id)) {
    current <- definition[definition$signature_id == signature_id, , drop = FALSE]
    present <- current$gene[current$gene %in% rownames(z)]
    minimum_present <- max(3L, ceiling(0.7 * nrow(current)))
    if (length(present) < minimum_present) next
    weights <- current$direction[match(present, current$gene)]
    score <- colSums(z[present, , drop = FALSE] * weights) / sum(abs(weights))
    output[[signature_id]] <- data.frame(
      pseudobulk = names(score), signature_id = signature_id,
      score = unname(score), orientation = current$orientation[1], role = current$role[1]
    )
  }
  do.call(rbind, output)
}

paired_difference <- function(scores, metadata, target) {
  merged <- merge(scores, metadata, by = "pseudobulk", all.x = TRUE)
  donors <- sort(unique(merged$donor))
  output <- lapply(donors, function(donor) {
    current <- merged[merged$donor == donor, , drop = FALSE]
    target_score <- current$score[current$compartment == target]
    other_score <- current$score[current$compartment != target]
    difference <- if (length(target_score) == 1L && length(other_score)) target_score - mean(other_score) else NA_real_
    data.frame(donor = donor, difference = difference)
  })
  do.call(rbind, output)
}

prepare_atlas_scores <- function(dataset_id) {
  input <- file.path(audit, "checkpoints", paste0("06_", dataset_id, "_compartments.rds"))
  obj <- readRDS(input)
  DefaultAssay(obj) <- "RNA"
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  meta <- obj@meta.data[colnames(counts), , drop = FALSE]
  eligible <- rownames(meta)[
    meta$tissue == "Tumor" & !is.na(meta$compartment) & meta$compartment != "Uncertain"
  ]
  counts <- counts[, eligible, drop = FALSE]
  meta <- meta[eligible, , drop = FALSE]
  meta$group <- paste(meta$donor, meta$compartment, sep = "||")
  cell_counts <- as.data.frame(table(meta$donor, meta$compartment), stringsAsFactors = FALSE)
  names(cell_counts) <- c("donor", "compartment", "cells")
  cell_counts$dataset <- dataset_id
  cell_counts$qc_status <- ifelse(cell_counts$cells >= minimum_cells, "qualified", "insufficient_cells")
  qualified <- with(cell_counts[cell_counts$cells >= minimum_cells, ], paste(donor, compartment, sep = "||"))
  cells <- rownames(meta)[meta$group %in% qualified]
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  groups <- sort(unique(meta$group))
  indicator <- sparseMatrix(
    i = seq_len(nrow(meta)), j = match(meta$group, groups), x = 1,
    dims = c(nrow(meta), length(groups)), dimnames = list(rownames(meta), groups)
  )
  pb_counts <- counts[, rownames(indicator), drop = FALSE] %*% indicator
  logcpm <- log1p(t(t(pb_counts) / Matrix::colSums(pb_counts) * 1e6))
  split_group <- do.call(rbind, strsplit(colnames(logcpm), "\\|\\|"))
  pb_meta <- data.frame(
    pseudobulk = colnames(logcpm), donor = split_group[, 1],
    compartment = split_group[, 2], dataset = dataset_id,
    stringsAsFactors = FALSE
  )
  list(scores = score_programs(logcpm, defs), metadata = pb_meta, cell_counts = cell_counts)
}

primary <- read.delim(primary_summary_path, check.names = FALSE)
atlas_ids <- c("GSE191288", "GSE281736")
atlas_objects <- lapply(atlas_ids, prepare_atlas_scores)
names(atlas_objects) <- atlas_ids
write.table(
  do.call(rbind, lapply(atlas_objects, `[[`, "cell_counts")),
  file.path(audit, "results", "07_external_atlas_donor_compartment_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

atlas_rows <- list()
donor_rows <- list()
for (dataset_id in atlas_ids) {
  atlas <- atlas_objects[[dataset_id]]
  for (row_index in seq_len(nrow(primary))) {
    signature_id <- primary$signature_id[row_index]
    target <- primary$target_compartment[row_index]
    current <- atlas$scores[atlas$scores$signature_id == signature_id, , drop = FALSE]
    differences <- paired_difference(current[, c("pseudobulk", "score")], atlas$metadata, target)
    finite <- differences$difference[is.finite(differences$difference)]
    compartment_means <- merge(current, atlas$metadata, by = "pseudobulk")
    compartment_means <- aggregate(score ~ compartment, compartment_means, mean)
    top_compartment <- compartment_means$compartment[which.max(compartment_means$score)]
    estimable <- length(finite) >= minimum_donors
    mean_difference <- if (estimable) mean(finite) else NA_real_
    atlas_rows[[paste(dataset_id, signature_id, sep = "||")]] <- data.frame(
      dataset = dataset_id, signature_id = signature_id,
      frozen_target_compartment = target,
      paired_donors = length(finite), estimable = estimable,
      mean_target_vs_other_difference = mean_difference,
      direction = ifelse(estimable, sign(mean_difference), NA),
      top_compartment = top_compartment,
      target_is_top_compartment = top_compartment == target
    )
    differences$dataset <- dataset_id
    differences$signature_id <- signature_id
    differences$target_compartment <- target
    donor_rows[[paste(dataset_id, signature_id, sep = "||")]] <- differences
  }
}
atlas_table <- do.call(rbind, atlas_rows)
donor_table <- do.call(rbind, donor_rows)
write.table(atlas_table, file.path(audit, "results", "07_external_atlas_preferences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(donor_table, file.path(audit, "results", "07_external_atlas_donor_differences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

gate_rows <- lapply(seq_len(nrow(primary)), function(row_index) {
  signature_id <- primary$signature_id[row_index]
  current <- atlas_table[atlas_table$signature_id == signature_id, , drop = FALSE]
  primary_direction <- sign(primary$full_mean_paired_preference[row_index])
  all_estimable <- nrow(current) == 2L && all(current$estimable)
  all_same_direction <- all_estimable && all(current$direction == primary_direction)
  pass <- isTRUE(primary$donor_gate_pass[row_index]) && all_same_direction
  data.frame(
    signature_id = signature_id,
    target_compartment = primary$target_compartment[row_index],
    GSE184362_preference = primary$full_mean_paired_preference[row_index],
    GSE184362_donor_gate = primary$donor_gate_pass[row_index],
    GSE191288_preference = current$mean_target_vs_other_difference[current$dataset == "GSE191288"],
    GSE191288_estimable = current$estimable[current$dataset == "GSE191288"],
    GSE281736_preference = current$mean_target_vs_other_difference[current$dataset == "GSE281736"],
    GSE281736_estimable = current$estimable[current$dataset == "GSE281736"],
    three_atlas_same_direction = all_same_direction,
    three_atlas_AND_gate = pass,
    figure_status = ifelse(pass, "eligible_for_main_text", "move_cellular_layer_to_supplement")
  )
})
gate <- do.call(rbind, gate_rows)
write.table(gate, file.path(audit, "results", "07_three_atlas_AND_gate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

plot_data <- rbind(
  data.frame(dataset = "GSE184362", signature_id = primary$signature_id, preference = primary$full_mean_paired_preference),
  data.frame(dataset = atlas_table$dataset, signature_id = atlas_table$signature_id, preference = atlas_table$mean_target_vs_other_difference)
)
pdf(file.path(audit, "figures", "07_three_atlas_direction_map.pdf"), width = 9, height = 7, useDingbats = FALSE)
print(
  ggplot(plot_data, aes(dataset, signature_id, fill = preference)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey80") +
    theme_classic() + labs(x = NULL, y = NULL, fill = "Target preference")
)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(audit, "manifests", "07_sessionInfo.txt"))
writeLines("Stage 07 complete", file.path(audit, "checkpoints", "07_COMPLETE"))
log_line("Stage 07 complete")
