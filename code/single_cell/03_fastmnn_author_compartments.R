#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 96 * 1024^3)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(scuttle)
  library(scran)
  library(batchelor)
  library(BiocParallel)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(future)
})
future::plan("sequential")

project <- "/home/dony/ThyroidCancer_Project"
audit <- file.path(project, "rair_audit")
input_path <- file.path(audit, "checkpoints", "02_author_qc_singlets.rds")
output_path <- file.path(audit, "checkpoints", "03_mnn_clustered_counts.rds")

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), ..., "\n", sep = "")
  flush.console()
}

if (!file.exists(input_path)) stop("Stage 02 checkpoint is missing: ", input_path)
log_line("Loading Stage 02 singlet object")
obj <- readRDS(input_path)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]])) != 1L || !"counts" %in% Layers(obj[["RNA"]])) {
  obj <- JoinLayers(obj, assay = "RNA")
}
counts <- LayerData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data[colnames(counts), , drop = FALSE]
log_line("Loaded ", ncol(counts), " cells across ", length(unique(meta$sample)), " samples")

sample_to_short <- c(
  GSM5585102_PTC1_T = "T1", GSM5585103_PTC1_P = "P1",
  GSM5585104_PTC2_T = "T2", GSM5585105_PTC2_P = "P2",
  GSM5585106_PTC2_LeftLN = "LN2l", GSM5585107_PTC3_T = "T3",
  GSM5585108_PTC3_P = "P3", GSM5585109_PTC3_LeftLN = "LN3l",
  GSM5585110_PTC3_RightLN = "LN3r", GSM5585111_PTC4_SC = "SC4",
  GSM5585112_PTC5_T = "T5", GSM5585113_PTC5_P = "P5",
  GSM5585114_PTC5_RightLN = "LN5r", GSM5585115_PTC6_RightLN = "LN6r",
  GSM5585116_PTC7_RightLN = "LN7r", GSM5585117_PTC8_T = "T8",
  GSM5585118_PTC8_P = "P8", GSM5585119_PTC9_T = "T9",
  GSM5585120_PTC9_P = "P9", GSM5585121_PTC10_T = "T10",
  GSM5585122_PTC10_RightLN = "LN10r", GSM5585123_PTC11_RightLN = "LN11r",
  GSM5585124_PTC11_SC = "SC11"
)
if (!setequal(unique(meta$sample), names(sample_to_short))) {
  stop("Server sample identifiers differ from the frozen 23-sample map.")
}
meta$author_short_id <- unname(sample_to_short[meta$sample])

author_batch_order <- c(
  "LN10r", "LN11r", "LN2l", "LN3l", "LN3r", "LN5r", "LN6r", "LN7r",
  "P1", "P2", "P3", "P5", "P8", "P9", "SC11", "SC4", "T10", "T1",
  "T2", "T3", "T5", "T8", "T9"
)

log_line("Building sample-level SingleCellExperiment objects")
sce_list <- vector("list", length(author_batch_order))
names(sce_list) <- author_batch_order
variance_list <- vector("list", length(author_batch_order))
names(variance_list) <- author_batch_order
for (short_id in author_batch_order) {
  cells <- rownames(meta)[meta$author_short_id == short_id]
  sample_counts <- counts[, cells, drop = FALSE]
  sce <- SingleCellExperiment(
    assays = list(counts = sample_counts),
    colData = DataFrame(meta[cells, , drop = FALSE])
  )
  sce <- logNormCounts(sce)
  variance_list[[short_id]] <- modelGeneVar(sce)
  sce_list[[short_id]] <- sce
  log_line(short_id, ": ", ncol(sce), " cells")
}

combined_variance <- combineVar(variance_list, equiweight = TRUE)
combined_variance$gene <- rownames(combined_variance)
combined_variance <- combined_variance[
  !grepl("^(MT-|RPL|RPS)", combined_variance$gene) & is.finite(combined_variance$bio),
  , drop = FALSE
]
combined_variance <- combined_variance[order(combined_variance$bio, decreasing = TRUE), , drop = FALSE]
hvg <- head(combined_variance$gene, 5000L)
write.table(
  combined_variance,
  file.path(audit, "results", "03_combined_gene_variance.tsv.gz"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(hvg, file.path(audit, "manifests", "03_fastmnn_hvg_5000.txt"))
log_line("Selected ", length(hvg), " label-blind HVGs")

log_line("Starting batchelor::fastMNN in the author-declared batch order")
mnn_args <- c(
  sce_list,
  list(
    subset.row = hvg,
    d = 50,
    k = 20,
    correct.all = FALSE,
    BPPARAM = SerialParam()
  )
)
mnn <- do.call(fastMNN, mnn_args)
corrected <- reducedDim(mnn, "corrected")
cell_order <- colnames(mnn)
rownames(corrected) <- cell_order
log_line("fastMNN complete: ", nrow(corrected), " cells x ", ncol(corrected), " dimensions")

obj <- obj[, cell_order]
obj <- AddMetaData(obj, meta[cell_order, , drop = FALSE])
obj[["mnn"]] <- CreateDimReducObject(embeddings = corrected, key = "MNN_", assay = "RNA")
obj <- RunUMAP(obj, reduction = "mnn", dims = 1:20, n.neighbors = 30, min.dist = 0.3, seed.use = 42, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "mnn", dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 1, random.seed = 42, verbose = FALSE)
log_line("Clustering complete: ", length(unique(Idents(obj))), " clusters")

cluster <- as.character(Idents(obj))
cluster_levels <- sort(unique(cluster), method = "radix")
indicator <- sparseMatrix(
  i = seq_along(cluster),
  j = match(cluster, cluster_levels),
  x = 1,
  dims = c(length(cluster), length(cluster_levels)),
  dimnames = list(colnames(obj), cluster_levels)
)
ordered_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, rownames(indicator), drop = FALSE]
cluster_counts <- ordered_counts %*% indicator
cluster_library <- Matrix::colSums(cluster_counts)
cluster_logcpm <- log1p(t(t(cluster_counts) / cluster_library * 1e6))

author_markers <- list(
  `T & NK cells` = c("IL7R", "CCR7", "CD3D", "CD3E", "CD3G", "GZMK", "GZMA", "CD8A", "NKG7"),
  `B cells` = c("CD79A", "CD74", "IGHM", "MS4A1", "CD19", "JCHAIN"),
  `Thyroid cells` = c("TG", "TSHR", "TPO", "SLC5A5", "EPCAM", "KRT19", "CITED1", "FOXE1", "PAX8"),
  `Myeloid cells` = c("LYZ", "CD14", "CD68", "FCGR3A", "C1QA", "C1QB", "LST1"),
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "SPARC", "ACTA2", "TAGLN", "LUM"),
  `Endothelial cells` = c("PECAM1", "VWF", "PLVAP", "SPARCL1", "MMRN1", "LYVE1", "KDR")
)

gene_z <- t(scale(t(as.matrix(cluster_logcpm))))
gene_z[!is.finite(gene_z)] <- 0
marker_scores <- vapply(
  author_markers,
  function(genes) {
    present <- intersect(genes, rownames(gene_z))
    if (length(present) < 2L) return(rep(NA_real_, ncol(gene_z)))
    colMeans(gene_z[present, , drop = FALSE])
  },
  numeric(ncol(gene_z))
)
marker_scores <- t(marker_scores)
colnames(marker_scores) <- colnames(gene_z)

log_line("Running cluster-level SingleR as an orthogonal annotation check")
hpca <- HumanPrimaryCellAtlasData()
single_r <- SingleR(
  test = cluster_logcpm,
  ref = hpca,
  labels = hpca$label.main,
  de.method = "classic"
)

single_r_broad <- function(label) {
  if (grepl("B_cell|Plasma", label, ignore.case = TRUE)) return("B cells")
  if (grepl("T_cell|NK_cell", label, ignore.case = TRUE)) return("T & NK cells")
  if (grepl("Monocyte|Macrophage|DC|Neutrophil", label, ignore.case = TRUE)) return("Myeloid cells")
  if (grepl("Fibroblast|Smooth_muscle|Mesenchymal", label, ignore.case = TRUE)) return("Fibroblasts")
  if (grepl("Endothelial", label, ignore.case = TRUE)) return("Endothelial cells")
  if (grepl("Epithelial|Keratinocyte|Tissue_stem", label, ignore.case = TRUE)) return("Thyroid cells")
  return("Unmapped")
}

annotation <- data.frame(
  cluster = cluster_levels,
  cells = as.integer(table(factor(cluster, levels = cluster_levels))),
  author_marker_top = apply(marker_scores[, cluster_levels, drop = FALSE], 2, function(x) rownames(marker_scores)[which.max(x)]),
  author_marker_top_score = apply(marker_scores[, cluster_levels, drop = FALSE], 2, max, na.rm = TRUE),
  author_marker_margin = apply(marker_scores[, cluster_levels, drop = FALSE], 2, function(x) diff(sort(x, decreasing = TRUE)[1:2]) * -1),
  singleR_label = single_r$labels[match(cluster_levels, rownames(single_r))],
  stringsAsFactors = FALSE
)
annotation$singleR_broad <- vapply(annotation$singleR_label, single_r_broad, character(1))
annotation$annotation_agreement <- annotation$author_marker_top == annotation$singleR_broad
annotation$final_compartment <- ifelse(
  annotation$annotation_agreement,
  annotation$author_marker_top,
  ifelse(annotation$author_marker_margin >= 0.5, annotation$author_marker_top, "Uncertain")
)
annotation$certainty <- ifelse(
  annotation$annotation_agreement,
  "author_markers_plus_SingleR_agree",
  ifelse(annotation$author_marker_margin >= 0.5, "author_marker_dominant", "uncertain_discordant")
)

score_table <- data.frame(
  cluster = rep(colnames(marker_scores), each = nrow(marker_scores)),
  compartment = rep(rownames(marker_scores), times = ncol(marker_scores)),
  marker_score = as.vector(marker_scores)
)
write.table(annotation, file.path(audit, "results", "03_cluster_compartment_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(score_table, file.path(audit, "results", "03_cluster_author_marker_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cluster_map <- setNames(annotation$final_compartment, annotation$cluster)
certainty_map <- setNames(annotation$certainty, annotation$cluster)
obj$author_reconstructed_cluster <- as.character(Idents(obj))
obj$compartment <- unname(cluster_map[obj$author_reconstructed_cluster])
obj$compartment_certainty <- unname(certainty_map[obj$author_reconstructed_cluster])

cell_metadata <- obj@meta.data
cell_metadata$cell <- rownames(cell_metadata)
write.table(
  cell_metadata,
  gzfile(file.path(audit, "results", "03_cell_metadata.tsv.gz"), "wt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

umap <- Embeddings(obj, "umap")
umap_df <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], obj@meta.data)
pdf(file.path(audit, "figures", "03_author_reconstructed_compartments.pdf"), width = 12, height = 5.5, useDingbats = FALSE)
print(
  ggplot(umap_df, aes(UMAP_1, UMAP_2, color = compartment)) +
    geom_point(size = 0.05, alpha = 0.35) +
    facet_wrap(~tissue) +
    theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    labs(color = NULL)
)
dev.off()

obj@graphs <- list()
obj@neighbors <- list()
obj <- DietSeurat(obj, assays = "RNA", dimreducs = c("mnn", "umap"), layers = "counts", misc = TRUE)
saveRDS(obj, output_path, compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(audit, "manifests", "03_sessionInfo.txt"))
log_line("Stage 03 complete: ", output_path)
