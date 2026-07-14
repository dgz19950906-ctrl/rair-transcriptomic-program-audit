#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 64 * 1024^3)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(harmony)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(future)
})
future::plan("sequential")

project <- "/home/dony/ThyroidCancer_Project"
audit <- file.path(project, "rair_audit")
external <- file.path(audit, "external_data")

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), ..., "\n", sep = "")
  flush.console()
}

author_markers <- list(
  `T & NK cells` = c("IL7R", "CCR7", "CD3D", "CD3E", "CD3G", "GZMK", "GZMA", "CD8A", "NKG7"),
  `B cells` = c("CD79A", "CD74", "IGHM", "MS4A1", "CD19", "JCHAIN"),
  `Thyroid cells` = c("TG", "TSHR", "TPO", "SLC5A5", "EPCAM", "KRT19", "CITED1", "FOXE1", "PAX8"),
  `Myeloid cells` = c("LYZ", "CD14", "CD68", "FCGR3A", "C1QA", "C1QB", "LST1"),
  Fibroblasts = c("COL1A1", "COL1A2", "DCN", "SPARC", "ACTA2", "TAGLN", "LUM"),
  `Endothelial cells` = c("PECAM1", "VWF", "PLVAP", "SPARCL1", "MMRN1", "LYVE1", "KDR")
)

single_r_broad <- function(label) {
  if (grepl("B_cell|Plasma", label, ignore.case = TRUE)) return("B cells")
  if (grepl("T_cell|NK_cell", label, ignore.case = TRUE)) return("T & NK cells")
  if (grepl("Monocyte|Macrophage|DC|Neutrophil", label, ignore.case = TRUE)) return("Myeloid cells")
  if (grepl("Fibroblast|Smooth_muscle|Mesenchymal", label, ignore.case = TRUE)) return("Fibroblasts")
  if (grepl("Endothelial", label, ignore.case = TRUE)) return("Endothelial cells")
  if (grepl("Epithelial|Keratinocyte|Tissue_stem", label, ignore.case = TRUE)) return("Thyroid cells")
  "Unmapped"
}

read_gse191288 <- function() {
  raw_dir <- file.path(external, "GSE191288", "raw")
  files <- sort(list.files(raw_dir, pattern = "\\.h5$", full.names = TRUE))
  map <- data.frame(
    gsm = paste0("GSM574302", 1:7),
    donor = c("P1", "P1", "P2", "P2", "P3", "P3", "P4"),
    tissue = c(rep("Tumor", 6), "Normal"),
    site = c("Left", "Right", "Left", "Right", "Left", "Right", "Normal"),
    stringsAsFactors = FALSE
  )
  objects <- list()
  for (file in files) {
    gsm <- sub("_.*", "", basename(file))
    info <- map[map$gsm == gsm, , drop = FALSE]
    if (nrow(info) != 1L) stop("Unmapped GSE191288 file: ", basename(file))
    counts <- Read10X_h5(file, use.names = TRUE, unique.features = TRUE)
    if (is.list(counts) && !inherits(counts, "dgCMatrix")) {
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }
    obj <- CreateSeuratObject(counts, project = gsm, min.cells = 3, min.features = 200)
    obj$sample <- gsm
    obj$donor <- info$donor
    obj$tissue <- info$tissue
    obj$site <- info$site
    objects[[gsm]] <- obj
    log_line("GSE191288 ", gsm, ": ", ncol(obj), " cells loaded")
  }
  objects
}

read_gse281736 <- function() {
  raw_dir <- file.path(external, "GSE281736", "raw")
  matrix_files <- sort(list.files(raw_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
  map <- data.frame(
    gsm = paste0("GSM86272", sprintf("%02d", 9:20)),
    donor = rep(c("ADULT1", paste0("CAYA", 1:5)), each = 2),
    tissue = rep(c("Tumor", "Paratumor"), 6),
    age_group = rep(c("Adult", rep("CAYA", 5)), each = 2),
    stringsAsFactors = FALSE
  )
  objects <- list()
  for (matrix_file in matrix_files) {
    gsm <- sub("_.*", "", basename(matrix_file))
    prefix <- sub("_matrix\\.mtx\\.gz$", "", matrix_file)
    barcode_file <- paste0(prefix, "_barcodes.tsv.gz")
    gene_file <- paste0(prefix, "_genes.tsv.gz")
    info <- map[map$gsm == gsm, , drop = FALSE]
    if (nrow(info) != 1L) stop("Unmapped GSE281736 file: ", basename(matrix_file))
    gene_preview <- read.delim(gene_file, header = FALSE, nrows = 5)
    feature_column <- if (ncol(gene_preview) >= 2L) 2L else 1L
    counts <- ReadMtx(
      mtx = matrix_file, cells = barcode_file, features = gene_file,
      feature.column = feature_column, unique.features = TRUE
    )
    obj <- CreateSeuratObject(counts, project = gsm, min.cells = 3, min.features = 200)
    obj$sample <- gsm
    obj$donor <- info$donor
    obj$tissue <- info$tissue
    obj$age_group <- info$age_group
    objects[[gsm]] <- obj
    log_line("GSE281736 ", gsm, ": ", ncol(obj), " cells loaded")
  }
  objects
}

annotate_dataset <- function(objects, dataset_id) {
  log_line(dataset_id, ": merging ", length(objects), " samples")
  obj <- merge(objects[[1]], y = objects[-1], project = dataset_id)
  obj <- JoinLayers(obj, assay = "RNA")
  raw_cells <- ncol(obj)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  keep <- with(obj@meta.data, nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
  obj <- subset(obj, cells = rownames(obj@meta.data)[keep])
  log_line(dataset_id, ": author-comparable QC retained ", ncol(obj), " / ", raw_cells, " cells")

  qc_by_sample <- aggregate(
    x = list(post_qc_cells = rep(1L, ncol(obj))),
    by = list(sample = obj$sample, donor = obj$donor, tissue = obj$tissue),
    FUN = sum
  )
  write.table(qc_by_sample, file.path(audit, "results", paste0("06_", dataset_id, "_qc_by_sample.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)
  selected <- setdiff(VariableFeatures(obj), grep("^(MT-|RPL|RPS)", VariableFeatures(obj), value = TRUE))
  obj <- ScaleData(obj, features = selected, verbose = FALSE)
  obj <- RunPCA(obj, features = selected, npcs = 30, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", theta = 2, max.iter.harmony = 20, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, n.neighbors = 30, min.dist = 0.3, seed.use = 42, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1, random.seed = 42, verbose = FALSE)
  log_line(dataset_id, ": ", length(unique(Idents(obj))), " clusters")

  cluster <- as.character(Idents(obj))
  cluster_levels <- sort(unique(cluster), method = "radix")
  indicator <- sparseMatrix(
    i = seq_along(cluster), j = match(cluster, cluster_levels), x = 1,
    dims = c(length(cluster), length(cluster_levels)),
    dimnames = list(colnames(obj), cluster_levels)
  )
  counts <- LayerData(obj, assay = "RNA", layer = "counts")[, rownames(indicator), drop = FALSE]
  cluster_counts <- counts %*% indicator
  cluster_logcpm <- log1p(t(t(cluster_counts) / Matrix::colSums(cluster_counts) * 1e6))
  gene_z <- t(scale(t(as.matrix(cluster_logcpm))))
  gene_z[!is.finite(gene_z)] <- 0
  marker_scores <- vapply(author_markers, function(genes) {
    present <- intersect(genes, rownames(gene_z))
    if (length(present) < 2L) return(rep(NA_real_, ncol(gene_z)))
    colMeans(gene_z[present, , drop = FALSE])
  }, numeric(ncol(gene_z)))
  marker_scores <- t(marker_scores)
  colnames(marker_scores) <- colnames(gene_z)

  hpca <- HumanPrimaryCellAtlasData()
  single_r <- SingleR(test = cluster_logcpm, ref = hpca, labels = hpca$label.main, de.method = "classic")
  annotation <- data.frame(
    cluster = cluster_levels,
    cells = as.integer(table(factor(cluster, levels = cluster_levels))),
    author_marker_top = apply(marker_scores[, cluster_levels, drop = FALSE], 2, function(x) rownames(marker_scores)[which.max(x)]),
    author_marker_top_score = apply(marker_scores[, cluster_levels, drop = FALSE], 2, max, na.rm = TRUE),
    author_marker_margin = apply(marker_scores[, cluster_levels, drop = FALSE], 2, function(x) -diff(sort(x, decreasing = TRUE)[1:2])),
    singleR_label = single_r$labels[match(cluster_levels, rownames(single_r))],
    stringsAsFactors = FALSE
  )
  annotation$singleR_broad <- vapply(annotation$singleR_label, single_r_broad, character(1))
  annotation$annotation_agreement <- annotation$author_marker_top == annotation$singleR_broad
  annotation$final_compartment <- ifelse(
    annotation$annotation_agreement, annotation$author_marker_top,
    ifelse(annotation$author_marker_margin >= 0.5, annotation$author_marker_top, "Uncertain")
  )
  annotation$certainty <- ifelse(
    annotation$annotation_agreement, "author_markers_plus_SingleR_agree",
    ifelse(annotation$author_marker_margin >= 0.5, "author_marker_dominant", "uncertain_discordant")
  )
  write.table(annotation, file.path(audit, "results", paste0("06_", dataset_id, "_cluster_annotation.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

  cluster_map <- setNames(annotation$final_compartment, annotation$cluster)
  certainty_map <- setNames(annotation$certainty, annotation$cluster)
  obj$reconstructed_cluster <- as.character(Idents(obj))
  obj$compartment <- unname(cluster_map[obj$reconstructed_cluster])
  obj$compartment_certainty <- unname(certainty_map[obj$reconstructed_cluster])

  umap <- Embeddings(obj, "umap")
  plot_data <- data.frame(UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], obj@meta.data)
  pdf(file.path(audit, "figures", paste0("06_", dataset_id, "_compartments.pdf")), width = 11, height = 5.5, useDingbats = FALSE)
  print(
    ggplot(plot_data, aes(UMAP_1, UMAP_2, color = compartment)) +
      geom_point(size = 0.08, alpha = 0.4) + facet_wrap(~tissue) + theme_classic() +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) + labs(color = NULL)
  )
  dev.off()

  obj@graphs <- list()
  obj@neighbors <- list()
  obj <- DietSeurat(obj, assays = "RNA", dimreducs = c("harmony", "umap"), layers = "counts", misc = TRUE)
  out <- file.path(audit, "checkpoints", paste0("06_", dataset_id, "_compartments.rds"))
  saveRDS(obj, out, compress = FALSE)
  log_line(dataset_id, ": saved ", out)
  rm(obj, counts, cluster_counts, cluster_logcpm, gene_z, marker_scores)
  gc()
}

log_line("Stage 06 starting GSE191288")
annotate_dataset(read_gse191288(), "GSE191288")
log_line("Stage 06 starting GSE281736")
annotate_dataset(read_gse281736(), "GSE281736")
writeLines(capture.output(sessionInfo()), file.path(audit, "manifests", "06_sessionInfo.txt"))
writeLines("Stage 06 complete", file.path(audit, "checkpoints", "06_COMPLETE"))
log_line("Stage 06 complete")
