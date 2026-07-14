#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 64 * 1024^3)
set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(DoubletFinder)
  library(Matrix)
  library(future)
})
future::plan("sequential")

project <- "/home/dony/ThyroidCancer_Project"
audit <- file.path(project, "rair_audit")
checkpoint_dir <- file.path(audit, "checkpoints", "02_doubletfinder")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(audit, "results"), recursive = TRUE, showWarnings = FALSE)

log_line <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), ..., "\n", sep = "")
  flush.console()
}

get_df_function <- function(current, legacy) {
  exports <- getNamespaceExports("DoubletFinder")
  if (current %in% exports) return(getExportedValue("DoubletFinder", current))
  if (legacy %in% exports) return(getExportedValue("DoubletFinder", legacy))
  stop("DoubletFinder function unavailable: ", current, " / ", legacy)
}

param_sweep <- get_df_function("paramSweep", "paramSweep_v3")
summarize_sweep <- get_df_function("summarizeSweep", "summarizeSweep")
find_pk <- get_df_function("find.pK", "find.pK")
model_homotypic <- get_df_function("modelHomotypic", "modelHomotypic")
doublet_finder <- get_df_function("doubletFinder", "doubletFinder_v3")

raw_path <- file.path(project, "data", "GSE184362_merged_raw.rds")
log_line("Loading raw merged object: ", raw_path)
raw <- readRDS(raw_path)
DefaultAssay(raw) <- "RNA"
log_line("Joining Seurat v5 count layers")
raw <- JoinLayers(raw, assay = "RNA")

raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")
qc_keep <- with(
  raw@meta.data,
  nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10
)
qc_cells <- rownames(raw@meta.data)[qc_keep]
qc_flow <- data.frame(
  stage = c("raw", "author_qc"),
  cells = c(ncol(raw), length(qc_cells)),
  removed = c(NA_integer_, ncol(raw) - length(qc_cells)),
  rule = c("CreateSeuratObject_min_features_200_min_cells_3", "nFeature_RNA_gt500_lt5000_percent_mt_lt10")
)
write.table(
  qc_flow, file.path(audit, "results", "02_qc_flow.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
log_line("Author QC retained ", length(qc_cells), " / ", ncol(raw), " cells")
raw <- subset(raw, cells = qc_cells)
gc()

samples <- sort(unique(raw@meta.data$sample))
all_calls <- vector("list", length(samples))
all_summaries <- vector("list", length(samples))
names(all_calls) <- names(all_summaries) <- samples

for (index in seq_along(samples)) {
  sample_id <- samples[index]
  checkpoint <- file.path(checkpoint_dir, paste0(sample_id, "_doublet_calls.tsv"))
  summary_checkpoint <- file.path(checkpoint_dir, paste0(sample_id, "_summary.tsv"))

  if (file.exists(checkpoint) && file.exists(summary_checkpoint)) {
    log_line("[", index, "/", length(samples), "] Reusing checkpoint for ", sample_id)
    all_calls[[sample_id]] <- read.delim(checkpoint, check.names = FALSE)
    all_summaries[[sample_id]] <- read.delim(summary_checkpoint, check.names = FALSE)
    next
  }

  sample_cells <- rownames(raw@meta.data)[raw@meta.data$sample == sample_id]
  log_line("[", index, "/", length(samples), "] Processing ", sample_id, " (", length(sample_cells), " post-QC cells)")
  counts <- LayerData(raw, assay = "RNA", layer = "counts")[, sample_cells, drop = FALSE]
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 0, min.features = 0)
  obj <- AddMetaData(obj, raw@meta.data[sample_cells, , drop = FALSE])
  rm(counts)

  obj <- NormalizeData(obj, verbose = FALSE)
  feature_number <- min(10000L, nrow(obj) - 1L)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = feature_number, verbose = FALSE)
  variable_features <- VariableFeatures(obj)
  obj <- ScaleData(obj, features = variable_features, verbose = FALSE)
  npcs <- min(50L, length(variable_features) - 1L, ncol(obj) - 1L)
  if (npcs < 10L) stop("Too few PCs available for ", sample_id)
  obj <- RunPCA(obj, features = variable_features, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(npcs), verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)

  doublet_rate <- 0.075
  t10_igkc_filter <- FALSE
  if (grepl("PTC10_T$", sample_id)) {
    doublet_rate <- 0.20
    if ("IGKC" %in% rownames(obj)) {
      igkc <- FetchData(obj, vars = "IGKC", layer = "data")[[1]]
      retained <- names(igkc)[is.finite(igkc) & igkc < 1]
      obj <- subset(obj, cells = retained)
      t10_igkc_filter <- TRUE
      log_line(sample_id, ": applied author IGKC normalized-expression <1 filter; retained ", ncol(obj), " cells")
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(feature_number, nrow(obj) - 1L), verbose = FALSE)
      variable_features <- VariableFeatures(obj)
      obj <- ScaleData(obj, features = variable_features, verbose = FALSE)
      npcs <- min(50L, length(variable_features) - 1L, ncol(obj) - 1L)
      obj <- RunPCA(obj, features = variable_features, npcs = npcs, verbose = FALSE)
      obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(npcs), verbose = FALSE)
      obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
    }
  }

  pcs <- seq_len(npcs)
  log_line(sample_id, ": starting DoubletFinder parameter sweep")
  sweep <- param_sweep(obj, PCs = pcs, sct = FALSE)
  sweep_stats <- summarize_sweep(sweep, GT = FALSE)
  bcmvn <- find_pk(sweep_stats)
  valid_pk <- bcmvn[is.finite(bcmvn$BCmetric) & !is.na(bcmvn$pK), , drop = FALSE]
  if (!nrow(valid_pk)) stop("No valid DoubletFinder pK for ", sample_id)
  selected_pk <- as.numeric(as.character(valid_pk$pK[which.max(valid_pk$BCmetric)]))

  homotypic <- model_homotypic(obj@meta.data$seurat_clusters)
  expected_raw <- round(doublet_rate * ncol(obj))
  expected_adjusted <- max(1L, round(expected_raw * (1 - homotypic)))
  log_line(
    sample_id, ": pK=", selected_pk, ", expected_raw=", expected_raw,
    ", expected_adjusted=", expected_adjusted
  )
  obj <- doublet_finder(
    obj, PCs = pcs, pN = 0.25, pK = selected_pk,
    nExp = expected_adjusted, reuse.pANN = NULL, sct = FALSE
  )
  class_column <- tail(grep("^DF.classifications", colnames(obj@meta.data), value = TRUE), 1)
  if (!length(class_column)) stop("DoubletFinder classification column missing for ", sample_id)

  classification <- as.character(obj@meta.data[[class_column]])
  calls <- data.frame(
    cell = rownames(obj@meta.data),
    sample = sample_id,
    patient = obj@meta.data$patient,
    tissue = obj@meta.data$tissue,
    doublet_call = classification,
    selected_pK = selected_pk,
    expected_doublets_raw = expected_raw,
    expected_doublets_adjusted = expected_adjusted,
    homotypic_proportion = homotypic,
    author_doublet_rate = doublet_rate,
    t10_igkc_filter_applied = t10_igkc_filter
  )
  summary <- data.frame(
    sample = sample_id,
    patient = unique(obj@meta.data$patient)[1],
    tissue = unique(obj@meta.data$tissue)[1],
    cells_entering_doubletfinder = nrow(calls),
    singlets = sum(classification == "Singlet"),
    doublets = sum(classification == "Doublet"),
    selected_pK = selected_pk,
    expected_doublets_raw = expected_raw,
    expected_doublets_adjusted = expected_adjusted,
    homotypic_proportion = homotypic,
    author_doublet_rate = doublet_rate,
    t10_igkc_filter_applied = t10_igkc_filter
  )
  write.table(calls, checkpoint, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary, summary_checkpoint, sep = "\t", quote = FALSE, row.names = FALSE)
  all_calls[[sample_id]] <- calls
  all_summaries[[sample_id]] <- summary
  rm(obj, sweep, sweep_stats, bcmvn, valid_pk)
  gc()
}

calls <- do.call(rbind, all_calls)
summaries <- do.call(rbind, all_summaries)
write.table(calls, file.path(audit, "results", "02_doublet_calls.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summaries, file.path(audit, "results", "02_doublet_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

singlets <- calls$cell[calls$doublet_call == "Singlet"]
log_line("Subsetting joined raw object to ", length(singlets), " singlets")
singlet_object <- subset(raw, cells = singlets)
singlet_object@graphs <- list()
singlet_object@neighbors <- list()
singlet_object@reductions <- list()
VariableFeatures(singlet_object) <- character(0)
saveRDS(
  singlet_object,
  file.path(audit, "checkpoints", "02_author_qc_singlets.rds"),
  compress = FALSE
)

final_flow <- rbind(
  qc_flow,
  data.frame(
    stage = "author_qc_plus_doubletfinder",
    cells = length(singlets),
    removed = length(qc_cells) - length(singlets),
    rule = "per_sample_DoubletFinder_author_rates_and_T10_IGKC_rule"
  )
)
write.table(final_flow, file.path(audit, "results", "02_qc_flow.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(audit, "manifests", "02_sessionInfo.txt"))
log_line("Stage 02 complete: ", length(singlets), " singlets saved")
