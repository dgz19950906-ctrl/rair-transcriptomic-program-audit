options(stringsAsFactors = FALSE)

project <- "/home/dony/ThyroidCancer_Project"
audit <- file.path(project, "rair_audit")
dirs <- file.path(audit, c("scripts", "logs", "checkpoints", "results", "figures", "manifests"))
for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)

cat("=== OPTION A PREFLIGHT ===\n")
cat("Timestamp:", format(Sys.time(), tz = "UTC"), "UTC\n")
cat("R:", R.version.string, "\n")

packages <- c(
  "Seurat", "SeuratObject", "Matrix", "DoubletFinder", "batchelor",
  "SingleCellExperiment", "scater", "scran", "scuttle", "BiocParallel",
  "SingleR", "celldex", "dplyr", "data.table", "ggplot2", "future"
)
package_table <- data.frame(
  package = packages,
  installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    packages,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    },
    character(1)
  )
)
print(package_table, row.names = FALSE)
write.table(
  package_table,
  file.path(audit, "manifests/package_preflight.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

assets <- c(
  merged_raw = file.path(project, "data/GSE184362_merged_raw.rds"),
  processed = file.path(project, "data/GSE184362_processed.rds"),
  sample_map = file.path(project, "data/GSE184362_sample_map.rds")
)
asset_table <- data.frame(
  asset = names(assets),
  path = unname(assets),
  exists = file.exists(assets),
  bytes = unname(file.info(assets)$size),
  modified = format(file.info(assets)$mtime, tz = "UTC")
)
print(asset_table, row.names = FALSE)
write.table(
  asset_table,
  file.path(audit, "manifests/asset_preflight.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (!all(asset_table$exists)) stop("Required server-side assets are missing.")

cat("\nLoading server-side raw merged object...\n")
raw <- readRDS(assets[["merged_raw"]])
cat("Class:", paste(class(raw), collapse = ", "), "\n")
cat("Dimensions:", nrow(raw), "genes x", ncol(raw), "cells\n")
cat("Metadata columns:\n")
print(colnames(raw@meta.data))
cat("Sample count:", length(unique(raw@meta.data$sample)), "\n")
cat("Donor count:", length(unique(raw@meta.data$patient)), "\n")
print(table(raw@meta.data$tissue, useNA = "ifany"))
print(table(raw@meta.data$patient, useNA = "ifany"))

object_summary <- data.frame(
  genes = nrow(raw),
  cells = ncol(raw),
  samples = length(unique(raw@meta.data$sample)),
  donors = length(unique(raw@meta.data$patient)),
  has_counts = "counts" %in% SeuratObject::Layers(raw[["RNA"]]),
  has_data = "data" %in% SeuratObject::Layers(raw[["RNA"]])
)
write.table(
  object_summary,
  file.path(audit, "manifests/raw_object_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(audit, "manifests/preflight_sessionInfo.txt"))
cat("Preflight complete. The 13 GB processed object was not loaded.\n")
