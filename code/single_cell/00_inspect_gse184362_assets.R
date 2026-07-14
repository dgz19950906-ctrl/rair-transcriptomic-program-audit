options(stringsAsFactors = FALSE)

project <- "/home/dony/ThyroidCancer_Project"

cat("=== R environment ===\n")
cat("R:", R.version.string, "\n")
packages <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "ggplot2", "data.table", "future")
versions <- vapply(
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
print(versions)

cat("\n=== Sample map ===\n")
sample_map <- readRDS(file.path(project, "data/GSE184362_sample_map.rds"))
cat("Entries:", length(sample_map), "\n")
str(sample_map, max.level = 2)

cat("\n=== Compact metadata ===\n")
compact_meta <- readRDS(file.path(project, "data/GSE184362_meta.rds"))
cat("Class:", paste(class(compact_meta), collapse = ", "), "\n")
str(compact_meta, max.level = 2)

cat("\n=== Existing small annotated object ===\n")
small_path <- file.path(project, "data/phase3_annotated_final.rds")
small <- readRDS(small_path)
cat("Class:", paste(class(small), collapse = ", "), "\n")
cat("Dimensions:", paste(dim(small), collapse = " x "), "\n")
cat("Metadata columns:\n")
print(colnames(small@meta.data))
for (field in intersect(c("sample", "patient", "tissue", "celltype", "seurat_clusters"), colnames(small@meta.data))) {
  cat("\n--", field, "--\n")
  print(table(small@meta.data[[field]], useNA = "ifany"))
}

cat("\n=== Large object file ===\n")
large_path <- file.path(project, "data/GSE184362_processed.rds")
info <- file.info(large_path)
print(info[, c("size", "mtime")])
cat("Large object was not loaded by this inspection script.\n")
