#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --name value pairs")
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("primary", "gmt", "out")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")
if (file.exists(args$out)) stop("Refusing to overwrite an existing registry")

primary <- fread(args$primary)
required_cols <- c("signature_id", "gene", "direction", "orientation", "source", "role", "lock_status")
if (length(setdiff(required_cols, names(primary)))) stop("Primary registry schema mismatch")

control_ids <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

gmt_lines <- strsplit(readLines(args$gmt, warn = FALSE), "\t", fixed = TRUE)
names(gmt_lines) <- vapply(gmt_lines, `[[`, character(1), 1L)
if (length(setdiff(control_ids, names(gmt_lines)))) stop("A frozen Hallmark control is absent from the GMT file")

controls <- rbindlist(lapply(control_ids, function(program) {
  genes <- unique(gmt_lines[[program]][-(1:2)])
  data.table(
    signature_id = program,
    gene = genes,
    direction = 1L,
    orientation = "higher_more_aggressive_program",
    source = "MSigDB Hallmark v2025.1.Hs",
    role = "aggressiveness_negative_control",
    lock_status = "locked"
  )
}))

registry <- rbindlist(list(primary[, ..required_cols], controls), use.names = TRUE)
if (uniqueN(registry$signature_id) != 9L) stop("Expected nine frozen programs")
if (registry[, anyDuplicated(gene), by = signature_id][, any(V1 > 0L)]) stop("Duplicate genes within a program")
if (any(!registry$direction %in% c(-1L, 1L)) || any(registry$lock_status != "locked")) stop("Registry lock failure")

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
fwrite(registry, args$out, sep = "\t")

summary <- registry[, .(
  genes_requested = .N,
  positive = sum(direction == 1L),
  negative = sum(direction == -1L),
  source = source[[1]],
  role = role[[1]]
), by = signature_id]
fwrite(summary, paste0(args$out, ".summary.tsv"), sep = "\t")

manifest <- list(
  registry = normalizePath(args$out),
  registry_sha256 = digest(args$out, algo = "sha256", file = TRUE, serialize = FALSE),
  primary_registry_sha256 = digest(args$primary, algo = "sha256", file = TRUE, serialize = FALSE),
  hallmark_gmt_sha256 = digest(args$gmt, algo = "sha256", file = TRUE, serialize = FALSE),
  programs = 9L,
  rows = nrow(registry),
  clinical_labels_loaded = FALSE
)
write_json(manifest, paste0(args$out, ".manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), "\n")
