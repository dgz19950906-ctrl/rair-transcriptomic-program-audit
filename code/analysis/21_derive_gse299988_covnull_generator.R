#!/usr/bin/env Rscript

# Derive the cohort-specific GSE299988 label-blind null generator from the
# independently frozen GSE151179 implementation using audited substitutions.

suppressPackageStartupMessages({ library(digest); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || args[[1]] != "--source" || args[[3]] != "--out") {
  stop("Usage: --source SOURCE_R --out OUTPUT_R")
}
source_path <- normalizePath(args[[2]])
out_path <- normalizePath(args[[4]], mustWork = FALSE)
if (file.exists(out_path)) stop("Refusing to overwrite derived generator")

x <- readLines(source_path, warn = FALSE)
original <- paste(x, collapse = "\n")
assert_count <- function(pattern, expected, fixed = TRUE) {
  observed <- lengths(regmatches(original, gregexpr(pattern, original, fixed = fixed)))
  if (observed != expected) stop(sprintf("Substitution precondition failed for %s: %d != %d", pattern, observed, expected))
}
assert_count("GSE151179_primary_preRAI_gene_expression.tsv.gz", 1L)
assert_count("if (ncol(expr_dt) != 18L", 1L)
assert_count("151179", 13L)
assert_count('cohort = "GSE151179"', 2L)

y <- gsub("GSE151179_primary_preRAI_gene_expression.tsv.gz",
          "GSE299988_tumor_gene_expression.tsv.gz", x, fixed = TRUE)
y <- gsub("if (ncol(expr_dt) != 18L", "if (ncol(expr_dt) != 11L", y, fixed = TRUE)
y <- gsub("151179", "299988", y, fixed = TRUE)

if (identical(x, y)) stop("No substitutions were made")
writeLines(y, out_path, useBytes = TRUE)
Sys.chmod(out_path, mode = "0755")

sha <- function(p) digest(p, algo = "sha256", file = TRUE, serialize = FALSE)
report <- list(
  source = source_path,
  source_sha256 = sha(source_path),
  output = normalizePath(out_path),
  output_sha256 = sha(out_path),
  substitutions = list(
    expression_basename = "GSE299988_tumor_gene_expression.tsv.gz",
    expression_schema_columns = 11,
    cohort = "GSE299988",
    seed_prefix = 299988
  ),
  clinical_labels_loaded = FALSE
)
write_json(report, paste0(out_path, ".derivation.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
