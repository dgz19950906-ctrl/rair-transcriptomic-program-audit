#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table); library(digest); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1]] != "--base") stop("Usage: --base DIR")
base <- normalizePath(args[[2]])
out <- file.path(base, "manifests")
master_path <- file.path(out, "MASTER_FROZEN_SHA256.tsv")
if (file.exists(master_path)) stop("Refusing to overwrite master freeze")

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
programs <- c(
  "TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
)

rows <- lapply(programs, function(program) {
  result <- file.path(base, "results", "GSE299988", program, "tolerance_0.05")
  validation <- file.path(base, "validation", program, "tolerance_0.05")
  tier <- fread(file.path(result, "tier_manifest.tsv"))
  cov <- fread(file.path(result, "covariance_diagnostic.tsv"))
  chain <- fread(file.path(result, "chain_diagnostics.tsv"))
  jac <- fread(file.path(result, "jaccard_summary.tsv"))
  label <- read_json(file.path(result, "label_blind_audit.json"), simplifyVector = TRUE)
  valid <- read_json(file.path(validation, "validation_summary.json"), simplifyVector = TRUE)
  frozen <- fread(file.path(result, "FROZEN_SHA256.tsv"))
  result_hash_ok <- all(vapply(file.path(result, frozen$file), sha256_file, character(1)) == frozen$sha256)
  validation_hash <- fread(file.path(validation, "VALIDATION_SHA256.tsv"))
  validation_hash_ok <- all(vapply(validation_hash$file, sha256_file, character(1)) == validation_hash$sha256)
  pc <- cov[feature == "pc1_variance_proportion"]
  mr <- cov[feature == "mean_pairwise_pearson_r"]
  data.table(
    program = program,
    status = tier$status[[1]],
    genes_requested = tier$genes_requested[[1]],
    genes_present = tier$genes_present[[1]],
    positive_directions = tier$positive_directions[[1]],
    negative_directions = tier$negative_directions[[1]],
    universe_genes = tier$universe_genes[[1]],
    chains_feasible = tier$chains_feasible[[1]],
    chains_contributing = tier$chains_contributing[[1]],
    unique_sets = tier$unique_sets[[1]],
    selected_sets = tier$selected_sets[[1]],
    target_pc1 = pc$target[[1]],
    target_mean_r = mr$target[[1]],
    maximum_absolute_pc1_deviation = pc$maximum_absolute_deviation[[1]],
    maximum_absolute_mean_r_deviation = mr$maximum_absolute_deviation[[1]],
    median_within_region_acceptance = median(chain$within_region_acceptance_rate, na.rm = TRUE),
    jaccard_median = jac$median[[1]],
    jaccard_maximum = jac$maximum[[1]],
    label_blind_passed = isTRUE(label$passed) && !isTRUE(label$clinical_label_files_loaded) &&
      !isTRUE(label$endpoint_scoring_performed) && label$n_expression_samples == 10L,
    independent_validation_passed = isTRUE(valid$passed),
    result_hashes_passed = result_hash_ok,
    validation_hashes_passed = validation_hash_ok,
    generator_script_sha256 = tier$script_sha256[[1]],
    selected_gene_manifest_sha256 = tier$gene_set_manifest_sha256[[1]]
  )
})

summary <- rbindlist(rows)
generator <- file.path(base, "scripts", "21_covariance_null_generator_gse299988.R")
expected_generator_hash <- sha256_file(generator)
all_pass <- nrow(summary) == 9L && all(summary$status == "estimable_at_0.05") &&
  all(summary$selected_sets == 1000L) && all(summary$chains_contributing >= 5L) &&
  all(summary$maximum_absolute_pc1_deviation <= 0.05 + 1e-12) &&
  all(summary$maximum_absolute_mean_r_deviation <= 0.05 + 1e-12) &&
  all(summary$label_blind_passed) && all(summary$independent_validation_passed) &&
  all(summary$result_hashes_passed) && all(summary$validation_hashes_passed) &&
  all(summary$generator_script_sha256 == expected_generator_hash)

summary_path <- file.path(out, "all_programs_tolerance_0.05_summary.tsv")
audit_path <- file.path(out, "master_audit_summary.json")
fwrite(summary, summary_path, sep = "\t")
audit <- list(
  passed = all_pass,
  cohort = "GSE299988",
  programs = nrow(summary),
  estimable_at_0.05 = sum(summary$status == "estimable_at_0.05"),
  selected_null_programs = sum(summary$selected_sets),
  all_label_blind = all(summary$label_blind_passed),
  clinical_label_files_loaded = FALSE,
  endpoint_scoring_performed = FALSE,
  expression_input_sha256 = sha256_file(file.path(base, "inputs", "GSE299988_tumor_gene_expression.tsv.gz")),
  complete_registry_sha256 = sha256_file(file.path(base, "inputs", "frozen_programs_all9.tsv")),
  generator_sha256 = expected_generator_hash,
  tolerance_ladder_used = 0.05,
  relaxed_tolerance_required = FALSE
)
write_json(audit, audit_path, pretty = TRUE, auto_unbox = TRUE)

core_files <- c(
  file.path(base, "raw", c("GSE299988_Processed_data_files.xlsx", "GSE299988_series_matrix.txt.gz",
                            "GSE299988_family.soft.gz", "SOURCE_SHA256.txt")),
  file.path(base, "inputs", c("GSE299988_all13_gene_expression.tsv.gz",
                               "GSE299988_tumor_gene_expression.tsv.gz", "frozen_programs_all9.tsv")),
  file.path(base, "scripts", c("20_prepare_gse299988_covnull.py",
                                 "21_derive_gse299988_covnull_generator.R",
                                 "21_covariance_null_generator_gse299988.R",
                                 "21_covariance_null_generator_gse299988.R.derivation.json",
                                 "12_validate_covariance_null_general.R",
                                 "22_run_gse299988_covnull_sequence_005.sh",
                                 "23_freeze_gse299988_covnull_master.R")),
  file.path(base, "protocol", "COVARIANCE_PRESERVING_NULL_PROTOCOL_FROZEN_2026-07-13.md"),
  file.path(base, "logs", c("20_prepare_gse299988_covnull.log", "21_derive_generator.log",
                             "sequence_0.05_progress.tsv", "sequence_0.05.complete")),
  file.path(out, c("20_prepare_gse299988_covnull_script.sha256", "GSE299988_matrix_rebuild_SHA256.tsv",
                   "GSE299988_matrix_rebuild_qc.json", "GSE299988_sample_accession_map.tsv",
                   "all_programs_tolerance_0.05_summary.tsv", "master_audit_summary.json"))
)
program_files <- unlist(lapply(programs, function(program) {
  result <- file.path(base, "results", "GSE299988", program, "tolerance_0.05")
  validation <- file.path(base, "validation", program, "tolerance_0.05")
  c(file.path(result, c("FROZEN_SHA256.tsv", "tier_manifest.tsv", "label_blind_audit.json",
                        "selected_null_gene_sets.tsv", "covariance_diagnostic.tsv")),
    file.path(validation, c("VALIDATION_SHA256.tsv", "validation_report.tsv", "validation_summary.json")))
}))
files <- c(core_files, program_files)
if (any(!file.exists(files))) stop("A master-freeze input is missing")
master <- data.table(file = substring(normalizePath(files), nchar(base) + 2L),
                     sha256 = vapply(files, sha256_file, character(1)))
fwrite(master, master_path, sep = "\t")
cat(toJSON(audit, pretty = TRUE, auto_unbox = TRUE), "\n")
if (!all_pass) quit(status = 2L)
