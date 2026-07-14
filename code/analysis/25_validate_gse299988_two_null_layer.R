#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table); library(digest); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || args[[1]] != "--analysis" || args[[3]] != "--out") {
  stop("Usage: --analysis DIR --out DIR")
}
analysis <- normalizePath(args[[2]])
out <- normalizePath(args[[4]], mustWork = FALSE)
if (file.exists(out)) stop("Refusing to overwrite validation output")
dir.create(out, recursive = TRUE)
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

tables <- file.path(analysis, "tables"); manifests <- file.path(analysis, "manifests")
freeze <- fread(file.path(manifests, "GSE299988_TWO_NULL_FROZEN_SHA256.tsv"))
freeze[, absolute_path := file.path(analysis, file)]
freeze[, observed_sha256 := vapply(absolute_path, sha256_file, character(1))]
freeze[, hash_pass := observed_sha256 == sha256]
if (!all(freeze$hash_pass)) stop("Frozen clinical-layer hashes failed")

run <- read_json(file.path(manifests, "GSE299988_two_null_run_manifest.json"), simplifyVector = TRUE)
if (!isTRUE(run$passed) || run$patients != 10L || run$programs != 9L || run$tests != 9L ||
    run$exact_label_allocations_per_test != 252L || run$covariance_null_sets_per_test != 1000L ||
    !isTRUE(run$challenge_nonidentifiable) || !isTRUE(run$clinical_labels_accessed_after_null_master_freeze)) {
  stop("Run-manifest gate failed")
}

results <- fread(file.path(tables, "GSE299988_two_null_challenge.tsv"))
scores <- fread(file.path(tables, "GSE299988_patient_program_scores.tsv"))
exact <- fread(file.path(tables, "GSE299988_exact_label_nulls.tsv.gz"))
cov <- fread(file.path(tables, "GSE299988_covariance_program_nulls.tsv.gz"))
lopo <- fread(file.path(tables, "GSE299988_lopo_estimates.tsv"))
if (nrow(results) != 9L || uniqueN(results$signature_id) != 9L ||
    nrow(scores) != 90L || uniqueN(scores$geo_accession) != 10L ||
    nrow(exact) != 9L * 252L || nrow(cov) != 9L * 1000L || nrow(lopo) != 90L) {
  stop("Saved table dimensions failed")
}

exact_check <- exact[, .(
  allocations = .N,
  extreme = sum(abs(permuted_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12),
  p = sum(abs(permuted_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12) / .N
), by = signature_id]
cov_check <- cov[, .(
  sets = .N,
  extreme = sum(abs(null_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12),
  p = (1 + sum(abs(null_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12)) / (1 + .N)
), by = signature_id]
check <- merge(results, exact_check, by = "signature_id", suffixes = c("", "_exact"))
check <- merge(check, cov_check, by = "signature_id", suffixes = c("", "_cov"))
check[, exact_pass := allocations == 252L & extreme == exact_label_extreme & abs(p - exact_label_p) < 1e-15]
check[, cov_pass := sets == 1000L & extreme_cov == covariance_extreme &
        abs(p_cov - covariance_program_p) < 1e-15]
if (!all(check$exact_pass & check$cov_pass)) stop("Null P-value reproduction failed")

exact_q <- p.adjust(results$exact_label_p, method = "BH")
cov_q <- p.adjust(results$covariance_program_p, method = "BH")
if (max(abs(exact_q - results$exact_label_q_bh9)) > 1e-14 ||
    max(abs(cov_q - results$covariance_program_q_bh9)) > 1e-14) stop("BH9 reproduction failed")

lopo_check <- lopo[, .(
  n = .N, n_estimable = sum(is.finite(adverse_aligned_leave_one_out_hedges_g)),
  same_direction_n = sum(sign(adverse_aligned_leave_one_out_hedges_g) ==
                           sign(full_adverse_aligned_hedges_g[[1]]), na.rm = TRUE),
  min_value = min(adverse_aligned_leave_one_out_hedges_g),
  max_value = max(adverse_aligned_leave_one_out_hedges_g)
), by = signature_id]
lopo_check <- merge(lopo_check, results, by = "signature_id")
lopo_check[, pass := n == 10L & n_estimable == lopo_estimable_n &
             same_direction_n == lopo_same_direction_n &
             abs(min_value - lopo_min) < 1e-14 & abs(max_value - lopo_max) < 1e-14]
if (!all(lopo_check$pass)) stop("LOPO reproduction failed")

if (max(results$max_absolute_score_reproduction_delta) > 1e-12 ||
    max(results$absolute_raw_g_reproduction_delta, results$absolute_aligned_g_reproduction_delta) > 1e-12) {
  stop("Legacy reproduction threshold failed")
}

validation <- list(
  passed = TRUE, cohort = "GSE299988", frozen_files_checked = nrow(freeze),
  patients = 10L, programs = 9L, tests = 9L,
  exact_allocations_per_test = 252L, covariance_null_sets_per_test = 1000L,
  exact_p_reproduced = TRUE, covariance_p_reproduced = TRUE,
  bh9_families_reproduced = TRUE, lopo_reproduced = TRUE,
  legacy_scores_and_effects_reproduced = TRUE,
  challenge_nonidentifiable = TRUE,
  master_manifest_sha256 = run$master_manifest_sha256
)
write_json(validation, file.path(out, "GSE299988_two_null_validation.json"),
           pretty = TRUE, auto_unbox = TRUE)
fwrite(check, file.path(out, "pvalue_reproduction.tsv"), sep = "\t")
fwrite(lopo_check, file.path(out, "lopo_reproduction.tsv"), sep = "\t")
fwrite(freeze[, .(file, expected_sha256 = sha256, observed_sha256, hash_pass)],
       file.path(out, "frozen_hash_validation.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
outputs <- list.files(out, full.names = TRUE)
fwrite(data.table(file = basename(outputs), sha256 = vapply(outputs, sha256_file, character(1))),
       file.path(out, "VALIDATION_SHA256.tsv"), sep = "\t")
cat(toJSON(validation, pretty = TRUE, auto_unbox = TRUE), "\n")
