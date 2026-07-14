#!/usr/bin/env Rscript

# Independent validation of the frozen clinical-label analysis layer.
# This script does not refit, regenerate, or select any null program.

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
required <- c("analysis", "out")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")
analysis <- normalizePath(args$analysis)
out <- normalizePath(args$out, mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

tables <- file.path(analysis, "tables")
manifests <- file.path(analysis, "manifests")
freeze_path <- file.path(manifests, "CLINICAL_LAYER_FROZEN_SHA256.tsv")
run_path <- file.path(manifests, "clinical_layer_run_manifest.json")
if (!file.exists(freeze_path) || !file.exists(run_path)) stop("Frozen clinical layer is incomplete")

freeze <- fread(freeze_path)
freeze[, absolute_path := file.path(analysis, file)]
freeze[, exists := file.exists(absolute_path)]
freeze[, observed_sha256 := ifelse(exists, vapply(absolute_path, sha256_file, character(1)), NA_character_)]
freeze[, hash_pass := exists & sha256 == observed_sha256]
if (!all(freeze$hash_pass)) stop("Frozen clinical-layer hash validation failed")

run <- read_json(run_path, simplifyVector = TRUE)
if (!isTRUE(run$passed) || run$patients != 17L || run$programs != 9L || run$tests != 18L ||
    !isTRUE(run$endpoint_scoring_performed)) stop("Clinical run manifest gate failed")
if (run$label_permutations$uptake_failure != 1716L ||
    run$label_permutations$response_failure_with_uptake != 330L) stop("Permutation counts changed")

results <- fread(file.path(tables, "clinical_endpoint_two_null_results.tsv"))
label_nulls <- fread(file.path(tables, "exact_label_permutation_nulls.tsv.gz"))
cov_nulls <- fread(file.path(tables, "covariance_program_identity_nulls.tsv.gz"))
lopo <- fread(file.path(tables, "lopo_estimates.tsv"))
lopo_summary <- fread(file.path(tables, "lopo_summary.tsv"))
patient_scores <- fread(file.path(tables, "patient_program_scores.tsv"))

if (nrow(results) != 18L || uniqueN(results$signature_id) != 9L || uniqueN(results$contrast) != 2L ||
    anyDuplicated(results[, .(signature_id, contrast)])) stop("18-test result family invalid")
if (nrow(patient_scores) != 9L * 17L || uniqueN(patient_scores$patient_id) != 17L) stop("Patient-score table invalid")
expected_counts <- c(RAI_avid_remission = 4L, RAI_avid_persistent = 7L, RAI_nonavid_persistent = 6L)
observed_counts <- table(unique(patient_scores[, .(patient_id, analysis_group)])$analysis_group)
if (!identical(as.integer(observed_counts[names(expected_counts)]), unname(expected_counts))) stop("Clinical counts changed")

label_check <- label_nulls[, .(
  permutations = .N,
  extreme = sum(abs(permuted_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12),
  recomputed_p = sum(abs(permuted_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12) / .N
), by = .(signature_id, contrast)]
cov_check <- cov_nulls[, .(
  null_sets = .N,
  extreme = sum(abs(null_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12),
  recomputed_p = (1 + sum(abs(null_raw_hedges_g) >= abs(observed_raw_hedges_g[[1]]) - 1e-12)) / (1 + .N)
), by = .(signature_id, contrast)]

check <- merge(results, label_check, by = c("signature_id", "contrast"), suffixes = c("", "_label"))
check <- merge(check, cov_check, by = c("signature_id", "contrast"), suffixes = c("", "_cov"))
check[, label_p_pass := abs(exact_label_p - recomputed_p) < 1e-15]
check[, label_count_pass := exact_label_permutations == permutations & exact_label_extreme == extreme]
check[, covariance_p_pass := abs(covariance_program_p - recomputed_p_cov) < 1e-15]
check[, covariance_count_pass := covariance_null_sets == null_sets & covariance_extreme == extreme_cov]

expected_label_n <- ifelse(check$contrast == "uptake_failure", 1716L, 330L)
if (any(check$permutations != expected_label_n) || any(check$null_sets != 1000L)) stop("Null distribution size mismatch")
if (!all(check$label_p_pass & check$label_count_pass & check$covariance_p_pass & check$covariance_count_pass)) {
  stop("Saved P values do not reproduce from complete null distributions")
}

label_q <- p.adjust(results$exact_label_p, method = "BH")
cov_q <- p.adjust(results$covariance_program_p, method = "BH")
# Text serialization can introduce sub-femtoscale differences in otherwise
# identical adjusted probabilities; the threshold remains far below reporting precision.
bh_pass <- max(abs(label_q - results$exact_label_q_bh18)) < 1e-14 &&
  max(abs(cov_q - results$covariance_program_q_bh18)) < 1e-14
if (!bh_pass) stop("BH correction does not reproduce across the 18-test families")

numeric_fields <- c("raw_hedges_g", "adverse_aligned_hedges_g", "bootstrap_ci_low", "bootstrap_ci_high",
                    "exact_label_p", "exact_label_q_bh18", "covariance_program_p", "covariance_program_q_bh18")
if (any(!is.finite(as.matrix(results[, ..numeric_fields])))) stop("Non-finite statistic in final results")
if (any(results$bootstrap_ci_low > results$bootstrap_ci_high)) stop("Bootstrap CI ordering invalid")

lopo_check <- lopo[, .(
  n = .N,
  n_finite = sum(is.finite(adverse_aligned_leave_one_out_hedges_g)),
  min_recomputed = min(adverse_aligned_leave_one_out_hedges_g),
  max_recomputed = max(adverse_aligned_leave_one_out_hedges_g),
  same_direction_recomputed = sum(sign(adverse_aligned_leave_one_out_hedges_g) == sign(full_adverse_aligned_hedges_g[[1]]))
), by = .(signature_id, contrast)]
lopo_check <- merge(lopo_check, lopo_summary, by = c("signature_id", "contrast"))
lopo_check[, denominator_pass := n == ifelse(contrast == "uptake_failure", 13L, 11L)]
lopo_check[, summary_pass := n_finite == n_estimable & abs(min_recomputed - lopo_min) < 1e-15 &
             abs(max_recomputed - lopo_max) < 1e-15 & same_direction_recomputed == same_direction_n]
if (!all(lopo_check$denominator_pass & lopo_check$summary_pass)) stop("LOPO validation failed")

validation <- list(
  passed = TRUE,
  analysis_directory = analysis,
  frozen_files_checked = nrow(freeze),
  frozen_hashes_passed = sum(freeze$hash_pass),
  patients = 17L,
  programs = 9L,
  tests = 18L,
  exact_label_allocations = list(uptake_failure = 1716L, response_failure_with_uptake = 330L),
  covariance_null_sets_per_test = 1000L,
  exact_p_reproduced = TRUE,
  covariance_p_reproduced = TRUE,
  bh18_families_reproduced = TRUE,
  lopo_reproduced = TRUE,
  clinical_input_sha256 = run$clinical_input_sha256,
  null_master_sha256 = run$null_master_sha256
)
write_json(validation, file.path(out, "clinical_layer_validation.json"), pretty = TRUE, auto_unbox = TRUE)
fwrite(check, file.path(out, "pvalue_reproduction.tsv"), sep = "\t")
fwrite(lopo_check, file.path(out, "lopo_validation.tsv"), sep = "\t")
fwrite(freeze[, .(file, expected_sha256 = sha256, observed_sha256, hash_pass)],
       file.path(out, "frozen_hash_validation.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

outputs <- sort(list.files(out, full.names = TRUE))
outputs <- outputs[basename(outputs) != "VALIDATION_SHA256.tsv"]
fwrite(data.table(file = basename(outputs), sha256 = vapply(outputs, sha256_file, character(1))),
       file.path(out, "VALIDATION_SHA256.tsv"), sep = "\t")

cat("Independent clinical-label validation passed.\n")
print(results[, .(program_label, contrast, adverse_aligned_hedges_g,
                  exact_label_p, exact_label_q_bh18,
                  covariance_program_p, covariance_program_q_bh18)])
