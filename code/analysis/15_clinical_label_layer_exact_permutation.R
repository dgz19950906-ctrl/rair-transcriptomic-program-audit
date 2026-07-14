#!/usr/bin/env Rscript

# Analysis: exact patient-label permutation and covariance-null endpoint scoring
# Date: 2026-07-13
# Random seed: 42
# Unit: one pretreatment primary tumor per unique patient

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(digest)
  library(jsonlite)
})
set.seed(42)

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --name value pairs")
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("null_base", "clinical", "out")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")
null_base <- normalizePath(args$null_base); clinical_path <- normalizePath(args$clinical)
out <- normalizePath(args$out, mustWork = FALSE)
tables <- file.path(out, "tables"); manifests <- file.path(out, "manifests")
dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(manifests, "CLINICAL_LAYER_FROZEN_SHA256.tsv"))) stop("Refusing to overwrite a frozen clinical layer")
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

programs <- c(
  "TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
)
program_labels <- c(
  TDS_16 = "TDS-16", IODIDE_HANDLING_11 = "Iodide-handling-11",
  CONDELLO_2025_SIX = "Condello-6",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_HYPOXIA = "Hypoxia", HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint", HALLMARK_E2F_TARGETS = "E2F targets",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
)

# Gate: the complete label-blind null layer must already be frozen and valid.
master_audit_path <- file.path(null_base, "manifests", "master_audit_summary.json")
master_hash_path <- file.path(null_base, "manifests", "MASTER_FROZEN_SHA256.tsv.sha256")
master_audit <- read_json(master_audit_path, simplifyVector = TRUE)
if (!isTRUE(master_audit$passed) || master_audit$programs != 9L || !isTRUE(master_audit$all_label_blind) ||
    isTRUE(master_audit$endpoint_scoring_performed)) stop("Label-blind master gate failed")
master_hash_line <- strsplit(readLines(master_hash_path, warn = FALSE)[1], "[[:space:]]+")[[1]][1]
if (sha256_file(file.path(null_base, "manifests", "MASTER_FROZEN_SHA256.tsv")) != master_hash_line) {
  stop("Master null manifest hash mismatch")
}

expr_path <- file.path(null_base, "inputs", "GSE151179_primary_preRAI_gene_expression.tsv.gz")
registry_path <- file.path(null_base, "inputs", "frozen_programs_all9.tsv")
expr_dt <- fread(expr_path, check.names = FALSE)
symbols <- as.character(expr_dt[[1]]); expr <- as.matrix(expr_dt[, -1L]); storage.mode(expr) <- "double"
rownames(expr) <- symbols; sample_ids <- colnames(expr)
gene_sd <- apply(expr, 1L, sd); eligible <- is.finite(gene_sd) & gene_sd > 0
z <- (expr[eligible, , drop = FALSE] - rowMeans(expr[eligible, , drop = FALSE])) / gene_sd[eligible]
registry <- fread(registry_path); registry[, direction := as.integer(direction)]
if (!setequal(unique(registry$signature_id), programs)) stop("Nine-program registry mismatch")

clinical <- fread(clinical_path)
required_clinical <- c("geo_accession", "patient_id", "eligible_primary_pre_rai", "analysis_group")
if (length(setdiff(required_clinical, names(clinical)))) stop("Clinical sample schema mismatch")
clinical <- clinical[match(sample_ids, geo_accession)]
if (nrow(clinical) != 17L || anyNA(clinical$geo_accession) || !identical(clinical$geo_accession, sample_ids)) {
  stop("Clinical-expression sample alignment failed")
}
if (uniqueN(clinical$patient_id) != 17L || any(!clinical$eligible_primary_pre_rai)) stop("Patient independence gate failed")
expected_groups <- c(RAI_avid_remission = 4L, RAI_avid_persistent = 7L, RAI_nonavid_persistent = 6L)
observed_groups <- table(clinical$analysis_group)
if (!identical(as.integer(observed_groups[names(expected_groups)]), unname(expected_groups))) stop("Clinical group counts changed")

contrasts <- list(
  uptake_failure = list(adverse = "RAI_nonavid_persistent", reference = "RAI_avid_persistent"),
  response_failure_with_uptake = list(adverse = "RAI_avid_persistent", reference = "RAI_avid_remission")
)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y); df <- nx + ny - 2L
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / df)
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
  (1 - 3 / (4 * df - 1)) * (mean(x) - mean(y)) / pooled_sd
}

hedges_rows <- function(score_matrix, x_index, y_index) {
  nx <- length(x_index); ny <- length(y_index); df <- nx + ny - 2L
  mx <- rowMeans(score_matrix[, x_index, drop = FALSE]); my <- rowMeans(score_matrix[, y_index, drop = FALSE])
  vx <- apply(score_matrix[, x_index, drop = FALSE], 1L, var)
  vy <- apply(score_matrix[, y_index, drop = FALSE], 1L, var)
  pooled_sd <- sqrt(((nx - 1) * vx + (ny - 1) * vy) / df)
  g <- (1 - 3 / (4 * df - 1)) * (mx - my) / pooled_sd
  g[!is.finite(g)] <- NA_real_; g
}

exact_permutation <- function(x, y) {
  pooled <- c(x, y); nx <- length(x)
  allocations <- combn(seq_along(pooled), nx)
  values <- apply(allocations, 2L, function(i) hedges_g(pooled[i], pooled[-i]))
  list(values = values, allocations = ncol(allocations))
}

bootstrap_ci <- function(x, y, B = 10000L) {
  values <- replicate(B, hedges_g(sample(x, length(x), replace = TRUE), sample(y, length(y), replace = TRUE)))
  unname(quantile(values, c(0.025, 0.975), na.rm = TRUE, type = 6))
}

orientation_multiplier <- function(orientation) {
  if (orientation %in% c("higher_more_differentiated", "higher_more_iodide_handling")) -1 else 1
}

observed_rows <- list(); patient_score_rows <- list(); label_null_rows <- list()
cov_null_rows <- list(); lopo_rows <- list(); lopo_summary_rows <- list()

for (p_idx in seq_along(programs)) {
  program <- programs[[p_idx]]
  definition_all <- registry[signature_id == program]
  definition <- definition_all[gene %in% rownames(z)]
  tier_path <- file.path(null_base, "results", "GSE151179", program, "tolerance_0.05", "tier_manifest.tsv")
  tier <- fread(tier_path)
  if (tier$status[[1]] != "estimable_at_0.05" || tier$selected_sets[[1]] != 1000L ||
      nrow(definition) != tier$genes_present[[1]]) stop("Program-level null gate failed for ", program)

  observed_score <- colSums(z[definition$gene, , drop = FALSE] * definition$direction) / nrow(definition)
  orientation <- definition$orientation[[1]]; role <- definition$role[[1]]
  multiplier <- orientation_multiplier(orientation)
  patient_score_rows[[program]] <- data.table(
    geo_accession = sample_ids, patient_id = clinical$patient_id, analysis_group = clinical$analysis_group,
    signature_id = program, program_label = unname(program_labels[[program]]),
    score = unname(observed_score), adverse_orienting_multiplier = multiplier
  )

  null_long_path <- file.path(null_base, "results", "GSE151179", program, "tolerance_0.05", "selected_null_gene_sets.tsv")
  null_long <- fread(null_long_path)
  if (uniqueN(null_long$set_id) != 1000L || any(!null_long$gene %in% rownames(z))) stop("Null gene-set integrity failed")
  set_levels <- unique(null_long$set_id)
  weights <- sparseMatrix(
    i = match(null_long$set_id, set_levels), j = match(null_long$gene, rownames(z)),
    x = null_long$direction, dims = c(length(set_levels), nrow(z))
  )
  denominators <- as.numeric(rowSums(abs(weights)))
  null_scores <- Diagonal(x = 1 / denominators) %*% weights %*% z
  null_scores <- as.matrix(null_scores); rownames(null_scores) <- set_levels; colnames(null_scores) <- sample_ids

  for (c_idx in seq_along(contrasts)) {
    contrast <- names(contrasts)[[c_idx]]; contrast_def <- contrasts[[contrast]]
    x_idx <- which(clinical$analysis_group == contrast_def$adverse)
    y_idx <- which(clinical$analysis_group == contrast_def$reference)
    x <- observed_score[x_idx]; y <- observed_score[y_idx]
    observed_raw_g <- hedges_g(x, y); observed_adverse_g <- multiplier * observed_raw_g

    exact <- exact_permutation(x, y)
    exact_extreme <- sum(abs(exact$values) >= abs(observed_raw_g) - 1e-12, na.rm = TRUE)
    exact_p <- exact_extreme / exact$allocations
    ci <- bootstrap_ci(x, y, B = 10000L) * multiplier
    ci <- sort(ci)

    null_raw_g <- hedges_rows(null_scores, x_idx, y_idx)
    if (sum(is.finite(null_raw_g)) != 1000L) stop("Non-finite covariance-null effects for ", program)
    cov_extreme <- sum(abs(null_raw_g) >= abs(observed_raw_g) - 1e-12)
    cov_p <- (1 + cov_extreme) / (1 + length(null_raw_g))

    observed_rows[[paste(program, contrast)]] <- data.table(
      signature_id = program, program_label = unname(program_labels[[program]]), role = role,
      orientation = orientation, contrast = contrast,
      adverse_group = contrast_def$adverse, reference_group = contrast_def$reference,
      n_adverse = length(x), n_reference = length(y), genes_requested = nrow(definition_all),
      genes_present = nrow(definition), raw_hedges_g = observed_raw_g,
      adverse_aligned_hedges_g = observed_adverse_g,
      bootstrap_ci_low = ci[[1]], bootstrap_ci_high = ci[[2]], bootstrap_resamples = 10000L,
      exact_label_permutations = exact$allocations, exact_label_extreme = exact_extreme,
      exact_label_p = exact_p, covariance_null_sets = length(null_raw_g),
      covariance_extreme = cov_extreme, covariance_program_p = cov_p
    )

    label_null_rows[[paste(program, contrast)]] <- data.table(
      signature_id = program, contrast = contrast, permutation_id = seq_along(exact$values),
      permuted_raw_hedges_g = exact$values, permuted_adverse_hedges_g = multiplier * exact$values,
      observed_raw_hedges_g = observed_raw_g
    )
    cov_null_rows[[paste(program, contrast)]] <- data.table(
      signature_id = program, contrast = contrast, set_id = set_levels,
      null_raw_hedges_g = null_raw_g, null_adverse_hedges_g = multiplier * null_raw_g,
      observed_raw_hedges_g = observed_raw_g, observed_adverse_hedges_g = observed_adverse_g
    )

    omit_ids <- c(x_idx, y_idx); omit_groups <- c(rep(contrast_def$adverse, length(x_idx)), rep(contrast_def$reference, length(y_idx)))
    lopo_values <- vapply(seq_along(omit_ids), function(k) {
      omit <- omit_ids[[k]]
      multiplier * hedges_g(observed_score[setdiff(x_idx, omit)], observed_score[setdiff(y_idx, omit)])
    }, numeric(1))
    lopo_rows[[paste(program, contrast)]] <- data.table(
      signature_id = program, contrast = contrast, omitted_geo_accession = sample_ids[omit_ids],
      omitted_group = omit_groups, adverse_aligned_leave_one_out_hedges_g = lopo_values,
      full_adverse_aligned_hedges_g = observed_adverse_g
    )
    lopo_summary_rows[[paste(program, contrast)]] <- data.table(
      signature_id = program, contrast = contrast, full_effect = observed_adverse_g,
      lopo_min = min(lopo_values), lopo_max = max(lopo_values), n_estimable = sum(is.finite(lopo_values)),
      same_direction_n = sum(sign(lopo_values) == sign(observed_adverse_g), na.rm = TRUE)
    )
  }
}

observed <- rbindlist(observed_rows)
if (nrow(observed) != 18L || uniqueN(observed$signature_id) != 9L || uniqueN(observed$contrast) != 2L) stop("18-test family gate failed")
observed[, exact_label_q_bh18 := p.adjust(exact_label_p, method = "BH")]
observed[, covariance_program_q_bh18 := p.adjust(covariance_program_p, method = "BH")]
observed <- observed[order(match(signature_id, programs), match(contrast, names(contrasts)))]

patient_scores <- rbindlist(patient_score_rows); label_nulls <- rbindlist(label_null_rows)
cov_nulls <- rbindlist(cov_null_rows); lopo <- rbindlist(lopo_rows); lopo_summary <- rbindlist(lopo_summary_rows)
lopo_summary <- merge(lopo_summary, observed[, .(signature_id, contrast, program_label, role, orientation,
                                                 exact_label_p, exact_label_q_bh18,
                                                 covariance_program_p, covariance_program_q_bh18)],
                      by = c("signature_id", "contrast"), all.x = TRUE)

fwrite(observed, file.path(tables, "clinical_endpoint_two_null_results.tsv"), sep = "\t")
fwrite(observed, file.path(tables, "clinical_endpoint_two_null_results.csv"))
fwrite(patient_scores, file.path(tables, "patient_program_scores.tsv"), sep = "\t")
fwrite(lopo, file.path(tables, "lopo_estimates.tsv"), sep = "\t")
fwrite(lopo_summary, file.path(tables, "lopo_summary.tsv"), sep = "\t")
fwrite(label_nulls, file.path(tables, "exact_label_permutation_nulls.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(cov_nulls, file.path(tables, "covariance_program_identity_nulls.tsv.gz"), sep = "\t", compress = "gzip")

access <- data.table(
  input_role = c("expression", "complete_program_registry", "public_deidentified_clinical_labels", "label_blind_master_manifest"),
  path = c(expr_path, registry_path, clinical_path, file.path(null_base, "manifests", "MASTER_FROZEN_SHA256.tsv")),
  sha256 = vapply(c(expr_path, registry_path, clinical_path, file.path(null_base, "manifests", "MASTER_FROZEN_SHA256.tsv")), sha256_file, character(1)),
  label_stage = c(FALSE, FALSE, TRUE, FALSE)
)
fwrite(access, file.path(manifests, "clinical_layer_input_access.tsv"), sep = "\t")

script_hash <- if (length(grep("^--file=", commandArgs(FALSE)))) {
  sha256_file(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
} else NA_character_
run_manifest <- list(
  passed = TRUE, cohort = "GSE151179", independent_unit = "patient",
  patients = 17L, group_counts = as.list(expected_groups), programs = 9L, tests = 18L,
  label_permutations = list(uptake_failure = choose(13, 6), response_failure_with_uptake = choose(11, 7)),
  covariance_null_sets_per_test = 1000L, bh_families = c("18 exact-label tests", "18 covariance-program tests"),
  seed = 42L, clinical_input_sha256 = sha256_file(clinical_path),
  null_master_sha256 = master_hash_line, analysis_script_sha256 = script_hash,
  endpoint_scoring_performed = TRUE
)
write_json(run_manifest, file.path(manifests, "clinical_layer_run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(manifests, "sessionInfo.txt"))

output_files <- sort(c(list.files(tables, full.names = TRUE),
                       file.path(manifests, c("clinical_layer_input_access.tsv", "clinical_layer_run_manifest.json", "sessionInfo.txt"))))
fwrite(data.table(file = substring(output_files, nchar(out) + 2L),
                  sha256 = vapply(output_files, sha256_file, character(1))),
       file.path(manifests, "CLINICAL_LAYER_FROZEN_SHA256.tsv"), sep = "\t")

analysis_outputs <- c(
  "# Analysis Outputs", "", "Generated: 2026-07-13",
  "Study type: exact patient-label permutation and covariance-matched program-identity calibration", "",
  "## Tables", "- `clinical_endpoint_two_null_results.tsv/.csv` -- 18 observed tests with effect sizes, 95% bootstrap CIs, exact-label P/q and covariance-program P/q.",
  "- `patient_program_scores.tsv` -- frozen patient-level program scores.",
  "- `lopo_estimates.tsv` and `lopo_summary.tsv` -- leave-one-patient-out robustness.",
  "- `exact_label_permutation_nulls.tsv.gz` -- complete exact allocation distributions.",
  "- `covariance_program_identity_nulls.tsv.gz` -- complete frozen program-identity null effects.", "",
  "## Integrity", "- Exact allocations include the observed allocation.",
  "- BH correction is performed separately across the frozen 18-test family for each null layer.",
  "- No null program was regenerated or selected after clinical labels were read."
)
writeLines(analysis_outputs, file.path(out, "_analysis_outputs.md"))

cat("Clinical label layer complete.\n")
print(observed[, .(program_label, contrast, adverse_aligned_hedges_g, bootstrap_ci_low, bootstrap_ci_high,
                   exact_label_p, exact_label_q_bh18, covariance_program_p, covariance_program_q_bh18)])
