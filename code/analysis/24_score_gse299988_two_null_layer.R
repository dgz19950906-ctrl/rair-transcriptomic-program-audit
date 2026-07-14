#!/usr/bin/env Rscript

# Post-freeze clinical-label layer for the non-identifiable GSE299988 challenge.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(Matrix)
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
required <- c("base", "labels", "legacy-scores", "legacy-challenge", "out")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")

base <- normalizePath(args$base)
labels_path <- normalizePath(args$labels)
legacy_scores_path <- normalizePath(args[["legacy-scores"]])
legacy_challenge_path <- normalizePath(args[["legacy-challenge"]])
out <- normalizePath(args$out, mustWork = FALSE)
if (file.exists(out)) stop("Refusing to overwrite clinical-label layer")
dir.create(file.path(out, "tables"), recursive = TRUE)
dir.create(file.path(out, "manifests"), recursive = TRUE)
dir.create(file.path(out, "logs"), recursive = TRUE)
tables <- file.path(out, "tables"); manifests <- file.path(out, "manifests")

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
master_path <- file.path(base, "manifests", "MASTER_FROZEN_SHA256.tsv")
master <- fread(master_path)
master[, absolute_path := file.path(base, file)]
master[, observed_sha256 := vapply(absolute_path, sha256_file, character(1))]
if (any(master$observed_sha256 != master$sha256)) stop("Pre-label master-freeze hash validation failed")

expr_path <- file.path(base, "inputs", "GSE299988_tumor_gene_expression.tsv.gz")
registry_path <- file.path(base, "inputs", "frozen_programs_all9.tsv")
expr_dt <- fread(expr_path, check.names = FALSE)
if (ncol(expr_dt) != 11L || names(expr_dt)[[1]] != "symbol") stop("Expression schema mismatch")
symbols <- expr_dt[[1]]; sample_ids <- names(expr_dt)[-1L]
if (anyDuplicated(symbols) || !all(grepl("^GSM[0-9]+$", sample_ids))) stop("Expression identity gate failed")
expr <- as.matrix(expr_dt[, -1L]); storage.mode(expr) <- "double"; rownames(expr) <- symbols
z <- t(scale(t(expr))); z <- z[apply(z, 1L, function(v) all(is.finite(v))), , drop = FALSE]
rm(expr_dt, expr)

registry <- fread(registry_path)
labels <- fread(labels_path)
if (nrow(labels) != 10L || !setequal(labels$geo_accession, sample_ids) ||
    !setequal(as.integer(table(labels$analysis_group)), c(5L, 5L))) stop("Frozen label schema mismatch")
setkey(labels, geo_accession)
labels <- labels[sample_ids]
if (!identical(labels$geo_accession, sample_ids)) stop("Label ordering failed")

programs <- c(
  "TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
)
if (!setequal(unique(registry$signature_id), programs)) stop("Registry family mismatch")
adverse_group <- "RAI_nonavid_LN_positive"; reference_group <- "RAI_avid_LN_negative"
x_idx <- which(labels$analysis_group == adverse_group); y_idx <- which(labels$analysis_group == reference_group)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y); df <- nx + ny - 2L
  pooled_sd <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / df)
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
  (1 - 3 / (4 * df - 1)) * (mean(x) - mean(y)) / pooled_sd
}
hedges_rows <- function(score_matrix, xi, yi) {
  nx <- length(xi); ny <- length(yi); df <- nx + ny - 2L
  mx <- rowMeans(score_matrix[, xi, drop = FALSE]); my <- rowMeans(score_matrix[, yi, drop = FALSE])
  vx <- apply(score_matrix[, xi, drop = FALSE], 1L, var)
  vy <- apply(score_matrix[, yi, drop = FALSE], 1L, var)
  pooled_sd <- sqrt(((nx - 1) * vx + (ny - 1) * vy) / df)
  g <- (1 - 3 / (4 * df - 1)) * (mx - my) / pooled_sd
  g[!is.finite(g)] <- NA_real_; g
}
orientation_multiplier <- function(orientation) {
  if (orientation %in% c("higher_more_differentiated", "higher_more_iodide_handling")) -1 else 1
}

legacy_scores <- fread(legacy_scores_path)
legacy_scores <- merge(legacy_scores, labels[, .(workbook_column, geo_accession)],
                       by.x = "sample", by.y = "workbook_column", all.x = TRUE, sort = FALSE)
legacy_challenge <- fread(legacy_challenge_path)
observed_rows <- list(); score_rows <- list(); exact_rows <- list(); cov_rows <- list(); lopo_rows <- list()

for (program in programs) {
  definition_all <- registry[signature_id == program]
  definition <- definition_all[gene %in% rownames(z)]
  tier_path <- file.path(base, "results", "GSE299988", program, "tolerance_0.05", "tier_manifest.tsv")
  tier <- fread(tier_path)
  if (tier$status[[1]] != "estimable_at_0.05" || tier$selected_sets[[1]] != 1000L ||
      nrow(definition) != tier$genes_present[[1]]) stop("Program-level null gate failed for ", program)
  score <- colSums(z[definition$gene, , drop = FALSE] * definition$direction) / nrow(definition)
  orientation <- definition$orientation[[1]]; multiplier <- orientation_multiplier(orientation)
  raw_g <- hedges_g(score[x_idx], score[y_idx]); aligned_g <- multiplier * raw_g

  legacy_s <- legacy_scores[signature_id == program][match(sample_ids, geo_accession)]
  score_delta <- max(abs(score - legacy_s$score))
  legacy_r <- legacy_challenge[signature_id == program]
  raw_delta <- abs(raw_g - legacy_r$raw_hedges_g[[1]])
  aligned_delta <- abs(aligned_g - legacy_r$adverse_aligned_hedges_g[[1]])
  if (!is.finite(score_delta) || score_delta > 1e-12 || raw_delta > 1e-12 || aligned_delta > 1e-12) {
    stop("Legacy numerical reproduction gate failed for ", program)
  }

  allocations <- combn(seq_along(score), length(x_idx))
  exact_values <- apply(allocations, 2L, function(i) hedges_g(score[i], score[-i]))
  exact_extreme <- sum(abs(exact_values) >= abs(raw_g) - 1e-12, na.rm = TRUE)
  exact_p <- exact_extreme / ncol(allocations)

  null_path <- file.path(base, "results", "GSE299988", program, "tolerance_0.05",
                         "selected_null_gene_sets.tsv")
  null_long <- fread(null_path)
  if (uniqueN(null_long$set_id) != 1000L || any(!null_long$gene %in% rownames(z))) stop("Null integrity failed")
  set_levels <- unique(null_long$set_id)
  weights <- sparseMatrix(i = match(null_long$set_id, set_levels),
                          j = match(null_long$gene, rownames(z)), x = null_long$direction,
                          dims = c(length(set_levels), nrow(z)))
  null_scores <- as.matrix(Diagonal(x = 1 / as.numeric(rowSums(abs(weights)))) %*% weights %*% z)
  rownames(null_scores) <- set_levels; colnames(null_scores) <- sample_ids
  null_g <- hedges_rows(null_scores, x_idx, y_idx)
  if (sum(is.finite(null_g)) != 1000L) stop("Non-finite covariance-null effects")
  cov_extreme <- sum(abs(null_g) >= abs(raw_g) - 1e-12)
  cov_p <- (1 + cov_extreme) / (1 + length(null_g))

  omit_ids <- c(x_idx, y_idx)
  lopo <- vapply(omit_ids, function(i) multiplier * hedges_g(score[setdiff(x_idx, i)], score[setdiff(y_idx, i)]), numeric(1))
  score_rows[[program]] <- data.table(
    geo_accession = sample_ids, workbook_column = labels$workbook_column,
    analysis_group = labels$analysis_group, signature_id = program, score = unname(score),
    adverse_orienting_multiplier = multiplier
  )
  exact_rows[[program]] <- data.table(
    signature_id = program, permutation_id = seq_along(exact_values),
    permuted_raw_hedges_g = exact_values, permuted_adverse_hedges_g = multiplier * exact_values,
    observed_raw_hedges_g = raw_g
  )
  cov_rows[[program]] <- data.table(
    signature_id = program, set_id = set_levels, null_raw_hedges_g = null_g,
    null_adverse_hedges_g = multiplier * null_g, observed_raw_hedges_g = raw_g,
    observed_adverse_hedges_g = aligned_g
  )
  lopo_rows[[program]] <- data.table(
    signature_id = program, omitted_geo_accession = sample_ids[omit_ids],
    omitted_group = labels$analysis_group[omit_ids], adverse_aligned_leave_one_out_hedges_g = lopo,
    full_adverse_aligned_hedges_g = aligned_g
  )
  observed_rows[[program]] <- data.table(
    signature_id = program, role = definition$role[[1]], orientation = orientation,
    contrast = "RAI_nonavid_LN_positive_minus_RAI_avid_LN_negative",
    n_adverse = length(x_idx), n_reference = length(y_idx), genes_requested = nrow(definition_all),
    genes_present = nrow(definition), raw_hedges_g = raw_g, adverse_aligned_hedges_g = aligned_g,
    frozen_bootstrap_ci_low = legacy_r$adverse_aligned_bootstrap_g_ci_low[[1]],
    frozen_bootstrap_ci_high = legacy_r$adverse_aligned_bootstrap_g_ci_high[[1]],
    exact_label_allocations = ncol(allocations), exact_label_extreme = exact_extreme,
    exact_label_p = exact_p, covariance_null_sets = length(null_g),
    covariance_extreme = cov_extreme, covariance_program_p = cov_p,
    lopo_same_direction_n = sum(sign(lopo) == sign(aligned_g), na.rm = TRUE),
    lopo_estimable_n = sum(is.finite(lopo)), lopo_min = min(lopo), lopo_max = max(lopo),
    max_absolute_score_reproduction_delta = score_delta,
    absolute_raw_g_reproduction_delta = raw_delta,
    absolute_aligned_g_reproduction_delta = aligned_delta
  )
}

observed <- rbindlist(observed_rows)
if (nrow(observed) != 9L) stop("Nine-test family gate failed")
observed[, exact_label_q_bh9 := p.adjust(exact_label_p, method = "BH")]
observed[, covariance_program_q_bh9 := p.adjust(covariance_program_p, method = "BH")]
observed[, program_order__ := match(signature_id, programs)]
setorder(observed, program_order__)
observed[, program_order__ := NULL]
patient_scores <- rbindlist(score_rows); exact_nulls <- rbindlist(exact_rows)
cov_nulls <- rbindlist(cov_rows); lopo <- rbindlist(lopo_rows)

fwrite(observed, file.path(tables, "GSE299988_two_null_challenge.tsv"), sep = "\t")
fwrite(observed, file.path(tables, "GSE299988_two_null_challenge.csv"))
fwrite(patient_scores, file.path(tables, "GSE299988_patient_program_scores.tsv"), sep = "\t")
fwrite(lopo, file.path(tables, "GSE299988_lopo_estimates.tsv"), sep = "\t")
fwrite(exact_nulls, file.path(tables, "GSE299988_exact_label_nulls.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(cov_nulls, file.path(tables, "GSE299988_covariance_program_nulls.tsv.gz"), sep = "\t", compress = "gzip")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
access <- data.table(
  input_role = c("expression", "complete_program_registry", "public_challenge_labels",
                 "prelabel_master_manifest", "legacy_frozen_scores", "legacy_frozen_effects", "scoring_script"),
  path = c(expr_path, registry_path, labels_path, master_path, legacy_scores_path, legacy_challenge_path, script_path)
)
access[, sha256 := vapply(path, sha256_file, character(1))]
fwrite(access, file.path(manifests, "input_access_manifest.tsv"), sep = "\t")

run_manifest <- list(
  passed = TRUE, cohort = "GSE299988", patients = 10L, programs = 9L, tests = 9L,
  challenge_nonidentifiable = TRUE, endpoint_scoring_performed = TRUE,
  exact_label_allocations_per_test = 252L, covariance_null_sets_per_test = 1000L,
  bh_families = c("9 exact-label tests", "9 covariance-program tests"),
  clinical_labels_accessed_after_null_master_freeze = TRUE,
  master_manifest_sha256 = sha256_file(master_path),
  labels_sha256 = sha256_file(labels_path),
  maximum_score_reproduction_delta = max(observed$max_absolute_score_reproduction_delta),
  maximum_effect_reproduction_delta = max(observed$absolute_raw_g_reproduction_delta,
                                           observed$absolute_aligned_g_reproduction_delta)
)
write_json(run_manifest, file.path(manifests, "GSE299988_two_null_run_manifest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 15)
capture.output(sessionInfo(), file = file.path(out, "logs", "sessionInfo.txt"))

freeze_files <- c(list.files(tables, full.names = TRUE),
                  file.path(manifests, c("input_access_manifest.tsv", "GSE299988_two_null_run_manifest.json")),
                  file.path(out, "logs", "sessionInfo.txt"))
freeze <- data.table(file = substring(normalizePath(freeze_files), nchar(out) + 2L),
                     sha256 = vapply(freeze_files, sha256_file, character(1)))
fwrite(freeze, file.path(manifests, "GSE299988_TWO_NULL_FROZEN_SHA256.tsv"), sep = "\t")

cat(toJSON(run_manifest, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
print(observed[, .(signature_id, adverse_aligned_hedges_g, exact_label_p, exact_label_q_bh9,
                   covariance_program_p, covariance_program_q_bh9)])
