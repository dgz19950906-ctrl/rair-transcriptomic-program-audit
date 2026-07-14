#!/usr/bin/env Rscript

# Label-blind program redundancy audit. Parallel analysis independently permutes
# each program column, preserving every program's marginal score distribution
# while breaking between-program covariance. Endpoint labels are not used to
# select the number of components.

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
base <- file.path(root, "phase1_cross_definition")
td <- file.path(base, "results", "tables")
fd <- file.path(base, "results", "figures")

scores151_long <- read.delim(file.path(td, "frozen_signature_sample_scores.tsv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
meta151 <- read.delim(file.path(base, "processed", "GSE151179_primary_preRAI_samples.tsv"),
                      stringsAsFactors = FALSE, check.names = FALSE)
scores299_long <- read.delim(file.path(td, "GSE299988_confounded_challenge_scores.tsv"),
                             stringsAsFactors = FALSE, check.names = FALSE)

to_matrix <- function(d, sample_col) {
  samples <- unique(d[[sample_col]])
  programs <- unique(d$signature_id)
  out <- matrix(NA_real_, nrow = length(samples), ncol = length(programs),
                dimnames = list(samples, programs))
  for (i in seq_len(nrow(d))) out[d[[sample_col]][i], d$signature_id[i]] <- d$score[i]
  stopifnot(!anyNA(out))
  out
}

m151 <- to_matrix(scores151_long, "geo_accession")
m299 <- to_matrix(scores299_long, "sample")
generic_ids <- grep("^HALLMARK_", colnames(m151), value = TRUE)
candidate_ids <- setdiff(colnames(m151), generic_ids)
stopifnot(setequal(colnames(m151), colnames(m299)))

parallel_analysis <- function(mat, cohort, B = 1000L, seed = 1L) {
  x <- scale(mat)
  obs <- eigen(cor(x), symmetric = TRUE, only.values = TRUE)$values
  set.seed(seed)
  perm <- replicate(B, {
    xp <- apply(x, 2, sample, replace = FALSE)
    eigen(cor(xp), symmetric = TRUE, only.values = TRUE)$values
  })
  q95 <- apply(perm, 1, quantile, probs = .95, names = FALSE, type = 6)
  raw_retained <- obs > q95
  retained_sequential <- as.logical(cumprod(raw_retained))
  data.frame(cohort = cohort, component = seq_along(obs), observed_eigenvalue = obs,
             permuted_95th_percentile = q95,
             exceeds_95th_percentile = raw_retained,
             retained_sequentially = retained_sequential, stringsAsFactors = FALSE)
}

pa151_generic <- parallel_analysis(m151[, generic_ids, drop = FALSE],
                                   "GSE151179_generic_controls", seed = 151179)
pa299_generic <- parallel_analysis(m299[, generic_ids, drop = FALSE],
                                   "GSE299988_generic_controls", seed = 299988)
pa151_all <- parallel_analysis(m151, "GSE151179_all_programs", seed = 151180)
pa299_all <- parallel_analysis(m299, "GSE299988_all_programs", seed = 299989)
pa <- rbind(pa151_generic, pa299_generic, pa151_all, pa299_all)
write.table(pa, file.path(td, "program_parallel_analysis.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

bootstrap_correlations <- function(mat, cohort, B = 1000L, seed = 1L) {
  obs <- cor(mat, method = "spearman")
  pairs <- which(upper.tri(obs), arr.ind = TRUE)
  set.seed(seed)
  rows <- lapply(seq_len(nrow(pairs)), function(k) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    vals <- replicate(B, {
      idx <- sample(seq_len(nrow(mat)), replace = TRUE)
      suppressWarnings(cor(mat[idx, i], mat[idx, j], method = "spearman"))
    })
    data.frame(cohort = cohort, program_1 = colnames(mat)[i], program_2 = colnames(mat)[j],
               spearman_rho = obs[i, j],
               bootstrap_ci_low = unname(quantile(vals, .025, na.rm = TRUE)),
               bootstrap_ci_high = unname(quantile(vals, .975, na.rm = TRUE)),
               bootstrap_sign_stability = mean(sign(vals) == sign(obs[i, j]), na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

corr151 <- bootstrap_correlations(m151, "GSE151179", seed = 7151179)
corr299 <- bootstrap_correlations(m299, "GSE299988", seed = 7299988)
cor_boot <- rbind(corr151, corr299)
write.table(cor_boot, file.path(td, "program_pair_correlations_bootstrap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y); df <- nx + ny - 2
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / df)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (1 - 3 / (4 * df - 1)) * (mean(x) - mean(y)) / sp
}
exact_p <- function(x, y) {
  pooled <- c(x, y); nx <- length(x); obs <- mean(x) - mean(y)
  idx <- combn(seq_along(pooled), nx)
  null <- apply(idx, 2, function(i) mean(pooled[i]) - mean(pooled[-i]))
  (sum(abs(null) >= abs(obs) - 1e-12) + 1) / (length(null) + 1)
}
orientation_multiplier <- function(sid) ifelse(sid %in% c("TDS_16", "IODIDE_HANDLING_11"), -1, 1)

generic_components <- function(mat, pa_table) {
  k <- sum(pa_table$retained_sequentially)
  if (k == 0) return(NULL)
  pc <- prcomp(mat[, generic_ids, drop = FALSE], center = TRUE, scale. = TRUE)
  scores <- pc$x[, seq_len(k), drop = FALSE]
  loadings <- pc$rotation[, seq_len(k), drop = FALSE]
  # Sign is arbitrary. Orient each component so its largest absolute loading is positive.
  for (j in seq_len(k)) {
    mult <- ifelse(loadings[which.max(abs(loadings[, j])), j] < 0, -1, 1)
    scores[, j] <- scores[, j] * mult
    loadings[, j] <- loadings[, j] * mult
  }
  list(scores = scores, loadings = loadings, k = k,
       variance_each = summary(pc)$importance[2, seq_len(k)],
       variance_cumulative = sum(summary(pc)$importance[2, seq_len(k)]))
}

pc151 <- generic_components(m151, pa151_generic)
pc299 <- generic_components(m299, pa299_generic)
loading_rows <- list()
long_loadings <- function(obj, cohort) do.call(rbind, lapply(seq_len(obj$k), function(j) {
  data.frame(cohort = cohort, component = paste0("PC", j), program = rownames(obj$loadings),
             loading = obj$loadings[, j], component_variance_explained = obj$variance_each[j],
             retained_components = obj$k,
             retained_variance_cumulative = obj$variance_cumulative, stringsAsFactors = FALSE)
}))
if (!is.null(pc151)) loading_rows[["GSE151179"]] <- long_loadings(pc151, "GSE151179")
if (!is.null(pc299)) loading_rows[["GSE299988"]] <- long_loadings(pc299, "GSE299988")
if (length(loading_rows)) write.table(do.call(rbind, loading_rows),
                                      file.path(td, "generic_dimension_loadings.tsv"),
                                      sep = "\t", quote = FALSE, row.names = FALSE)

# With two retained components, individual PC axes may rotate. Assess stability of
# the retained subspace rather than over-interpreting any single loading vector.
bootstrap_subspace <- function(mat, obj, cohort, B = 1000L, seed = 1L) {
  if (is.null(obj)) return(NULL)
  ref <- obj$loadings
  set.seed(seed)
  vals <- lapply(seq_len(B), function(b) {
    idx <- sample(seq_len(nrow(mat)), replace = TRUE)
    xb <- scale(mat[idx, generic_ids, drop = FALSE])
    if (any(!is.finite(xb))) return(c(min_canonical_correlation = NA,
                                      mean_canonical_correlation = NA,
                                      projection_distance = NA,
                                      cumulative_variance = NA))
    pb <- prcomp(xb, center = FALSE, scale. = FALSE)
    vb <- pb$rotation[, seq_len(obj$k), drop = FALSE]
    cc <- svd(t(ref) %*% vb, nu = 0, nv = 0)$d
    proj_dist <- sqrt(sum((ref %*% t(ref) - vb %*% t(vb))^2))
    c(min_canonical_correlation = min(cc), mean_canonical_correlation = mean(cc),
      projection_distance = proj_dist,
      cumulative_variance = sum(summary(pb)$importance[2, seq_len(obj$k)]))
  })
  d <- as.data.frame(do.call(rbind, vals))
  d$cohort <- cohort; d$iteration <- seq_len(nrow(d))
  d
}
sub151 <- bootstrap_subspace(m151, pc151, "GSE151179", seed = 8151179)
sub299 <- bootstrap_subspace(m299, pc299, "GSE299988", seed = 8299988)
subspace <- do.call(rbind, Filter(Negate(is.null), list(sub151, sub299)))
if (!is.null(subspace)) {
  sub_con <- gzfile(file.path(td, "generic_subspace_bootstrap_iterations.tsv.gz"), "wt")
  write.table(subspace, sub_con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(sub_con)
  subspace_summary <- do.call(rbind, lapply(split(subspace, subspace$cohort), function(d) {
    data.frame(cohort = d$cohort[1], retained_components = ifelse(d$cohort[1] == "GSE151179", pc151$k, pc299$k),
               min_canonical_correlation_median = median(d$min_canonical_correlation, na.rm = TRUE),
               min_canonical_correlation_q025 = unname(quantile(d$min_canonical_correlation, .025, na.rm = TRUE)),
               min_canonical_correlation_q975 = unname(quantile(d$min_canonical_correlation, .975, na.rm = TRUE)),
               projection_distance_median = median(d$projection_distance, na.rm = TRUE),
               cumulative_variance_median = median(d$cumulative_variance, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  write.table(subspace_summary, file.path(td, "generic_subspace_bootstrap_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

residual_results <- list()
if (!is.null(pc151)) {
  groups <- setNames(meta151$analysis_group[match(rownames(m151), meta151$geo_accession)],
                     rownames(m151))
  contrasts <- list(
    uptake_failure = c("RAI_nonavid_persistent", "RAI_avid_persistent"),
    response_failure_with_uptake = c("RAI_avid_persistent", "RAI_avid_remission"))
  for (sid in candidate_ids) {
    fit <- lm(m151[, sid] ~ ., data = as.data.frame(pc151$scores))
    residual <- resid(fit)
    r2 <- summary(fit)$r.squared
    for (cn in names(contrasts)) {
      gd <- contrasts[[cn]]; x <- residual[groups == gd[1]]; y <- residual[groups == gd[2]]
      mult <- orientation_multiplier(sid)
      lopo <- c(vapply(seq_along(x), function(i) hedges_g(x[-i], y), numeric(1)),
                vapply(seq_along(y), function(i) hedges_g(x, y[-i]), numeric(1)))
      raw_g <- hedges_g(x, y)
      residual_results[[paste("151", sid, cn)]] <- data.frame(
        cohort = "GSE151179", signature_id = sid, contrast = cn,
        generic_components_retained = pc151$k,
        generic_variance_cumulative = pc151$variance_cumulative,
        variance_removed_from_candidate = r2, residual_raw_hedges_g = raw_g,
        residual_adverse_aligned_g = mult * raw_g, residual_exact_permutation_p = exact_p(x, y),
        residual_LOPO_direction_count = sum(sign(lopo) == sign(raw_g), na.rm = TRUE),
        residual_LOPO_total = length(lopo), stringsAsFactors = FALSE)
    }
  }
}
if (!is.null(pc299)) {
  groups <- setNames(scores299_long$group[match(rownames(m299), scores299_long$sample)], rownames(m299))
  for (sid in candidate_ids) {
    fit <- lm(m299[, sid] ~ ., data = as.data.frame(pc299$scores))
    residual <- resid(fit)
    r2 <- summary(fit)$r.squared
    x <- residual[groups == "RAI_nonavid_LN_positive"]
    y <- residual[groups == "RAI_avid_LN_negative"]
    mult <- orientation_multiplier(sid); raw_g <- hedges_g(x, y)
    lopo <- c(vapply(seq_along(x), function(i) hedges_g(x[-i], y), numeric(1)),
              vapply(seq_along(y), function(i) hedges_g(x, y[-i]), numeric(1)))
    residual_results[[paste("299", sid)]] <- data.frame(
      cohort = "GSE299988", signature_id = sid, contrast = "confounded_challenge",
      generic_components_retained = pc299$k,
      generic_variance_cumulative = pc299$variance_cumulative,
      variance_removed_from_candidate = r2, residual_raw_hedges_g = raw_g,
      residual_adverse_aligned_g = mult * raw_g, residual_exact_permutation_p = exact_p(x, y),
      residual_LOPO_direction_count = sum(sign(lopo) == sign(raw_g), na.rm = TRUE),
      residual_LOPO_total = length(lopo), stringsAsFactors = FALSE)
  }
}
if (length(residual_results)) {
  write.table(do.call(rbind, residual_results), file.path(td, "residual_specificity_sensitivity.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  write.table(data.frame(note = "No generic component passed the sequential parallel-analysis threshold; residual alignment was not estimated."),
              file.path(td, "residual_specificity_sensitivity.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}

# Condello correlation-neighborhood shift. Independent patient bootstraps quantify
# uncertainty; this is descriptive because the two cohorts differ in endpoint,
# platform and confounding structure.
set.seed(88001)
condello_rows <- lapply(generic_ids, function(gid) {
  r151 <- suppressWarnings(cor(m151[, "CONDELLO_2025_SIX"], m151[, gid], method = "spearman"))
  r299 <- suppressWarnings(cor(m299[, "CONDELLO_2025_SIX"], m299[, gid], method = "spearman"))
  diffs <- replicate(1000L, {
    i151 <- sample(seq_len(nrow(m151)), replace = TRUE)
    i299 <- sample(seq_len(nrow(m299)), replace = TRUE)
    a <- suppressWarnings(cor(m151[i151, "CONDELLO_2025_SIX"], m151[i151, gid], method = "spearman"))
    b <- suppressWarnings(cor(m299[i299, "CONDELLO_2025_SIX"], m299[i299, gid], method = "spearman"))
    b - a
  })
  data.frame(generic_program = gid, rho_GSE151179 = r151, rho_GSE299988 = r299,
             rho_difference_299_minus_151 = r299 - r151,
             bootstrap_difference_ci_low = unname(quantile(diffs, .025, na.rm = TRUE)),
             bootstrap_difference_ci_high = unname(quantile(diffs, .975, na.rm = TRUE)),
             stringsAsFactors = FALSE)
})
condello_shift <- do.call(rbind, condello_rows)
write.table(condello_shift, file.path(td, "condello_correlation_neighborhood_shift.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_corr <- function(mat, main) {
  cm <- cor(mat, method = "spearman")
  ord <- hclust(as.dist(1 - cm), method = "average")$order
  cm <- cm[ord, ord]
  image(seq_len(nrow(cm)), seq_len(ncol(cm)), t(cm[nrow(cm):1, ]),
        col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
        zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "", main = main)
  axis(1, at = seq_len(ncol(cm)), labels = colnames(cm), las = 2, cex.axis = .48)
  axis(2, at = seq_len(nrow(cm)), labels = rev(rownames(cm)), las = 2, cex.axis = .48)
}

pdf(file.path(fd, "program_redundancy_parallel_analysis.pdf"), width = 12, height = 9,
    useDingbats = FALSE)
par(mfrow = c(2, 2), mar = c(8, 8, 3, 1))
plot_corr(m151, "A  GSE151179 program correlations")
plot_corr(m299, "B  GSE299988 program correlations")
par(mar = c(5, 5, 3, 1))
plot(pa151_generic$component, pa151_generic$observed_eigenvalue, type = "b", pch = 19,
     xlab = "Component", ylab = "Eigenvalue", main = "C  GSE151179 parallel analysis")
lines(pa151_generic$component, pa151_generic$permuted_95th_percentile, type = "b",
      pch = 1, lty = 2, col = "#B2182B")
legend("topright", c("Observed", "Permutation 95th"), lty = c(1, 2), pch = c(19, 1),
       col = c("black", "#B2182B"), bty = "n", cex = .75)
plot(pa299_generic$component, pa299_generic$observed_eigenvalue, type = "b", pch = 19,
     xlab = "Component", ylab = "Eigenvalue", main = "D  GSE299988 parallel analysis")
lines(pa299_generic$component, pa299_generic$permuted_95th_percentile, type = "b",
      pch = 1, lty = 2, col = "#B2182B")
dev.off()

pdf(file.path(fd, "condello_correlation_neighborhood_shift.pdf"), width = 8, height = 5.5,
    useDingbats = FALSE)
yy <- seq_len(nrow(condello_shift))
plot(condello_shift$rho_difference_299_minus_151, yy,
     xlim = range(c(condello_shift$bootstrap_difference_ci_low,
                    condello_shift$bootstrap_difference_ci_high, 0), na.rm = TRUE),
     ylim = c(.5, nrow(condello_shift) + .5), yaxt = "n",
     xlab = "Spearman rho difference (GSE299988 - GSE151179)", ylab = "", pch = 19)
segments(condello_shift$bootstrap_difference_ci_low, yy,
         condello_shift$bootstrap_difference_ci_high, yy, lwd = 2)
abline(v = 0, lty = 3, col = "grey50")
axis(2, at = yy, labels = condello_shift$generic_program, las = 2, cex.axis = .7)
dev.off()

cat("Generic PCs retained in GSE151179:", sum(pa151_generic$retained_sequentially), "\n")
cat("Generic PCs retained in GSE299988:", sum(pa299_generic$retained_sequentially), "\n")
if (!is.null(pc151)) cat("GSE151179 retained variance:", pc151$variance_cumulative, "\n")
if (!is.null(pc299)) cat("GSE299988 retained variance:", pc299$variance_cumulative, "\n")
