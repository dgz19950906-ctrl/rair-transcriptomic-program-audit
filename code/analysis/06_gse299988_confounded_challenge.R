#!/usr/bin/env Rscript

# Adversarial challenge only: RAI status and LN status are perfectly co-linear.
# No adjusted or RAIR-specific effect is identifiable in this dataset.

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
base <- file.path(root, "phase1_cross_definition")

expr_all <- read.delim(file.path(base, "processed", "GSE299988_gene_expression.tsv.gz"),
                       row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(base, "processed", "GSE299988_samples.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
defs <- read.delim(file.path(base, "config", "frozen_signatures.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)

gmt <- strsplit(readLines(file.path(base, "raw", "h.all.v2025.1.Hs.symbols.gmt")), "\t")
names(gmt) <- vapply(gmt, `[`, character(1), 1)
control_ids <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
                 "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
                 "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE")
control_defs <- do.call(rbind, lapply(control_ids, function(sid) {
  data.frame(signature_id = sid, gene = gmt[[sid]][-(1:2)], direction = 1,
             orientation = "higher_more_aggressive_program",
             source = "MSigDB Hallmark v2025.1.Hs", role = "aggressiveness_negative_control",
             lock_status = "locked", stringsAsFactors = FALSE)
}))
defs <- rbind(defs, control_defs)

tumor <- meta$sample[meta$analysis_group != "adjacent_normal"]
stopifnot(length(tumor) == 10L)
expr <- expr_all[, tumor, drop = FALSE]
z <- t(scale(t(expr)))
group <- setNames(meta$analysis_group[match(tumor, meta$sample)], tumor)

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
multiplier <- function(x) ifelse(x %in% c("higher_more_differentiated",
                                          "higher_more_iodide_handling"), -1, 1)

gene_mean <- rowMeans(expr); gene_sd <- apply(expr, 1, sd)
valid <- names(gene_mean)[is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0]
decile <- function(x) pmin(10L, ceiling(rank(x, ties.method = "average") / length(x) * 10))
nf <- data.frame(gene = valid, mb = decile(gene_mean[valid]), sb = decile(gene_sd[valid]))
nf$key <- paste(nf$mb, nf$sb, sep = "_")
exclude <- unique(defs$gene)
pools <- split(nf$gene[!nf$gene %in% exclude], nf$key[!nf$gene %in% exclude])
rownames(nf) <- nf$gene
sample_null <- function(genes, weights) {
  keys <- nf[genes, "key"]; picked <- character(length(genes))
  for (key in unique(keys)) {
    ii <- which(keys == key); picked[ii] <- sample(pools[[key]], length(ii), replace = FALSE)
  }
  colSums(z[picked, , drop = FALSE] * weights) / sum(abs(weights))
}

set.seed(299988)
result <- list(); lopo_out <- list(); null_out <- list(); scores_out <- list()
for (sid in unique(defs$signature_id)) {
  d <- defs[defs$signature_id == sid, , drop = FALSE]
  present <- d$gene[d$gene %in% rownames(z)]
  if (length(present) < max(3L, ceiling(0.7 * nrow(d)))) next
  weights <- d$direction[match(present, d$gene)]
  score <- colSums(z[present, , drop = FALSE] * weights) / sum(abs(weights))
  x <- score[group == "RAI_nonavid_LN_positive"]
  y <- score[group == "RAI_avid_LN_negative"]
  raw_g <- hedges_g(x, y); mult <- multiplier(d$orientation[1]); aligned_g <- mult * raw_g
  boot <- replicate(10000L, hedges_g(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  raw_ci <- unname(quantile(boot, c(.025, .975), na.rm = TRUE))
  aligned_ci <- sort(mult * raw_ci)
  lopo <- c(vapply(seq_along(x), function(i) hedges_g(x[-i], y), numeric(1)),
            vapply(seq_along(y), function(i) hedges_g(x, y[-i]), numeric(1)))
  ng <- present[present %in% valid]; nw <- weights[match(ng, present)]
  null_g <- replicate(1000L, {
    ns <- sample_null(ng, nw)
    hedges_g(ns[group == "RAI_nonavid_LN_positive"], ns[group == "RAI_avid_LN_negative"])
  })
  result[[sid]] <- data.frame(
    signature_id = sid, role = d$role[1], orientation = d$orientation[1],
    contrast = "RAI_nonavid_LN_positive_minus_RAI_avid_LN_negative",
    n_each_group = 5, raw_hedges_g = raw_g, adverse_aligned_hedges_g = aligned_g,
    raw_bootstrap_g_ci_low = raw_ci[1], raw_bootstrap_g_ci_high = raw_ci[2],
    adverse_aligned_bootstrap_g_ci_low = aligned_ci[1],
    adverse_aligned_bootstrap_g_ci_high = aligned_ci[2],
    exact_permutation_p = exact_p(x, y),
    lopo_direction_stability = mean(sign(lopo) == sign(raw_g), na.rm = TRUE),
    lopo_aligned_low = min(mult * lopo, na.rm = TRUE),
    lopo_aligned_high = max(mult * lopo, na.rm = TRUE),
    genes_requested = nrow(d), genes_present = length(present),
    empirical_random_p = (sum(abs(null_g) >= abs(raw_g), na.rm = TRUE) + 1) /
                         (sum(is.finite(null_g)) + 1),
    null_abs_95th_percentile = unname(quantile(abs(null_g), .95, na.rm = TRUE)),
    stringsAsFactors = FALSE)
  lopo_out[[sid]] <- data.frame(signature_id = sid, omission = seq_along(lopo),
                                adverse_aligned_leave_one_out_g = mult * lopo)
  null_out[[sid]] <- data.frame(signature_id = sid, iteration = seq_along(null_g),
                                null_hedges_g = null_g)
  scores_out[[sid]] <- data.frame(sample = names(score), signature_id = sid,
                                  score = unname(score), group = group[names(score)])
}

results <- do.call(rbind, result)
results$null_evidence_log10_ratio <- log10(results$empirical_random_p /
                                           results$exact_permutation_p)
generic_rows <- results$role == "aggressiveness_negative_control"
generic_max <- max(results$adverse_aligned_hedges_g[generic_rows], na.rm = TRUE)
best_generic_id <- results$signature_id[generic_rows][which.max(
  results$adverse_aligned_hedges_g[generic_rows])]
results$best_observed_generic_program <- best_generic_id
results$specificity_margin_vs_strongest_generic <- ifelse(
  results$role == "aggressiveness_negative_control", NA,
  results$adverse_aligned_hedges_g - generic_max)

# Paired group bootstrap for the candidate-minus-best-generic-family increment.
# The maximum generic control is reselected within each bootstrap replicate.
scores_by_sig <- lapply(scores_out, function(d) setNames(d$score, d$sample))
orient <- setNames(results$orientation, results$signature_id)
generic_ids <- results$signature_id[generic_rows]
candidate_ids <- results$signature_id[!generic_rows]
idx_x <- which(group == "RAI_nonavid_LN_positive")
idx_y <- which(group == "RAI_avid_LN_negative")
set.seed(299989)
increment_rows <- lapply(candidate_ids, function(sid) {
  vals <- replicate(10000L, {
    bx <- sample(idx_x, replace = TRUE); by <- sample(idx_y, replace = TRUE)
    candidate_g <- multiplier(orient[[sid]]) * hedges_g(scores_by_sig[[sid]][bx],
                                                         scores_by_sig[[sid]][by])
    generic_g <- vapply(generic_ids, function(gid) {
      multiplier(orient[[gid]]) * hedges_g(scores_by_sig[[gid]][bx], scores_by_sig[[gid]][by])
    }, numeric(1))
    candidate_g - max(generic_g, na.rm = TRUE)
  })
  ci <- unname(quantile(vals, c(.025, .975), na.rm = TRUE))
  data.frame(signature_id = sid, specificity_increment_bootstrap_low = ci[1],
             specificity_increment_bootstrap_high = ci[2],
             bootstrap_probability_increment_le_zero = mean(vals <= 0, na.rm = TRUE),
             stringsAsFactors = FALSE)
})
increment_table <- do.call(rbind, increment_rows)
results <- merge(results, increment_table, by = "signature_id", all.x = TRUE, sort = FALSE)
results$incremental_evidence_status <- ifelse(
  results$role == "aggressiveness_negative_control", "generic_control",
  ifelse(results$specificity_increment_bootstrap_low > 0 & results$empirical_random_p < .05,
         "positive_increment_and_random_null_pass",
         "no_robust_incremental_evidence"))

td <- file.path(base, "results", "tables"); fd <- file.path(base, "results", "figures")

# Explicit identifiability audit. With five observations in each diagonal cell and
# zero off-diagonal observations, RAI status and LN status are deterministically
# co-linear. Fisher's exact test is reported because expected cell counts are <5.
rai_ln_table <- matrix(c(5, 0, 0, 5), nrow = 2, byrow = TRUE,
                       dimnames = list(RAI = c("avid", "nonavid"), LN = c("negative", "positive")))
rai_ln_fisher_p <- fisher.test(rai_ln_table)$p.value
rai_ln_phi <- 1
write.table(cbind(RAI = rownames(rai_ln_table), as.data.frame(rai_ln_table)),
            file.path(td, "GSE299988_RAI_LN_contingency.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(test = "Fisher_exact_two_sided", p_value = rai_ln_fisher_p,
                       phi = rai_ln_phi, interpretation = "deterministic_co-linearity"),
            file.path(td, "GSE299988_RAI_LN_identifiability.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(results, file.path(td, "GSE299988_confounded_challenge.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, lopo_out), file.path(td, "GSE299988_confounded_challenge_LOPO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
score_table <- do.call(rbind, scores_out)
write.table(score_table, file.path(td, "GSE299988_confounded_challenge_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
con <- gzfile(file.path(td, "GSE299988_confounded_challenge_nulls.tsv.gz"), "wt")
write.table(do.call(rbind, null_out), con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)

ord <- order(results$adverse_aligned_hedges_g)
f <- results[ord, ]; yy <- seq_len(nrow(f))
pdf(file.path(fd, "GSE299988_confounded_challenge_forest.pdf"), width = 9, height = 6.8,
    useDingbats = FALSE)
plot(f$adverse_aligned_hedges_g, yy,
     xlim = range(c(f$lopo_aligned_low, f$lopo_aligned_high, 0), na.rm = TRUE),
     ylim = c(.5, nrow(f) + .5), yaxt = "n", xlab = "Adverse-oriented Hedges g (RAI-nonavid/LN+ minus RAI-avid/LN-)",
     ylab = "", pch = 21, bg = ifelse(f$role == "aggressiveness_negative_control", "#7C3AED", "#0F766E"))
segments(f$lopo_aligned_low, yy, f$lopo_aligned_high, yy, lwd = 2,
         col = ifelse(f$role == "aggressiveness_negative_control", "#7C3AED", "#0F766E"))
abline(v = 0, lty = 3, col = "grey50")
axis(2, at = yy, labels = f$signature_id, las = 2, cex.axis = .68)
mtext("RAI status and LN status are perfectly co-linear; this is a stress test, not RAIR validation",
      side = 3, cex = .8)
dev.off()

pdf(file.path(fd, "GSE299988_confounded_challenge_specificity.pdf"), width = 15, height = 6.8,
    useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(5, 6, 3, 1))
image(x = 1:2, y = 1:2, z = t(rai_ln_table[nrow(rai_ln_table):1, ]),
      col = colorRampPalette(c("white", "#DC2626"))(6), axes = FALSE,
      xlab = "LN status", ylab = "RAI status", main = "A  Non-separable design")
axis(1, at = 1:2, labels = c("LN-", "LN+"))
axis(2, at = 1:2, labels = c("nonavid", "avid"), las = 2)
text(rep(1:2, each = 2), rep(1:2, 2), labels = c(0, 5, 5, 0), cex = 1.5)
mtext("phi = 1 | DESIGN-LEVEL NON-IDENTIFIABILITY", side = 3, line = -1.05,
      cex = .75, col = "#B91C1C", font = 2)
mtext(sprintf("Fisher exact P = %.5f", rai_ln_fisher_p), side = 3,
      line = -2.0, cex = .58, col = "grey25")
mtext("Phi=1 denotes deterministic correspondence; P tests independence and cannot separate RAI from LN effects",
      side = 1, line = 3.7, cex = .47, col = "grey35")
plot(f$adverse_aligned_hedges_g, yy,
     xlim = range(c(f$lopo_aligned_low, f$lopo_aligned_high, 0), na.rm = TRUE),
     ylim = c(.5, nrow(f) + .5), yaxt = "n",
     xlab = "Adverse-oriented Hedges g", ylab = "", pch = 21,
     bg = ifelse(f$role == "aggressiveness_negative_control", "#7C3AED", "#0F766E"),
     main = "B  Confounded challenge effects")
segments(f$lopo_aligned_low, yy, f$lopo_aligned_high, yy, lwd = 2,
         col = ifelse(f$role == "aggressiveness_negative_control", "#7C3AED", "#0F766E"))
abline(v = c(0, generic_max), lty = c(3, 2), col = c("grey50", "#7C3AED"))
axis(2, at = yy, labels = f$signature_id, las = 2, cex.axis = .62)

inc <- results[results$role != "aggressiveness_negative_control", ]
iy <- seq_len(nrow(inc))
plot(inc$specificity_margin_vs_strongest_generic, iy,
     xlim = range(c(inc$specificity_increment_bootstrap_low,
                    inc$specificity_increment_bootstrap_high, 0), na.rm = TRUE),
     ylim = c(.5, nrow(inc) + .5), yaxt = "n", xlab = "g increment vs bootstrap-wise best generic control",
     ylab = "", pch = 21, bg = "#0F766E", main = "C  Specificity increment")
segments(inc$specificity_increment_bootstrap_low, iy,
         inc$specificity_increment_bootstrap_high, iy, lwd = 2, col = "#0F766E")
abline(v = 0, lty = 3, col = "grey50")
axis(2, at = iy, labels = inc$signature_id, las = 2, cex.axis = .72)
dev.off()

print(results, row.names = FALSE)
