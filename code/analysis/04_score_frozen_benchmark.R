#!/usr/bin/env Rscript

# Pilot endpoint-alignment benchmark for signatures locked independently of the
# 17-patient GSE151179 endpoint analysis.

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
root <- normalizePath(file.path(dirname(script_path), "..", ".."))
base <- file.path(root, "phase1_cross_definition")

expr <- read.delim(file.path(base, "processed", "GSE151179_primary_preRAI_gene_expression.tsv.gz"),
                   row.names = 1, check.names = FALSE)
samples <- read.delim(file.path(base, "processed", "GSE151179_primary_preRAI_samples.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE)
defs <- read.delim(file.path(base, "config", "frozen_signatures.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expr), samples$geo_accession))

# Prespecified aggressiveness controls from the official MSigDB Hallmark file.
gmt_path <- file.path(base, "raw", "h.all.v2025.1.Hs.symbols.gmt")
control_ids <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
                 "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
                 "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE")
if (file.exists(gmt_path)) {
  gmt <- strsplit(readLines(gmt_path), "\t", fixed = FALSE)
  names(gmt) <- vapply(gmt, `[`, character(1), 1)
  gmt <- gmt[intersect(control_ids, names(gmt))]
  control_defs <- do.call(rbind, lapply(names(gmt), function(sid) {
    data.frame(signature_id = sid, gene = gmt[[sid]][-(1:2)], direction = 1,
               orientation = "higher_more_aggressive_program",
               source = "MSigDB Hallmark v2025.1.Hs",
               role = "aggressiveness_negative_control", lock_status = "locked",
               stringsAsFactors = FALSE)
  }))
  defs <- rbind(defs, control_defs)
}

table_dir <- file.path(base, "results", "tables")
figure_dir <- file.path(base, "results", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

z <- t(scale(t(expr)))

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

bootstrap_ci <- function(x, y, B = 10000L) {
  vals <- replicate(B, mean(sample(x, replace = TRUE)) - mean(sample(y, replace = TRUE)))
  unname(quantile(vals, c(0.025, 0.975), names = FALSE, type = 6))
}

get_group <- function(score, group) score[samples$analysis_group == group]
contrast_defs <- list(
  uptake_failure = c("RAI_nonavid_persistent", "RAI_avid_persistent"),
  response_failure_with_uptake = c("RAI_avid_persistent", "RAI_avid_remission")
)

orientation_multiplier <- function(x) {
  ifelse(x %in% c("higher_more_differentiated", "higher_more_iodide_handling"), -1, 1)
}

# Null programs match each real gene on deciles of cohort mean expression and
# expression SD. Direction counts are preserved. Co-expression is deliberately
# not matched and is treated as a limitation, not as evidence against a set.
gene_mean <- rowMeans(expr, na.rm = TRUE)
gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
valid_null <- names(gene_mean)[is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0]
decile <- function(x) pmin(10L, ceiling(rank(x, ties.method = "average") / length(x) * 10))
null_frame <- data.frame(gene = valid_null,
                         mean_bin = decile(gene_mean[valid_null]),
                         sd_bin = decile(gene_sd[valid_null]), stringsAsFactors = FALSE)
rownames(null_frame) <- null_frame$gene

null_frame$key <- paste(null_frame$mean_bin, null_frame$sd_bin, sep = "_")
excluded_real_genes <- unique(defs$gene)
pool_by_key <- split(null_frame$gene[!null_frame$gene %in% excluded_real_genes],
                     null_frame$key[!null_frame$gene %in% excluded_real_genes])

sample_matched_set <- function(target_genes, target_weights) {
  target <- null_frame[target_genes, c("mean_bin", "sd_bin"), drop = FALSE]
  target$key <- paste(target$mean_bin, target$sd_bin, sep = "_")
  picked <- character(length(target_genes))
  for (key in unique(target$key)) {
    ii <- which(target$key == key)
    pool <- pool_by_key[[key]]
    if (length(pool) < length(ii)) stop("Insufficient matched null pool for bin ", key)
    picked[ii] <- sample(pool, length(ii), replace = FALSE)
  }
  list(genes = picked, weights = target_weights)
}

score_rows <- list(); result_rows <- list(); lopo_rows <- list(); null_rows <- list()
set.seed(151179)
for (sid in unique(defs$signature_id)) {
  d <- defs[defs$signature_id == sid, , drop = FALSE]
  present <- d$gene[d$gene %in% rownames(z)]
  missing <- setdiff(d$gene, present)
  if (length(present) < max(3L, ceiling(0.7 * nrow(d)))) next
  weights <- d$direction[match(present, d$gene)]
  score <- colSums(z[present, , drop = FALSE] * weights) / sum(abs(weights))
  score_rows[[sid]] <- data.frame(
    geo_accession = names(score), signature_id = sid, score = unname(score),
    genes_requested = nrow(d), genes_present = length(present),
    missing_genes = ifelse(length(missing), paste(missing, collapse = ";"), "none"),
    orientation = d$orientation[1], role = d$role[1], stringsAsFactors = FALSE
  )
  for (cn in names(contrast_defs)) {
    groups <- contrast_defs[[cn]]
    x <- get_group(score, groups[1]); y <- get_group(score, groups[2])
    ci <- bootstrap_ci(x, y)
    full_diff <- mean(x) - mean(y)
    lx <- vapply(seq_along(x), function(i) hedges_g(x[-i], y), numeric(1))
    ly <- vapply(seq_along(y), function(i) hedges_g(x, y[-i]), numeric(1))
    lopo <- c(lx, ly); full_g <- hedges_g(x, y)
    mult <- orientation_multiplier(d$orientation[1])
    result_rows[[paste(sid, cn)]] <- data.frame(
      signature_id = sid, orientation = d$orientation[1], role = d$role[1],
      contrast = cn, adverse_group = groups[1], reference_group = groups[2],
      n_adverse = length(x), n_reference = length(y),
      mean_difference = full_diff, bootstrap_ci_low = ci[1], bootstrap_ci_high = ci[2],
      hedges_g = full_g, adverse_aligned_hedges_g = mult * full_g,
      exact_permutation_p = exact_p(x, y),
      lopo_min_hedges_g = min(lopo, na.rm = TRUE), lopo_max_hedges_g = max(lopo, na.rm = TRUE),
      lopo_direction_stability = mean(sign(lopo) == sign(full_g), na.rm = TRUE),
      genes_requested = nrow(d), genes_present = length(present),
      stringsAsFactors = FALSE
    )
    omitted <- c(names(x), names(y))
    omitted_group <- c(rep(groups[1], length(x)), rep(groups[2], length(y)))
    lopo_rows[[paste(sid, cn)]] <- data.frame(
      signature_id = sid, contrast = cn, omitted_sample = omitted,
      omitted_group = omitted_group,
      leave_one_out_hedges_g = lopo, full_hedges_g = full_g,
      adverse_aligned_leave_one_out_hedges_g = mult * lopo,
      stringsAsFactors = FALSE
    )

    # Empirical null distribution: 1,000 matched programs per signature/contrast.
    ng <- present[present %in% valid_null]
    nw <- weights[match(ng, present)]
    if (length(ng) >= 3L) {
      null_g <- replicate(1000L, {
        ns <- sample_matched_set(ng, nw)
        null_score <- colSums(z[ns$genes, , drop = FALSE] * ns$weights) / sum(abs(ns$weights))
        hedges_g(get_group(null_score, groups[1]), get_group(null_score, groups[2]))
      })
      null_rows[[paste(sid, cn)]] <- data.frame(
        signature_id = sid, contrast = cn, null_iteration = seq_along(null_g),
        null_hedges_g = null_g, observed_hedges_g = full_g,
        empirical_two_sided_p = (sum(abs(null_g) >= abs(full_g), na.rm = TRUE) + 1) /
                                (sum(is.finite(null_g)) + 1), stringsAsFactors = FALSE)
    }
  }
}

scores <- do.call(rbind, score_rows)
results <- do.call(rbind, result_rows)
lopo <- do.call(rbind, lopo_rows)
nulls <- do.call(rbind, null_rows)
write.table(scores, file.path(table_dir, "frozen_signature_sample_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lopo, file.path(table_dir, "frozen_signature_LOPO_estimates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
null_con <- gzfile(file.path(table_dir, "matched_random_program_nulls.tsv.gz"), "wt")
write.table(nulls, null_con, sep = "\t", quote = FALSE, row.names = FALSE)
close(null_con)
null_summary <- do.call(rbind, lapply(split(nulls, interaction(nulls$signature_id, nulls$contrast,
                                                               drop = TRUE)), function(d) {
  data.frame(signature_id = d$signature_id[1], contrast = d$contrast[1],
             observed_hedges_g = d$observed_hedges_g[1],
             empirical_two_sided_p = d$empirical_two_sided_p[1],
             null_abs_95th_percentile = unname(quantile(abs(d$null_hedges_g), 0.95, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}))
write.table(null_summary, file.path(table_dir, "matched_random_program_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Specificity is incremental alignment beyond the strongest prespecified generic
# aggressiveness control on the same endpoint. In this small cohort the margin is
# descriptive; no causal or predictive meaning is attached.
results <- merge(results,
                 null_summary[, c("signature_id", "contrast", "empirical_two_sided_p",
                                  "null_abs_95th_percentile")],
                 by = c("signature_id", "contrast"), all.x = TRUE, sort = FALSE)
results$null_evidence_log10_ratio <- log10(results$empirical_two_sided_p /
                                           results$exact_permutation_p)
best_generic <- aggregate(adverse_aligned_hedges_g ~ contrast,
                          results[results$role == "aggressiveness_negative_control", ], max)
names(best_generic)[2] <- "best_generic_adverse_aligned_g"
results <- merge(results, best_generic, by = "contrast", all.x = TRUE, sort = FALSE)
results$specificity_increment_vs_best_generic <- ifelse(
  results$role == "aggressiveness_negative_control", NA_real_,
  results$adverse_aligned_hedges_g - results$best_generic_adverse_aligned_g)
write.table(results, file.path(table_dir, "frozen_signature_endpoint_alignment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

up <- results[results$contrast == "uptake_failure", c("signature_id", "adverse_aligned_hedges_g")]
names(up)[2] <- "uptake_g"
rp <- results[results$contrast == "response_failure_with_uptake", c("signature_id", "adverse_aligned_hedges_g")]
names(rp)[2] <- "response_g"
plotdat <- merge(up, rp, by = "signature_id")
role_map <- unique(results[, c("signature_id", "role")])
plotdat <- merge(plotdat, role_map, by = "signature_id", all.x = TRUE)
cols <- c(biological_reference = "#0F766E", independent_ex_vivo_avidity_candidate = "#DC2626",
          aggressiveness_negative_control = "#7C3AED")
pdf(file.path(figure_dir, "two_axis_endpoint_alignment_adverse_oriented.pdf"), width = 8.2, height = 7.0,
    useDingbats = FALSE)
par(mar = c(5, 5, 1.2, 1.2))
lim <- range(c(plotdat$uptake_g, plotdat$response_g, -0.2, 0.2))
lim <- max(abs(lim)) * c(-1, 1)
plot(plotdat$uptake_g, plotdat$response_g, xlim = lim, ylim = lim,
     xlab = "Alignment with uptake failure (adverse-oriented Hedges g)",
     ylab = "Alignment with response failure despite uptake (adverse-oriented Hedges g)",
     pch = 21, cex = 1.6, bg = cols[plotdat$role], col = "white", lwd = 1.2)
abline(h = 0, v = 0, lty = 3, col = "grey60")
text(plotdat$uptake_g, plotdat$response_g, labels = plotdat$signature_id,
     pos = 3, cex = 0.75)
dev.off()

# Main-text LOPO forest: interval is the range of leave-one-patient-out Hedges g.
forest <- merge(results[, c("signature_id", "contrast", "adverse_aligned_hedges_g",
                            "lopo_direction_stability", "bootstrap_ci_low", "bootstrap_ci_high", "orientation")],
                aggregate(adverse_aligned_leave_one_out_hedges_g ~ signature_id + contrast,
                          lopo, function(v) c(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))))
forest$lopo_low <- forest$adverse_aligned_leave_one_out_hedges_g[, "min"]
forest$lopo_high <- forest$adverse_aligned_leave_one_out_hedges_g[, "max"]
forest$class <- ifelse(forest$lopo_direction_stability >= 0.8, "direction-stable; estimate near zero allowed",
                       "direction-unstable")
write.table(forest[, setdiff(names(forest), "adverse_aligned_leave_one_out_hedges_g")],
            file.path(table_dir, "LOPO_forest_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(figure_dir, "LOPO_stability_forest.pdf"), width = 9, height = 8, useDingbats = FALSE)
ord <- order(forest$contrast, forest$adverse_aligned_hedges_g)
f <- forest[ord, ]; yy <- seq_len(nrow(f))
plot(f$adverse_aligned_hedges_g, yy, xlim = range(c(f$lopo_low, f$lopo_high, 0), na.rm = TRUE),
     ylim = c(0.5, nrow(f) + 0.5), yaxt = "n", xlab = "Adverse-oriented Hedges g", ylab = "",
     pch = 21, bg = ifelse(f$lopo_direction_stability >= 0.8, "#0F766E", "#9CA3AF"))
segments(f$lopo_low, yy, f$lopo_high, yy, col = ifelse(f$lopo_direction_stability >= 0.8, "#0F766E", "#9CA3AF"), lwd = 2)
abline(v = 0, lty = 3, col = "grey50")
axis(2, at = yy, labels = paste(f$signature_id, f$contrast, sep = " | "), las = 2, cex.axis = 0.58)
mtext("Iodide-handling uptake axis: 12/13 LOPO directions preserved; the sole reversal followed omission of GSM4567934 from the n=6 nonavid/persistent group.",
      side = 3, cex = 0.62, line = 0.2)
dev.off()

print(results, row.names = FALSE)
