#!/usr/bin/env Rscript

# Figure 3: patient-deletion stability and two-null calibration
# Date: 2026-07-13
# Backend: R only

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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
required <- c("analysis", "validation", "out")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")
analysis <- normalizePath(args$analysis)
validation <- normalizePath(args$validation)
out <- normalizePath(args$out, mustWork = FALSE)
figures <- file.path(out, "figures"); tables <- file.path(out, "tables")
manuscript <- file.path(out, "manuscript"); manifests <- file.path(out, "manifests")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript, recursive = TRUE, showWarnings = FALSE)
dir.create(manifests, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(manifests, "FIGURE3_TWO_NULL_FROZEN_SHA256.tsv"))) stop("Refusing to overwrite frozen Figure 3")
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

# Input integrity gates.
validation_json <- read_json(file.path(validation, "clinical_layer_validation.json"), simplifyVector = TRUE)
if (!isTRUE(validation_json$passed) || !isTRUE(validation_json$exact_p_reproduced) ||
    !isTRUE(validation_json$covariance_p_reproduced) || !isTRUE(validation_json$bh18_families_reproduced) ||
    !isTRUE(validation_json$lopo_reproduced)) stop("Independent validation gate failed")
freeze <- fread(file.path(analysis, "manifests", "CLINICAL_LAYER_FROZEN_SHA256.tsv"))
freeze[, absolute_path := file.path(analysis, file)]
if (any(!file.exists(freeze$absolute_path)) ||
    any(vapply(freeze$absolute_path, sha256_file, character(1)) != freeze$sha256)) stop("Clinical layer changed after validation")

results <- fread(file.path(analysis, "tables", "clinical_endpoint_two_null_results.tsv"))
lopo_summary <- fread(file.path(analysis, "tables", "lopo_summary.tsv"))
if (nrow(results) != 18L || uniqueN(results$signature_id) != 9L || uniqueN(results$contrast) != 2L) stop("18-test family invalid")
if (any(lopo_summary$n_estimable != ifelse(lopo_summary$contrast == "uptake_failure", 13L, 11L))) stop("LOPO denominator changed")

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
contrast_labels <- c(
  uptake_failure = "Uptake failure",
  response_failure_with_uptake = "Response failure despite uptake"
)
role_labels <- c(
  biological_reference = "Biological reference",
  independent_ex_vivo_avidity_candidate = "Ex vivo candidate",
  aggressiveness_negative_control = "Generic control"
)

font_family <- "Nimbus Sans"
width_mm <- 183; height_mm <- 196
width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
palette <- c(
  ink = "#202124", neutral = "#73777C", neutral_light = "#D5D9DD",
  biological = "#3977A8", condello = "#C56F3E", generic = "#777B80",
  endpoint_null = "#356D9D", identity_null = "#D0802F",
  strip = "#F1F3F4", grid = "#E5E7E9"
)
role_colors <- c(
  "Biological reference" = palette[["biological"]],
  "Ex vivo candidate" = palette[["condello"]],
  "Generic control" = palette[["generic"]]
)
null_colors <- c(
  "Exact patient-label null" = palette[["endpoint_null"]],
  "Covariance-matched program null" = palette[["identity_null"]]
)

theme_nature <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black", size = base_size - 0.2),
      axis.title = element_text(color = "black", size = base_size),
      strip.background = element_rect(fill = palette[["strip"]], color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.1, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.35),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.15, color = palette[["neutral"]], hjust = 0),
      plot.tag = element_text(face = "bold", size = 9.3),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      plot.margin = margin(3.0, 3.5, 3.0, 3.5, unit = "mm")
    )
}

# Panel a: patient deletion stability for all frozen programs.
forest <- copy(lopo_summary)
forest[, program := factor(program_labels[signature_id], levels = rev(unname(program_labels[programs])))]
forest[, endpoint := factor(contrast_labels[contrast], levels = unname(contrast_labels))]
forest[, role_label := role_labels[role]]
forest[, direction_label := paste0(same_direction_n, "/", n_estimable)]
forest[, label_x := 2.32]
forest <- forest[order(match(signature_id, programs), match(contrast, names(contrast_labels)))]

p_a <- ggplot(forest, aes(x = full_effect, y = program, color = role_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.38, color = palette[["neutral"]]) +
  geom_errorbar(aes(xmin = lopo_min, xmax = lopo_max), orientation = "y", width = 0.18, linewidth = 0.62) +
  geom_point(size = 2.25) +
  geom_text(aes(x = label_x, label = direction_label), hjust = 1, color = palette[["ink"]],
            family = font_family, size = 2.35) +
  facet_wrap(~endpoint, nrow = 1) +
  scale_color_manual(values = role_colors, guide = guide_legend(nrow = 1, override.aes = list(size = 2.4))) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  coord_cartesian(xlim = c(-1.62, 2.39), clip = "off") +
  labs(
    x = "Adverse-oriented Hedges' g (point: full cohort; line: LOPO range)", y = NULL,
    title = "Patient-deletion sensitivity",
    subtitle = "Fractions show sign-preserving deletions; LOPO ranges are not confidence intervals",
    tag = "a", color = NULL
  ) +
  theme_nature() +
  theme(
    axis.text.y = element_text(size = 6.55),
    strip.text = element_text(size = 6.9),
    legend.position = "bottom",
    panel.spacing.x = unit(0.38, "cm")
  )

# Panel b: two separately adjusted 18-test null families.
q_long <- rbindlist(list(
  results[, .(signature_id, program_label, role, contrast,
              null_type = "Exact patient-label null", p_value = exact_label_p, q_value = exact_label_q_bh18)],
  results[, .(signature_id, program_label, role, contrast,
              null_type = "Covariance-matched program null", p_value = covariance_program_p, q_value = covariance_program_q_bh18)]
))
q_long[, program := factor(program_labels[signature_id], levels = rev(unname(program_labels[programs])))]
q_long[, endpoint := factor(contrast_labels[contrast], levels = unname(contrast_labels))]
q_long[, null_type := factor(null_type, levels = names(null_colors))]
q_long[, null_column := factor(
  ifelse(null_type == "Exact patient-label null", "Patient-label null", "Program-identity null"),
  levels = c("Patient-label null", "Program-identity null")
)]
q_long[, evidence_strength := pmin(-log10(q_value) / -log10(0.009), 1)]
q_long[, fill_alpha := 0.16 + 0.84 * evidence_strength]
q_long[, q_label := sprintf("q = %.3f", q_value)]
q_long <- q_long[order(match(signature_id, programs), match(contrast, names(contrast_labels)), null_type)]

p_b <- ggplot(q_long, aes(x = null_column, y = program)) +
  geom_tile(aes(fill = null_type, alpha = fill_alpha), width = 0.92, height = 0.86,
            color = "white", linewidth = 0.65) +
  geom_tile(data = q_long[q_value < 0.05], fill = NA, color = palette[["ink"]],
            width = 0.92, height = 0.86, linewidth = 0.72) +
  geom_text(data = q_long[evidence_strength <= 0.72], aes(label = q_label),
            color = palette[["ink"]], family = font_family, fontface = "bold", size = 2.35) +
  geom_text(data = q_long[evidence_strength > 0.72], aes(label = q_label),
            color = "white", family = font_family, fontface = "bold", size = 2.35) +
  facet_wrap(~endpoint, nrow = 1) +
  scale_fill_manual(values = null_colors, guide = "none") +
  scale_alpha_identity() +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Two-null evidence matrix",
    subtitle = "BH-adjusted q is printed in each cell; darker fill indicates smaller q; black border marks q < 0.05",
    tag = "b"
  ) +
  theme_nature() +
  theme(
    axis.text.y = element_text(size = 6.55),
    axis.text.x.top = element_text(size = 6.55, face = "bold", margin = margin(b = 2.0)),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 6.9),
    legend.position = "none",
    panel.spacing.x = unit(0.38, "cm")
  )

figure3 <- p_a / p_b + plot_layout(heights = c(1.05, 1.00)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

base_file <- file.path(figures, "Figure3_lopo_two_null_calibration")
svglite::svglite(paste0(base_file, ".svg"), width = width_in, height = height_in)
print(figure3); dev.off()
grDevices::cairo_pdf(paste0(base_file, ".pdf"), width = width_in, height = height_in, family = font_family)
print(figure3); dev.off()
ragg::agg_tiff(paste0(base_file, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, compression = "lzw")
print(figure3); dev.off()
ragg::agg_png(paste0(base_file, ".png"), width = width_in, height = height_in,
              units = "in", res = 600)
print(figure3); dev.off()

# Figure source data and audit summary.
panel_a_source <- forest[, .(signature_id, program = as.character(program), role, role_label,
                             contrast, endpoint = as.character(endpoint), full_effect,
                             lopo_min, lopo_max, same_direction_n, n_estimable, direction_label,
                             exact_label_p, exact_label_q_bh18,
                             covariance_program_p, covariance_program_q_bh18)]
panel_b_source <- q_long[, .(signature_id, program = as.character(program), role,
                             contrast, endpoint = as.character(endpoint), null_type = as.character(null_type),
                             p_value, q_value, evidence_strength, fill_alpha,
                             threshold_pass = q_value < 0.05)]
fwrite(panel_a_source, file.path(tables, "Figure3_panel_a_lopo_source.tsv"), sep = "\t")
fwrite(panel_b_source, file.path(tables, "Figure3_panel_b_two_null_bh_source.tsv"), sep = "\t")
fwrite(results, file.path(tables, "Figure3_all18_results.tsv"), sep = "\t")

fmt <- function(x) sprintf("%.3f", x)
get_result <- function(sig, endpoint) results[signature_id == sig & contrast == endpoint]
primary <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
primary_q_range_exact <- range(results[signature_id %in% primary, exact_label_q_bh18])
primary_q_range_cov <- range(results[signature_id %in% primary, covariance_program_q_bh18])
hyp_u <- get_result("HALLMARK_HYPOXIA", "uptake_failure")
hyp_r <- get_result("HALLMARK_HYPOXIA", "response_failure_with_uptake")
infl_r <- get_result("HALLMARK_INFLAMMATORY_RESPONSE", "response_failure_with_uptake")

legend_text <- paste0(
  "## Figure 3 legend\n\n",
  "**Fig. 3 | Patient-deletion stability and two-null calibration of endpoint alignment.** ",
  "**a,** Leave-one-patient-out (LOPO) sensitivity for all nine frozen programs in GSE151179. Points show the full-cohort adverse-oriented Hedges' g; horizontal lines span the minimum and maximum estimate after deleting each eligible patient in turn. Fractions report sign-preserving deletions over estimable deletions (13 for uptake failure and 11 for response failure despite retained uptake). LOPO ranges are sensitivity summaries, not confidence intervals. ",
  "**b,** Evidence matrix for two prespecified null models. The exact patient-label null enumerated all 1,716 and 330 allocations for the uptake- and response-failure contrasts, respectively. The covariance-matched program-identity null compared each observed program with 1,000 label-blind programs matched on size, direction, marginal expression and co-expression summaries. Benjamini-Hochberg-adjusted q values are printed in the cells; adjustment was performed separately across the same 18 program-by-endpoint tests for each null family. Darker fill indicates smaller q, and a black border marks q < 0.05. The patient-label null asks whether a fixed program is associated with an endpoint under exchangeability; the program-identity null asks whether its observed effect is more extreme than effects from structurally comparable programs. Neither null establishes RAIR specificity, prediction or mechanism. Source data are provided as a Source Data file.\n"
)

methods_text <- paste0(
  "## Methods: patient-deletion stability and two-null calibration\n\n",
  "Patient-level sensitivity was assessed separately for each frozen program and endpoint by deleting one patient at a time and recalculating adverse-oriented Hedges' g without changing the program. The uptake-failure contrast yielded 13 LOPO estimates (six nonavid-persistent versus seven avid-persistent patients), and the response-failure-with-retained-uptake contrast yielded 11 (seven avid-persistent versus four avid-remission patients). Directional stability was the number of estimable deletions retaining the sign of the full-cohort effect. The LOPO minimum-to-maximum interval was treated as a patient-deletion sensitivity range rather than a confidence interval.\n\n",
  "Two prespecified null models were evaluated separately. For the endpoint-label null, group labels were enumerated exactly at the patient level while preserving group sizes, yielding choose(13,6) = 1,716 allocations for uptake failure and choose(11,7) = 330 for response failure despite retained uptake. Two-sided exact probabilities were the proportion of allocations with an absolute Hedges' g at least as large as the observed value; the observed allocation was included and no finite-sampling plus-one correction was applied. For the program-identity null, 1,000 unique label-blind programs were generated for each tested program before clinical labels were accessed. Null programs preserved gene number and direction counts, matched each gene position within the frozen joint expression-mean and expression-SD stratum, and were constrained at absolute tolerance 0.05 for both the PC1 variance proportion and mean pairwise Pearson correlation. The two-sided empirical probability was (1 + the number of null absolute effects at least as large as the observed absolute effect)/(1 + 1,000). Benjamini-Hochberg correction was applied separately across the same prespecified 18 program-by-endpoint comparisons for each null family.\n"
)

results_text <- paste0(
  "## Patient stability and program-identity extremeness do not establish endpoint association\n\n",
  "LOPO analysis distinguished sensitivity to individual patients from statistical evidence (Fig. 3a). For the three primary programs, sign-preserving deletion counts were 11/13 for TDS-16, 12/13 for iodide handling-11 and 10/13 for Condello-6 on the uptake-failure axis, and 10/11, 7/11 and 10/11, respectively, on the response-failure-with-retained-uptake axis. Thus, a small full-cohort effect could retain its direction under most deletions without showing evidence against either null model.\n\n",
  "No program was significant after exact patient-label permutation and correction across 18 tests (minimum BH-adjusted q = ",
  fmt(min(results$exact_label_q_bh18)), "; Fig. 3b). The three primary programs likewise showed no calibrated evidence under either the exact-label null (q range ", fmt(primary_q_range_exact[1]), "-", fmt(primary_q_range_exact[2]), ") or the covariance-matched program-identity null (q range ", fmt(primary_q_range_cov[1]), "-", fmt(primary_q_range_cov[2]), "). Under the latter null, Hypoxia on the uptake-failure axis was unusually extreme (P = ", fmt(hyp_u$covariance_program_p), ", q = ", fmt(hyp_u$covariance_program_q_bh18), "), whereas its response-axis result did not cross the adjusted threshold (P = ", fmt(hyp_r$covariance_program_p), ", q = ", fmt(hyp_r$covariance_program_q_bh18), "). Inflammatory response was also extreme on the response-failure axis (P = ", fmt(infl_r$covariance_program_p), ", q = ", fmt(infl_r$covariance_program_q_bh18), "). These program-identity results indicate coherence relative to structurally matched gene sets; because the corresponding exact patient-label tests were not significant after multiplicity correction, they do not establish robust endpoint association or RAIR specificity in this cohort.\n"
)
writeLines(c(legend_text, "", methods_text, "", results_text),
           file.path(manuscript, "Figure3_two_null_methods_results_legend.md"))

contract <- c(
  "# Figure 3 contract and freeze",
  "",
  "- Core conclusion: patient-deletion direction stability, endpoint-label association and program-identity extremeness are distinct properties.",
  "- Panel a: all nine frozen programs, two endpoint contrasts, full effect, LOPO range and sign-preserving counts.",
  "- Panel b: a four-column evidence matrix displaying exact patient-label and covariance-matched program-identity q values for both endpoints, with separate BH correction across 18 tests per null family.",
  "- Patient-label permutations: 1,716 for uptake failure and 330 for response failure despite uptake.",
  "- Program-identity nulls: 1,000 per program and endpoint; generated label-blind and frozen before endpoint scoring.",
  "- Boundary: neither null establishes a gold-standard RAIR construct, prediction, mechanism or causality.",
  "- Backend: R only.",
  paste0("- Final size: ", width_mm, " x ", height_mm, " mm."),
  "- Exports: editable SVG, vector PDF, 600 dpi TIFF and 600 dpi PNG."
)
writeLines(contract, file.path(figures, "Figure3_two_null_figure_contract.md"))

capture.output(sessionInfo(), file = file.path(manifests, "sessionInfo.txt"))
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
run_manifest <- list(
  passed = TRUE, backend = "R", width_mm = width_mm, height_mm = height_mm,
  programs = 9L, endpoints = 2L, tests_per_null_family = 18L,
  label_permutations = list(uptake_failure = 1716L, response_failure_with_uptake = 330L),
  covariance_null_sets_per_test = 1000L,
  exact_label_min_q = min(results$exact_label_q_bh18),
  covariance_significant_tests = q_long[null_type == "Covariance-matched program null" & q_value < 0.05,
                                       .(signature_id, contrast, q_value)],
  analysis_freeze_sha256 = sha256_file(file.path(analysis, "manifests", "CLINICAL_LAYER_FROZEN_SHA256.tsv")),
  validation_sha256 = sha256_file(file.path(validation, "clinical_layer_validation.json")),
  figure_script_sha256 = sha256_file(script_path)
)
write_json(run_manifest, file.path(manifests, "figure3_two_null_run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)

outputs <- sort(c(list.files(figures, full.names = TRUE), list.files(tables, full.names = TRUE),
                  list.files(manuscript, full.names = TRUE),
                  file.path(manifests, c("sessionInfo.txt", "figure3_two_null_run_manifest.json"))))
fwrite(data.table(file = substring(outputs, nchar(out) + 2L),
                  sha256 = vapply(outputs, sha256_file, character(1))),
       file.path(manifests, "FIGURE3_TWO_NULL_FROZEN_SHA256.tsv"), sep = "\t")
writeLines(paste0(sha256_file(file.path(manifests, "FIGURE3_TWO_NULL_FROZEN_SHA256.tsv")),
                  "  FIGURE3_TWO_NULL_FROZEN_SHA256.tsv"),
           file.path(manifests, "FIGURE3_TWO_NULL_FROZEN_SHA256.tsv.sha256"))

cat("Figure 3 two-null rebuild complete.\n")
print(results[, .(program_label, contrast, adverse_aligned_hedges_g,
                  exact_label_p, exact_label_q_bh18,
                  covariance_program_p, covariance_program_q_bh18)])
