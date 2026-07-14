#!/usr/bin/env Rscript

# Analysis: GSE151179 three-group, two-endpoint program-alignment figure
# Date: 2026-07-14
# Random seed: 42
# Unit: patient (one pretreatment primary tumor per patient)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(scales)
  library(grid)
})

set.seed(42)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript build_figure2_nature.R <rair_audit_dir>")
root <- normalizePath(args[[1]], mustWork = TRUE)
input_dir <- file.path(root, "bulk_audit")
result_dir <- file.path(root, "bootstrap_unified_v2", "figure2", "source_data")
figure_dir <- file.path(root, "figures", "figure2_unified_bootstrap_v2")
manuscript_dir <- file.path(root, "bootstrap_unified_v2", "figure2", "manuscript_text")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Nimbus Sans"
width_mm <- 183
height_mm <- 155
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

palette <- c(
  ink = "#222222",
  neutral = "#6B7280",
  neutral_light = "#D9DEE4",
  remission = "#4C78A8",
  avid_persistent = "#E3A33B",
  nonavid_persistent = "#C65D57",
  biological = "#377EB8",
  condello = "#D97745",
  generic = "#7A7F87",
  blue_light = "#E3EEF6",
  orange_light = "#F8EAD2",
  red_light = "#F3DEDC",
  grey_fill = "#F5F6F7"
)

theme_nature <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.text = element_text(color = "black", size = base_size - 0.3),
      axis.title = element_text(color = "black", size = base_size),
      strip.background = element_rect(fill = "#F1F2F3", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.1, color = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      plot.title = element_text(face = "bold", size = base_size + 0.7, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.1, color = palette[["neutral"]], hjust = 0),
      plot.tag = element_text(face = "bold", size = 9),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

samples <- read_tsv(file.path(input_dir, "GSE151179_primary_preRAI_samples.tsv"), show_col_types = FALSE)
scores <- read_tsv(file.path(input_dir, "frozen_signature_sample_scores.tsv"), show_col_types = FALSE)
frozen_effects <- read_tsv(file.path(input_dir, "frozen_signature_endpoint_alignment.tsv"), show_col_types = FALSE)
canonical_effects <- read_tsv(
  file.path(root, "covariance_null_v2", "clinical_label_layer_v1", "tables",
            "clinical_endpoint_two_null_results.tsv"),
  show_col_types = FALSE
)

stopifnot(nrow(samples) == 17L)
stopifnot(n_distinct(samples$patient_id) == 17L)
stopifnot(all(samples$eligible_primary_pre_rai))
stopifnot(all(samples$tissue_type == "Primary tumor"))
stopifnot(all(samples$collection_before_after_rai == "Before"))

group_levels <- c("RAI_avid_remission", "RAI_avid_persistent", "RAI_nonavid_persistent")
group_labels <- c(
  RAI_avid_remission = "RAI-avid / remission",
  RAI_avid_persistent = "RAI-avid / persistent",
  RAI_nonavid_persistent = "RAI-nonavid / persistent"
)
group_short <- c(
  RAI_avid_remission = "Avid\nremission\n(n = 4)",
  RAI_avid_persistent = "Avid\npersistent\n(n = 7)",
  RAI_nonavid_persistent = "Nonavid\npersistent\n(n = 6)"
)
group_colors <- c(
  RAI_avid_remission = palette[["remission"]],
  RAI_avid_persistent = palette[["avid_persistent"]],
  RAI_nonavid_persistent = palette[["nonavid_persistent"]]
)

observed_counts <- table(factor(samples$analysis_group, levels = group_levels))
stopifnot(identical(as.integer(observed_counts), c(4L, 7L, 6L)))

primary_programs <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
program_labels <- c(
  TDS_16 = "TDS-16",
  IODIDE_HANDLING_11 = "Iodide handling-11",
  CONDELLO_2025_SIX = "Condello-6",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_G2M_CHECKPOINT = "G2M",
  HALLMARK_E2F_TARGETS = "E2F",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory"
)

orientation_multiplier <- function(orientation) {
  ifelse(orientation %in% c("higher_more_differentiated", "higher_more_iodide_handling"), -1, 1)
}

hedges_g <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  df <- nx + ny - 2
  if (nx < 2 || ny < 2 || df <= 0) return(NA_real_)
  pooled_sd <- sqrt(((nx - 1) * stats::var(x) + (ny - 1) * stats::var(y)) / df)
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
  correction <- 1 - 3 / (4 * df - 1)
  correction * (mean(x) - mean(y)) / pooled_sd
}

contrast_definitions <- tibble(
  contrast = c("uptake_failure", "response_failure_with_uptake"),
  adverse_group = c("RAI_nonavid_persistent", "RAI_avid_persistent"),
  reference_group = c("RAI_avid_persistent", "RAI_avid_remission"),
  axis_label = c(
    "Uptake failure\n(nonavid persistent - avid persistent)",
    "Response failure despite uptake\n(avid persistent - avid remission)"
  )
)

score_meta <- scores %>%
  left_join(samples %>% select(geo_accession, patient_id, analysis_group), by = "geo_accession") %>%
  mutate(
    multiplier = orientation_multiplier(orientation),
    adverse_score = multiplier * score,
    program = unname(program_labels[signature_id])
  )

if (any(is.na(score_meta$analysis_group))) stop("Sample-score join produced missing clinical groups")

effects <- canonical_effects %>%
  transmute(
    signature_id, contrast, adverse_group, reference_group,
    adverse_aligned_hedges_g,
    ci_low = bootstrap_ci_low,
    ci_high = bootstrap_ci_high,
    bootstrap_finite = bootstrap_resamples,
    n_adverse, n_reference, role, orientation,
    exact_permutation_p = exact_label_p,
    p_bh_18 = exact_label_q_bh18
  ) %>%
  left_join(
    frozen_effects %>% select(
      signature_id, contrast,
      frozen_adverse_aligned_hedges_g = adverse_aligned_hedges_g
    ),
    by = c("signature_id", "contrast")
  ) %>%
  mutate(
    point_difference_from_frozen = adverse_aligned_hedges_g - frozen_adverse_aligned_hedges_g,
    program = unname(program_labels[signature_id])
  )

if (nrow(effects) != 18L) stop("Canonical two-null table did not contain 18 endpoint rows")
if (any(abs(effects$point_difference_from_frozen) > 1e-10)) {
  stop("Canonical Hedges g does not reproduce the frozen point estimates")
}

write_tsv(effects, file.path(result_dir, "10_figure2_bootstrap_hedges_g.tsv"))
write_tsv(
  score_meta %>%
    select(geo_accession, patient_id, analysis_group, signature_id, program, role,
           orientation, score, adverse_score),
  file.path(result_dir, "10_figure2_patient_program_scores.tsv")
)
write_tsv(
  samples %>% count(analysis_group, name = "n") %>%
    mutate(group_label = unname(group_labels[analysis_group])),
  file.path(result_dir, "10_figure2_group_counts.tsv")
)

# Panel a: three mutually exclusive groups and the two prespecified contrasts.
group_boxes <- tibble(
  group = factor(group_levels, levels = group_levels),
  x = c(18, 50, 82),
  xmin = c(6, 38, 70),
  xmax = c(30, 62, 94),
  ymin = 9,
  ymax = 27,
  label = unname(group_short[group_levels]),
  fill = c(palette[["blue_light"]], palette[["orange_light"]], palette[["red_light"]]),
  outline = unname(group_colors[group_levels])
)

p_a <- ggplot() +
  geom_rect(
    data = group_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = outline),
    linewidth = 0.7
  ) +
  scale_fill_identity() + scale_color_identity() +
  geom_text(
    data = group_boxes,
    aes(x = x, y = 18, label = label),
    family = font_family, size = 3.0, fontface = "bold", lineheight = 0.95
  ) +
  geom_segment(
    aes(x = 30.5, xend = 37.5, y = 18, yend = 18),
    color = palette[["neutral"]], linewidth = 0.55,
    arrow = arrow(length = unit(2.2, "mm"), type = "closed")
  ) +
  geom_segment(
    aes(x = 62.5, xend = 69.5, y = 18, yend = 18),
    color = palette[["neutral"]], linewidth = 0.55,
    arrow = arrow(length = unit(2.2, "mm"), type = "closed")
  ) +
  annotate(
    "label", x = 34, y = 31,
    label = "Response failure despite uptake\n7 persistent vs 4 remission",
    family = font_family, size = 2.35, lineheight = 0.95,
    fill = "white", color = palette[["neutral"]], linewidth = 0.25
  ) +
  geom_segment(
    aes(x = 34, xend = 34, y = 27.8, yend = 19.8),
    color = palette[["neutral"]], linewidth = 0.35
  ) +
  annotate(
    "label", x = 66, y = 31,
    label = "Uptake failure\n6 nonavid vs 7 avid persistent",
    family = font_family, size = 2.35, lineheight = 0.95,
    fill = "white", color = palette[["neutral"]], linewidth = 0.25
  ) +
  geom_segment(
    aes(x = 66, xend = 66, y = 27.8, yend = 19.8),
    color = palette[["neutral"]], linewidth = 0.35
  ) +
  annotate(
    "text", x = 50, y = 3.5,
    label = "17 disjoint patients | pretreatment primary tumors | three mutually exclusive clinical groups",
    family = font_family, size = 2.25, color = palette[["neutral"]]
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 36), clip = "off") +
  labs(title = "Three clinical groups define two controlled endpoint contrasts", tag = "a") +
  theme_void(base_size = 7.2, base_family = font_family) +
  theme(
    plot.title = element_text(face = "bold", size = 7.9),
    plot.tag = element_text(face = "bold", size = 9),
    plot.tag.position = c(0, 1),
    plot.margin = margin(3, 3, 3, 3, unit = "mm")
  )

# Panel b: patient-level primary-program distributions in adverse orientation.
plot_scores <- score_meta %>%
  filter(signature_id %in% primary_programs) %>%
  group_by(signature_id) %>%
  mutate(adverse_score_z = as.numeric(scale(adverse_score))) %>%
  ungroup() %>%
  mutate(
    analysis_group = factor(analysis_group, levels = group_levels),
    program = factor(program, levels = unname(program_labels[primary_programs]))
  )

p_b <- ggplot(plot_scores, aes(x = analysis_group, y = adverse_score_z, color = analysis_group)) +
  geom_hline(yintercept = 0, color = "#D4D7DB", linewidth = 0.35) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, linewidth = 0.4,
    color = "#55585C", fill = "white"
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 42),
    size = 1.45, alpha = 0.85
  ) +
  facet_wrap(~program, ncol = 1, scales = "fixed") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_x_discrete(labels = c("Avid\nremission", "Avid\npersistent", "Nonavid\npersistent")) +
  labs(
    x = NULL,
    y = "Adverse-oriented standardized program score",
    title = "Patient-level primary-program distributions",
    subtitle = "Higher values indicate the program direction aligned with an adverse state",
    tag = "b"
  ) +
  theme_nature(7.2) +
  theme(
    axis.text.x = element_text(size = 6.2, lineheight = 0.9),
    strip.text = element_text(size = 6.6),
    panel.spacing.y = unit(0.15, "cm")
  )

# Panel c: two-dimensional adverse-oriented effect map with bootstrap CIs.
effect_wide <- effects %>%
  select(signature_id, program, role, contrast, adverse_aligned_hedges_g, ci_low, ci_high,
         exact_permutation_p, p_bh_18) %>%
  pivot_wider(
    names_from = contrast,
    values_from = c(adverse_aligned_hedges_g, ci_low, ci_high, exact_permutation_p, p_bh_18)
  ) %>%
  mutate(
    role_label = case_when(
      role == "biological_reference" ~ "Biological reference",
      role == "independent_ex_vivo_avidity_candidate" ~ "Ex vivo candidate",
      TRUE ~ "Generic control"
    ),
    point_color = case_when(
      role == "biological_reference" ~ palette[["biological"]],
      role == "independent_ex_vivo_avidity_candidate" ~ palette[["condello"]],
      TRUE ~ palette[["generic"]]
    ),
    label_dx = case_when(
      signature_id == "CONDELLO_2025_SIX" ~ -0.32,
      signature_id == "TDS_16" ~ -0.18,
      signature_id == "IODIDE_HANDLING_11" ~ 0.22,
      signature_id == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ 0.30,
      signature_id == "HALLMARK_HYPOXIA" ~ 0.04,
      signature_id == "HALLMARK_ANGIOGENESIS" ~ 0.33,
      signature_id == "HALLMARK_G2M_CHECKPOINT" ~ -0.27,
      signature_id == "HALLMARK_E2F_TARGETS" ~ -0.30,
      TRUE ~ -0.28
    ),
    label_dy = case_when(
      signature_id == "CONDELLO_2025_SIX" ~ 0.25,
      signature_id == "TDS_16" ~ -0.23,
      signature_id == "IODIDE_HANDLING_11" ~ 0.22,
      signature_id == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" ~ 0.20,
      signature_id == "HALLMARK_HYPOXIA" ~ -0.28,
      signature_id == "HALLMARK_ANGIOGENESIS" ~ 0.10,
      signature_id == "HALLMARK_G2M_CHECKPOINT" ~ -0.20,
      signature_id == "HALLMARK_E2F_TARGETS" ~ 0.20,
      TRUE ~ -0.25
    ),
    label_x = adverse_aligned_hedges_g_uptake_failure + label_dx,
    label_y = adverse_aligned_hedges_g_response_failure_with_uptake + label_dy,
    label_show = role != "aggressiveness_negative_control" |
      signature_id %in% c("HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS")
  )

write_tsv(effect_wide, file.path(result_dir, "10_figure2_two_axis_source.tsv"))

x_range <- range(c(effect_wide$ci_low_uptake_failure, effect_wide$ci_high_uptake_failure,
                   effect_wide$label_x, 0), na.rm = TRUE)
y_range <- range(c(effect_wide$ci_low_response_failure_with_uptake,
                   effect_wide$ci_high_response_failure_with_uptake,
                   effect_wide$label_y, 0), na.rm = TRUE)
x_pad <- max(0.18, diff(x_range) * 0.07)
y_pad <- max(0.18, diff(y_range) * 0.07)

p_c <- ggplot(effect_wide, aes(
  x = adverse_aligned_hedges_g_uptake_failure,
  y = adverse_aligned_hedges_g_response_failure_with_uptake
)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#F8F3EA", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#EEF4F8", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#F5F1F7", alpha = 0.45) +
  geom_hline(yintercept = 0, linetype = "dashed", color = palette[["neutral"]], linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = palette[["neutral"]], linewidth = 0.4) +
  geom_segment(
    aes(x = ci_low_uptake_failure, xend = ci_high_uptake_failure,
        y = adverse_aligned_hedges_g_response_failure_with_uptake,
        yend = adverse_aligned_hedges_g_response_failure_with_uptake,
        color = role_label),
    linewidth = 0.55, alpha = 0.80
  ) +
  geom_segment(
    aes(x = adverse_aligned_hedges_g_uptake_failure,
        xend = adverse_aligned_hedges_g_uptake_failure,
        y = ci_low_response_failure_with_uptake,
        yend = ci_high_response_failure_with_uptake,
        color = role_label),
    linewidth = 0.55, alpha = 0.80
  ) +
  geom_segment(
    data = effect_wide %>% filter(label_show),
    aes(xend = label_x, yend = label_y),
    color = "#B8BCC2", linewidth = 0.28
  ) +
  geom_point(aes(color = role_label), size = 2.3, stroke = 0.3) +
  geom_text(
    data = effect_wide %>% filter(label_show),
    aes(x = label_x, y = label_y, label = program, color = role_label),
    family = font_family, size = 2.35, fontface = "bold"
  ) +
  scale_color_manual(values = c(
    "Biological reference" = palette[["biological"]],
    "Ex vivo candidate" = palette[["condello"]],
    "Generic control" = palette[["generic"]]
  )) +
  coord_cartesian(
    xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
    ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
    clip = "off"
  ) +
  labs(
    x = "Uptake-failure alignment\nAdverse-oriented Hedges' g",
    y = "Response-failure-with-uptake alignment\nAdverse-oriented Hedges' g",
    title = "Programs align differently with the two endpoint axes",
    subtitle = "Crosshairs denote patient-bootstrap 95% confidence intervals",
    color = NULL,
    tag = "c"
  ) +
  theme_nature(7.2) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(0.45, "cm"),
    plot.subtitle = element_text(size = 6.7)
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.2, linewidth = 0.8), nrow = 1))

figure2 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.31, 0.69), widths = c(0.82, 1.45)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

base_file <- file.path(figure_dir, "Figure2_three_group_endpoint_alignment_unified_bootstrap_v2")

svglite::svglite(paste0(base_file, ".svg"), width = width_in, height = height_in)
print(figure2)
dev.off()

grDevices::cairo_pdf(paste0(base_file, ".pdf"), width = width_in, height = height_in, family = font_family)
print(figure2)
dev.off()

ragg::agg_tiff(
  paste0(base_file, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure2)
dev.off()

ragg::agg_png(
  paste0(base_file, ".png"), width = width_in, height = height_in,
  units = "in", res = 600
)
print(figure2)
dev.off()

fmt_effect <- function(sig, contrast) {
  z <- effects %>% filter(signature_id == sig, .data$contrast == !!contrast)
  sprintf("%.2f (95%% bootstrap CI %.2f to %.2f; exact P = %.3f)",
          z$adverse_aligned_hedges_g, z$ci_low, z$ci_high, z$exact_permutation_p)
}

legend <- paste0(
  "## Figure 2 legend\n\n",
  "**Fig. 2 | Three clinical groups reveal endpoint-dependent alignment of radioiodine-related transcriptomic programs.** ",
  "**a,** Seventeen patients with pretreatment primary tumors from GSE151179 were partitioned into three mutually exclusive groups: RAI-avid disease with remission (n = 4), RAI-avid disease with persistence (n = 7), and RAI-nonavid disease with persistence (n = 6). The response-failure-with-uptake contrast compares the two RAI-avid groups, whereas the uptake-failure contrast compares persistent tumors with and without uptake. ",
  "**b,** Adverse-oriented standardized patient-level scores for TDS-16, iodide handling-11, and Condello-6. Points denote patients; boxes show the median and interquartile range, and whiskers extend to 1.5 times the interquartile range. Scores were sign-oriented so that higher values indicate alignment with an adverse state. ",
  "**c,** Two-dimensional endpoint-alignment map for the three primary programs and six prespecified generic aggressiveness controls. Coordinates are adverse-oriented Hedges' g values; horizontal and vertical crosshairs are percentile 95% confidence intervals from 10,000 within-group patient-bootstrap resamples. Positive values indicate alignment with the adverse group on the corresponding axis. The three primary programs, Hypoxia and Angiogenesis are directly labelled for legibility; all program identities, exact permutation P values and Benjamini-Hochberg-adjusted values for the 18 comparisons are provided in Source Data. The small cohort supports an alignment audit rather than predictive or causal inference. Source data are provided as a Source Data file.\n\n"
)

methods <- paste0(
  "## Methods: three-group endpoint-alignment analysis\n\n",
  "The GSE151179 analysis set comprised 17 disjoint patients represented by one pretreatment primary-tumor transcriptome each. Clinical annotations defined RAI-avid disease with remission (n = 4), RAI-avid disease with persistence (n = 7), and RAI-nonavid disease with persistence (n = 6). Two contrasts were prespecified: uptake failure, defined as RAI-nonavid persistent minus RAI-avid persistent disease, and response failure despite uptake, defined as RAI-avid persistent minus RAI-avid remission. Thus, persistence was held constant in the uptake contrast and uptake was held constant in the response contrast.\n\n",
  "Frozen program scores were sign-oriented so that positive effects indicated alignment with the adverse group: differentiation and iodide-handling programs were multiplied by -1, whereas Condello-6 and generic aggressiveness programs retained their reported direction. Standardized mean differences were estimated as small-sample-corrected Hedges' g. Percentile 95% confidence intervals were obtained from 10,000 within-group patient-bootstrap resamples using a fixed random seed. Two-sided exact permutation P values were inherited from the frozen endpoint analysis and adjusted across the 18 program-by-endpoint comparisons using the Benjamini-Hochberg method. Inference emphasized effect direction and uncertainty; no classifier performance or causal effect was estimated.\n\n"
)

results_text <- paste0(
  "## Radioiodine-related programs show endpoint-dependent alignment in the three-group clinical frame\n\n",
  "The three-group GSE151179 design separated failure of clinical uptake from failure of structural response despite retained uptake (Fig. 2a). The uptake-failure comparison included six RAI-nonavid persistent and seven RAI-avid persistent tumors, whereas the response-failure-with-uptake comparison included seven RAI-avid persistent and four RAI-avid remission tumors. Patient-level distributions showed substantial within-group overlap for all three primary programs, consistent with the limited cohort size and arguing against a predictive interpretation (Fig. 2b).\n\n",
  "The primary programs occupied different regions of the two-axis effect map (Fig. 2c). TDS-16 showed weak positive alignment with uptake failure [", fmt_effect("TDS_16", "uptake_failure"), "] but weak inverse alignment with response failure despite uptake [", fmt_effect("TDS_16", "response_failure_with_uptake"), "]. Iodide handling-11 showed similarly weak uptake-axis alignment [", fmt_effect("IODIDE_HANDLING_11", "uptake_failure"), "] and was near zero on the response axis [", fmt_effect("IODIDE_HANDLING_11", "response_failure_with_uptake"), "]. Condello-6, which was trained against ex vivo iodine concentration, did not positively align with clinical uptake failure [", fmt_effect("CONDELLO_2025_SIX", "uptake_failure"), "] but showed modest positive alignment with failure despite retained uptake [", fmt_effect("CONDELLO_2025_SIX", "response_failure_with_uptake"), "].\n\n",
  "Generic aggressiveness controls also differed between axes. Hypoxia and angiogenesis showed the largest positive uptake-axis effects, whereas several generic programs were inversely aligned with the response-failure-with-uptake axis. Confidence intervals were wide and no point estimate was interpreted as a validated predictive effect. Instead, the divergent coordinates show that uptake failure and treatment failure despite retained uptake capture different transcriptomic associations within this cohort.\n\n"
)

writeLines(c(legend, methods, results_text), file.path(manuscript_dir, "Figure2_methods_results_legend.md"))

results_path <- file.path(manuscript_dir, "RESULTS_DRAFT.md")
if (file.exists(results_path)) {
  existing <- paste(readLines(results_path, warn = FALSE), collapse = "\n")
  block <- paste0(
    "<!-- FIGURE2_RESULTS_START -->\n", results_text,
    "<!-- FIGURE2_RESULTS_END -->\n\n",
    "<!-- Results sections corresponding to Figs. 3-4 will be inserted here after final figure freeze. -->"
  )
  if (grepl("<!-- FIGURE2_RESULTS_START -->", existing, fixed = TRUE)) {
    existing <- sub(
      "(?s)<!-- FIGURE2_RESULTS_START -->.*?<!-- Results sections corresponding to Figs\\. 3-4 will be inserted here after final figure freeze\\. -->",
      block, existing, perl = TRUE
    )
  } else {
    existing <- sub(
      "<!-- Results sections corresponding to Figs\\. 2-4 will be inserted here after final figure freeze\\. -->",
      block, existing, perl = TRUE
    )
  }
  writeLines(existing, results_path)
}

contract <- c(
  "# Figure 2 contract and manifest",
  "",
  "- Core conclusion: The same RAIR-associated programs show different alignment with uptake failure and response failure despite retained uptake; the two endpoint axes are not interchangeable.",
  "- Archetype: asymmetric quantitative figure.",
  "- Panel a: three mutually exclusive clinical groups and two prespecified contrasts.",
  "- Panel b: patient-level primary-program distributions.",
  "- Panel c: hero two-axis adverse-oriented Hedges' g map with patient-bootstrap 95% CIs.",
  "- Analysis unit: patient; one pretreatment primary tumor per patient.",
  "- Sample sizes: avid/remission n=4, avid/persistent n=7, nonavid/persistent n=6.",
  "- Statistics: small-sample-corrected Hedges' g, 10,000 within-group bootstrap resamples, exact permutation P, BH correction across 18 comparisons.",
  "- Boundary: alignment audit only; no prediction, causality, equivalence or endpoint specificity is claimed.",
  "- Backend: R only.",
  paste0("- Final size: ", width_mm, " x ", height_mm, " mm."),
  "- Exports: editable SVG, vector PDF, 600 dpi TIFF and 600 dpi PNG.",
  "- Source data: results/10_figure2_*.tsv.",
  "- Figure 3 retains LOPO stability and matched-random null distributions; these are not duplicated here."
)
writeLines(contract, file.path(figure_dir, "Figure2_figure_contract.md"))

analysis_manifest <- paste(
  "# Analysis outputs",
  "",
  "Generated: 2026-07-14",
  "Study type: three-group observational transcriptomic endpoint-alignment audit",
  "",
  "## Tables",
  "- `10_figure2_bootstrap_hedges_g.tsv` -- adverse-oriented Hedges' g with bootstrap 95% CIs and exact/BH-adjusted P values.",
  "- `10_figure2_two_axis_source.tsv` -- wide-format source data for the two-axis panel.",
  "- `10_figure2_patient_program_scores.tsv` -- patient-level frozen scores and adverse orientation.",
  "- `10_figure2_group_counts.tsv` -- reconciled mutually exclusive group counts.",
  "",
  "## Figures",
  "- `Figure2_three_group_endpoint_alignment_unified_bootstrap_v2.svg/.pdf/.tiff/.png` -- manuscript Figure 2 candidate using the canonical bootstrap source.",
  "",
  "## Reproducibility",
  "- Random seed: 42.",
  "- Bootstrap resamples: 10,000 within clinical groups.",
  "- Bootstrap source: canonical two-null clinical-label layer; the same confidence intervals are reused in Figures 2 and 4.",
  "- Point estimates asserted against the previously frozen endpoint-alignment table.",
  sep = "\n"
)
writeLines(analysis_manifest, file.path(result_dir, "_figure2_analysis_outputs.md"))

cat("Figure 2 generated.\n")
cat("Groups:", paste(names(observed_counts), observed_counts, collapse = "; "), "\n")
cat("Primary effects:\n")
print(effects %>% filter(signature_id %in% primary_programs) %>%
        select(signature_id, contrast, adverse_aligned_hedges_g, ci_low, ci_high,
               exact_permutation_p, p_bh_18, n_adverse, n_reference))
