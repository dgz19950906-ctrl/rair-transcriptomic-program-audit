#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(svglite)
  library(ragg)
})

project_dir <- "/home/dony/ThyroidCancer_Project/rair_audit"
input_dir <- file.path(project_dir, "bulk_audit", "gse299988")
output_dir <- file.path(project_dir, "figures", "final")
results_dir <- file.path(project_dir, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

challenge <- read_tsv(file.path(input_dir, "GSE299988_confounded_challenge.tsv"), show_col_types = FALSE)
contingency <- read_tsv(file.path(input_dir, "GSE299988_RAI_LN_contingency.tsv"), show_col_types = FALSE)
identifiability <- read_tsv(file.path(input_dir, "GSE299988_RAI_LN_identifiability.tsv"), show_col_types = FALSE)
gse151179 <- read_tsv(file.path(results_dir, "10_figure2_bootstrap_hedges_g.tsv"), show_col_types = FALSE)

program_labels <- c(
  TDS_16 = "TDS-16",
  IODIDE_HANDLING_11 = "Iodide handling-11",
  CONDELLO_2025_SIX = "Condello-6",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
)

role_labels <- c(
  biological_reference = "Biological reference",
  independent_ex_vivo_avidity_candidate = "Ex vivo candidate",
  aggressiveness_negative_control = "Generic control"
)

palette <- c(
  "Biological reference" = "#4C78A8",
  "Ex vivo candidate" = "#D9824B",
  "Generic control" = "#858A91",
  "GSE151179\nuptake axis" = "#4C78A8",
  "GSE299988\nRAI/LN challenge" = "#C65353"
)

theme_nature <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#252525"),
      axis.title = element_text(size = base_size, colour = "#252525"),
      axis.text = element_text(size = base_size - 0.4, colour = "#252525"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = "#252525"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold", colour = "#111111"),
      plot.tag.position = c(0, 1),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_set(theme_nature())

# Panel a: deterministic RAI/LN correspondence.
ct_long <- contingency %>%
  pivot_longer(c(negative, positive), names_to = "LN", values_to = "n") %>%
  mutate(
    RAI = recode(RAI, avid = "RAI-avid", nonavid = "RAI-nonavid"),
    LN = recode(LN, negative = "LN-negative", positive = "LN-positive"),
    RAI = factor(RAI, levels = c("RAI-avid", "RAI-nonavid")),
    LN = factor(LN, levels = c("LN-negative", "LN-positive"))
  )

p_a <- ggplot(ct_long, aes(LN, RAI, fill = n)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = n), size = 4.1, family = "Arial", fontface = "bold", colour = "#202020") +
  scale_fill_gradient(low = "#F2F2F2", high = "#A9C4D9", limits = c(0, 5), guide = "none") +
  coord_equal(clip = "off") +
  labs(
    title = "Deterministic cohort structure",
    subtitle = sprintf("phi = %.0f; Fisher exact P = %.5f", identifiability$phi[1], identifiability$p_value[1]),
    x = NULL, y = NULL
  ) +
  annotate("text", x = 1.5, y = 2.73, label = "Design-level non-identifiability",
           colour = "#B23A3A", fontface = "bold", size = 2.9, family = "Arial") +
  annotate("text", x = 1.5, y = 2.99,
           label = "P tests independence; it cannot separate RAI and LN contributions.",
           colour = "#555555", size = 2.2, family = "Arial") +
  scale_y_discrete(expand = expansion(add = c(0.35, 1.25))) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.ticks = element_blank(), axis.line = element_blank(),
    plot.margin = margin(4, 7, 4, 5)
  )

# Panel b: challenge-cohort apparent effects.
effect_df <- challenge %>%
  mutate(
    program = unname(program_labels[signature_id]),
    role_label = unname(role_labels[role]),
    program = factor(program, levels = rev(unname(program_labels)))
  )

p_b <- ggplot(effect_df, aes(adverse_aligned_hedges_g, program, colour = role_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = adverse_aligned_bootstrap_g_ci_low,
                    xmax = adverse_aligned_bootstrap_g_ci_high),
                orientation = "y", width = 0, linewidth = 0.55, alpha = 0.9) +
  geom_point(size = 2.15) +
  scale_colour_manual(values = palette[c("Biological reference", "Ex vivo candidate", "Generic control")]) +
  scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4.1, 3.3)) +
  labs(
    title = "Apparent effects in the collinear challenge",
    subtitle = "RAI-nonavid/LN-positive minus RAI-avid/LN-negative (n = 5 per group)",
    x = "Adverse-aligned Hedges' g (bootstrap 95% CI)", y = NULL, colour = NULL
  ) +
  theme(legend.position = "top", legend.justification = "left")

# Panel c: prespecified increment relative to strongest observed generic control.
increment_df <- effect_df %>%
  filter(role != "aggressiveness_negative_control") %>%
  mutate(program = factor(as.character(program), levels = rev(c("TDS-16", "Iodide handling-11", "Condello-6"))))

p_c <- ggplot(increment_df, aes(specificity_margin_vs_strongest_generic, program, colour = role_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = specificity_increment_bootstrap_low,
                    xmax = specificity_increment_bootstrap_high),
                orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.15) +
  scale_colour_manual(values = palette[c("Biological reference", "Ex vivo candidate")], guide = "none") +
  scale_x_continuous(limits = c(-5.3, 2.8), breaks = seq(-5, 2, 1)) +
  labs(
    title = "No robust specificity increment",
    subtitle = "Primary program minus strongest observed generic control",
    x = "Specificity increment (bootstrap 95% CI)", y = NULL
  )

# Panel d: same programs across a clinically separated and a collinear structure.
primary_ids <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
clean_df <- gse151179 %>%
  filter(signature_id %in% primary_ids, contrast == "uptake_failure") %>%
  transmute(
    signature_id, program,
    dataset = "GSE151179\nuptake axis",
    estimate = adverse_aligned_hedges_g, lower = ci_low, upper = ci_high
  )

challenge_primary <- challenge %>%
  filter(signature_id %in% primary_ids) %>%
  transmute(
    signature_id,
    program = unname(program_labels[signature_id]),
    dataset = "GSE299988\nRAI/LN challenge",
    estimate = adverse_aligned_hedges_g,
    lower = adverse_aligned_bootstrap_g_ci_low,
    upper = adverse_aligned_bootstrap_g_ci_high
  )

shift_df <- bind_rows(clean_df, challenge_primary) %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE151179\nuptake axis", "GSE299988\nRAI/LN challenge")),
    program = factor(program, levels = c("TDS-16", "Iodide handling-11", "Condello-6"))
  )

p_d <- ggplot(shift_df, aes(dataset, estimate, group = program)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_line(colour = "#A8A8A8", linewidth = 0.55) +
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = dataset), width = 0.055, linewidth = 0.5) +
  geom_point(aes(colour = dataset), size = 2.1) +
  geom_text(
    data = shift_df %>% filter(program == "Condello-6"),
    aes(label = sprintf("%.2f", estimate)),
    nudge_x = c(-0.09, 0.09), nudge_y = c(-0.22, 0.22),
    family = "Arial", size = 2.25, fontface = "bold", show.legend = FALSE
  ) +
  facet_wrap(~ program, nrow = 1) +
  scale_colour_manual(values = palette[c("GSE151179\nuptake axis", "GSE299988\nRAI/LN challenge")], guide = "none") +
  scale_y_continuous(limits = c(-3.25, 3.5), breaks = seq(-3, 3, 1)) +
  labs(
    title = "Dataset-structure-dependent shifts",
    subtitle = "Identical frozen program definitions; adverse-aligned effects shown numerically for Condello-6",
    x = NULL, y = "Adverse-aligned Hedges' g\n(bootstrap 95% CI)"
  ) +
  theme(
    axis.text.x = element_text(size = 5.2, lineheight = 0.9),
    strip.text = element_text(size = 6),
    panel.spacing.x = unit(4, "mm")
  )

figure <- ((p_a | p_b) + plot_layout(widths = c(0.82, 1.68))) /
  ((p_c | p_d) + plot_layout(widths = c(0.9, 1.6))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

base <- file.path(output_dir, "Figure4_GSE299988_confounded_challenge")
width_in <- 183 / 25.4
height_in <- 142 / 25.4

svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()

cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(figure)
dev.off()

agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw")
print(figure)
dev.off()

agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure)
dev.off()

write_tsv(ct_long, file.path(results_dir, "11_figure4_contingency_source.tsv"))
write_tsv(effect_df %>% select(signature_id, program, role_label, adverse_aligned_hedges_g,
                               adverse_aligned_bootstrap_g_ci_low, adverse_aligned_bootstrap_g_ci_high,
                               exact_permutation_p, empirical_random_p, lopo_direction_stability),
          file.path(results_dir, "11_figure4_challenge_effects_source.tsv"))
write_tsv(increment_df %>% select(signature_id, program, specificity_margin_vs_strongest_generic,
                                  specificity_increment_bootstrap_low, specificity_increment_bootstrap_high,
                                  bootstrap_probability_increment_le_zero, incremental_evidence_status),
          file.path(results_dir, "11_figure4_specificity_increment_source.tsv"))
write_tsv(shift_df, file.path(results_dir, "11_figure4_cross_cohort_shift_source.tsv"))

message("Figure 4 written to: ", base)
