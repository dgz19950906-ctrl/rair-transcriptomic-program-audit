#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(readr); library(scales); library(svglite); library(ragg); library(grid)
})

project_dir <- "/home/dony/ThyroidCancer_Project/rair_audit"
old_dir <- file.path(project_dir, "bulk_audit", "gse299988")
new_dir <- file.path(project_dir, "covariance_null_gse299988_v1", "clinical_label_layer_v1", "tables")
results_dir <- file.path(project_dir, "results")
out_dir <- file.path(project_dir, "figures", "figure4_two_null_v3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

challenge_old <- read_tsv(file.path(old_dir, "GSE299988_confounded_challenge.tsv"), show_col_types = FALSE)
two_null <- read_tsv(file.path(new_dir, "GSE299988_two_null_challenge.tsv"), show_col_types = FALSE)
contingency <- read_tsv(file.path(old_dir, "GSE299988_RAI_LN_contingency.tsv"), show_col_types = FALSE)
identifiability <- read_tsv(file.path(old_dir, "GSE299988_RAI_LN_identifiability.tsv"), show_col_types = FALSE)
gse151179 <- read_tsv(file.path(results_dir, "10_figure2_bootstrap_hedges_g.tsv"), show_col_types = FALSE)

program_ids <- c(
  "TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
)
program_labels <- c(
  TDS_16 = "TDS-16", IODIDE_HANDLING_11 = "Iodide handling-11",
  CONDELLO_2025_SIX = "Condello-6",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_HYPOXIA = "Hypoxia", HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint", HALLMARK_E2F_TARGETS = "E2F targets",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
)
display_order <- unname(program_labels[program_ids])
role_labels <- c(
  biological_reference = "Biological reference",
  independent_ex_vivo_avidity_candidate = "Ex vivo candidate",
  aggressiveness_negative_control = "Generic control"
)
palette <- c(
  "Biological reference" = "#4C78A8", "Ex vivo candidate" = "#D9824B",
  "Generic control" = "#858A91", "GSE151179\nuptake axis" = "#4C78A8",
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
      plot.tag.position = c(0, 1), legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5), panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_set(theme_nature())

# a, deterministic RAI/LN correspondence.
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
  labs(title = "Deterministic cohort structure",
       subtitle = sprintf("phi = %.0f; Fisher exact P = %.5f", identifiability$phi[1], identifiability$p_value[1]),
       x = NULL, y = NULL) +
  annotate("text", x = 1.5, y = 2.72, label = "Design-level non-identifiability",
           colour = "#B23A3A", fontface = "bold", size = 2.8, family = "Arial") +
  annotate("text", x = 1.5, y = 2.99,
           label = "Independent RAI and LN contributions cannot be estimated.",
           colour = "#555555", size = 2.15, family = "Arial") +
  scale_y_discrete(expand = expansion(add = c(0.35, 1.25))) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.ticks = element_blank(),
        axis.line = element_blank(), plot.margin = margin(4, 7, 4, 5))

# b, frozen observed-effect layer with updated two-null columns joined for traceability.
effect_df <- challenge_old %>%
  select(signature_id, role, adverse_aligned_hedges_g,
         adverse_aligned_bootstrap_g_ci_low, adverse_aligned_bootstrap_g_ci_high) %>%
  left_join(two_null %>% select(signature_id, exact_label_p, exact_label_q_bh9,
                                covariance_program_p, covariance_program_q_bh9), by = "signature_id") %>%
  mutate(program = factor(unname(program_labels[signature_id]), levels = rev(display_order)),
         role_label = unname(role_labels[role]))
p_b <- ggplot(effect_df, aes(adverse_aligned_hedges_g, program, colour = role_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = adverse_aligned_bootstrap_g_ci_low,
                    xmax = adverse_aligned_bootstrap_g_ci_high),
                orientation = "y", width = 0, linewidth = 0.55, alpha = 0.9) +
  geom_point(size = 2.15) +
  scale_colour_manual(values = palette[c("Biological reference", "Ex vivo candidate", "Generic control")]) +
  scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4.1, 3.3)) +
  labs(title = "Apparent effects in the collinear challenge",
       subtitle = "RAI-nonavid/LN-positive minus RAI-avid/LN-negative (n = 5 per group)",
       x = "Adverse-aligned Hedges' g (frozen bootstrap 95% CI)", y = NULL, colour = NULL) +
  theme(legend.position = "top", legend.justification = "left")

# c, two distinct null models; q values are descriptive under non-identifiability.
fmt_q <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
q_long <- two_null %>%
  transmute(signature_id, program = unname(program_labels[signature_id]),
            `Exact patient labels` = exact_label_q_bh9,
            `Covariance-matched programs` = covariance_program_q_bh9) %>%
  pivot_longer(c(`Exact patient labels`, `Covariance-matched programs`),
               names_to = "null_model", values_to = "q") %>%
  mutate(program = factor(program, levels = rev(display_order)),
         null_model = factor(null_model, levels = c("Exact patient labels", "Covariance-matched programs")),
         q_label = fmt_q(q), signal = -log10(pmax(q, 1e-4)))
p_c <- ggplot(q_long, aes(null_model, program)) +
  geom_tile(aes(fill = signal), colour = "white", linewidth = 0.8) +
  geom_tile(data = q_long %>% filter(q < 0.05), fill = NA, colour = "#111111", linewidth = 0.65) +
  geom_text(aes(label = q_label), family = "Arial", size = 2.25, colour = "#202020") +
  scale_fill_gradient(low = "#F1F3F4", high = "#8EB4CF", limits = c(0, 1.35),
                      oob = squish, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(title = "Two-null calibration",
       subtitle = "BH9 q; descriptive under non-identifiability",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 5.7, face = "bold", lineheight = 0.9),
        axis.text.y = element_text(size = 5.8), axis.ticks = element_blank(), axis.line = element_blank())

# d, identical primary definitions across a separated and a collinear structure.
primary_ids <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
clean_df <- gse151179 %>%
  filter(signature_id %in% primary_ids, contrast == "uptake_failure") %>%
  transmute(signature_id, program, dataset = "GSE151179\nuptake",
            estimate = adverse_aligned_hedges_g, lower = ci_low, upper = ci_high)
challenge_primary <- two_null %>%
  filter(signature_id %in% primary_ids) %>%
  transmute(signature_id, program = unname(program_labels[signature_id]),
            dataset = "GSE299988\nRAI/LN", estimate = adverse_aligned_hedges_g,
            lower = frozen_bootstrap_ci_low, upper = frozen_bootstrap_ci_high)
shift_df <- bind_rows(clean_df, challenge_primary) %>%
  mutate(dataset = factor(dataset, levels = c("GSE151179\nuptake", "GSE299988\nRAI/LN")),
         program = factor(program, levels = c("TDS-16", "Iodide handling-11", "Condello-6")))
p_d <- ggplot(shift_df, aes(dataset, estimate, group = program)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_line(colour = "#A8A8A8", linewidth = 0.55) +
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = dataset), width = 0.055, linewidth = 0.5) +
  geom_point(aes(colour = dataset), size = 2.1) +
  geom_text(data = shift_df %>% filter(program == "Condello-6"),
            aes(label = sprintf("%.2f", estimate)), nudge_x = c(-0.09, 0.09),
            nudge_y = c(-0.22, 0.22), family = "Arial", size = 2.25,
            fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ program, nrow = 1) +
  scale_colour_manual(values = c("GSE151179\nuptake" = "#4C78A8", "GSE299988\nRAI/LN" = "#C65353"),
                      guide = "none") +
  scale_y_continuous(limits = c(-3.25, 3.5), breaks = seq(-3, 3, 1)) +
  labs(title = "Context-dependent shifts",
       subtitle = "Same definitions; no cross-cohort P comparison",
       x = NULL, y = "Adverse-aligned Hedges' g\n(bootstrap 95% CI)") +
  theme(axis.text.x = element_text(size = 5.2, lineheight = 0.9), strip.text = element_text(size = 6),
        panel.spacing.x = unit(6, "mm"), plot.margin = margin(5, 5, 5, 10))

figure <- ((p_a | p_b) + plot_layout(widths = c(0.82, 1.68))) /
  ((p_c | p_d) + plot_layout(widths = c(0.98, 1.52))) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 8, face = "bold"))

base <- file.path(out_dir, "Figure4_GSE299988_two_null_challenge_v3")
width_in <- 183 / 25.4; height_in <- 142 / 25.4
svglite(paste0(base, ".svg"), width = width_in, height = height_in); print(figure); dev.off()
cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Arial"); print(figure); dev.off()
agg_tiff(paste0(base, ".tiff"), width = width_in, height = height_in, units = "in", res = 600,
         compression = "lzw"); print(figure); dev.off()
agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300);
print(figure); dev.off()

write_tsv(ct_long, file.path(out_dir, "Figure4_panel_a_contingency_source.tsv"))
write_tsv(effect_df, file.path(out_dir, "Figure4_panel_b_effects_source.tsv"))
write_tsv(q_long, file.path(out_dir, "Figure4_panel_c_two_null_q_source.tsv"))
write_tsv(shift_df, file.path(out_dir, "Figure4_panel_d_shift_source.tsv"))
contract <- c(
  "Core conclusion: deterministic RAI/LN correspondence makes independent effects non-identifiable; large apparent shifts do not produce calibrated evidence under either declared null.",
  "Figure archetype: quantitative grid.", "Backend: R only.",
  "Target/output: double-column 183 x 142 mm; editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG.",
  "Panels: a design identifiability; b frozen apparent effects; c two-null BH9 calibration; d cross-cohort point-estimate shifts.",
  "Reviewer risk: q values are descriptive within a non-identifiable challenge and must not be interpreted as RAI-specific evidence."
)
writeLines(contract, file.path(out_dir, "Figure4_two_null_figure_contract.txt"))
message("Figure 4 v2 written to: ", base)
