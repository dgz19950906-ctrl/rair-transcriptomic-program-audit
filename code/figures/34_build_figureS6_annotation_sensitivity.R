#!/usr/bin/env Rscript

# Figure: annotation-circularity sensitivity for Figure 5
# Date: 2026-07-14
# Random seed: 42
# Backend: R only

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(42)
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 34_build_figureS6_annotation_sensitivity.R <rair_audit_dir>")
audit <- normalizePath(args[[1]], mustWork = TRUE)
sensitivity_dir <- file.path(audit, "singlecell_annotation_sensitivity_v1")
out_dir <- file.path(audit, "figures", "figureS6_annotation_sensitivity_v2")
if (dir.exists(out_dir)) stop("Refusing to overwrite existing output directory: ", out_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

primary_ids <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
program_labels <- c(
  TDS_16 = "TDS-16",
  IODIDE_HANDLING_11 = "Iodide handling-11",
  CONDELLO_2025_SIX = "Condello-6"
)
method_labels <- c(
  primary_hybrid = "Primary hybrid",
  leave_program_genes_out_hybrid = "Leave-program-genes-out",
  SingleR_only = "SingleR only"
)
method_colors <- c(
  `Primary hybrid` = "#30343B",
  `Leave-program-genes-out` = "#4C78A8",
  `SingleR only` = "#D9824B"
)
atlas_order <- c("GSE184362", "GSE191288", "GSE281736")

theme_paper <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Nimbus Sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222"),
      strip.background = element_rect(fill = "#F2F3F4", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = 8),
      plot.tag.position = c(0, 1),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

original_effects <- read_tsv(
  file.path(audit, "results", "08_three_atlas_bootstrap_effects.tsv"),
  show_col_types = FALSE
) %>%
  filter(signature_id %in% primary_ids) %>%
  transmute(
    dataset, signature_id, annotation_method = "primary_hybrid",
    target_compartment, paired_donors = n_donors,
    estimate, bootstrap_ci_low = conf_low, bootstrap_ci_high = conf_high
  )

external_top <- read_tsv(
  file.path(audit, "results", "07_external_atlas_preferences.tsv"),
  show_col_types = FALSE
) %>%
  filter(signature_id %in% primary_ids) %>%
  select(dataset, signature_id, target_is_top_compartment)
original_top <- bind_rows(
  tibble(dataset = "GSE184362", signature_id = primary_ids, target_is_top_compartment = TRUE),
  external_top
)
original_effects <- original_effects %>%
  left_join(original_top, by = c("dataset", "signature_id"))

sensitivity_effects <- read_tsv(
  file.path(sensitivity_dir, "annotation_sensitivity_effects.tsv"),
  show_col_types = FALSE
) %>%
  filter(signature_id %in% primary_ids) %>%
  select(
    dataset, signature_id, annotation_method, target_compartment,
    paired_donors, estimate, bootstrap_ci_low, bootstrap_ci_high,
    target_is_top_compartment
  )

effects <- bind_rows(original_effects, sensitivity_effects) %>%
  mutate(
    program = factor(unname(program_labels[signature_id]), levels = unname(program_labels[primary_ids])),
    method = factor(unname(method_labels[annotation_method]), levels = unname(method_labels)),
    dataset = factor(dataset, levels = rev(atlas_order)),
    top_status = ifelse(target_is_top_compartment, "Frozen target was top compartment", "Frozen target was not top compartment")
  )
write_tsv(effects, file.path(out_dir, "FigureS6_panel_a_effects_source.tsv"))

p_a <- ggplot(effects, aes(estimate, dataset, colour = method, shape = top_status)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_errorbar(
    aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
    orientation = "y", width = 0.10, linewidth = 0.45,
    position = position_dodge(width = 0.48)
  ) +
  geom_point(size = 2.0, stroke = 0.7, position = position_dodge(width = 0.48)) +
  facet_wrap(~program, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = method_colors) +
  scale_shape_manual(values = c(
    `Frozen target was top compartment` = 16,
    `Frozen target was not top compartment` = 1
  )) +
  labs(
    x = "Mean paired donor difference\n(frozen target - other qualified compartments)",
    y = NULL,
    subtitle = "Target-compartment preference under alternative annotation rules",
    tag = "a"
  ) +
  theme_paper() +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    panel.spacing.x = unit(4, "mm"),
    plot.subtitle = element_text(face = "bold")
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1)
  )

sensitivity_gate <- read_tsv(
  file.path(sensitivity_dir, "annotation_sensitivity_AND_gate.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    signature_id, annotation_method,
    lodo_same = GSE184362_lodo_same_direction,
    lodo_total = GSE184362_lodo_total,
    three_atlas = three_atlas_same_frozen_direction,
    gate_pass = sensitivity_AND_gate
  )
primary_cross <- read_tsv(file.path(audit, "results", "07_three_atlas_AND_gate.tsv"), show_col_types = FALSE) %>%
  filter(signature_id %in% primary_ids) %>%
  select(signature_id, three_atlas_same_direction, three_atlas_AND_gate)
primary_lodo <- read_tsv(file.path(audit, "results", "04_lodo_summary.tsv"), show_col_types = FALSE) %>%
  filter(signature_id %in% primary_ids) %>%
  left_join(primary_cross, by = "signature_id") %>%
  transmute(
    signature_id, annotation_method = "primary_hybrid",
    lodo_same = same_direction_estimable_iterations,
    lodo_total = total_iterations,
    three_atlas = three_atlas_same_direction,
    gate_pass = three_atlas_AND_gate
  )
gate <- bind_rows(primary_lodo, sensitivity_gate) %>%
  mutate(
    program = factor(unname(program_labels[signature_id]), levels = rev(unname(program_labels[primary_ids]))),
    method = factor(unname(method_labels[annotation_method]), levels = unname(method_labels)),
    label = sprintf("%d/%d\n3-atlas %s", lodo_same, lodo_total, ifelse(three_atlas, "+", "-"))
  )
write_tsv(gate, file.path(out_dir, "FigureS6_panel_b_gate_source.tsv"))

p_b <- ggplot(gate, aes(method, program, fill = gate_pass)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = label), family = "Nimbus Sans", size = 2.45, lineheight = 0.9) +
  scale_fill_manual(values = c(`TRUE` = "#D9EAD3", `FALSE` = "#F4CCCC"), guide = "none") +
  scale_x_discrete(labels = c(
    `Primary hybrid` = "Primary\nhybrid",
    `Leave-program-genes-out` = "LPGO\nhybrid",
    `SingleR only` = "SingleR\nonly"
  )) +
  labs(
    x = NULL, y = NULL,
    subtitle = "Frozen LODO and three-atlas AND gate",
    tag = "b"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1, size = 6.3),
    axis.ticks = element_blank(), axis.line = element_blank(),
    plot.subtitle = element_text(face = "bold")
  )

gse191 <- read_tsv(
  file.path(sensitivity_dir, "GSE191288_donor_differences_primary_and_sensitivity.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    program = factor(unname(program_labels[signature_id]), levels = unname(program_labels[primary_ids])),
    method = factor(unname(method_labels[annotation_method]), levels = unname(method_labels)),
    donor = factor(donor, levels = c("P1", "P2", "P3"))
  )
write_tsv(gse191, file.path(out_dir, "FigureS6_panel_c_GSE191288_donor_source.tsv"))

p_c <- ggplot(gse191, aes(method, difference, colour = donor, group = donor)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#777777") +
  geom_line(linewidth = 0.4, alpha = 0.65) +
  geom_point(size = 1.9) +
  facet_wrap(~program, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c(P1 = "#4C78A8", P2 = "#E3A33B", P3 = "#59A14F")) +
  scale_x_discrete(labels = c(
    `Primary hybrid` = "Primary\nhybrid",
    `Leave-program-genes-out` = "LPGO\nhybrid",
    `SingleR only` = "SingleR\nonly"
  )) +
  labs(
    x = NULL, y = "Paired donor difference",
    subtitle = "GSE191288: all three donors retained the frozen direction",
    caption = "n = 3 donors. LODO is not estimable:\neach deletion leaves two donors (<3).",
    tag = "c"
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(angle = 22, hjust = 1, size = 6.1),
    legend.position = "bottom", plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 6.2),
    plot.margin = margin(3, 3, 8, 3, unit = "mm")
  ) +
  guides(colour = guide_legend(nrow = 1))

figure <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.05, 1.0), widths = c(0.78, 1.45)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

base <- file.path(out_dir, "FigureS6_annotation_circularity_sensitivity")
width_in <- 183 / 25.4
height_in <- 158 / 25.4
svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(figure)
dev.off()
cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "Nimbus Sans")
print(figure)
dev.off()
agg_tiff(
  paste0(base, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure)
dev.off()
agg_png(paste0(base, ".png"), width = width_in, height = height_in, units = "in", res = 300)
print(figure)
dev.off()

contract <- c(
  "Core conclusion: frozen target-compartment preferences remain directionally positive under leave-program-genes-out and SingleR-only annotation, with Condello-6 retaining T/NK top-compartment identity in all three atlases.",
  "Figure archetype: quantitative grid.",
  "Backend: R only.",
  "Target/output: 183 x 158 mm; editable SVG/PDF, 600-dpi TIFF and 300-dpi PNG.",
  "Panel a: effect sizes and donor-bootstrap intervals under three annotation rules.",
  "Panel b: frozen GSE184362 LODO plus three-atlas AND gate.",
  "Panel c: descriptive GSE191288 donor directions; LODO is explicitly not estimable at n=3.",
  "Reviewer risk: positive target-vs-other mean is not equivalent to top-compartment identity; the open point flags the SingleR-only TDS-16 exception in GSE184362."
)
writeLines(contract, file.path(out_dir, "FigureS6_figure_contract.txt"))

legend <- paste0(
  "**Supplementary Fig. S6 | Annotation-circularity sensitivity of cross-atlas compartment preferences.** ",
  "**a,** Mean paired donor differences between each frozen target compartment and the same donor's other QC-qualified compartments under the primary hybrid annotation, a leave-program-genes-out hybrid annotation, and SingleR-only annotation. Program genes were excluded from every marker panel and from the shared SingleR test/reference feature space in the leave-program-genes-out analysis. Error bars are percentile 95% confidence intervals from 10,000 donor-level bootstrap resamples; filled points indicate that the frozen target was also the highest-scoring compartment, whereas open points indicate a positive target-versus-other mean without top-compartment identity. ",
  "**b,** Pre-specified GSE184362 leave-one-donor-out and three-atlas AND gate. All six program-by-sensitivity combinations retained 11/11 leave-one-donor-out directions and the same positive direction across the three atlases. ",
  "**c,** Individual GSE191288 donor differences. All three donors retained the frozen direction for every primary program under all annotation rules. Because deletion of one donor leaves only two donors, leave-one-donor-out effects were not estimable under the pre-specified three-donor minimum. These sensitivity analyses address broad-compartment annotation circularity and do not establish RAIR specificity, patient-level provenance or causality."
)
writeLines(legend, file.path(out_dir, "FigureS6_legend.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Figure S6 written to: ", base)
