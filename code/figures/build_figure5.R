#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

set.seed(42)
options(stringsAsFactors = FALSE)

root <- normalizePath(".")
out_dir <- file.path(root, "figures", "final")
res_dir <- file.path(root, "results")
text_dir <- file.path(root, "manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)

primary <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
program_label <- c(
  TDS_16 = "TDS-16",
  IODIDE_HANDLING_11 = "Iodide handling-11",
  CONDELLO_2025_SIX = "Condello-6",
  HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_E2F_TARGETS = "E2F targets",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_G2M_CHECKPOINT = "G2M checkpoint",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
)
target_map <- c(
  TDS_16 = "Thyroid cells",
  IODIDE_HANDLING_11 = "Thyroid cells",
  CONDELLO_2025_SIX = "T & NK cells",
  HALLMARK_ANGIOGENESIS = "Fibroblasts",
  HALLMARK_E2F_TARGETS = "B cells",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "Fibroblasts",
  HALLMARK_G2M_CHECKPOINT = "B cells",
  HALLMARK_HYPOXIA = "Fibroblasts",
  HALLMARK_INFLAMMATORY_RESPONSE = "Myeloid cells"
)
compartment_order <- c(
  "Thyroid cells", "T & NK cells", "B cells", "Myeloid cells",
  "Fibroblasts", "Endothelial cells"
)
compartment_short <- c(
  "Thyroid cells" = "Thyroid",
  "T & NK cells" = "T/NK",
  "B cells" = "B",
  "Myeloid cells" = "Myeloid",
  "Fibroblasts" = "Fibroblast",
  "Endothelial cells" = "Endothelial"
)
compartment_colors <- c(
  "Thyroid cells" = "#0072B2",
  "T & NK cells" = "#E69F00",
  "B cells" = "#56B4E9",
  "Myeloid cells" = "#D55E00",
  "Fibroblasts" = "#009E73",
  "Endothelial cells" = "#CC79A7"
)
atlas_order <- c("GSE184362", "GSE191288", "GSE281736")
atlas_n <- c(GSE184362 = 11L, GSE191288 = 3L, GSE281736 = 6L)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Nimbus Sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "#F2F2F2", color = NA),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_blank(),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = c(0, 1)
    )
}

# Reconstruct frozen target-versus-other donor differences for GSE184362.
pb <- read_tsv(file.path(res_dir, "04_pseudobulk_program_scores.tsv"), show_col_types = FALSE) %>%
  filter(compartment %in% compartment_order)

gse184_diff <- pb %>%
  mutate(target_compartment = unname(target_map[signature_id])) %>%
  group_by(donor, signature_id, target_compartment) %>%
  summarise(
    target_score = score[compartment == target_compartment][1],
    other_mean = mean(score[compartment != target_compartment], na.rm = TRUE),
    difference = target_score - other_mean,
    .groups = "drop"
  ) %>%
  filter(is.finite(difference)) %>%
  mutate(dataset = "GSE184362") %>%
  select(donor, difference, dataset, signature_id, target_compartment)

external_diff <- read_tsv(
  file.path(res_dir, "07_external_atlas_donor_differences.tsv"),
  show_col_types = FALSE
)
all_diff <- bind_rows(gse184_diff, external_diff) %>%
  mutate(
    dataset = factor(dataset, levels = atlas_order),
    program = unname(program_label[signature_id]),
    target_color = unname(compartment_colors[target_compartment])
  )

boot_mean <- function(x, B = 10000L, seed = 42L) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) {
    return(tibble(estimate = mean(x), conf_low = NA_real_, conf_high = NA_real_, n_donors = length(x)))
  }
  set.seed(seed)
  vals <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  tibble(
    estimate = mean(x),
    conf_low = unname(quantile(vals, 0.025, type = 7)),
    conf_high = unname(quantile(vals, 0.975, type = 7)),
    n_donors = length(x)
  )
}

boot_summary <- all_diff %>%
  group_by(dataset, signature_id, program, target_compartment) %>%
  group_modify(~boot_mean(.x$difference, B = 10000L, seed = 42L + as.integer(.y$dataset))) %>%
  ungroup() %>%
  mutate(
    dataset = as.character(dataset),
    direction = case_when(
      is.na(conf_low) ~ "not_estimable",
      estimate > 0 ~ "target_preference",
      TRUE ~ "direction_reversal"
    )
  )

write_tsv(boot_summary, file.path(res_dir, "08_three_atlas_bootstrap_effects.tsv"))

# Panel A: structural UMAP, deliberately independent of program scoring.
obj <- readRDS(file.path(root, "checkpoints", "03_mnn_clustered_counts.rds"))
umap <- as.data.frame(Embeddings(obj, "umap"))
meta <- obj[[]]
umap$cell <- rownames(umap)
umap$compartment <- meta[umap$cell, "compartment"]
set.seed(42)
umap <- umap %>%
  filter(compartment %in% compartment_order) %>%
  group_by(compartment) %>%
  group_modify(~slice_sample(.x, n = min(3000L, nrow(.x)))) %>%
  ungroup() %>%
  mutate(compartment = factor(compartment, levels = compartment_order))

pA <- ggplot(umap, aes(x = umap_1, y = umap_2, color = compartment)) +
  geom_point(size = 0.18, alpha = 0.55, stroke = 0) +
  scale_color_manual(values = compartment_colors, labels = compartment_short) +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 2)) +
  labs(x = "UMAP 1", y = "UMAP 2", subtitle = "GSE184362 compartment structure", tag = "A") +
  theme_paper(8.5) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.5),
    legend.key.width = unit(0.38, "cm"),
    plot.subtitle = element_text(face = "bold", size = 8.5)
  )

# Panel B: donor-level discovery-atlas pseudobulk evidence.
pb_primary <- pb %>%
  filter(signature_id %in% primary) %>%
  mutate(
    program = factor(unname(program_label[signature_id]), levels = unname(program_label[primary])),
    compartment = factor(compartment, levels = compartment_order)
  )

pB <- ggplot(pb_primary, aes(x = compartment, y = score, color = compartment)) +
  geom_boxplot(width = 0.60, outlier.shape = NA, linewidth = 0.35, color = "#595959", fill = "white") +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 42), size = 0.85, alpha = 0.72) +
  facet_wrap(~program, ncol = 1, scales = "free_y") +
  scale_color_manual(values = compartment_colors, guide = "none") +
  scale_x_discrete(labels = compartment_short) +
  labs(x = NULL, y = "Donor pseudobulk program score", subtitle = "GSE184362 donor-level scores", tag = "B") +
  theme_paper(8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7),
    strip.text = element_text(size = 7.8),
    plot.subtitle = element_text(face = "bold", size = 8.5),
    panel.spacing.y = unit(0.08, "cm")
  )

# Panel C: pre-specified target-compartment preference across atlases.
forest <- boot_summary %>%
  filter(signature_id %in% primary) %>%
  mutate(
    program = factor(program, levels = unname(program_label[primary])),
    dataset = factor(dataset, levels = rev(atlas_order)),
    label = sprintf("%.2f [%.2f, %.2f]", estimate, conf_low, conf_high)
  )

x_lim <- range(c(0, forest$conf_low, forest$conf_high), na.rm = TRUE)
x_pad <- diff(x_lim) * 0.08

pC <- ggplot(forest, aes(x = estimate, y = dataset, color = target_compartment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777", linewidth = 0.4) +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), orientation = "y", width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.3) +
  facet_wrap(~program, nrow = 1) +
  scale_color_manual(values = compartment_colors, labels = compartment_short) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    x = "Mean paired donor difference\n(frozen target compartment - other compartments)",
    y = NULL,
    subtitle = "Target-compartment preference across independent PTC atlases",
    tag = "C"
  ) +
  theme_paper(9) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8.5),
    plot.subtitle = element_text(face = "bold", size = 9),
    panel.spacing.x = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(nrow = 1))

condello_note <- ggplot() +
  annotate(
    "label", x = 0, y = 0,
    label = "Condello-6: T/NK was the highest-scoring compartment in all three atlases",
    hjust = 0, vjust = 0.5, size = 3.0, family = "Nimbus Sans",
    fill = alpha(compartment_colors[["T & NK cells"]], 0.14),
    color = "#4D4D4D", linewidth = 0.25
  ) +
  xlim(0, 1) + ylim(-0.5, 0.5) + theme_void()

main_fig <- ((pA + pB) + plot_layout(widths = c(0.95, 1.35))) /
  pC /
  condello_note +
  plot_layout(heights = c(1.35, 1.05, 0.14)) +
  plot_annotation() &
  theme(plot.background = element_rect(fill = "white", color = NA))

pdf_file <- file.path(out_dir, "Figure5_cellular_compartment_validation.pdf")
png_file <- file.path(out_dir, "Figure5_cellular_compartment_validation.png")
ggsave(pdf_file, main_fig, width = 7.0, height = 8.4, units = "in", device = cairo_pdf)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(png_file, main_fig, width = 7.0, height = 8.4, units = "in", dpi = 600, device = ragg::agg_png)
} else {
  ggsave(png_file, main_fig, width = 7.0, height = 8.4, units = "in", dpi = 600, device = "png", type = "cairo")
}

# Supplementary Figure S5: all frozen programs and explicit non-estimable cells.
gate <- read_tsv(file.path(res_dir, "07_three_atlas_AND_gate.tsv"), show_col_types = FALSE)
supp <- gate %>%
  select(signature_id, target_compartment,
         GSE184362_preference, GSE191288_preference, GSE281736_preference,
         GSE184362_donor_gate, GSE191288_estimable, GSE281736_estimable) %>%
  pivot_longer(
    cols = c(GSE184362_preference, GSE191288_preference, GSE281736_preference),
    names_to = "dataset", values_to = "estimate"
  ) %>%
  mutate(dataset = sub("_preference$", "", dataset)) %>%
  left_join(
    tibble(
      signature_id = rep(gate$signature_id, each = 3),
      dataset = rep(atlas_order, times = nrow(gate)),
      estimable = c(rbind(gate$GSE184362_donor_gate, gate$GSE191288_estimable, gate$GSE281736_estimable))
    ),
    by = c("signature_id", "dataset")
  ) %>%
  mutate(
    program = unname(program_label[signature_id]),
    program = factor(program, levels = rev(unname(program_label[gate$signature_id]))),
    dataset = factor(dataset, levels = atlas_order),
    cell_label = ifelse(estimable & is.finite(estimate), sprintf("%.2f", estimate), "NE"),
    plot_value = ifelse(estimable, estimate, NA_real_)
  )

pS <- ggplot(supp, aes(x = dataset, y = program, fill = plot_value)) +
  geom_tile(color = "white", linewidth = 1.1) +
  geom_text(aes(label = cell_label), size = 3.2, family = "Nimbus Sans", color = "black") +
  scale_fill_gradient2(
    low = "#F2F2F2", mid = "#9ECAE1", high = "#08519C",
    midpoint = 1.0, na.value = "#D9D9D9",
    name = "Mean paired\ndifference"
  ) +
  labs(
    x = NULL, y = NULL,
    subtitle = "Frozen target-compartment preference for all audited programs",
    caption = "NE, not estimable under the pre-specified donor rule (GSE191288 B-cell pseudobulk: 2 qualified donors)."
  ) +
  theme_minimal(base_size = 9, base_family = "Nimbus Sans") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(hjust = 1),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, size = 8),
    legend.position = "right"
  )

supp_pdf <- file.path(out_dir, "FigureS5_all_program_compartment_gate.pdf")
supp_png <- file.path(out_dir, "FigureS5_all_program_compartment_gate.png")
ggsave(supp_pdf, pS, width = 7.0, height = 5.2, units = "in", device = cairo_pdf)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(supp_png, pS, width = 7.0, height = 5.2, units = "in", dpi = 600, device = ragg::agg_png)
} else {
  ggsave(supp_png, pS, width = 7.0, height = 5.2, units = "in", dpi = 600, device = "png", type = "cairo")
}

# Manuscript-ready text uses the exact remotely calculated estimates and CIs.
fmt_effect <- function(sig, ds) {
  z <- boot_summary %>% filter(signature_id == sig, dataset == ds)
  sprintf("%.2f (95%% bootstrap CI %.2f to %.2f)", z$estimate, z$conf_low, z$conf_high)
}

caption <- paste0(
  "## Figure 5 legend\n\n",
  "**Figure 5 | Cross-atlas cellular-compartment preferences of radioiodine-related transcriptomic programs.** ",
  "**A,** UMAP of GSE184362 cells colored by the frozen broad-compartment annotation; a deterministic, compartment-stratified display sample (maximum 3,000 cells per compartment) is shown for legibility. ",
  "**B,** Donor-level pseudobulk scores for TDS-16, iodide handling-11, and Condello-6 in GSE184362. Points denote donors and boxes summarize the donor distributions. ",
  "**C,** Mean paired donor differences between each program's pre-specified target compartment and the mean of the same donor's other QC-qualified compartments across three independent papillary thyroid cancer atlases. Error bars are percentile 95% confidence intervals from 10,000 donor-level bootstrap resamples. ",
  "Positive values indicate preference for the frozen target compartment. TDS-16 and iodide handling-11 preferentially mapped to thyroid cells, whereas Condello-6 preferentially mapped to T/NK cells in all three atlases. ",
  "Pseudobulk units required at least 20 cells per donor-by-cell-type combination. In leave-one-donor-out analyses, cell types represented by fewer than three qualified donors in an iteration were classified as not estimable rather than as direction reversals. ",
  "These analyses quantify cross-dataset, population-level cell-type expression preferences and do not establish patient-level malignant-cell origin, RAIR specificity, or causality.\n\n"
)

methods <- paste0(
  "## Methods: cross-atlas cellular-compartment preference analysis\n\n",
  "We evaluated cellular-compartment preferences in three independent PTC single-cell atlases (GSE184362, GSE191288, and GSE281736) using analysis decisions frozen before cross-atlas evaluation. Gene expression was aggregated into donor-by-cell-type pseudobulk units; a unit was QC-qualified when it contained at least 20 cells. Program scores were calculated at pseudobulk level. For each donor and program, the compartment-preference effect was defined as the score in the pre-specified target compartment minus the mean score across the same donor's other QC-qualified broad compartments. Target compartments were thyroid cells for TDS-16 and iodide handling-11 and T/NK cells for Condello-6. Uncertainty was estimated by 10,000 percentile bootstrap resamples of donors (seed 42). Cross-donor consistency in GSE184362 was assessed by leave-one-donor-out analysis. A cell type represented by fewer than three QC-qualified donors after donor removal was excluded from that iteration and recorded as not estimable; such iterations were not counted as direction reversals. Main-text eligibility required both at least 9 of 11 leave-one-donor-out iterations in the frozen direction and the same direction in all three atlases. No thresholds, annotations, target compartments, or eligibility rules were modified after inspecting cross-atlas results.\n\n"
)

results <- paste0(
  "## Results: distinct lineage and immune-compartment preferences replicate across PTC atlases\n\n",
  "All three primary programs satisfied the pre-specified cross-donor and three-atlas AND gate. TDS-16 preferentially mapped to thyroid cells in GSE184362 [", fmt_effect("TDS_16", "GSE184362"), "], GSE191288 [", fmt_effect("TDS_16", "GSE191288"), "], and GSE281736 [", fmt_effect("TDS_16", "GSE281736"), "]. Iodide handling-11 showed the same thyroid-cell preference in GSE184362 [", fmt_effect("IODIDE_HANDLING_11", "GSE184362"), "], GSE191288 [", fmt_effect("IODIDE_HANDLING_11", "GSE191288"), "], and GSE281736 [", fmt_effect("IODIDE_HANDLING_11", "GSE281736"), "]. In contrast, Condello-6 preferentially mapped to T/NK cells in all three atlases: GSE184362 [", fmt_effect("CONDELLO_2025_SIX", "GSE184362"), "], GSE191288 [", fmt_effect("CONDELLO_2025_SIX", "GSE191288"), "], and GSE281736 [", fmt_effect("CONDELLO_2025_SIX", "GSE281736"), "]. T/NK cells were the highest-scoring broad compartment for Condello-6 in each atlas. This reproducible immune-compartment preference challenges a thyroid-lineage-intrinsic interpretation of Condello-6 and indicates sensitivity to immune-cell composition; because the original bulk cohort and the single-cell atlases were unpaired, it does not establish that immune infiltration caused the original bulk association.\n\n"
)

supp_caption <- paste0(
  "## Supplementary Figure S5 legend\n\n",
  "**Supplementary Figure S5 | Cross-atlas compartment-preference gate for all audited programs.** Cells report the mean paired donor difference between each program's frozen target compartment and the same donor's other QC-qualified compartments. NE denotes not estimable under the pre-specified donor rule. E2F targets and G2M checkpoint were not estimable in GSE191288 because only two donors contributed QC-qualified B-cell pseudobulk units; these entries were retained as honest nulls and were not interpreted as direction reversals.\n"
)

writeLines(c(caption, methods, results, supp_caption), file.path(text_dir, "Figure5_methods_results_legend.md"))

manifest <- c(
  "# Figure manifest",
  "",
  "- Figure: Figure 5 and Supplementary Figure S5",
  "- Key message: Condello-6 reproducibly prefers the T/NK compartment across three independent PTC atlases, whereas TDS-16 and iodide handling-11 prefer thyroid cells.",
  "- Analysis unit: donor-level QC-qualified pseudobulk; never cells.",
  "- Effect: frozen target-compartment score minus the same donor's mean score in other qualified compartments.",
  "- Uncertainty: 10,000 donor bootstrap resamples, percentile 95% CI, seed 42.",
  "- Frozen eligibility: GSE184362 LODO >=9/11 in the same direction AND same direction in GSE184362, GSE191288, and GSE281736.",
  "- Pseudobulk threshold: >=20 cells per donor-by-cell-type unit.",
  "- Non-estimable rule: <3 qualified donors after donor removal; not counted as reversal.",
  "- Main figure size: 7.0 x 8.4 inches; PDF vector and PNG 600 dpi.",
  "- Supplement size: 7.0 x 5.2 inches; PDF vector and PNG 600 dpi.",
  "- Typeface: Nimbus Sans; colorblind-safe Wong-derived compartment palette.",
  "- UMAP role: structural context only; program evidence is donor-pseudobulk based.",
  "- Boundary: population-level cell-type preference, not patient-level provenance, RAIR specificity, or causality.",
  "- Meta-analysis-only forest requirements (pooled diamond, I2, tau2): not applicable; this is a donor-level comparative effect display without pooled meta-analysis."
)
writeLines(manifest, file.path(out_dir, "_figure_manifest.md"))

cat("Generated:\n")
cat(pdf_file, "\n", png_file, "\n", supp_pdf, "\n", supp_png, "\n", sep = "")
cat(file.path(res_dir, "08_three_atlas_bootstrap_effects.tsv"), "\n")
cat(file.path(text_dir, "Figure5_methods_results_legend.md"), "\n")
