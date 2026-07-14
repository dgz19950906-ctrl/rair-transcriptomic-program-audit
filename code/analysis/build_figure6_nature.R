#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(svglite)
  library(ragg)
  library(grid)
})

project_dir <- "/home/dony/ThyroidCancer_Project/rair_audit"
results_dir <- file.path(project_dir, "results")
fig_dir <- file.path(project_dir, "figures", "final")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fig3 <- read_tsv(file.path(results_dir, "11_figure3_display_source.tsv"), show_col_types = FALSE)
fig2_effect <- read_tsv(file.path(results_dir, "10_figure2_bootstrap_hedges_g.tsv"), show_col_types = FALSE)
fig4_effect <- read_tsv(file.path(results_dir, "11_figure4_challenge_effects_source.tsv"), show_col_types = FALSE)
fig4_inc <- read_tsv(file.path(results_dir, "11_figure4_specificity_increment_source.tsv"), show_col_types = FALSE)
atlas <- read_tsv(file.path(results_dir, "07_three_atlas_AND_gate.tsv"), show_col_types = FALSE)

ids <- c("TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX")
program_order <- c("TDS-16", "Iodide handling-11", "Condello-6")

fmt <- function(x) sprintf("%+.2f", x)
fmt_p <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))

clinical <- fig3 %>%
  filter(signature_id %in% ids) %>%
  mutate(axis = if_else(contrast == "uptake_failure", "uptake", "response")) %>%
  select(signature_id, axis, full_effect, same_direction_n, n_estimable, matched_random_p) %>%
  pivot_wider(names_from = axis,
              values_from = c(full_effect, same_direction_n, n_estimable, matched_random_p))

challenge <- fig2_effect %>%
  filter(signature_id %in% ids, contrast == "uptake_failure") %>%
  transmute(signature_id, gse151179_uptake = adverse_aligned_hedges_g) %>%
  left_join(
    fig4_effect %>%
      filter(signature_id %in% ids) %>%
      transmute(signature_id, gse299988_challenge = adverse_aligned_hedges_g),
    by = "signature_id"
  ) %>%
  left_join(fig4_inc %>% select(signature_id, specificity_margin_vs_strongest_generic,
                                specificity_increment_bootstrap_low,
                                specificity_increment_bootstrap_high,
                                incremental_evidence_status), by = "signature_id")

passport <- tibble(
  signature_id = ids,
  program = program_order,
  construct = c("Thyroid\ndifferentiation", "Iodide-handling\nphysiology", "Ex vivo iodine\ncontent"),
  recovery = c("Exact 16-gene\nbiological proxy", "Exact 11-gene\nprespecified reference", "Partial 6-gene\nsubset"),
  classification = c("Lineage proxy;\nnot RAIR-\nspecific", "Physiology proxy;\nnot RAIR-\nspecific", "Immune-composition-\nsensitive\nex vivo candidate")
) %>%
  left_join(clinical, by = "signature_id") %>%
  left_join(challenge, by = "signature_id") %>%
  left_join(atlas %>% select(signature_id, target_compartment, three_atlas_AND_gate), by = "signature_id") %>%
  mutate(
    endpoint = paste0("U ", fmt(full_effect_uptake), "\nR ", fmt(full_effect_response)),
    lopo = paste0(same_direction_n_uptake, "/", n_estimable_uptake,
                  " U\n", same_direction_n_response, "/", n_estimable_response, " R"),
    random = paste0("Not extreme\nP ", fmt_p(matched_random_p_uptake),
                    " / ", fmt_p(matched_random_p_response)),
    challenge_shift = paste0("U ", fmt(gse151179_uptake), " to ",
                             fmt(gse299988_challenge),
                             "\nNo robust\nincrement"),
    compartment = paste0(str_replace_all(target_compartment, " & ", "/"),
                         "\n3/3 atlases;\nAND pass"),
    recovery = if_else(signature_id == "IODIDE_HANDLING_11",
                       "Exact 11-gene\nprespecified\nreference", recovery)
  )

long <- passport %>%
  select(program, construct, recovery, endpoint, lopo, random, challenge_shift, compartment, classification) %>%
  pivot_longer(-program, names_to = "dimension", values_to = "display") %>%
  mutate(
    dimension = factor(dimension,
      levels = c("construct", "recovery", "endpoint", "lopo", "random",
                 "challenge_shift", "compartment", "classification"),
      labels = c("Intended\nconstruct", "Recovery", "Clinical alignment\nHedges' g",
                 "LOPO direction\nU / R", "Matched-random\ncalibration",
                 "RAI/LN challenge", "Cell compartment", "Allowed\nclassification")),
    program = factor(program, levels = rev(program_order)),
    domain = case_when(
      dimension %in% c("Intended\nconstruct", "Recovery") ~ "Provenance",
      dimension %in% c("Clinical alignment\nHedges' g", "LOPO direction\nU / R", "Matched-random\ncalibration") ~ "Endpoint/robustness",
      dimension == "RAI/LN challenge" ~ "Confounding",
      dimension == "Cell compartment" ~ "Compartment",
      TRUE ~ "Interpretation"
    )
  )

domain_cols <- c(
  Provenance = "#E9EDF2",
  `Endpoint/robustness` = "#E6EEF7",
  Confounding = "#F8E6DF",
  Compartment = "#E4F1EA",
  Interpretation = "#EEEAF5"
)

theme_nature <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(colour = "#252525"),
      plot.title = element_text(size = 7.1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.1, colour = "#555555", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_nature())

# Panel a: independent audit domains.
domains <- tibble(
  x = 1:5,
  title = c("Recoverability", "Endpoint alignment", "Robustness", "Confounding", "Compartment"),
  sub = c("genes + direction\n+ formula", "uptake vs response", "LOPO + matched\nrandom sets", "non-identifiable\nchallenge", "donor-aware\nthree-atlas gate"),
  fill = c("#E9EDF2", "#E6EEF7", "#E6EEF7", "#F8E6DF", "#E4F1EA")
)

p_a <- ggplot(domains, aes(x, 1)) +
  geom_segment(data = domains %>% filter(x < 5),
               aes(x = x + 0.44, xend = x + 0.56, y = 1, yend = 1),
               colour = "#A6A6A6", linewidth = 0.55,
               arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_tile(aes(fill = fill), width = 0.86, height = 0.64, colour = "white", linewidth = 0.8) +
  geom_text(aes(label = title), y = 1.10, family = "Arial", fontface = "bold", size = 2.6) +
  geom_text(aes(label = sub), y = 0.89, family = "Arial", colour = "#555555", size = 2.15, lineheight = 0.92) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.5, 5.5), ylim = c(0.62, 1.38), clip = "off") +
  labs(title = "Independent audit dimensions",
       subtitle = "Passing one layer does not imply validity in another") +
  theme_void(base_family = "Arial", base_size = 6.4) +
  theme(plot.title = element_text(size = 7.1, face = "bold"),
        plot.subtitle = element_text(size = 6.1, colour = "#555555"),
        plot.margin = margin(3, 5, 2, 5))

# Panel b: hero passport matrix.
p_b <- ggplot(long, aes(dimension, program)) +
  geom_tile(aes(fill = domain), colour = "white", linewidth = 1.0) +
  geom_text(aes(label = display), family = "Arial", size = 1.88, lineheight = 0.88, colour = "#222222") +
  scale_fill_manual(values = domain_cols, guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Frozen program passports",
       subtitle = "U, uptake failure; R, response failure despite retained uptake; P, expression-matched random-program P value") +
  theme(
    axis.text.x = element_text(size = 5.65, face = "bold", lineheight = 0.9, margin = margin(b = 4)),
    axis.text.y = element_text(size = 6.2, face = "bold", margin = margin(r = 5)),
    plot.margin = margin(3, 5, 3, 5)
  )

# Panel c: evidence-calibrated closing boundary.
boundary <- tibble(
  side = rep(c("Supported by this audit", "Not established"), each = 3),
  y = rep(3:1, 2),
  text = c(
    "Frozen audit programs were projected without outcome-guided refitting",
    "Endpoint alignment and stability can be audited separately",
    "Population-level compartment preference replicates across atlases",
    "A portable RAIR-specific clinical signature",
    "Patient-level malignant-cell origin or immune causality",
    "Independent RAI effects in the RAI/LN-collinear challenge"
  )
) %>%
  mutate(
    side = factor(side, levels = c("Supported by this audit", "Not established")),
    text = str_wrap(text, width = 47)
  )

p_c <- ggplot(boundary, aes(y = y)) +
  geom_rect(data = tibble(side = factor(c("Supported by this audit", "Not established"),
                                        levels = c("Supported by this audit", "Not established")),
                          xmin = 0.02, xmax = 0.98, ymin = 0.45, ymax = 3.55,
                          fill = c("#E4F1EA", "#F3ECE9")),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, colour = NA) +
  geom_point(aes(x = 0.09, colour = side), size = 2.0) +
  geom_text(aes(x = 0.15, label = text), hjust = 0, family = "Arial", size = 2.15,
            lineheight = 0.92, colour = "#252525") +
  facet_wrap(~ side, nrow = 1) +
  scale_fill_identity() +
  scale_colour_manual(values = c("Supported by this audit" = "#3A7D5B", "Not established" = "#A85B4A"), guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.35, 3.65), clip = "off") +
  labs(title = "Shared claim boundary") +
  theme_void(base_family = "Arial", base_size = 6.4) +
  theme(
    strip.text = element_text(size = 6.4, face = "bold"),
    plot.title = element_text(size = 7.1, face = "bold"),
    panel.spacing.x = unit(4, "mm"),
    plot.margin = margin(3, 5, 4, 5)
  )

figure <- p_a / p_b / p_c +
  plot_layout(heights = c(0.75, 2.1, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

base <- file.path(fig_dir, "Figure6_program_passport_synthesis")
width_in <- 183 / 25.4
height_in <- 145 / 25.4

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

sanitize_tsv <- function(df) {
  df %>% mutate(across(
    everything(),
    ~ if (is.character(.x) || is.factor(.x)) {
      str_replace_all(as.character(.x), "[\\r\\n]+", " | ")
    } else {
      .x
    }
  ))
}
write_tsv(sanitize_tsv(passport), file.path(results_dir, "12_figure6_program_passport_source.tsv"))
write_tsv(sanitize_tsv(long), file.path(results_dir, "12_figure6_passport_matrix_source.tsv"))
write_tsv(sanitize_tsv(boundary), file.path(results_dir, "12_figure6_claim_boundary_source.tsv"))

message("Figure 6 written to: ", base)
