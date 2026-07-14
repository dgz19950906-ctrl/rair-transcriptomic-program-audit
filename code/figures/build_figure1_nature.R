#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(42)

root <- normalizePath(".")
audit_dir <- file.path(root, "literature_audit")
out_dir <- file.path(root, "figures", "final")
res_dir <- file.path(root, "results")
manuscript_dir <- file.path(root, "manuscript")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Nimbus Sans"
width_mm <- 183
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

palette <- c(
  ink = "#222222",
  neutral = "#6B7280",
  neutral_light = "#E5E7EB",
  blue = "#377EB8",
  blue_mid = "#8FBAD2",
  blue_light = "#DCEAF4",
  teal = "#5AAE9F",
  teal_light = "#DCEFEA",
  orange = "#E39A3B",
  orange_light = "#F7E7CE",
  red = "#C95A54",
  red_light = "#F4DAD8",
  grey_fill = "#F5F5F3"
)

theme_nature_schematic <- function(base_size = 7.2) {
  theme_void(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0, color = palette[["ink"]]),
      plot.subtitle = element_text(size = base_size - 0.2, color = palette[["neutral"]], hjust = 0),
      plot.tag = element_text(size = 9, face = "bold", color = "black"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

flow <- fromJSON(file.path(audit_dir, "literature_flow_counts.json"))
audit <- read_tsv(file.path(audit_dir, "literature_tool_audit.tsv"), show_col_types = FALSE)
stopifnot(nrow(audit) == 22L)

# Strict, sequential full-satisfaction rules. Partial evidence remains in the
# study passport but is not counted as a pass in the visual funnel.
audit_flags <- audit %>%
  mutate(
    level1_full = !is.na(claim_present) & claim_present != "",
    level2_full = !grepl("^(no_|partial)", clinical_endpoint_recoverable) &
      !is.na(clinical_endpoint_recoverable),
    level3_full = grepl(
      "^exact_|thresholds?_recoverable|miR-200c-3p_direction_and_threshold_recoverable|^TDS16_recoverable$|^eTDS_recoverable",
      genes_directions_formula_recoverable
    ),
    level4_full = grepl("^yes($|_)", coalesce(paper_code_endpoint_consistent, "")),
    level5_full = patient_tissue_timepoint_auditable == "yes",
    level6_full = grepl("^yes($|_)", directly_portable_without_refitting),
    pass1 = level1_full,
    pass2 = pass1 & level2_full,
    pass3 = pass2 & level3_full,
    pass4 = pass3 & level4_full,
    pass5 = pass4 & level5_full,
    pass6 = pass5 & level6_full
  )

funnel_counts <- c(
  sum(audit_flags$pass1), sum(audit_flags$pass2), sum(audit_flags$pass3),
  sum(audit_flags$pass4), sum(audit_flags$pass5), sum(audit_flags$pass6)
)
expected_counts <- c(22L, 15L, 5L, 1L, 0L, 0L)
if (!identical(as.integer(funnel_counts), expected_counts)) {
  stop("Frozen funnel counts changed: observed ", paste(funnel_counts, collapse = ", "),
       "; expected ", paste(expected_counts, collapse = ", "))
}

flow_source <- tibble(
  stage = c(
    "Database retrieval", "Database deduplication", "Citation expansion",
    "Title/abstract screening", "Excluded before full-text/code audit",
    "Full-text/code audit", "GEO datasets reported separately"
  ),
  count = c(
    flow$database_records_raw_total,
    flow$database_records_unique_before_citation_search,
    flow$citation_search_candidates,
    flow$unique_records_after_automated_deduplication,
    flow$unique_records_after_automated_deduplication - flow$total_full_text_or_code_audited,
    flow$total_full_text_or_code_audited,
    flow$geo_dataset_records_reported_separately
  ),
  note = c(
    "PubMed 702; Europe PMC 157; Semantic Scholar top 100 of 242",
    "PMID, DOI and title deduplication",
    "Eight seed studies; rate-limited routes declared",
    "Frozen master index",
    "Title/abstract exclusions",
    "21 database records plus one supplementary-source study",
    "Datasets, not article records"
  )
)

funnel_source <- tibble(
  level = 1:6,
  criterion = c(
    "Radioiodine-related expression claim present",
    "Clinical or biological endpoint recoverable",
    "Genes, directions and formula recoverable",
    "Paper and public code endpoint-consistent",
    "Patient, tissue site and timepoint fully auditable",
    "Directly portable without refitting"
  ),
  strict_sequential_full_count = funnel_counts,
  denominator = nrow(audit),
  counting_rule = "Full satisfaction at this and every preceding level; partial evidence not counted as a pass"
)

construct_source <- tibble(
  construct_1 = c("Ex vivo iodine content", "Clinical lesion uptake", "Ex vivo iodine content"),
  construct_2 = c("Clinical lesion uptake", "Structural treatment response", "Structural treatment response"),
  evidence_type = c("observed in current audit", "external clinical context", "untested"),
  annotation = c(
    "Condello-6 did not positively align with clinical uptake failure (adverse-oriented g = -0.11)",
    "Uptake and structural response are distinct clinical measurements (Boucai context)",
    "No paired public transcriptomic dataset tested this bridge"
  ),
  causal_interpretation = "none"
)

write_tsv(flow_source, file.path(res_dir, "09_figure1_search_flow_source.tsv"))
write_tsv(funnel_source, file.path(res_dir, "09_figure1_funnel_source.tsv"))
write_tsv(construct_source, file.path(res_dir, "09_figure1_construct_edges_source.tsv"))
write_tsv(audit_flags, file.path(res_dir, "09_figure1_study_level_flags.tsv"))

# Panel a: PRISMA-like search flow.
flow_boxes <- tibble(
  xmin = c(2, 27, 57, 80), xmax = c(21, 46, 76, 98),
  ymin = c(16, 16, 16, 16), ymax = c(31, 31, 31, 31),
  label = c(
    "Article databases\n959 raw records\nPubMed 702\nEurope PMC 157\nSemantic Scholar 100/242",
    "Deduplicated\n797 unique records",
    "Title/abstract screened\n892 unique records",
    "Full-text / code audit\n22 studies\n21 database +\n1 supplementary"
  ),
  fill = c(palette[["blue_light"]], palette[["blue_light"]], palette[["teal_light"]], palette[["orange_light"]]),
  text_size = c(2.25, 2.55, 2.55, 2.40)
)

p_a <- ggplot() +
  geom_rect(
    data = flow_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = palette[["ink"]], linewidth = 0.35
  ) +
  scale_fill_identity() +
  geom_text(
    data = flow_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label, size = text_size),
    family = font_family, lineheight = 0.94, color = palette[["ink"]]
  ) +
  scale_size_identity() +
  geom_segment(
    data = tibble(x = c(21, 46, 76), xend = c(27, 57, 80), y = 23.5, yend = 23.5),
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45, color = palette[["neutral"]],
    arrow = arrow(length = unit(2.2, "mm"), type = "closed")
  ) +
  annotate(
    "rect", xmin = 47.5, xmax = 56, ymin = 33, ymax = 40,
    fill = palette[["grey_fill"]], color = palette[["neutral"]], linewidth = 0.3
  ) +
  annotate(
    "text", x = 51.75, y = 36.5, label = "Citation expansion\n+95 candidates",
    family = font_family, size = 2.35, lineheight = 0.95
  ) +
  geom_segment(
    aes(x = 51.75, xend = 59, y = 33, yend = 31),
    linewidth = 0.4, color = palette[["neutral"]],
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  annotate(
    "rect", xmin = 57, xmax = 76, ymin = 3, ymax = 12,
    fill = "white", color = palette[["red"]], linewidth = 0.4
  ) +
  annotate(
    "text", x = 66.5, y = 7.5, label = "870 excluded before\nfull-text / code audit",
    family = font_family, size = 2.45, lineheight = 0.95, color = palette[["red"]]
  ) +
  geom_segment(
    aes(x = 66.5, xend = 66.5, y = 16, yend = 12),
    linewidth = 0.4, color = palette[["red"]],
    arrow = arrow(length = unit(2, "mm"), type = "closed")
  ) +
  annotate(
    "label", x = 12, y = 7.5,
    label = "27 GEO datasets\nreported separately",
    family = font_family, size = 2.35, lineheight = 0.95,
    fill = "white", color = palette[["neutral"]], linewidth = 0.25
  ) +
  annotate(
    "text", x = 12, y = 1.7, label = "Datasets, not article records",
    family = font_family, size = 2.05, color = palette[["neutral"]]
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 42), clip = "off") +
  labs(
    title = "Frozen systematic search and code-level audit",
    subtitle = "Search frozen on 12 July 2026; Semantic Scholar was a declared non-exhaustive sensitivity route",
    tag = "a"
  ) +
  theme_nature_schematic(7.2)

# Panel b: six-level strict sequential reproducibility funnel.
funnel_plot <- funnel_source %>%
  mutate(
    y = rev(level),
    half_width = 0.9 + 4.4 * sqrt(strict_sequential_full_count / max(denominator)),
    next_half_width = lead(half_width, default = half_width[n()]),
    fill = c("#DCEAF4", "#C4DCEB", "#9FC5D9", "#76ADC8", "#F2D4AE", "#E5A455")
  )

funnel_poly <- bind_rows(lapply(seq_len(nrow(funnel_plot)), function(i) {
  r <- funnel_plot[i, ]
  tibble(
    level = r$level,
    x = c(-r$half_width, r$half_width, r$next_half_width, -r$next_half_width),
    y = c(r$y + 0.43, r$y + 0.43, r$y - 0.43, r$y - 0.43),
    fill = r$fill
  )
}))

p_b <- ggplot() +
  geom_polygon(
    data = funnel_poly,
    aes(x = x, y = y, group = level, fill = fill),
    color = "white", linewidth = 0.55
  ) +
  scale_fill_identity() +
  geom_text(
    data = funnel_plot,
    aes(x = 0, y = y, label = paste0("L", level, "   ", strict_sequential_full_count, "/", denominator)),
    family = font_family, size = 2.75, fontface = "bold", color = palette[["ink"]]
  ) +
  geom_segment(
    data = funnel_plot,
    aes(x = half_width + 0.15, xend = 5.25, y = y, yend = y),
    linewidth = 0.3, color = palette[["neutral_light"]]
  ) +
  geom_text(
    data = funnel_plot %>% mutate(short_criterion = c(
      "Expression claim present",
      "Endpoint recoverable",
      "Genes + directions + formula",
      "Paper-code endpoint concordance",
      "Sampling context fully auditable",
      "Portable without refitting"
    )),
    aes(x = 5.45, y = y, label = short_criterion),
    hjust = 0, family = font_family, size = 2.25, color = palette[["ink"]]
  ) +
  annotate(
    "text", x = 2.55, y = 0.27,
    label = "Strict sequential counts: full satisfaction at the current and every preceding level.\nPartial evidence remains in the study passport and is not counted as a pass.",
    family = font_family, size = 2.05, lineheight = 0.96, color = palette[["neutral"]]
  ) +
  annotate(
    "label", x = 6.7, y = 0.72,
    label = "0 tools passed all six levels",
    family = font_family, fontface = "bold", size = 2.35,
    fill = palette[["red_light"]], color = palette[["red"]], linewidth = 0.3
  ) +
  coord_cartesian(xlim = c(-5.7, 11.5), ylim = c(0, 6.7), clip = "off") +
  labs(
    title = "Six-level reproducibility funnel",
    subtitle = "Endpoint specificity and direct portability progressively narrow the eligible tool set",
    tag = "b"
  ) +
  theme_nature_schematic(7.2)

# Panel c: non-transitive endpoint-construct triangle, no causal arrows.
nodes <- tibble(
  x = c(0.5, 0.16, 0.84), y = c(0.82, 0.18, 0.18),
  label = c(
    "Ex vivo iodine\ncontent / handling",
    "Clinical lesion\nRAI uptake",
    "Structural treatment\nresponse / persistence"
  ),
  fill = c(palette[["teal_light"]], palette[["blue_light"]], palette[["orange_light"]]),
  outline = c(palette[["teal"]], palette[["blue"]], palette[["orange"]])
)

p_c <- ggplot() +
  geom_curve(
    aes(x = 0.47, y = 0.76, xend = 0.20, yend = 0.26),
    curvature = 0.02, linetype = "dashed", linewidth = 0.7, color = palette[["blue"]]
  ) +
  geom_curve(
    aes(x = 0.25, y = 0.18, xend = 0.75, yend = 0.18),
    curvature = -0.04, linetype = "dotted", linewidth = 0.8, color = palette[["neutral"]]
  ) +
  geom_curve(
    aes(x = 0.54, y = 0.76, xend = 0.80, yend = 0.26),
    curvature = -0.02, linetype = "longdash", linewidth = 0.6, color = palette[["neutral_light"]]
  ) +
  geom_point(
    data = nodes, aes(x = x, y = y, fill = fill, color = outline),
    shape = 21, size = 31, stroke = 0.8
  ) +
  scale_fill_identity() + scale_color_identity() +
  geom_text(
    data = nodes, aes(x = x, y = y, label = label),
    family = font_family, size = 2.35, fontface = "bold", lineheight = 0.95, color = palette[["ink"]]
  ) +
  annotate(
    "label", x = 0.22, y = 0.57,
    label = "Condello-6\ng = -0.11\nno positive alignment",
    family = font_family, size = 2.05, lineheight = 0.93,
    fill = "white", color = palette[["blue"]], linewidth = 0.25
  ) +
  annotate(
    "label", x = 0.50, y = 0.055,
    label = "Uptake ≠ response\nclinical context only",
    family = font_family, size = 2.0, lineheight = 0.93,
    fill = "white", color = palette[["neutral"]], linewidth = 0.2
  ) +
  annotate(
    "label", x = 0.79, y = 0.57,
    label = "?\nUntested in paired\npublic transcriptomes",
    family = font_family, size = 2.0, lineheight = 0.93,
    fill = "white", color = palette[["neutral"]], linewidth = 0.2
  ) +
  annotate(
    "text", x = 0.5, y = -0.06,
    label = "Parallel measurement constructs - not biological stages or a causal sequence",
    family = font_family, size = 2.05, color = palette[["neutral"]]
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1), clip = "off") +
  labs(
    title = "Three measurement constructs form an open triangle",
    subtitle = "Their transcriptomic associations cannot be assumed interchangeable",
    tag = "c"
  ) +
  theme_nature_schematic(7.2)

figure1 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.38, 0.62), widths = c(1.36, 1.0)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

base_file <- file.path(out_dir, "Figure1_reproducibility_endpoint_architecture")

svglite::svglite(paste0(base_file, ".svg"), width = width_in, height = height_in)
print(figure1)
dev.off()

grDevices::cairo_pdf(paste0(base_file, ".pdf"), width = width_in, height = height_in, family = font_family)
print(figure1)
dev.off()

ragg::agg_tiff(
  paste0(base_file, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, compression = "lzw"
)
print(figure1)
dev.off()

ragg::agg_png(
  paste0(base_file, ".png"), width = width_in, height = height_in,
  units = "in", res = 600
)
print(figure1)
dev.off()

figure_contract <- c(
  "# Figure 1 contract and manifest",
  "",
  "- Core conclusion: The frozen literature contains RAIR-associated transcriptomic candidates, but full endpoint definition, formula recovery, code concordance, sampling provenance and direct portability progressively collapse; ex vivo iodine content, clinical uptake and structural response remain distinct measurement constructs.",
  "- Archetype: schematic-led composite.",
  "- Hero evidence: six-level strict sequential reproducibility funnel.",
  "- Coverage evidence: PRISMA-like frozen search flow.",
  "- Conceptual synthesis: open non-transitive three-construct triangle.",
  "- Backend: R only (ggplot2 and patchwork).",
  paste0("- Final size: ", width_mm, " x ", height_mm, " mm."),
  "- Exports: editable SVG, vector PDF, 600 dpi TIFF and 600 dpi PNG.",
  "- Source data: results/09_figure1_*.tsv.",
  "- Counting boundary: strict full satisfaction at the current and every preceding level; partial evidence remains in passports.",
  "- Interpretive boundary: zero tools passing all six levels does not show that candidate markers or biological associations are false.",
  "- Endpoint boundary: triangle edges are non-causal; Boucai is clinical context, not validation.",
  "- Search boundary: Semantic Scholar retrieved the top 100 of 242 relevance-ranked records and was not treated as exhaustive."
)
writeLines(figure_contract, file.path(out_dir, "Figure1_figure_contract.md"))

legend <- paste0(
  "## Figure 1 legend\n\n",
  "**Fig. 1 | Systematic recovery audit and endpoint architecture of RAIR-associated transcriptomic tools.** ",
  "**a,** PRISMA-like flow of the frozen systematic evidence-and-reproducibility search. PubMed, Europe PMC and the retrieved Semantic Scholar sensitivity set yielded 959 raw article records and 797 unique database records after deduplication. Citation expansion added 95 candidates, producing 892 records for title/abstract screening; 22 studies underwent full-text and/or code-level audit. GEO datasets (n = 27) are reported separately because they are datasets rather than article records. ",
  "**b,** Six-level strict sequential reproducibility funnel. Counts require full satisfaction of the displayed level and every preceding level; partially recoverable evidence remains in the study passports but is not counted as a pass. Of 22 audited studies, 15 had a fully recoverable clinical or biological endpoint, five additionally supplied fully recoverable genes, directions and formula, one additionally showed paper-code endpoint concordance, and none had fully auditable patient, tissue-site and treatment-timepoint context or passed all six levels. ",
  "**c,** Open non-transitive triangle of three measurement constructs frequently conflated under RAIR terminology. Condello-6 did not positively align with the clinical uptake-failure axis (adverse-oriented Hedges' g = -0.11); the uptake-response edge is clinical literature context, and the ex vivo iodine-response bridge remains untested in paired public transcriptomic data. Lines do not denote causality. The observed gaps may reflect methodological limitations, biological heterogeneity, limited measurement precision or a combination of these factors. Source data are provided as a Source Data file.\n\n"
)

methods <- paste0(
  "## Methods: systematic evidence and reproducibility audit\n\n",
  "The literature search was frozen on 12 July 2026 and was designed as a systematic evidence-and-reproducibility audit rather than an intervention-effect systematic review. Human thyroid-cancer RNA or gene-expression studies were eligible when they explicitly related a gene, directed program, score or expression-derived subgroup to ex vivo iodine content, clinical lesion uptake, structural response or a composite RAIR label. PubMed, Europe PMC and a relevance-ranked Semantic Scholar sensitivity search were combined with citation expansion; PMID, DOI and normalized title were used for deduplication. Semantic Scholar retrieval was limited to the top 100 of 242 results because of API rate limiting and was not treated as exhaustive. GEO dataset records were audited separately from article records.\n\n",
  "Each full-text/code-audited study was evaluated at six prespecified levels: presence of a radioiodine-related expression claim; recoverability of the clinical or biological endpoint; recoverability of genes, directions and formula; paper-code endpoint concordance; auditability of patient overlap, tissue site and treatment timepoint; and direct projection into a matching external data type without refitting. Figure 1b reports strict sequential full-satisfaction counts: a study was retained at a level only if it fully satisfied that level and every preceding level. Partial evidence and non-applicable code fields were retained in the study-level passport but were not counted as passes. The three-construct diagram was an analytical framework rather than a causal model. The Condello-6 edge used the frozen adverse-oriented GSE151179 effect; the uptake-response edge used external clinical context without claiming validation; and the ex vivo iodine-response edge was classified as untested.\n\n"
)

lit_results <- paste0(
  "## A systematic recovery audit reveals progressive loss of endpoint fidelity and portability\n\n",
  "The frozen search identified 959 records across PubMed, Europe PMC and the retrieved Semantic Scholar sensitivity set. Deduplication produced 797 unique database records, and citation expansion added 95 candidates, yielding 892 records for title/abstract screening (Fig. 1a). Twenty-two studies underwent full-text and/or code-level audit, including one 2026 online-first lncRNA study captured only by the supplementary search. The eligible literature contained multiple radioiodine-related candidate markers and programs, but these represented heterogeneous specimens and constructs, including circulating RNA markers of lesion uptake, single tissue markers, targeted expression panels, prognostic scores in RAI-treated or established-RAIR cohorts, lineage/NIS proxies, paired redifferentiation studies and an ex vivo iodine-content program.\n\n",
  "Strict sequential audit showed a steep loss of complete recoverability (Fig. 1b). All 22 studies contained a radioiodine-related expression claim, 15 additionally provided a fully recoverable clinical or biological endpoint, five also supplied recoverable genes, directions and formula, and only one additionally demonstrated paper-code endpoint concordance. None fully documented patient overlap, tissue site and treatment timepoint, and no study satisfied all six levels required for direct endpoint-matched tissue-transcriptome portability without refitting. This result does not indicate an absence of candidate markers. It indicates that no recovered tool simultaneously supplied endpoint specificity, formula recovery, paper-code concordance, auditable sampling context and direct external portability.\n\n",
  "The audit also separated three measurements frequently grouped under RAIR terminology: ex vivo iodine content or iodide-handling capacity, clinical lesion uptake and structural response despite retained uptake (Fig. 1c). Condello-6, trained against ex vivo iodine concentration, did not positively align with the clinical uptake-failure contrast in GSE151179 (adverse-oriented Hedges' g = -0.11). The distinction between uptake and structural response was retained as clinical context rather than external validation, and no paired public transcriptomic dataset tested the bridge between ex vivo iodine content and structural treatment response. Together, these findings indicate that the three measurements capture distinct transcriptomic associations that cannot be assumed interchangeable. The observed gaps may reflect methodological limitations, genuine biological heterogeneity, limited measurement precision or a combination of these factors; current public data cannot adjudicate among these possibilities. These observations do not establish that the original biological conclusions are false or imply research misconduct.\n\n"
)

fig5_results <- paste0(
  "## Radioiodine-related programs show reproducible lineage and immune-compartment preferences across PTC atlases\n\n",
  "All three primary programs satisfied the prespecified cross-donor and three-atlas AND gate (Fig. 5). TDS-16 preferentially mapped to thyroid cells in GSE184362 [mean paired donor difference 1.61, 95% bootstrap confidence interval (CI) 1.29 to 1.95], GSE191288 (1.82, 95% CI 1.42 to 2.15) and GSE281736 (1.93, 95% CI 1.49 to 2.48). Iodide handling-11 showed the same thyroid-cell preference in GSE184362 (1.80, 95% CI 1.48 to 2.12), GSE191288 (2.21, 95% CI 1.91 to 2.50) and GSE281736 (2.23, 95% CI 1.83 to 2.77).\n\n",
  "In contrast, Condello-6 preferentially mapped to T/NK cells in GSE184362 (0.56, 95% CI 0.40 to 0.74), GSE191288 (1.00, 95% CI 0.84 to 1.14) and GSE281736 (0.89, 95% CI 0.75 to 1.02). T/NK cells were the highest-scoring broad compartment for Condello-6 in each atlas. This reproducible immune-compartment preference challenges a thyroid-lineage-intrinsic interpretation of Condello-6 and indicates sensitivity to immune-cell composition. Because the original bulk cohort and the single-cell atlases were unpaired, these analyses do not establish that immune infiltration caused the original bulk association. E2F targets and G2M checkpoint were not estimable in GSE191288 because only two donors contributed QC-qualified B-cell pseudobulk units; these entries were retained as not estimable rather than interpreted as direction reversals (Supplementary Fig. S5).\n"
)

writeLines(c(legend, methods, lit_results), file.path(manuscript_dir, "Figure1_methods_results_legend.md"))
writeLines(
  c(
    "# Results", "",
    lit_results,
    "<!-- Results sections corresponding to Figs. 2-4 will be inserted here after final figure freeze. -->", "",
    fig5_results
  ),
  file.path(manuscript_dir, "RESULTS_DRAFT.md")
)

cat("Figure 1 generated and Results draft assembled.\n")
cat("Strict funnel counts:", paste(funnel_counts, collapse = " -> "), "\n")
cat("Outputs:\n", paste0(base_file, c(".svg", ".pdf", ".tiff", ".png"), collapse = "\n"), "\n")
