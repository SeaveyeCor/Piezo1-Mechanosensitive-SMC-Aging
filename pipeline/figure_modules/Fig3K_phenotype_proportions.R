# ==============================================================================
# PAPER PANEL: Fig 3B (Table S4c)   [internal module ID: Fig3K]
# PLOTS:  % of visceral SMCs classified contractile vs synthetic in each group
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Young vs Aged CTRL (aging); Aged CTRL vs Aged P1KO (rescue); Young CTRL vs Aged P1KO
# TEST:   chi-square test on contractile/synthetic counts
# ==============================================================================
# Stacked bar plot of contractile vs synthetic SMC percentages
# Style: theme_minimal, percentage labels, panel border
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate Fig3K
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/statistics.R"))
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))

  group_labels <- get_group_labels()

  # Colors
  phenotype_colors <- c("Contractile" = "#E31A1C", "Synthetic" = "#1F78B4")
  group_order <- c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  phenotype <- context$phenotype
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Merge phenotype with metadata
  pheno_with_meta <- merge(
    phenotype,
    metadata[, c("group"), drop = FALSE],
    by.x = "cell_barcode",
    by.y = "row.names"
  )

  # Calculate proportions by group
  proportions <- pheno_with_meta %>%
    group_by(group, phenotype) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(
      total = sum(n),
      pct = 100 * n / total,
      sem = sqrt(pct * (100 - pct) / total)  # binomial SEM approximation
    ) %>%
    ungroup() %>%
    mutate(group = factor(group, levels = group_order))

  # Add genotype for coloring border
  proportions$genotype <- ifelse(grepl("KO", proportions$group), "KO", "Ctrl")

  # Statistical tests
  # Compare Aged Ctrl vs Young Ctrl (contractile proportion)
  young_ctrl_data <- pheno_with_meta %>% filter(group == "Young_Ctrl")
  aged_ctrl_data <- pheno_with_meta %>% filter(group == "Aged_Ctrl")
  aged_ko_data <- pheno_with_meta %>% filter(group == "Aged_KO")

  # Chi-square tests for proportion differences
  # Test 1: Young Ctrl vs Aged Ctrl (aging effect)
  test_aging <- tryCatch({
    table_aging <- table(
      c(rep("Young", nrow(young_ctrl_data)), rep("Aged", nrow(aged_ctrl_data))),
      c(young_ctrl_data$phenotype, aged_ctrl_data$phenotype)
    )
    chisq.test(table_aging)
  }, error = function(e) list(p.value = NA))

  # Test 2: Aged Ctrl vs Aged KO (rescue effect)
  test_rescue <- tryCatch({
    table_rescue <- table(
      c(rep("Ctrl", nrow(aged_ctrl_data)), rep("KO", nrow(aged_ko_data))),
      c(aged_ctrl_data$phenotype, aged_ko_data$phenotype)
    )
    chisq.test(table_rescue)
  }, error = function(e) list(p.value = NA))

  # Test 3: Young Ctrl vs Aged KO
  test_yc_ako <- tryCatch({
    table_yc_ako <- table(
      c(rep("YoungCtrl", nrow(young_ctrl_data)), rep("AgedKO", nrow(aged_ko_data))),
      c(young_ctrl_data$phenotype, aged_ko_data$phenotype)
    )
    chisq.test(table_yc_ako)
  }, error = function(e) list(p.value = NA))

  # Create plot
  p <- ggplot(proportions, aes(x = group, y = pct, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.5) +
    # Percentage labels inside bars
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              position = position_stack(vjust = 0.5),
              color = "white", fontface = "bold", size = 6) +
    scale_fill_manual(values = phenotype_colors) +
    scale_x_discrete(labels = c(
      "Young_Ctrl" = "3",
      "Aged_Ctrl"  = "24",
      "Young_KO"   = "3",
      "Aged_KO"    = "24"
    )) +
    labs(
      title = "SMC Phenotype Composition by Group",
      x = NULL,
      y = "Percentage of Cells (%)",
      fill = "Phenotype"
    ) +
    # Theme (publication sizing)
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(face = "bold", size = 18),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 125))

  # Bracket comparisons with significance
  # Helper to draw a bracket: horizontal line with short vertical ticks at each end
  add_bracket <- function(p, x1, x2, y, label, size = 6, tick_height = 1.5) {
    p <- p +
      # Horizontal bar
      annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = 0.5) +
      # Left tick
      annotate("segment", x = x1, xend = x1, y = y - tick_height, yend = y, linewidth = 0.5) +
      # Right tick
      annotate("segment", x = x2, xend = x2, y = y - tick_height, yend = y, linewidth = 0.5) +
      # Label above bracket
      annotate("text", x = (x1 + x2) / 2, y = y + 1.5, label = label, size = size, fontface = "bold")
    return(p)
  }

  # Bracket 1: 3 CTRL vs 24 CTRL (positions 1 and 2) — aging
  sig_aging <- format_significance(test_aging$p.value)
  if (!is.na(test_aging$p.value)) {
    p <- add_bracket(p, 1, 2, 104, sig_aging)
  }

  # Bracket 2: 24 CTRL vs 24 KO (positions 2 and 4) — rescue
  sig_rescue <- format_significance(test_rescue$p.value)
  if (!is.na(test_rescue$p.value)) {
    p <- add_bracket(p, 2, 4, 112, sig_rescue)
  }

  # Bracket 3: 3 CTRL vs 24 KO (positions 1 and 4) — overall
  sig_yc_ako <- format_significance(test_yc_ako$p.value)
  if (!is.na(test_yc_ako$p.value)) {
    p <- add_bracket(p, 1, 4, 120, sig_yc_ako)
  }

  # Create output directory
  fig_dir <- file.path(context$output_dir, "Fig3K")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data for export
  export_data <- list(
    proportions = proportions,
    statistical_tests = data.frame(
      comparison = c("Aging (Aged Ctrl vs Young Ctrl)", "Rescue (Aged KO vs Aged Ctrl)", "Young Ctrl vs Aged KO"),
      test = c("Chi-square", "Chi-square", "Chi-square"),
      p_value = c(test_aging$p.value, test_rescue$p.value, test_yc_ako$p.value),
      significance = c(format_significance(test_aging$p.value),
                       format_significance(test_rescue$p.value),
                       format_significance(test_yc_ako$p.value))
    ),
    cell_counts = pheno_with_meta %>%
      group_by(group, phenotype) %>%
      summarize(n = n(), .groups = "drop")
  )

  # Methods text
  methods_text <- paste(
    "SMC phenotype proportions (Contractile vs Synthetic) were calculated for each",
    "experimental group. Phenotype classification was based on the ratio of",
    "contractile to synthetic marker expression (threshold > 1.5 for contractile).",
    "Statistical significance was assessed using chi-square tests comparing:",
    "(1) Aged Ctrl vs Young Ctrl for aging effects,",
    "(2) Aged KO vs Aged Ctrl for Piezo1 KO rescue effects, and",
    "(3) Aged KO vs Young Ctrl for overall rescue assessment."
  )

  technical_legend <- paste(
    "Stacked bar plot showing the percentage of SMCs classified as Contractile (red)",
    "or Synthetic (blue) in each experimental group.",
    "X-axis shows age in months (3 = young, 24 = aged).",
    "Statistical significance indicated by: * p<0.05, ** p<0.01, *** p<0.001."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = "Fig3K",
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 7,
    height = 7,
    code_path = "pipeline/figure_modules/Fig3K_phenotype_proportions.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    proportions = proportions,
    tests = export_data$statistical_tests
  ))
}
