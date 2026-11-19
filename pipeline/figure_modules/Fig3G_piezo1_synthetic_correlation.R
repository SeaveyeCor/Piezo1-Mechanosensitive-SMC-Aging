# ==============================================================================
# PAPER PANEL: Fig 3D (printed panel drawn by 04_cell_level/71_score_panels.R)   [internal module ID: Fig3G]
# PLOTS:  Piezo1 expression vs synthetic phenotype score (Vim, Col1a2, Dcn, Igfbp6, Efemp1) per cell
# CELLS:  CTRL visceral SMCs only (Young CTRL n = 469, Aged CTRL n = 477)
# GROUPS COMPARED: Young CTRL and Aged CTRL, each tested separately
# TEST:   Spearman correlation
# ==============================================================================
# Scatter plot with regression lines for Young Ctrl and Aged Ctrl
# Uses synthetic_phenotype gene set: Vim, Col1a2, Dcn, Igfbp6, Efemp1
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate Fig3G
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "utils/statistics.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))

  colors <- get_colors()
  gene_sets <- get_gene_sets()

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Get expression data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Filter to control SMCs
  ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$genotype == "Ctrl"]]

  # Get Piezo1 expression
  piezo1_expr <- if ("Piezo1" %in% rownames(expr)) {
    as.numeric(expr["Piezo1", ctrl_cells])
  } else {
    rep(NA, length(ctrl_cells))
  }

  # Calculate synthetic phenotype score
  synthetic_genes <- gene_sets$synthetic_phenotype
  synthetic_avail <- synthetic_genes[synthetic_genes %in% rownames(expr)]
  synthetic_score <- colMeans(as.matrix(expr[synthetic_avail, ctrl_cells, drop = FALSE]))

  cat(sprintf("  Synthetic phenotype genes: %d/%d available (%s)\n",
              length(synthetic_avail), length(synthetic_genes),
              paste(synthetic_avail, collapse = ", ")))

  # Create data frame
  plot_data <- data.frame(
    cell_barcode = ctrl_cells,
    piezo1 = piezo1_expr,
    synthetic_score = synthetic_score,
    age = metadata[ctrl_cells, "age"],
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(piezo1) & !is.na(synthetic_score))

  # Calculate correlations by age
  cor_young <- spearman_cor(
    plot_data$synthetic_score[plot_data$age == "Young"],
    plot_data$piezo1[plot_data$age == "Young"]
  )

  cor_aged <- spearman_cor(
    plot_data$synthetic_score[plot_data$age == "Aged"],
    plot_data$piezo1[plot_data$age == "Aged"]
  )

  cat(sprintf("  Young: rho = +%.4f, p = %.2e (n=%d)\n",
              cor_young$rho, cor_young$p_value, cor_young$n))
  cat(sprintf("  Aged:  rho = +%.4f, p = %.2e (n=%d)\n",
              cor_aged$rho, cor_aged$p_value, cor_aged$n))

  # Create plot - LINE OF BEST FIT ONLY (matching Fig3F convention)
  p <- ggplot(plot_data, aes(x = synthetic_score, y = piezo1, color = age)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.5) +
    scale_color_manual(
      values = c("Young" = "#4DAF4A", "Aged" = "#E41A1C"),
      labels = c("Young" = sprintf("3mo (rho=%.3f)", cor_young$rho),
                 "Aged"  = sprintf("24mo (rho=%.3f)", cor_aged$rho)),
      name = "Group"
    ) +
    labs(
      title = "Piezo1 vs Synthetic Phenotype Score",
      x = "Synthetic Phenotype Score",
      y = "Piezo1 Expression"
    ) +
    theme_publication() +
    theme(
      legend.position = "none"
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, "Fig3G")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data for export
  export_data <- list(
    correlation_summary = data.frame(
      age = c("Young", "Aged"),
      rho = c(cor_young$rho, cor_aged$rho),
      p_value = c(cor_young$p_value, cor_aged$p_value),
      n = c(cor_young$n, cor_aged$n)
    ),
    cell_level_data = plot_data,
    genes_used = data.frame(
      gene_set = "synthetic_phenotype",
      gene = synthetic_avail
    )
  )

  # Methods text
  methods_text <- paste(
    "Spearman correlation was calculated between Piezo1 expression and synthetic",
    "phenotype score in control SMCs, stratified by age (Young = 3 months, Aged = 24 months).",
    "The synthetic phenotype score was calculated as the mean normalized expression of",
    "five fibroblast-like markers (Vim, Col1a2, Dcn, Igfbp6, Efemp1) capturing the",
    "contractile-to-synthetic phenotype switch. This gene set has zero overlap with the",
    "NFAT Activity score, enabling non-autocorrelated correlation analyses. All five genes",
    "show significant aging increase and Piezo1-KO rescue."
  )

  technical_legend <- paste(
    "Scatter plot with linear regression lines showing the correlation between",
    "Piezo1 expression (y-axis) and synthetic phenotype score (x-axis) in control SMCs.",
    "Green = Young (3 mo), Red = Aged (24 mo).",
    "Shaded regions represent 95% confidence intervals.",
    "Spearman rho and p-values shown as text annotations.",
    "Synthetic phenotype score: Vim, Col1a2, Dcn, Igfbp6, Efemp1."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = "Fig3G",
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 7,
    height = 6,
    code_path = "pipeline/figure_modules/Fig3G_piezo1_synthetic_correlation.R"
  )

  # Save legend as separate PNG (sized for manuscript: legend ~1/3-1/6 of plot)
  p_for_legend <- p +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm")) +
    guides(color = guide_legend(override.aes = list(linewidth = 4)))
  legend_grob <- cowplot::get_legend(p_for_legend)
  ggsave(file.path(fig_dir, "Fig3G_legend.png"), legend_grob,
         width = 6, height = 3, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig3G_legend.pdf"), legend_grob,
         width = 6, height = 3, bg = "white")
  cat("  Saved: Fig3G_legend.png/pdf\n")

  # Save stats as separate PNG
  format_p_local <- function(pval) {
    if (is.na(pval)) return("NA")
    if (pval < 0.001) return("p<0.001")
    if (pval < 0.01) return(sprintf("p=%.3f", pval))
    return(sprintf("p=%.2f", pval))
  }
  stats_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.7,
             label = sprintf("3mo:  rho = +%.3f  %s", cor_young$rho, cor_young$significance),
             color = "#4DAF4A", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = 0.5, y = 0.3,
             label = sprintf("24mo:  rho = +%.3f  %s", cor_aged$rho, cor_aged$significance),
             color = "#E41A1C", size = 5, fontface = "bold", hjust = 0.5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(fig_dir, "Fig3G_stats.png"), stats_plot,
         width = 5, height = 1.5, dpi = 300, bg = "white")
  cat("  Saved: Fig3G_stats.png\n")

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    correlations = export_data$correlation_summary
  ))
}
