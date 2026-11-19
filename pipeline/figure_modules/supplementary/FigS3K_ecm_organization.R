# ==============================================================================
# PAPER PANEL: Fig 3M (bar version of the score; the printed violin uses the same per-cell score, drawn by 04_cell_level/75_score_superplot_inputs.R + 02_superplots.R)   [internal module ID: FigS3K]
# PLOTS:  ECM score per cell (Bgn, Col6a1, Dcn, Sparc, Timp2)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Wilcoxon (Mann-Whitney), Bonferroni x4
# ==============================================================================
# Bar plot of the unified ECM Structural score across groups.
#
# ACTIVE GENE SET:
#   Bgn, Col6a1, Dcn, Sparc, Timp2  (5 genes)
#
# Structural ECM genes only. The panel shares no gene with the NFAT-associated remodeling
# score (Wnt5a, Mmp2, Mmp9, Vcam1, Timp1), so the ECM and NFAT scores are not double-counted.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate FigS3K
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
  group_order <- get_group_order()
  group_labels <- get_group_labels()

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Get expression data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Calculate ECM Structural score (unified 5-gene ECM panel)
  ecm_genes <- gene_sets$ecm_structural
  ecm_avail <- ecm_genes[ecm_genes %in% rownames(expr)]

  cat("  ECM Score genes:", paste(ecm_avail, collapse = ", "), "\n")

  ecm_score <- colMeans(as.matrix(expr[ecm_avail, smc_cells, drop = FALSE]))

  # Create data frame
  plot_data <- data.frame(
    cell_barcode = smc_cells,
    ecm_score = ecm_score,
    group = metadata[smc_cells, "group"],
    stringsAsFactors = FALSE
  )

  # Calculate summary statistics
  summary_stats <- plot_data %>%
    group_by(group) %>%
    summarize(
      n = n(),
      mean = mean(ecm_score, na.rm = TRUE),
      sd = sd(ecm_score, na.rm = TRUE),
      sem = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      group = factor(group, levels = group_order),
      genotype = ifelse(grepl("KO", group), "KO", "Ctrl")
    )

  # Pairwise Mann-Whitney — 4 valid comparisons, Bonferroni-adjusted
  pairwise_tests <- pairwise_group_tests(plot_data$ecm_score, plot_data$group)

  # Create plot
  p <- ggplot(summary_stats, aes(x = group, y = mean, fill = genotype)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = colors$genotype) +
    scale_x_discrete(labels = group_labels) +
    labs(
      title = "ECM Score",
      x = NULL,
      y = "ECM Score"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # Stack pairwise brackets above the bars
  y_top <- max(summary_stats$mean + summary_stats$sem)
  p <- p + expand_limits(y = y_top * 1.6)
  p <- add_pairwise_brackets(p, pairwise_tests, group_order, y_top)

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data for export
  export_data <- list(
    summary_statistics = summary_stats,
    statistical_tests = pairwise_tests,
    genes_used = data.frame(gene = ecm_avail)
  )

  # Methods text
  methods_text <- paste(
    "ECM Score was computed as the mean normalized expression of 5 structural ECM",
    "genes (Bgn, Col6a1, Dcn, Sparc, Timp2). All five genes individually meet the",
    "data-validated criteria of >20% increase with aging and >20% decrease with",
    "Piezo1 KO (both p<0.05): Bgn (+50% aging, -21% rescue), Col6a1 (+29%, -31%),",
    "Dcn (+54%, -28%), Sparc (+26%, -42%), Timp2 (+23%, -40%). This unified panel",
    "is non-overlapping with the NFAT Activity Score gene set (Wnt5a, Mmp2, Mmp9,",
    "Vcam1, Fas, Timp1) — Mmp2/Mmp9/Timp1 live exclusively under NFAT Activity —",
    "to avoid statistical double-counting between the two headline scores.",
    "Group differences were tested by two-sided Mann-Whitney U with Bonferroni",
    "correction (m = 4 valid pairwise comparisons)."
  )

  technical_legend <- paste(
    "Bar plot showing mean ECM Score (+/- SEM) in bona fide SMCs across",
    "experimental groups. Score = mean normalized expression of 5 structural ECM",
    "genes: Bgn, Col6a1, Dcn, Sparc, Timp2. X-axis: age in months (3 = young,",
    "24 = aged). Grey bars = Control, Pink bars = Piezo1 KO. Significance",
    "brackets: Bonferroni-adjusted (m=4) Mann-Whitney U; non-significant",
    "comparisons are not shown. * p<0.05, ** p<0.01, *** p<0.001."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 6,
    height = 5,
    code_path = "pipeline/figure_modules/supplementary/FigS3K_ecm_organization.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    summary = summary_stats
  ))
}
