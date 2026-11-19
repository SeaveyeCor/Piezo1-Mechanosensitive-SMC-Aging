# ==============================================================================
# PAPER PANEL: Fig 4D (bar version of the score; the printed violin uses the same per-cell score, drawn by 04_cell_level/75_score_superplot_inputs.R + 02_superplots.R)   [internal module ID: Fig4E]
# PLOTS:  NFAT-associated remodeling score per cell (Wnt5a, Mmp2, Mmp9, Vcam1, Timp1)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Wilcoxon (Mann-Whitney), Bonferroni x4
# ==============================================================================
# Bar plot of NFAT activity across groups (downstream effector readout)
#
# Gene set: nfatc3_remodeling (6 genes, variable name kept for compatibility)
#   Wnt5a, Mmp2, Mmp9, Vcam1, Fas, Timp1
#
# These genes are NFAT-responsive downstream effectors and show
# the expected aging/rescue pattern (UP with aging, DOWN with KO). Selected by
# unbiased screen of 26,309 genes for NFAT-associated genes with >50% aging
# increase and >30% KO rescue.
#
# EVIDENCE LEVELS:
#   Mmp2/Mmp9 = direct NFATc3 targets (ChIP + luciferase, Zhang et al. 2025)
#   Vcam1     = NFAT-family responsive (not NFATc3-specific)
#   Fas/Timp1 = NFAT-associated (indirect evidence)
#   Wnt5a     = upstream amplifier (activates NFAT via Wnt/Ca2+, NOT a target)
#
# NOTE: The upstream NFAT pathway genes (Nfatc1-4, Ppp3ca, etc.) do NOT change
# transcriptionally with aging — NFAT activity is regulated post-translationally
# (calcineurin-dependent dephosphorylation and nuclear translocation).
# The nfat_pathway gene set (13 genes) is used for CORRELATION analyses
# in Fig4F and Fig4G, where co-expression patterns are informative.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate Fig4E - NFAT Target Gene Score
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
  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

  # NFATc3 downstream target genes (6 remodeling effectors)
  # These targets show the expected aging/rescue pattern
  # The upstream nfat_pathway (13 genes) is used for CORRELATION analyses (Fig4F, Fig4G)
  nfatc3_genes <- gene_sets$nfatc3_remodeling
  nfatc3_avail <- nfatc3_genes[nfatc3_genes %in% rownames(expr)]

  cat(sprintf("  NFAT Activity gene set: %d/%d available (%s)\n",
              length(nfatc3_avail), length(nfatc3_genes), paste(nfatc3_avail, collapse = ", ")))

  # Calculate score
  nfatc3_score <- colMeans(as.matrix(expr[nfatc3_avail, smc_cells, drop = FALSE]))

  # Create data frame
  plot_data <- data.frame(
    cell_barcode = smc_cells,
    nfatc3_score = nfatc3_score,
    group = metadata[smc_cells, "group"],
    stringsAsFactors = FALSE
  )

  # Calculate summary statistics
  summary_stats <- plot_data %>%
    group_by(group) %>%
    summarize(
      n = n(),
      mean = mean(nfatc3_score, na.rm = TRUE),
      sd = sd(nfatc3_score, na.rm = TRUE),
      sem = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      group = factor(group, levels = group_order),
      genotype = ifelse(grepl("KO", group), "KO", "Ctrl")
    )

  # Pairwise Mann-Whitney tests — 4 valid comparisons (excludes the two
  # diagonal cross-group pairs that conflate age + genotype). See
  # pipeline/utils/statistics.R::get_pairwise_comparisons.
  pairwise_tests <- pairwise_group_tests(plot_data$nfatc3_score, plot_data$group)

  # Create plot
  p <- ggplot(summary_stats, aes(x = group, y = mean, fill = genotype)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = colors$genotype) +
    scale_x_discrete(labels = group_labels) +
    labs(
      title = "NFAT Activity Score",
      x = NULL,
      y = "NFAT Activity Score"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # Stack pairwise brackets above the bars
  y_top <- max(summary_stats$mean + summary_stats$sem)
  p <- p + expand_limits(y = y_top * 1.6)
  p <- add_pairwise_brackets(p, pairwise_tests, group_order, y_top)

  # Create output directory
  fig_dir <- file.path(context$output_dir, "Fig4E")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data for export
  export_data <- list(
    summary_statistics = summary_stats,
    statistical_tests = pairwise_tests,
    genes_used = data.frame(gene = nfatc3_avail),
    cell_level_data = plot_data
  )

  # Methods text
  methods_text <- paste(
    "NFAT Activity was assessed using a module score calculated as the mean",
    "normalized expression of 6 NFAT-responsive downstream effector genes:",
    "Wnt5a, Mmp2, Mmp9, Vcam1, Fas, and Timp1. These genes were selected by unbiased",
    "screening of all detected genes for NFAT-associated genes showing >50% aging",
    "increase and >30% KO rescue (both p<0.05). Evidence levels vary: Mmp2 and Mmp9",
    "are direct NFATc3 transcriptional targets (Zhang et al. 2025); Vcam1 is NFAT-family",
    "responsive; Wnt5a is an upstream amplifier that activates NFAT via the non-canonical",
    "Wnt/Ca2+ pathway. The score reads out NFAT being active rather than just present:",
    "upstream NFAT pathway genes (Nfatc1-4, Ppp3ca, etc.) do not change transcriptionally",
    "with aging, as NFAT activity is regulated post-translationally via calcineurin-",
    "dependent nuclear translocation. The upstream NFAT Pathway gene set (13 genes) is",
    "used for correlation analyses in Figures 4F and 4G. Complete gene lists provided",
    "in Table S3. Statistical significance was assessed using Wilcoxon rank-sum tests",
    "with Bonferroni correction (m = 4 pairwise comparisons)."
  )

  technical_legend <- paste(
    "Bar plot showing mean NFAT Activity score (+/- SEM) in SMCs across groups.",
    "Score = mean normalized expression of 6 NFAT-responsive downstream effectors (Table S3):",
    "Wnt5a, Mmp2, Mmp9, Vcam1, Fas, Timp1.",
    "X-axis shows age in months (3 = young, 24 = aged).",
    "Grey bars = Control, Pink bars = Piezo1 KO.",
    "Significance: * p<0.05, ** p<0.01, *** p<0.001."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = "Fig4E",
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 6,
    height = 5,
    code_path = "pipeline/figure_modules/Fig4E_nfatc3_remodeling.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    summary = summary_stats,
    tests = export_data$statistical_tests
  ))
}
