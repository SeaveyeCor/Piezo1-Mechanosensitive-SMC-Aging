# ==============================================================================
# PAPER PANEL: Fig 3N (bar version of the score; the printed violin uses the same per-cell score, drawn by 04_cell_level/75_score_superplot_inputs.R + 02_superplots.R)   [internal module ID: FigS4H]
# PLOTS:  cytoskeletal remodeling score per cell (Ezr, Syne2, Lmna, Vim)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Wilcoxon (Mann-Whitney), Bonferroni x4
# ==============================================================================
# 4-group bar plot of the Cytoskeletal Remodeling score across the same SMC
# population used by Fig4E (NFAT-Resp Remodeling) and FigS3K (ECM Structural).
#
# Gene set: cyto_remodeling (4 genes; pipeline/core/gene_sets.R)
#   Ezr   — Ezrin, actin-membrane crosslinker, mesenchymal
#   Syne2 — Nesprin-2, nuclear envelope linker, cytoskeletal coupling
#   Lmna  — Lamin A/C, nuclear lamina, mechanotransduction
#   Vim   — Vimentin, intermediate filament, mesenchymal/synthetic marker
#
# These 4 genes were validated as the cytoskeletal step in the Fig4I pseudotime
# cascade (cascade_cyto_remodeling). They sit between the contractile state and
# the synthetic state on the dedifferentiation trajectory and are upregulated
# with aging and rescued by Piezo1 KO in the cluster-based SMC population.
#
# CAVEAT: Vim overlaps with the synthetic_phenotype score (Vim, Col1a2, Dcn,
# Igfbp6, Efemp1). This is intentional — Vim is a mesenchymal marker that drives
# both phenotypes — but worth noting if comparing absolute scores across the
# cyto and synthetic panels.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate FigS4H — Cytoskeletal Remodeling Score
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

  expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

  # Calculate Cytoskeletal Remodeling score
  cyto_genes <- gene_sets$cyto_remodeling
  cyto_avail <- cyto_genes[cyto_genes %in% rownames(expr)]

  cat("  Cyto Remodeling genes:", paste(cyto_avail, collapse = ", "), "\n")

  cyto_score <- colMeans(as.matrix(expr[cyto_avail, smc_cells, drop = FALSE]))

  plot_data <- data.frame(
    cell_barcode = smc_cells,
    cyto_score = cyto_score,
    group = metadata[smc_cells, "group"],
    stringsAsFactors = FALSE
  )

  # Summary statistics
  summary_stats <- plot_data %>%
    group_by(group) %>%
    summarize(
      n = n(),
      mean = mean(cyto_score, na.rm = TRUE),
      sd = sd(cyto_score, na.rm = TRUE),
      sem = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      group = factor(group, levels = group_order),
      genotype = ifelse(grepl("KO", group), "KO", "Ctrl")
    )

  # Pairwise Mann-Whitney — 4 valid comparisons, Bonferroni-adjusted
  pairwise_tests <- pairwise_group_tests(plot_data$cyto_score, plot_data$group)

  # Per-gene aging-rescue summary for the methods text
  per_gene <- data.frame()
  for (g in cyto_avail) {
    vals <- as.numeric(expr[g, smc_cells])
    yc <- vals[plot_data$group == "Young_Ctrl"]
    ac <- vals[plot_data$group == "Aged_Ctrl"]
    ak <- vals[plot_data$group == "Aged_KO"]
    aging_pct <- 100 * (mean(ac) - mean(yc)) / max(mean(yc), 1e-9)
    rescue_pct <- 100 * (mean(ak) - mean(ac)) / max(mean(ac), 1e-9)
    per_gene <- rbind(per_gene, data.frame(
      gene = g,
      mean_yc = mean(yc), mean_ac = mean(ac), mean_ak = mean(ak),
      aging_pct = aging_pct, rescue_pct = rescue_pct
    ))
  }

  # Create plot
  p <- ggplot(summary_stats, aes(x = group, y = mean, fill = genotype)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = colors$genotype) +
    scale_x_discrete(labels = group_labels) +
    labs(
      title = "Cytoskeletal Remodeling Score",
      x = NULL,
      y = "Cyto Remodeling Score"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # Stack pairwise brackets above the bars (ns suppressed by helper)
  y_top <- max(summary_stats$mean + summary_stats$sem)
  p <- p + expand_limits(y = y_top * 1.6)
  p <- add_pairwise_brackets(p, pairwise_tests, group_order, y_top)

  # Output
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  export_data <- list(
    summary_statistics = summary_stats,
    statistical_tests = pairwise_tests,
    per_gene_summary = per_gene,
    genes_used = data.frame(gene = cyto_avail),
    cell_level_data = plot_data
  )

  methods_text <- paste(
    "Cytoskeletal Remodeling score was computed as the mean normalized",
    "expression of 4 genes (Ezr, Syne2, Lmna, Vim) across bona fide SMCs.",
    "These genes form the cytoskeletal step of the contractile-to-synthetic",
    "pseudotime cascade (Fig 4F) and are upregulated with aging and rescued by",
    "Piezo1 KO in the cluster-based SMC population. Vim overlaps with the",
    "synthetic phenotype score by design (mesenchymal marker driving both",
    "scores). Pairwise group differences were tested by two-sided Mann-Whitney",
    "U tests with Bonferroni correction (m = 4)."
  )

  technical_legend <- paste(
    "Bar plot showing mean Cytoskeletal Remodeling score (+/- SEM) in bona fide",
    "SMCs across experimental groups. Score = mean normalized expression of 4",
    "genes: Ezr, Syne2, Lmna, Vim. X-axis shows age (3 = young, 24 = aged).",
    "Grey = Control, Pink = Piezo1 KO. Significance brackets show",
    "Bonferroni-adjusted Mann-Whitney U p-values for the 4 valid pairwise",
    "comparisons (Aging WT, Aging KO, KO@young, KO rescue@aged).",
    "Non-significant comparisons are not shown.",
    "* p<0.05, ** p<0.01, *** p<0.001."
  )

  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 6,
    height = 5,
    code_path = "pipeline/figure_modules/supplementary/FigS4H_cyto_remodeling.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    summary = summary_stats
  ))
}
