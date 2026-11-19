# ==============================================================================
# PAPER PANEL: Fig 3E (bar version of the score; the printed violin uses the same per-cell score, drawn by 04_cell_level/75_score_superplot_inputs.R + 02_superplots.R)   [internal module ID: Fig3H]
# PLOTS:  contractile score per cell (mean of Myh11, Acta2, Tagln, Cnn1, Mylk, Myocd, Smtn)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Wilcoxon (Mann-Whitney), Bonferroni x4
# ==============================================================================
# Mean expression of 7 canonical SMC contractile markers (contractile_core):
#   Myh11, Acta2, Tagln, Cnn1, Mylk, Myocd, Smtn
# Normalized to Young Control mean (= 100%)
# Bar plot with mean ± SEM, Kruskal-Wallis + Wilcoxon with Bonferroni correction
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

#' Generate Fig3H - Contractile Score
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/statistics.R"))
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))

  colors <- get_colors()
  group_order <- get_group_order()
  group_labels <- get_group_labels()

  gene_sets <- get_gene_sets()

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  fig_dir <- file.path(context$output_dir, "Fig3H")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  Computing Contractile Score (contractile_core, 7 genes)...\n")

  # Use contractile_core gene set
  contractile_genes <- gene_sets$contractile_core
  avail <- contractile_genes[contractile_genes %in% rownames(expr)]
  cat(sprintf("  Genes: %d/%d available (%s)\n",
              length(avail), length(contractile_genes), paste(avail, collapse = ", ")))

  # Calculate mean score per cell
  score <- colMeans(as.matrix(expr[avail, smc_cells, drop = FALSE]))

  # Get cells by group
  groups <- c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")
  cells_by_group <- lapply(groups, function(g) {
    smc_cells[smc_cells %in% rownames(metadata)[metadata$group == g]]
  })
  names(cells_by_group) <- groups

  # Normalize to Young Control mean = 100%
  yc_mean <- mean(score[cells_by_group[["Young_Ctrl"]]])
  score_pct <- (score / yc_mean) * 100

  cat(sprintf("  Young Ctrl mean (raw): %.4f\n", yc_mean))

  # Create data frame
  score_data <- data.frame(
    cell = smc_cells,
    score = score_pct,
    group = metadata[smc_cells, "group"],
    genotype = metadata[smc_cells, "genotype"],
    stringsAsFactors = FALSE
  )
  score_data$group <- factor(score_data$group, levels = groups)

  # Summary statistics
  summary_stats <- score_data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean = mean(score, na.rm = TRUE),
      sd = sd(score, na.rm = TRUE),
      sem = sd / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(genotype = ifelse(grepl("KO", group), "KO", "Ctrl"))

  cat("  Summary by group:\n")
  for (i in 1:nrow(summary_stats)) {
    cat(sprintf("    %s: %.1f ± %.1f (n=%d)\n",
                summary_stats$group[i], summary_stats$mean[i],
                summary_stats$sem[i], summary_stats$n[i]))
  }

  # Omnibus test (kept for reporting)
  kw_test <- kruskal.test(score ~ group, data = score_data)
  cat(sprintf("  Kruskal-Wallis: chi-squared = %.2f, p = %.2e\n",
              kw_test$statistic, kw_test$p.value))

  # Pairwise Mann-Whitney — 4 valid comparisons, Bonferroni-adjusted
  pairwise_tests <- pairwise_group_tests(score_data$score, score_data$group)
  cat("  Pairwise Wilcoxon tests (Bonferroni-adjusted, m=4):\n")
  for (i in seq_len(nrow(pairwise_tests))) {
    cat(sprintf("    %-18s %s vs %s: p_adj = %.2e (%s)\n",
                pairwise_tests$comparison[i],
                pairwise_tests$a[i], pairwise_tests$b[i],
                pairwise_tests$p_adjusted[i],
                pairwise_tests$significance[i]))
  }

  # Create bar plot (matching Fig4E style)
  p <- ggplot(summary_stats, aes(x = group, y = mean, fill = genotype)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem),
                  width = 0.2, linewidth = 0.5) +
    scale_fill_manual(values = colors$genotype) +
    scale_x_discrete(labels = group_labels) +
    labs(
      title = "Contractile Score",
      x = NULL,
      y = "Contractile Score\n(% of Young Control)"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # Stack pairwise brackets above the bars
  y_top <- max(summary_stats$mean + summary_stats$sem)
  p <- p + expand_limits(y = y_top * 1.6)
  p <- add_pairwise_brackets(p, pairwise_tests, groups, y_top)

  # Save
  ggsave(file.path(fig_dir, "Fig3H.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig3H.pdf"), p, width = 7, height = 6, bg = "white")
  cat("  Saved: Fig3H.png\n")
  cat("  Saved: Fig3H.pdf\n")

  # Export data
  export_data <- list(
    summary_statistics = summary_stats,
    kruskal_wallis = data.frame(
      test = "Kruskal-Wallis", statistic = kw_test$statistic,
      df = kw_test$parameter, p_value = kw_test$p.value
    ),
    wilcoxon_tests = pairwise_tests,
    gene_panel = data.frame(gene = avail, gene_set = "contractile_core"),
    cell_level_data = score_data
  )

  methods_text <- paste(
    "Contractile score was calculated as the mean normalized expression of seven",
    "canonical SMC contractile markers (Myh11, Acta2, Tagln, Cnn1, Mylk, Myocd, Smtn;",
    "Table S2). For each cell, the score was normalized to the Young Control mean",
    "(= 100%). Group differences were tested by Kruskal-Wallis, followed by pairwise",
    "Wilcoxon rank-sum tests with Bonferroni correction (m = 4)."
  )

  technical_legend <- paste(
    "Bar plot showing contractile score (mean +/- SEM) across experimental groups.",
    "Score = mean normalized expression of 7 canonical contractile markers as",
    "% of Young Control.",
    "X-axis shows age in months (3 = young, 24 = aged).",
    "Grey bars = Control, Pink bars = Piezo1 KO.",
    "Significance: * p<0.05, ** p<0.01, *** p<0.001 (Bonferroni-corrected)."
  )

  save_figure_outputs(
    figure_id = "Fig3H",
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 7,
    height = 6,
    code_path = "pipeline/figure_modules/Fig3H_contractile_capacity_score.R"
  )

  # Save legend as separate PNG (sized for manuscript: legend ~1/3-1/6 of plot)
  p_for_legend <- p +
    scale_fill_manual(values = colors$genotype, name = "Genotype") +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm"))
  legend_grob <- cowplot::get_legend(p_for_legend)
  ggsave(file.path(fig_dir, "Fig3H_legend.png"), legend_grob,
         width = 6, height = 3, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, "Fig3H_legend.pdf"), legend_grob,
         width = 6, height = 3, bg = "white")
  cat("  Saved: Fig3H_legend.png/pdf\n")

  cat("  Done!\n")

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    gene_panel = avail,
    summary = summary_stats,
    kruskal_wallis = kw_test,
    wilcoxon_tests = pairwise_tests
  ))
}
