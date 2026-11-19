# ==============================================================================
# PAPER PANEL: Fig 3H (Table S6)   [internal module ID: Fig3J]
# PLOTS:  log2 fold change of individual synthetic genes: grey = Aged vs Young CTRL (aging), pink = Aged P1KO vs Aged CTRL (rescue)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Aged CTRL vs Young CTRL; Aged P1KO vs Aged CTRL
# TEST:   Wilcoxon per gene, Bonferroni
# ==============================================================================
# Shows aging increase (positive) mirrored with KO rescue (negative)
# Style: theme_minimal, stacked bars, panel border
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(Seurat)
})

#' Generate Fig3J
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "core/output_manager.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))
  source(file.path(pipeline_root, "utils/statistics.R"))

  gene_sets <- get_gene_sets()

  # Colors
  color_ctrl <- "#7F7F7F"  # Grey for Ctrl/aging increase
  color_ko <- "#ff7fb2"    # Pink for KO rescue

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Get expression data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Get cells by group
  cells_young_ctrl <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Young_Ctrl"]]
  cells_aged_ctrl <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_Ctrl"]]
  cells_aged_ko <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_KO"]]

  # Combine synthetic and ECM genes for screening
  candidate_genes <- unique(c(
    gene_sets$synthetic_extended,
    gene_sets$ecm_glycoproteins,
    gene_sets$ecm_collagens
  ))
  available_genes <- candidate_genes[candidate_genes %in% rownames(expr)]

  # Calculate mean expression for each gene in each group
  results <- data.frame(
    gene = available_genes,
    young_ctrl = rowMeans(as.matrix(expr[available_genes, cells_young_ctrl, drop = FALSE])),
    aged_ctrl = rowMeans(as.matrix(expr[available_genes, cells_aged_ctrl, drop = FALSE])),
    aged_ko = rowMeans(as.matrix(expr[available_genes, cells_aged_ko, drop = FALSE]))
  )

  # Calculate log2FC
  pseudo <- 0.01
  results$log2fc_aging <- log2((results$aged_ctrl + pseudo) / (results$young_ctrl + pseudo))
  results$log2fc_rescue <- log2((results$aged_ko + pseudo) / (results$aged_ctrl + pseudo))

  # Select genes that increase with aging AND are rescued (decreased) by KO
  # Criteria: positive log2fc_aging AND negative log2fc_rescue
  rescue_genes <- results %>%
    filter(log2fc_aging > 0.1 & log2fc_rescue < -0.05) %>%
    mutate(rescue_magnitude = log2fc_aging - log2fc_rescue) %>%
    arrange(desc(rescue_magnitude)) %>%
    head(8)

  # Order genes by aging increase magnitude (highest first, left to right means decreasing)
  gene_order <- rescue_genes %>%
    arrange(desc(log2fc_aging)) %>%
    pull(gene)

  # --------------------------------------------------------------------------
  # STATISTICAL TESTS and BOOTSTRAP CI - Wilcoxon for each gene
  # --------------------------------------------------------------------------
  stat_results <- data.frame()

  # Bootstrap function for log2FC confidence intervals
  bootstrap_log2fc <- function(expr_g1, expr_g2, n_boot = 1000, pseudo = 0.01) {
    boot_fc <- replicate(n_boot, {
      s1 <- sample(expr_g1, replace = TRUE)
      s2 <- sample(expr_g2, replace = TRUE)
      log2((mean(s2) + pseudo) / (mean(s1) + pseudo))
    })
    c(se = sd(boot_fc), ci_low = quantile(boot_fc, 0.025), ci_high = quantile(boot_fc, 0.975))
  }

  for (gene in rescue_genes$gene) {
    # Get expression values for this gene
    gene_expr <- as.numeric(expr[gene, ])
    names(gene_expr) <- colnames(expr)

    # Create data frame for this gene
    gene_data <- data.frame(
      cell = smc_cells,
      expr = gene_expr[smc_cells],
      group = metadata[smc_cells, "group"]
    )

    # Pairwise Mann-Whitney — all 4 valid comparisons (Bonferroni-adjusted),
    # so the export carries the full pairwise table per gene. The mirror bars
    # display only aging-WT and KO-rescue-aged.
    pw <- pairwise_group_tests(gene_data$expr, gene_data$group)
    test_aging  <- pw[pw$comparison == "Aging (WT)",       ]
    test_rescue <- pw[pw$comparison == "KO rescue (aged)", ]

    # Bootstrap CI for aging and rescue log2FC
    expr_yc <- gene_expr[cells_young_ctrl]
    expr_ac <- gene_expr[cells_aged_ctrl]
    expr_ak <- gene_expr[cells_aged_ko]

    boot_aging <- bootstrap_log2fc(expr_yc, expr_ac)
    boot_rescue <- bootstrap_log2fc(expr_ac, expr_ak)

    stat_results <- rbind(stat_results, data.frame(
      gene = gene,
      aging_pvalue = test_aging$p_value,
      aging_signif = test_aging$significance,
      aging_se = boot_aging["se"],
      rescue_pvalue = test_rescue$p_value,
      rescue_signif = test_rescue$significance,
      rescue_se = boot_rescue["se"]
    ))

    pw$gene <- gene
    if (!exists("pairwise_full_table")) pairwise_full_table <- list()
    pairwise_full_table[[gene]] <- pw
  }
  pairwise_full_table <- do.call(rbind, pairwise_full_table)

  # Merge stats with rescue_genes
  rescue_genes <- merge(rescue_genes, stat_results, by = "gene")
  rescue_genes <- rescue_genes[match(gene_order, rescue_genes$gene), ]

  # Prepare data for stacked mirror plot
  plot_data <- rescue_genes %>%
    select(gene, log2fc_aging, log2fc_rescue, aging_signif, rescue_signif, aging_se, rescue_se) %>%
    pivot_longer(cols = c(log2fc_aging, log2fc_rescue),
                 names_to = "comparison",
                 values_to = "log2fc") %>%
    mutate(
      change_label = factor(comparison,
                            levels = c("log2fc_aging", "log2fc_rescue"),
                            labels = c("Age Increase\n(Aged Ctrl - Young Ctrl)",
                                      "KO Rescue\n(Aged KO - Aged Ctrl)")),
      gene = factor(gene, levels = gene_order),
      signif = ifelse(comparison == "log2fc_aging", aging_signif, rescue_signif),
      se = ifelse(comparison == "log2fc_aging", aging_se, rescue_se)
    )

  # Prepare asterisk annotation data - SEPARATE for aging and rescue
  # Aging asterisks: above the grey (positive) bars
  aging_asterisk_data <- rescue_genes %>%
    filter(aging_signif != "ns") %>%
    mutate(
      y_pos = log2fc_aging + aging_se + 0.08,  # Above the positive bar
      label = aging_signif
    ) %>%
    select(gene, y_pos, label)

  # Rescue asterisks: below the pink (negative) bars
  rescue_asterisk_data <- rescue_genes %>%
    filter(rescue_signif != "ns") %>%
    mutate(
      y_pos = log2fc_rescue - rescue_se - 0.08,  # Below the negative bar
      label = rescue_signif
    ) %>%
    select(gene, y_pos, label)

  # Create error bar data (separate for aging and rescue)
  errorbar_data <- plot_data %>%
    mutate(
      ymin = log2fc - se,
      ymax = log2fc + se
    )

  # Create stacked mirror bar plot
  # Bars are stacked using position = "identity" so aging (positive) goes up
  # and rescue (negative) goes down from the same baseline
  p <- ggplot(plot_data, aes(x = gene, y = log2fc, fill = change_label)) +
    # Stacked bars - identity position stacks them at same x position
    geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3,
             position = "identity") +
    # Error bars at the end of each bar segment
    geom_errorbar(data = errorbar_data,
                  aes(ymin = log2fc - se, ymax = log2fc + se),
                  width = 0.2, linewidth = 0.5, color = "black",
                  position = "identity") +
    # Zero line
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    # Add significance asterisks for AGING effect (above grey bars)
    geom_text(data = aging_asterisk_data, aes(x = gene, y = y_pos, label = label),
              inherit.aes = FALSE, size = 7, fontface = "bold", color = "black") +
    # Add significance asterisks for RESCUE effect (below pink bars)
    geom_text(data = rescue_asterisk_data, aes(x = gene, y = y_pos, label = label),
              inherit.aes = FALSE, size = 7, fontface = "bold", color = "black") +
    # Colors
    scale_fill_manual(values = c(
      "Age Increase\n(Aged Ctrl - Young Ctrl)" = color_ctrl,
      "KO Rescue\n(Aged KO - Aged Ctrl)" = color_ko
    )) +
    # Labels
    labs(
      title = "Synthetic Genes: Age Increase and Piezo1 KO Rescue",
      x = NULL,
      y = expression(paste("Log"[2], " Fold Change")),
      fill = NULL
    ) +
    # Theme (publication sizing)
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.text.x = element_text(face = "bold.italic", size = 22, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 22),
      axis.title.y = element_text(face = "bold", size = 20),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, "Fig3J")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Calculate rescue efficiency for export
  gene_summary <- rescue_genes %>%
    mutate(
      net_effect = log2fc_aging + log2fc_rescue,
      rescue_pct = (abs(log2fc_rescue) / log2fc_aging) * 100
    ) %>%
    arrange(desc(log2fc_aging))

  # Prepare data for export (multiple sheets)
  export_data <- list(
    summary = gene_summary,
    statistical_tests = stat_results,
    pairwise_full = pairwise_full_table,  # all 4 valid comparisons per gene, Bonferroni-adjusted
    all_screened_genes = results,
    plot_data = as.data.frame(plot_data),
    selection_criteria = data.frame(
      criterion = c("log2fc_aging_threshold", "log2fc_rescue_threshold", "n_genes_selected"),
      value = c("0.1", "-0.05", nrow(rescue_genes))
    ),
    methods = data.frame(
      parameter = c("Gene set source", "Selection criteria aging", "Selection criteria rescue",
                   "Pseudocount for log2FC", "Statistical test"),
      value = c("synthetic_extended + ECM genes", "log2FC > 0.1 (Aged Ctrl vs Young Ctrl)",
               "log2FC < -0.05 (Aged KO vs Aged Ctrl)", "0.01", "Wilcoxon rank-sum test")
    )
  )

  # Methods text
  methods_text <- paste(
    "Synthetic and ECM genes were screened for aging-associated increase and",
    "Piezo1 KO-mediated rescue. Genes were selected if they showed positive log2",
    "fold change in Aged Ctrl vs Young Ctrl (indicating aging increase) AND negative",
    "log2 fold change in Aged KO vs Aged Ctrl (indicating rescue/reduction).",
    "The top genes with the largest combined rescue magnitude are shown.",
    "Grey bars (positive) show aging increase; pink bars (negative) show KO rescue.",
    "Statistical significance was assessed using Wilcoxon rank-sum tests comparing",
    "expression levels between groups."
  )

  technical_legend <- paste(
    "Stacked mirror bar plot showing log2 fold changes for selected synthetic/ECM genes.",
    "Grey bars (pointing up) = Aged Ctrl vs Young Ctrl (aging increase).",
    "Pink bars (pointing down) = Aged KO vs Aged Ctrl (KO rescue).",
    "Genes ordered by aging increase magnitude. Asterisks above grey bars indicate significance",
    "of aging effect (Wilcoxon: Aged_Ctrl vs Young_Ctrl). Asterisks below pink bars indicate",
    "significance of rescue effect (Wilcoxon: Aged_KO vs Aged_Ctrl). * p<0.05, ** p<0.01, *** p<0.001."
  )

  # Save separate legend
  legend_plot <- ggplot(data.frame(
    x = 1:2, y = 1:2,
    Change = factor(c("Age Increase\n(Aged Ctrl - Young Ctrl)",
                       "KO Rescue\n(Aged KO - Aged Ctrl)"),
                     levels = c("Age Increase\n(Aged Ctrl - Young Ctrl)",
                                "KO Rescue\n(Aged KO - Aged Ctrl)"))
  ), aes(x = x, y = y, fill = Change)) +
    geom_bar(stat = "identity", width = 0.5, color = "black") +
    scale_fill_manual(values = c(
      "Age Increase\n(Aged Ctrl - Young Ctrl)" = color_ctrl,
      "KO Rescue\n(Aged KO - Aged Ctrl)" = color_ko
    )) +
    theme_classic(base_size = 28) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm"))
  legend_grob <- ggplot_gtable(ggplot_build(legend_plot))
  legend_idx <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend_only <- legend_grob$grobs[[legend_idx]]
    png(file.path(fig_dir, "Fig3J_legend.png"), width = 2400, height = 1200, res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, "Fig3J_legend.pdf"), width = 8, height = 4, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend: Fig3J_legend.png/pdf\n")
  }

  # Save outputs
  save_figure_outputs(
    figure_id = "Fig3J",
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 8,
    height = 7,
    code_path = "pipeline/figure_modules/Fig3J_synthetic_rescue.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    selected_genes = rescue_genes
  ))
}
