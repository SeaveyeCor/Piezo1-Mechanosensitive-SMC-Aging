# ==============================================================================
# PAPER PANEL: Fig S3J (internal module ID FigS4A)   [internal module ID: FigS3L]
# PLOTS:  heatmap of log2 fold change for ECM genes (collagens, fibronectins, laminins, MMPs, TIMPs, other ECM)
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Young P1KO vs Young CTRL; Aged P1KO vs Aged CTRL; Aged CTRL vs Young CTRL
# TEST:   none (fold changes only)
# ==============================================================================
# Heatmap showing Log2FC for ECM genes across comparisons
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

#' Generate FigS3L
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

  # Get cells by group
  cells_young_ctrl <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Young_Ctrl"]]
  cells_young_ko <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Young_KO"]]
  cells_aged_ctrl <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_Ctrl"]]
  cells_aged_ko <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_KO"]]

  # Define ECM gene categories
  ecm_categories <- list(
    "Collagens" = c("Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col5a1", "Col6a1"),
    "Fibronectins" = c("Fn1", "Fbln1", "Fbln2", "Fbln5"),
    "Laminins" = c("Lama2", "Lama4", "Lamb1", "Lamc1"),
    "MMPs" = c("Mmp2", "Mmp3", "Mmp9", "Mmp14"),
    "TIMPs" = c("Timp1", "Timp2", "Timp3"),
    "Other ECM" = c("Sparc", "Bgn", "Dcn", "Vcan", "Tnc", "Postn")
  )

  # Calculate mean expression and log2FC for each gene
  results <- list()

  for (category in names(ecm_categories)) {
    genes <- ecm_categories[[category]]

    for (gene in genes) {
      if (!gene %in% rownames(expr)) next

      mean_young_ctrl <- mean(as.numeric(expr[gene, cells_young_ctrl]))
      mean_young_ko <- mean(as.numeric(expr[gene, cells_young_ko]))
      mean_aged_ctrl <- mean(as.numeric(expr[gene, cells_aged_ctrl]))
      mean_aged_ko <- mean(as.numeric(expr[gene, cells_aged_ko]))

      pseudo <- 0.01

      results[[gene]] <- data.frame(
        gene = gene,
        category = category,
        young_ko_vs_ctrl = log2((mean_young_ko + pseudo) / (mean_young_ctrl + pseudo)),
        aged_ko_vs_ctrl = log2((mean_aged_ko + pseudo) / (mean_aged_ctrl + pseudo)),
        aged_vs_young_ctrl = log2((mean_aged_ctrl + pseudo) / (mean_young_ctrl + pseudo)),
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine results
  fc_df <- bind_rows(results)

  # Reshape for heatmap
  plot_data <- fc_df %>%
    pivot_longer(
      cols = c(young_ko_vs_ctrl, aged_ko_vs_ctrl, aged_vs_young_ctrl),
      names_to = "comparison",
      values_to = "log2fc"
    ) %>%
    mutate(
      comparison = factor(comparison,
                          levels = c("young_ko_vs_ctrl", "aged_ko_vs_ctrl", "aged_vs_young_ctrl"),
                          labels = c("P1KO 3\nvs\nCTRL 3", "P1KO 24\nvs\nCTRL 24", "CTRL 24\nvs\nCTRL 3")),  # comparison labels: 3 lines, horizontal
      category = factor(category, levels = names(ecm_categories))
    )

  # Order genes within categories
  gene_order <- fc_df %>%
    arrange(category, desc(aged_ko_vs_ctrl)) %>%
    pull(gene)

  plot_data$gene <- factor(plot_data$gene, levels = rev(gene_order))

  # Create heatmap
  p <- ggplot(plot_data, aes(x = comparison, y = gene, fill = log2fc)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = colors$heatmap["low"],
      mid = colors$heatmap["mid"],
      high = colors$heatmap["high"],
      midpoint = 0,
      limits = c(-2, 2),
      oob = scales::squish,
      name = expression(Log[2]~FC)
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "ECM Pathway: Comprehensive Comparison",
      x = "",
      y = ""
    ) +
    theme_heatmap() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 16, face = "bold", lineheight = 0.9),
      axis.text.y = element_text(face = "bold.italic", size = 16),
      strip.text.y = element_text(angle = 0, face = "bold", size = 18),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "none"   # legend exported separately so it can be resized independently
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Separate legend PNG/PDF — the heatmap fill scale extracted as a standalone
  # asset so the legend can be sized/positioned independently of the main plot
  # in the manuscript figure assembly.
  p_for_legend <- p +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm"),
          legend.key.height = unit(2.5, "cm"))
  legend_grob <- cowplot::get_legend(p_for_legend)
  if (!is.null(legend_grob)) {
    ggsave(file.path(fig_dir, paste0(context$fig_id, "_legend.png")),
           legend_grob, width = 4, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(fig_dir, paste0(context$fig_id, "_legend.pdf")),
           legend_grob, width = 4, height = 6, bg = "white")
    cat("  Saved separate legend\n")
  }

  # Prepare data for export
  export_data <- list(
    log2fc_matrix = fc_df,
    plot_data = plot_data,
    gene_categories = data.frame(
      gene = unlist(ecm_categories),
      category = rep(names(ecm_categories), sapply(ecm_categories, length))
    )
  )

  # Methods text
  methods_text <- paste(
    "Log2 fold changes were calculated for ECM genes across three comparisons:",
    "(1) Young KO vs Ctrl, (2) Aged KO vs Ctrl, and (3) Aged vs Young Ctrl.",
    "Genes are grouped by ECM category: Collagens, Fibronectins, Laminins,",
    "MMPs, TIMPs, and Other ECM components."
  )

  technical_legend <- paste(
    "Heatmap showing log2 fold change for ECM genes.",
    "Blue = negative FC (decreased), Red = positive FC (increased).",
    "Genes grouped by category (y-axis facets).",
    "Comparisons shown on x-axis."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 8,
    height = 10,
    code_path = "pipeline/figure_modules/supplementary/FigS3L_ecm_heatmap.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    log2fc = fc_df
  ))
}
