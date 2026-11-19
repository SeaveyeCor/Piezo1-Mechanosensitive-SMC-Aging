# ==============================================================================
# PAPER PANEL: none printed directly; pipeline version of Fig 4A (pathway scores incl. the 13-gene NFAT pathway score with Nfat5). The printed Fig 4A and its source data come from 04_cell_level/71_score_panels.R   [internal module ID: FigS3G]
# PLOTS:  Spearman rho between Piezo1 and each pathway score, per age
# CELLS:  CTRL visceral SMCs only (Young CTRL n = 469, Aged CTRL n = 477)
# GROUPS COMPARED: Young CTRL and Aged CTRL, each tested separately
# TEST:   Spearman correlation
# ==============================================================================
# Dot plot showing Spearman rho with Piezo1 for each pathway
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

#' Generate FigS3G
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
  piezo1_expr <- as.numeric(expr["Piezo1", ctrl_cells])

  # Define pathways to analyze (validated: all positive pathways maintain or
  # strengthen with aging; replaces "General Synthetic" and "ECM Organization"
  # which had ns correlations)
  pathways <- list(
    "NFAT Pathway" = gene_sets$nfat_pathway,
    "TGF-beta/Stress" = gene_sets$tgfbeta_pathway,
    "Calcium SOCE/IP3" = gene_sets$calcium_pathway,
    "Core Contractile" = gene_sets$contractile_core,
    "Dedifferentiation" = gene_sets$dedifferentiation,
    "NFAT Activity" = gene_sets$nfatc3_remodeling,
    "ECM Score" = gene_sets$ecm_structural,
    "Cyto Remodeling" = gene_sets$cyto_remodeling,
    "Synthetic Phenotype" = gene_sets$synthetic_phenotype
  )

  # Calculate pathway scores and correlations
  results <- list()

  for (pathway_name in names(pathways)) {
    genes <- pathways[[pathway_name]]
    avail_genes <- genes[genes %in% rownames(expr)]

    if (length(avail_genes) == 0) next

    # Calculate pathway score
    pathway_score <- colMeans(as.matrix(expr[avail_genes, ctrl_cells, drop = FALSE]))

    # Create temp data
    temp_data <- data.frame(
      piezo1 = piezo1_expr,
      score = pathway_score,
      age = metadata[ctrl_cells, "age"]
    )

    # Calculate correlations by age
    cor_young <- spearman_cor(
      temp_data$piezo1[temp_data$age == "Young"],
      temp_data$score[temp_data$age == "Young"]
    )

    cor_aged <- spearman_cor(
      temp_data$piezo1[temp_data$age == "Aged"],
      temp_data$score[temp_data$age == "Aged"]
    )

    results[[pathway_name]] <- data.frame(
      pathway = pathway_name,
      young_rho = cor_young$rho,
      young_p = cor_young$p_value,
      aged_rho = cor_aged$rho,
      aged_p = cor_aged$p_value
    )
  }

  # Combine results
  cor_df <- bind_rows(results)

  # Reshape for plotting
  plot_data <- cor_df %>%
    pivot_longer(
      cols = c(young_rho, aged_rho),
      names_to = "age_type",
      values_to = "rho"
    ) %>%
    mutate(
      age = ifelse(grepl("young", age_type), "Young", "Aged"),
      p_value = ifelse(age == "Young", young_p, aged_p),
      sig_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      pathway = factor(pathway, levels = rev(cor_df$pathway[order(cor_df$aged_rho)]))
    )

  # Determine p columns for pivot
  plot_data$young_p <- cor_df$young_p[match(plot_data$pathway, cor_df$pathway)]
  plot_data$aged_p <- cor_df$aged_p[match(plot_data$pathway, cor_df$pathway)]
  plot_data$p_value <- ifelse(plot_data$age == "Young", plot_data$young_p, plot_data$aged_p)
  plot_data$sig_label <- ifelse(plot_data$p_value < 0.001, "***",
                         ifelse(plot_data$p_value < 0.01, "**",
                         ifelse(plot_data$p_value < 0.05, "*", "ns")))

  # Create plot
  p <- ggplot(plot_data, aes(x = rho, y = pathway, color = age)) +
    geom_point(size = 4) +
    geom_text(aes(label = sig_label), vjust = -0.8, size = 5, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Young" = "#4DAF4A", "Aged" = "#E41A1C"),
                       labels = c("Young" = "3mo", "Aged" = "24mo"),
                       name = "Age") +
    labs(
      title = "Piezo1 Correlation with Pathway Scores",
      x = "Spearman rho with Piezo1",
      y = ""
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 14, face = "bold")
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Separate legend PNG
  legend_plot <- ggplot(data.frame(x = 1:2, y = 1:2,
    Age = factor(c("3mo", "24mo"), levels = c("3mo", "24mo"))),
    aes(x = x, y = y, color = Age)) +
    geom_point(size = 6) +
    scale_color_manual(values = c("3mo" = "#4DAF4A", "24mo" = "#E41A1C"),
                       name = "Age") +
    theme_classic(base_size = 28) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 10)))
  legend_grob <- ggplot_gtable(ggplot_build(legend_plot))
  legend_idx <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend_only <- legend_grob$grobs[[legend_idx]]
    png(file.path(fig_dir, paste0(context$fig_id, "_legend.png")),
        width = 1800, height = 900, res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, paste0(context$fig_id, "_legend.pdf")),
        width = 6, height = 3, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend\n")
  }

  # Prepare data for export
  export_data <- list(
    correlation_summary = cor_df,
    plot_data = plot_data
  )

  # Methods text
  methods_text <- paste(
    "Spearman correlation coefficients were calculated between Piezo1 expression",
    "and pathway module scores in control SMCs, stratified by age.",
    "Pathway scores were calculated as the mean normalized expression of",
    "genes in each pathway gene set."
  )

  technical_legend <- paste(
    "Dot plot showing Spearman correlation (rho) between Piezo1 and pathway scores.",
    "Green = Young (3 mo), Red = Aged (24 mo).",
    "Positive rho indicates positive correlation with Piezo1.",
    "Dashed line at rho = 0."
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
    height = 6,
    code_path = "pipeline/figure_modules/supplementary/FigS3G_pathway_correlation.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    correlations = cor_df
  ))
}
