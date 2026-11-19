# ==============================================================================
# PAPER PANEL: Fig 3A   [internal module ID: Fig3E]
# PLOTS:  UMAP of visceral SMCs coloured by phenotype (contractile vs synthetic), one panel per group, zoomed to the SMC region
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Young CTRL (3 mo), Aged CTRL (24 mo), Young P1KO (3 mo), Aged P1KO (24 mo) (shown side by side; no statistical test)
# TEST:   none
# ==============================================================================
# Four panels showing SMC phenotype (contractile/synthetic) for each group
# Zoomed into SMC coordinates with outlier trimming
# Separate legend PNG exported alongside panels
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate Fig3E
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "core/output_manager.R"))

  # Define colors
  color_contractile <- "#E41A1C"  # Red
  color_synthetic <- "#377EB8"    # Blue
  color_other <- "#D9D9D9"       # Light grey for non-SMC

  group_order <- c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  phenotype <- context$phenotype
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Get UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, "umap")
  umap_df <- data.frame(
    cell_barcode = rownames(umap_coords),
    UMAP_1 = umap_coords[, 1],
    UMAP_2 = umap_coords[, 2],
    stringsAsFactors = FALSE
  )

  # Merge with metadata
  plot_data <- merge(umap_df, metadata[, c("cell_barcode", "group")],
                     by = "cell_barcode", all.x = TRUE)

  # Assign phenotype only to cells in smc_cells
  phenotype_lookup <- setNames(phenotype$phenotype, phenotype$cell_barcode)
  plot_data$is_smc <- plot_data$cell_barcode %in% smc_cells
  plot_data$phenotype <- ifelse(
    plot_data$cell_barcode %in% names(phenotype_lookup),
    phenotype_lookup[plot_data$cell_barcode],
    "Other"
  )
  plot_data$phenotype[!plot_data$is_smc] <- "Other"

  # Verify counts
  cat("  Verifying SMC phenotype counts by group:\n")
  for (grp in group_order) {
    smc_in_grp <- sum(plot_data$is_smc & plot_data$group == grp)
    contr_in_grp <- sum(plot_data$is_smc & plot_data$group == grp &
                        plot_data$phenotype == "Contractile")
    synth_in_grp <- sum(plot_data$is_smc & plot_data$group == grp &
                        plot_data$phenotype == "Synthetic")
    cat(sprintf("    %s: %d SMCs (%d contractile, %d synthetic)\n",
                grp, smc_in_grp, contr_in_grp, synth_in_grp))
  }

  # ==========================================================================
  # COMPUTE ZOOM WINDOW: 1st-99th percentile of SMC UMAP coordinates + padding
  # ==========================================================================
  smc_umap <- plot_data %>% filter(is_smc)
  xlim <- quantile(smc_umap$UMAP_1, probs = c(0.01, 0.99))
  ylim <- quantile(smc_umap$UMAP_2, probs = c(0.01, 0.99))
  # Add 5% padding
  x_pad <- diff(xlim) * 0.05
  y_pad <- diff(ylim) * 0.05
  xlim <- c(xlim[1] - x_pad, xlim[2] + x_pad)
  ylim <- c(ylim[1] - y_pad, ylim[2] + y_pad)
  cat(sprintf("  Zoom window: UMAP1 [%.1f, %.1f], UMAP2 [%.1f, %.1f]\n",
              xlim[1], xlim[2], ylim[1], ylim[2]))

  # Group labels for titles — paper-consistent nomenclature
  group_titles <- c(
    "Young_Ctrl" = "CTRL 3",
    "Aged_Ctrl"  = "CTRL 24",
    "Young_KO"   = "P1KO 3",
    "Aged_KO"    = "P1KO 24"
  )

  # Create individual plots
  plots <- list()

  for (i in seq_along(group_order)) {
    grp <- group_order[i]

    grp_data <- plot_data %>%
      filter(group == grp) %>%
      mutate(
        highlight = case_when(
          is_smc & group == grp & phenotype == "Contractile" ~ "Contractile",
          is_smc & group == grp & phenotype == "Synthetic" ~ "Synthetic",
          TRUE ~ "Other"
        )
      )

    # Order: Other first, then Contractile, then Synthetic ON TOP
    grp_data$highlight <- factor(grp_data$highlight,
                                  levels = c("Other", "Contractile", "Synthetic"))
    grp_data <- grp_data[order(grp_data$highlight), ]

    p <- ggplot() +
      # Layer 1: Background (grey) cells - small, low alpha
      geom_point(data = filter(grp_data, highlight == "Other"),
                 aes(x = UMAP_1, y = UMAP_2),
                 color = color_other, size = 0.3, alpha = 0.4) +
      # Layer 2: Contractile (red) - bigger, semi-transparent
      geom_point(data = filter(grp_data, highlight == "Contractile"),
                 aes(x = UMAP_1, y = UMAP_2),
                 color = color_contractile, size = 1.5, alpha = 0.6) +
      # Layer 3: Synthetic (blue) ON TOP - bigger, semi-transparent
      geom_point(data = filter(grp_data, highlight == "Synthetic"),
                 aes(x = UMAP_1, y = UMAP_2),
                 color = color_synthetic, size = 1.5, alpha = 0.7) +
      # Zoom to SMC region
      coord_cartesian(xlim = xlim, ylim = ylim) +
      labs(
        title = group_titles[grp],
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      theme_classic(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        legend.position = "none",
        aspect.ratio = 1
      )

    plots[[as.character(i)]] <- p
  }

  # Create output directory
  fig_dir <- file.path(context$output_dir, "Fig3E")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # ==========================================================================
  # SEPARATE LEGEND PNG
  # ==========================================================================
  legend_plot <- ggplot(data.frame(
    x = c(1, 2), y = c(1, 2),
    Phenotype = factor(c("Contractile", "Synthetic"), levels = c("Contractile", "Synthetic"))
  ), aes(x = x, y = y, color = Phenotype)) +
    geom_point(size = 6) +
    scale_color_manual(
      values = c("Contractile" = color_contractile, "Synthetic" = color_synthetic),
      name = "SMC Phenotype"
    ) +
    theme_classic(base_size = 28) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 32),
      legend.text = element_text(size = 28),
      legend.key.size = unit(2, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 10, alpha = 1)))

  # Build the grob and extract the legend
  legend_grob <- ggplot_gtable(ggplot_build(legend_plot))
  legend_idx <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend_only <- legend_grob$grobs[[legend_idx]]
    png(file.path(fig_dir, "Fig3E_legend.png"), width = 1800, height = 900, res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, "Fig3E_legend.pdf"), width = 6, height = 3, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend: Fig3E_legend.png/pdf\n")
  } else {
    cat("  WARNING: Could not extract legend grob\n")
  }

  # Summary statistics
  phenotype_summary <- plot_data %>%
    filter(is_smc) %>%
    group_by(group, phenotype) %>%
    summarize(n = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()

  # Prepare data for export
  export_data <- list(
    phenotype_summary = phenotype_summary,
    phenotype_by_cell = phenotype,
    smc_cells_total = length(smc_cells),
    zoom_window = data.frame(
      axis = c("UMAP_1", "UMAP_2"),
      min = c(xlim[1], ylim[1]),
      max = c(xlim[2], ylim[2])
    ),
    smc_by_group = plot_data %>%
      filter(is_smc) %>%
      group_by(group) %>%
      summarize(n_smc = n(), .groups = "drop"),
    classification_parameters = data.frame(
      parameter = c("ratio_threshold", "contractile_genes", "synthetic_genes"),
      value = c("1.5", "Myh11,Acta2,Tagln,Cnn1,Mylk,Myocd,Smtn",
                "Col1a1,Fn1,Spp1,Klf4,Vim")
    )
  )

  # Methods text
  methods_text <- paste(
    "SMC phenotypic state was classified at the single-cell level based on the",
    "relative expression of established contractile and synthetic marker genes.",
    "ONLY cells from the SMC selection (cluster-based pure SMCs) were classified.",
    "For each SMC, phenotype scores were calculated as the mean normalized",
    "expression of contractile markers (Myh11, Acta2, Tagln, Cnn1, Mylk, Myocd, Smtn)",
    "and synthetic markers (Col1a1, Fn1, Spp1, Klf4, Vim).",
    "Cells with a contractile-to-synthetic score ratio greater than 1.5 were",
    "classified as having a contractile phenotype; all other SMCs were classified",
    "as synthetic. UMAP panels are zoomed to the SMC region (1st-99th percentile",
    "of SMC UMAP coordinates). Non-SMC cells within the zoom window shown in grey."
  )

  technical_legend <- paste(
    "UMAP visualization showing SMC phenotype classification for each experimental group.",
    "Fig3E1: Young Control, Fig3E2: Aged Control, Fig3E3: Young KO, Fig3E4: Aged KO.",
    "Red = Contractile SMCs, Blue = Synthetic SMCs, Grey = Non-SMC cells.",
    "View is zoomed to SMC UMAP coordinates (outliers trimmed at 1st/99th percentile).",
    sprintf("Total SMCs analyzed: %d", length(smc_cells))
  )

  # Save outputs
  save_multi_panel_figure(
    figure_id = "Fig3E",
    plots = plots,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 6,
    height = 5,
    code_path = "pipeline/figure_modules/Fig3E_umap_phenotype.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    phenotype_summary = phenotype_summary
  ))
}
