# ==============================================================================
# PAPER PANEL: Fig S3A   [internal module ID: FigS3D]
# PLOTS:  UMAP of all muscularis cells coloured by annotated cell population (Cl6 and Cl13 labelled separately)
# CELLS:  all muscularis cells after removal of non-muscularis populations (n = 11,738), all 6 mice
# GROUPS COMPARED: none (descriptive)
# TEST:   none
# ==============================================================================
# UMAP showing all cells colored by cell type
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate FigS3D
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))

  colors <- get_colors()

  # Get data
  seurat_obj <- context$seurat_obj
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
  plot_data <- merge(umap_df, metadata, by.x = "cell_barcode", by.y = "row.names")

  # Cell type labels: one labeled category per cell_type, matching Table S1 nomenclature
  cell_type_labels <- c(
    "SMC_Contractile"        = "Visceral SMC",
    "Fibroblast"             = "PDGFRα+ FLC",
    "Glial"                  = "Enteric Glia",
    "ICC"                    = "ICC",
    "Immune_Macrophage"      = "Macrophage",
    "Immune_TCell"           = "T Cell",
    "Immune_BCell"           = "B Cell",
    "Endothelial"            = "Vascular Endothelial",
    "Endothelial_Lymphatic"  = "Lymphatic Endothelial",
    "Pericyte"               = "Arterial SMC / Pericyte",
    "Neuronal"               = "Enteric Neuron",
    "Mesothelial"            = "Mesothelial",
    "Erythrocyte"            = "Erythrocyte",
    "Epithelial"             = "Epithelial",
    "Goblet"                 = "Goblet"
  )

  plot_data$cell_type_label <- cell_type_labels[plot_data$cell_type]
  # SMC-annotated clusters that are not visceral SMCs (Table S1): Cl6 and Cl13 get their own labels
  cl <- as.character(plot_data$seurat_clusters)
  plot_data$cell_type_label[plot_data$cell_type == "SMC_Contractile" & cl %in% "6"]  <- "Lymphatic vessel-associated (Cl6)"
  plot_data$cell_type_label[plot_data$cell_type == "SMC_Contractile" & cl %in% "13"] <- "Mixed glia/SMC (Cl13)"
  # Drop unmapped (should be zero post-filter)
  n_unmapped <- sum(is.na(plot_data$cell_type_label))
  if (n_unmapped > 0) {
    cat(sprintf("  WARNING: dropping %d cells with unmapped cell_type\n", n_unmapped))
    plot_data <- plot_data[!is.na(plot_data$cell_type_label), ]
  }

  # Distinct color palette for 15 cell types (S1-aligned)
  type_colors <- c(
    "Visceral SMC"            = "#E41A1C",
    "Lymphatic vessel-associated (Cl6)" = "#FB9A99",
    "Mixed glia/SMC (Cl13)"   = "#FFD92F",
    "PDGFRα+ FLC"             = "#377EB8",
    "Enteric Glia"            = "#FF7F00",
    "ICC"                     = "#4DAF4A",
    "Macrophage"              = "#F781BF",
    "T Cell"                  = "#984EA3",
    "B Cell"                  = "#A65628",
    "Vascular Endothelial"    = "#FFFF33",
    "Lymphatic Endothelial"   = "#00CED1",
    "Arterial SMC / Pericyte" = "#FF1493",
    "Enteric Neuron"          = "#1B9E77",
    "Mesothelial"             = "#66C2A5",
    "Erythrocyte"             = "#7F7F7F",
    "Epithelial"              = "#9999FF",
    "Goblet"                  = "#D95F02"
  )

  # Create plot
  p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = cell_type_label)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = type_colors, name = "Cell Type") +
    labs(
      title = "UMAP of All Cell Types",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_umap() +
    theme(legend.position = "none")

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Cell type summary
  cell_summary <- plot_data %>%
    group_by(cell_type_label) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(desc(n))

  # Save separate cell-type legend — only show cell types ACTUALLY present in the UMAP
  present <- as.character(unique(plot_data$cell_type_label))
  present <- intersect(names(type_colors), present)  # preserve color order
  legend_plot <- ggplot(data.frame(
    x = seq_along(present),
    y = seq_along(present),
    `Cell Type` = factor(present, levels = present),
    check.names = FALSE
  ), aes(x = x, y = y, color = `Cell Type`)) +
    geom_point(size = 5) +
    scale_color_manual(values = type_colors, name = "Cell Type") +
    theme_classic(base_size = 28) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 10, alpha = 1)))
  legend_grob <- ggplot_gtable(ggplot_build(legend_plot))
  legend_idx <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend_only <- legend_grob$grobs[[legend_idx]]
    # size the image to the legend itself so no label is cut off
    lw <- grid::convertWidth(sum(legend_only$widths), "in", valueOnly = TRUE) + 0.4
    lh <- grid::convertHeight(sum(legend_only$heights), "in", valueOnly = TRUE) + 0.4
    png(file.path(fig_dir, paste0(context$fig_id, "_legend.png")),
        width = lw, height = lh, units = "in", res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, paste0(context$fig_id, "_legend.pdf")),
        width = lw, height = lh, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend\n")
  }

  # Prepare data for export
  export_data <- list(
    umap_coordinates = umap_df,
    cell_type_counts = cell_summary
  )

  # Methods text
  methods_text <- paste(
    "UMAP visualization of the integrated scRNA-seq dataset showing all cell types.",
    "Cell types shown include: SMC (smooth muscle cells), PDGFRa+ (fibroblasts),",
    "ICC (interstitial cells of Cajal), Glia (enteric glia), Epithelial,",
    "Endothelial, and Macrophage."
  )

  technical_legend <- paste(
    "UMAP plot with cells colored by cell type annotation.",
    "Each dot represents a single cell."
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
    height = 7,
    code_path = "pipeline/figure_modules/supplementary/FigS3D_umap_celltypes.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    cell_summary = cell_summary
  ))
}
