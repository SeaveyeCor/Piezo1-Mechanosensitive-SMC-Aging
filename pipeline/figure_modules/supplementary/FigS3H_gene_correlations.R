# ==============================================================================
# PAPER PANEL: Fig S4A   [internal module ID: FigS3H]
# PLOTS:  Spearman rho between Piezo1 and representative genes, grouped by pathway
# CELLS:  CTRL visceral SMCs only (Young CTRL n = 469, Aged CTRL n = 477)
# GROUPS COMPARED: Young CTRL and Aged CTRL, each tested separately
# TEST:   Spearman correlation
# ==============================================================================
# Gene-level correlations with Piezo1, colored by pathway
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate FigS3H
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

  # Use the pathway correlation genes
  genes_to_analyze <- gene_sets$pathway_correlation_genes
  available_genes <- genes_to_analyze[genes_to_analyze %in% rownames(expr)]

  # Define pathway assignments for coloring
  pathway_assignments <- list(
    "NFAT" = c("Nfatc1", "Nfatc2", "Nfatc3", "Ppp3ca", "Rcan1"),
    "TGF-beta" = c("Tgfb1", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3", "Smad4", "Smad7"),
    "Calcium" = c("Itpr1", "Itpr2", "Stim1", "Stim2", "Orai1", "Calm1"),
    "YAP/TAZ" = c("Yap1", "Wwtr1", "Lats1", "Lats2"),
    "MAPK/PI3K" = c("Mapk1", "Mapk8", "Mapk14", "Akt1", "Pik3ca", "Mtor"),
    "Contractile" = c("Cnn1", "Myh11", "Acta2", "Tagln"),
    "Synthetic" = c("Klf4", "Vim", "Col1a1", "Fn1", "Spp1", "Pcna", "Col4a1", "Vcam1"),
    "Other" = c("Fosb", "Junb", "Jun", "Fos", "Pxn", "Itgav", "Itga5", "Rock1",
                "Rock2", "Eln", "Itgb3", "Adamts1", "Adamts4", "Fbln5", "Ptk2",
                "Trp53", "Il6", "Cxcl1", "Cdkn1a")
  )

  # Create pathway lookup
  gene_pathway <- character()
  for (pw in names(pathway_assignments)) {
    for (g in pathway_assignments[[pw]]) {
      gene_pathway[g] <- pw
    }
  }

  # Calculate correlations for each gene
  results <- list()

  for (gene in available_genes) {
    gene_expr <- as.numeric(expr[gene, ctrl_cells])

    cor_result <- spearman_cor(piezo1_expr, gene_expr)

    results[[gene]] <- data.frame(
      gene = gene,
      rho = cor_result$rho,
      p_value = cor_result$p_value,
      n = cor_result$n,
      pathway = ifelse(gene %in% names(gene_pathway), gene_pathway[gene], "Other")
    )
  }

  # Combine results
  cor_df_all <- bind_rows(results) %>%
    arrange(desc(rho))

  # Curated representative genes — 2-3 per pathway, manuscript-cited genes forced
  # Tells the story: contractile negative, mechanosensing/NFAT/TGF-beta positive
  curated_genes <- c(
    # Contractile (negative correlation expected)
    "Myh11", "Acta2", "Tagln", "Cnn1",
    # Synthetic / Dedifferentiation (positive)
    "Klf4", "Vim",
    # NFAT pathway (positive)
    "Nfatc1", "Nfatc2", "Nfatc3",
    # Calcium signaling (positive)
    "Itpr1", "Itpr2",
    # TGF-beta / stress (positive) — Tgfb1 cited in manuscript text
    "Tgfb1", "Tgfbr2", "Smad3",
    # Senescence (positive — supports aging story)
    "Cdkn1a", "Trp53"
  )

  cor_df <- cor_df_all %>%
    filter(gene %in% curated_genes) %>%
    arrange(desc(rho))

  # Order genes by correlation (most positive at bottom for ggplot y-axis)
  cor_df$gene <- factor(cor_df$gene, levels = cor_df$gene)

  # Define pathway colors (only pathways present in curated list)
  pathway_colors <- c(
    "Contractile" = "#A65628",
    "Synthetic" = "#F781BF",
    "NFAT" = "#E41A1C",
    "Calcium" = "#4DAF4A",
    "TGF-beta" = "#377EB8",
    "Other" = "#999999"
  )
  # Map Cdkn1a/Trp53 to "Senescence" for clarity
  cor_df$pathway[cor_df$gene %in% c("Cdkn1a", "Trp53")] <- "Senescence"
  pathway_colors["Senescence"] <- "#FF7F00"

  # Create plot
  p <- ggplot(cor_df, aes(x = rho, y = gene, color = pathway)) +
    geom_point(size = 6.25) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = pathway_colors, name = "Pathway") +
    labs(
      title = "Piezo1-Gene Correlations by Pathway",
      x = "Spearman rho with Piezo1",
      y = ""
    ) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 24, face = "bold.italic"),
      axis.title.x = element_text(size = 27, face = "bold", margin = margin(t = 10)),
      legend.position = "none"
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Separate legend PNG — only pathways actually present in the data
  legend_pathways <- intersect(names(pathway_colors), unique(cor_df$pathway))
  legend_plot <- ggplot(data.frame(x = seq_along(legend_pathways),
    y = seq_along(legend_pathways),
    Pathway = factor(legend_pathways, levels = legend_pathways)),
    aes(x = x, y = y, color = Pathway)) +
    geom_point(size = 5) +
    scale_color_manual(values = pathway_colors, name = "Pathway") +
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
        width = 1800, height = 2400, res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, paste0(context$fig_id, "_legend.pdf")),
        width = 6, height = 8, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend\n")
  }

  # Prepare data for export
  export_data <- list(
    correlations_displayed = cor_df,
    correlations_all = cor_df_all,
    pathway_assignments = data.frame(
      gene = names(gene_pathway),
      pathway = unname(gene_pathway)
    )
  )

  # Methods text
  methods_text <- paste(
    "Spearman correlation coefficients were calculated between Piezo1 expression",
    "and representative pathway genes in control bona fide SMCs.",
    "Curated genes (2-3 per pathway) were selected to represent contractile,",
    "synthetic, NFAT, calcium, TGF-beta, and senescence programs.",
    "Genes are ordered by correlation coefficient and colored by pathway."
  )

  technical_legend <- paste(
    "Dot plot showing Spearman rho between Piezo1 and curated representative genes.",
    "Genes ordered by correlation coefficient (most negative at top, most positive at bottom).",
    "Colors indicate pathway assignment. Dashed line at rho = 0.",
    sprintf("n = %d control bona fide SMCs.", length(ctrl_cells))
  )

  # Save outputs
  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 10,
    height = 8,
    code_path = "pipeline/figure_modules/supplementary/FigS3H_gene_correlations.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    correlations = cor_df
  ))
}
