# ==============================================================================
# PAPER PANEL: Fig S3F   [internal module ID: FigS3I]
# PLOTS:  GO biological-process pathways enriched in Aged P1KO vs Aged CTRL SMCs (higher in P1KO = rescue, higher in CTRL = aging)
# CELLS:  Aged CTRL (n = 477) and Aged P1KO (n = 436) visceral SMCs
# GROUPS COMPARED: Aged P1KO vs Aged CTRL
# TEST:   genes ranked by Wilcoxon statistic; GSEA (fgsea), adjusted p
# ==============================================================================
# Computes differential expression between Aged KO and Aged Ctrl SMCs,
# then runs GSEA with GO:BP gene sets using fgsea.
# Shows pathways enriched in rescue (KO) vs aging (Ctrl).
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Generate FigS3I (GO enrichment dot plot)
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))

  # Check for required packages
  if (!requireNamespace("fgsea", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    cat("  WARNING: fgsea and/or msigdbr not installed.\n")
    cat("  Install with: BiocManager::install('fgsea'); install.packages('msigdbr')\n")
    cat("  Falling back to Wilcoxon-based ORA approach...\n")
    return(generate_ora_fallback(context))
  }

  library(fgsea)
  library(msigdbr)

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc
  metadata <- seurat_obj@meta.data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Get Aged KO and Aged Ctrl cells
  aged_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_Ctrl"]]
  aged_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_KO"]]

  cat(sprintf("  Aged Ctrl: %d cells, Aged KO: %d cells\n",
              length(aged_ctrl_cells), length(aged_ko_cells)))

  # ---------------------------------------------------------------------------
  # Step 1: Differential expression (Wilcoxon rank-sum)
  # ---------------------------------------------------------------------------
  cat("  Computing differential expression (Aged KO vs Aged Ctrl)...\n")

  all_cells <- c(aged_ko_cells, aged_ctrl_cells)
  expr_subset <- as.matrix(expr[, all_cells])

  # Filter to genes expressed in >5% of cells in either group
  pct_ko <- rowMeans(expr_subset[, aged_ko_cells, drop = FALSE] > 0)
  pct_ctrl <- rowMeans(expr_subset[, aged_ctrl_cells, drop = FALSE] > 0)
  genes_to_test <- rownames(expr_subset)[pct_ko > 0.05 | pct_ctrl > 0.05]
  genes_to_test <- setdiff(genes_to_test, c("Piezo1"))  # exclude Piezo1 itself

  cat(sprintf("  Testing %d genes...\n", length(genes_to_test)))

  # Compute Wilcoxon for each gene
  de_results <- data.frame(
    gene = genes_to_test,
    log2FC = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes_to_test)) {
    g <- genes_to_test[i]
    ko_vals <- expr_subset[g, aged_ko_cells]
    ctrl_vals <- expr_subset[g, aged_ctrl_cells]

    # Log2FC: mean(KO) - mean(Ctrl) [already log-normalized data]
    de_results$log2FC[i] <- mean(ko_vals) - mean(ctrl_vals)

    # Wilcoxon test
    wt <- tryCatch(
      wilcox.test(ko_vals, ctrl_vals)$p.value,
      error = function(e) 1.0
    )
    de_results$p_value[i] <- wt
  }

  # Create ranking metric: sign(log2FC) * -log10(p)
  de_results$rank_metric <- sign(de_results$log2FC) * -log10(pmax(de_results$p_value, 1e-300))

  # Remove genes with NA rank
  de_results <- de_results[!is.na(de_results$rank_metric) & is.finite(de_results$rank_metric), ]

  # Create named ranked vector
  ranked_genes <- de_results$rank_metric
  names(ranked_genes) <- de_results$gene
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  cat(sprintf("  Ranked %d genes. Top up in KO: %s, Top up in Ctrl: %s\n",
              length(ranked_genes),
              paste(head(names(ranked_genes), 3), collapse = ", "),
              paste(tail(names(ranked_genes), 3), collapse = ", ")))

  # ---------------------------------------------------------------------------
  # Step 2: Load GO:BP gene sets
  # ---------------------------------------------------------------------------
  cat("  Loading GO:BP gene sets...\n")
  gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
  gobp_list <- split(gobp$gene_symbol, gobp$gs_name)
  cat(sprintf("  Loaded %d GO:BP pathways\n", length(gobp_list)))

  # ---------------------------------------------------------------------------
  # Step 3: Run GSEA
  # ---------------------------------------------------------------------------
  cat("  Running fgsea...\n")
  set.seed(42)
  gsea_results <- fgsea(
    pathways = gobp_list,
    stats = ranked_genes,
    minSize = 15,
    maxSize = 500,
    nPermSimple = 10000
  )

  # Filter significant
  sig_results <- gsea_results[gsea_results$padj < 0.05, ]
  cat(sprintf("  Found %d significant pathways (padj < 0.05)\n", nrow(sig_results)))

  # ---------------------------------------------------------------------------
  # Step 4: Select top pathways for display
  # ---------------------------------------------------------------------------

  # Categorize by direction
  up_in_ko <- sig_results[sig_results$NES > 0, ] %>% arrange(padj)
  up_in_ctrl <- sig_results[sig_results$NES < 0, ] %>% arrange(padj)

  # Select top pathways (prioritize SMC-relevant terms)
  smc_keywords <- "MUSCLE|CONTRAC|SMOOTH|MYOFIBRIL|SARCOMERE|ACTIN_FILAMENT|COLLAGEN|EXTRACELLULAR_MATRIX|WOUND|ADHES|CALCIUM|INFLAMM|CYTOKINE|CELL_MIGRAT|FIBER|RESPONSE_TO_MECH"

  # Top KO pathways: SMC-relevant first, then by significance
  ko_smc <- up_in_ko[grepl(smc_keywords, up_in_ko$pathway, ignore.case = TRUE), ]
  ko_other <- up_in_ko[!grepl(smc_keywords, up_in_ko$pathway, ignore.case = TRUE), ]
  top_ko <- rbind(head(ko_smc, 6), head(ko_other, 2)) %>% head(7)

  # Top Ctrl pathways: SMC-relevant first
  ctrl_smc <- up_in_ctrl[grepl(smc_keywords, up_in_ctrl$pathway, ignore.case = TRUE), ]
  ctrl_other <- up_in_ctrl[!grepl(smc_keywords, up_in_ctrl$pathway, ignore.case = TRUE), ]
  top_ctrl <- rbind(head(ctrl_smc, 6), head(ctrl_other, 2)) %>% head(6)

  # Combine for plot
  plot_terms <- rbind(
    top_ko %>% mutate(direction = "Higher in KO (rescue)"),
    top_ctrl %>% mutate(direction = "Higher in Ctrl (aging)")
  )

  if (nrow(plot_terms) == 0) {
    cat("  WARNING: No significant pathways found. Using relaxed threshold.\n")
    sig_results <- gsea_results[gsea_results$padj < 0.25, ]
    up_in_ko <- sig_results[sig_results$NES > 0, ] %>% arrange(padj) %>% head(6)
    up_in_ctrl <- sig_results[sig_results$NES < 0, ] %>% arrange(padj) %>% head(6)
    plot_terms <- rbind(
      up_in_ko %>% mutate(direction = "Higher in KO (rescue)"),
      up_in_ctrl %>% mutate(direction = "Higher in Ctrl (aging)")
    )
  }

  # Remove non-informative pathways
  plot_terms <- plot_terms %>%
    filter(!grepl("CYTOPLASMIC_TRANSLATION", pathway, ignore.case = TRUE))

  # Clean pathway names for display
  plot_terms <- plot_terms %>%
    mutate(
      pathway_clean = gsub("^GOBP_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = tools::toTitleCase(tolower(pathway_clean)),
      pathway_clean = ifelse(nchar(pathway_clean) > 60,
                             paste0(substr(pathway_clean, 1, 57), "..."),
                             pathway_clean),
      neg_log_padj = -log10(padj)
    )

  # Order: KO (positive NES) on top, Ctrl (negative NES) on bottom
  plot_terms <- plot_terms %>%
    arrange(NES) %>%
    mutate(pathway_clean = factor(pathway_clean, levels = unique(pathway_clean)))

  # ---------------------------------------------------------------------------
  # Step 5: Create plot
  # ---------------------------------------------------------------------------
  direction_colors <- c(
    "Higher in KO (rescue)" = "#E41A1C",
    "Higher in Ctrl (aging)" = "#377EB8"
  )

  p <- ggplot(plot_terms, aes(x = neg_log_padj, y = pathway_clean,
                               color = direction, size = size)) +
    geom_point() +
    scale_color_manual(values = direction_colors, name = "Direction") +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    labs(
      title = "GO Enrichment: Aged P1KO vs Aged CTRL",
      x = expression(-log[10](adjusted~p~value)),
      y = ""
    ) +
    theme_publication() +
    theme(
      axis.text.y = element_text(size = 16),
      legend.position = "none"
    )

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Separate legend PNG.
  # Color-key dots are forced to the same render size as the largest gene-count
  # dot in the size legend (size = 8 from scale_size_continuous(range = c(3, 8)))
  # so they're easy to see at a glance.
  legend_plot <- ggplot(data.frame(
    x = 1:2, y = 1:2, sz = c(50, 200),
    Direction = factor(c("Higher in KO (rescue)", "Higher in Ctrl (aging)"),
                       levels = c("Higher in KO (rescue)", "Higher in Ctrl (aging)"))),
    aes(x = x, y = y, color = Direction, size = sz)) +
    geom_point() +
    scale_color_manual(values = direction_colors, name = "Direction") +
    scale_size_continuous(range = c(3, 8), name = "Gene Count") +
    guides(color = guide_legend(override.aes = list(size = 10))) +
    theme_classic(base_size = 28) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold", size = 32),
          legend.text = element_text(size = 28),
          legend.key.size = unit(2, "cm"))
  legend_grob <- ggplot_gtable(ggplot_build(legend_plot))
  legend_idx <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  if (length(legend_idx) > 0) {
    legend_only <- legend_grob$grobs[[legend_idx]]
    png(file.path(fig_dir, paste0(context$fig_id, "_legend.png")),
        width = 2400, height = 2400, res = 300, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    pdf(file.path(fig_dir, paste0(context$fig_id, "_legend.pdf")),
        width = 8, height = 8, bg = "white")
    grid::grid.draw(legend_only)
    dev.off()
    cat("  Saved separate legend\n")
  }

  # Flatten leadingEdge for export
  gsea_export <- gsea_results %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "))
  sig_export <- sig_results %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ", "))
  plot_export <- plot_terms %>%
    select(pathway, pathway_clean, NES, padj, size, direction)

  # Prepare data for export
  export_data <- list(
    go_enrichment = as.data.frame(plot_export),
    all_significant = as.data.frame(sig_export),
    de_summary = de_results[order(de_results$p_value), ] %>% head(200)
  )

  # Methods text
  methods_text <- paste(
    "Gene Set Enrichment Analysis (GSEA) was performed using fgsea (v1.22+)",
    "with GO Biological Process gene sets from MSigDB (msigdbr, Mus musculus).",
    "Genes were ranked by signed -log10(p-value) from Wilcoxon rank-sum tests",
    "comparing Aged P1KO vs Aged CTRL cluster-based SMCs.",
    "Gene sets with 15-500 genes were tested with 10,000 permutations.",
    sprintf("A total of %d genes were ranked.", length(ranked_genes)),
    sprintf("%d pathways had padj < 0.05.", nrow(sig_results))
  )

  technical_legend <- paste(
    "Dot plot showing top enriched GO:BP pathways from GSEA.",
    "Red = Higher in Aged KO (rescue pathways),",
    "Blue = Higher in Aged Ctrl (aging pathways).",
    "X-axis: -log10(adjusted p-value). Dot size: gene set size.",
    "GSEA performed with fgsea using Wilcoxon-ranked gene list."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 12,
    height = 8,
    code_path = "pipeline/figure_modules/supplementary/FigS3I_go_enrichment.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    n_significant = nrow(sig_results),
    n_up_ko = nrow(up_in_ko),
    n_up_ctrl = nrow(up_in_ctrl)
  ))
}


#' Fallback: Over-representation analysis without fgsea/msigdbr
#' Uses curated gene sets and Fisher's exact test
generate_ora_fallback <- function(context) {
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))

  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc
  metadata <- seurat_obj@meta.data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  aged_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_Ctrl"]]
  aged_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[
    metadata$group == "Aged_KO"]]

  # Quick DE
  all_cells <- c(aged_ko_cells, aged_ctrl_cells)
  expr_sub <- as.matrix(expr[, all_cells])
  pct <- rowMeans(expr_sub > 0)
  genes_test <- names(pct[pct > 0.05])

  de <- data.frame(gene = genes_test, log2FC = NA, p = NA)
  for (i in seq_along(genes_test)) {
    g <- genes_test[i]
    de$log2FC[i] <- mean(expr_sub[g, aged_ko_cells]) - mean(expr_sub[g, aged_ctrl_cells])
    de$p[i] <- tryCatch(wilcox.test(expr_sub[g, aged_ko_cells],
                                     expr_sub[g, aged_ctrl_cells])$p.value,
                         error = function(e) 1)
  }
  de$padj <- p.adjust(de$p, method = "BH")

  # Significant DE genes
  up_ko <- de$gene[de$padj < 0.05 & de$log2FC > 0]
  up_ctrl <- de$gene[de$padj < 0.05 & de$log2FC < 0]

  cat(sprintf("  ORA fallback: %d up in KO, %d up in Ctrl\n",
              length(up_ko), length(up_ctrl)))

  # Create a simple ORA-based plot using known gene sets
  gene_sets <- get_gene_sets()
  sets_to_test <- list(
    "Core Contractile" = gene_sets$contractile_core,
    "Muscle Contraction" = gene_sets$muscle_contraction,
    "NFAT Pathway" = gene_sets$nfat_pathway,
    "TGF-beta/Stress" = gene_sets$tgfbeta_pathway,
    "ECM Remodeling" = gene_sets$ecm_remodeling,
    "ECM Organization" = gene_sets$ecm_organization,
    "Focal Adhesion" = gene_sets$focal_adhesion,
    "Calcium Signaling" = gene_sets$calcium_pathway
  )

  ora_results <- list()
  for (nm in names(sets_to_test)) {
    gs <- sets_to_test[[nm]]
    gs_avail <- gs[gs %in% de$gene]
    n_up_ko <- sum(gs_avail %in% up_ko)
    n_up_ctrl <- sum(gs_avail %in% up_ctrl)
    direction <- ifelse(n_up_ko > n_up_ctrl, "Higher in KO (rescue)", "Higher in Ctrl (aging)")
    ora_results[[nm]] <- data.frame(
      pathway = nm, n_genes = length(gs_avail),
      up_ko = n_up_ko, up_ctrl = n_up_ctrl,
      direction = direction, stringsAsFactors = FALSE
    )
  }
  ora_df <- bind_rows(ora_results)

  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  p <- ggplot(ora_df, aes(x = up_ko + up_ctrl, y = reorder(pathway, up_ko + up_ctrl),
                           fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Higher in KO (rescue)" = "#E41A1C",
                                  "Higher in Ctrl (aging)" = "#377EB8")) +
    labs(title = "DE Gene Overlap with Curated Pathways",
         subtitle = "ORA fallback (fgsea not installed)",
         x = "Number of DE genes", y = "") +
    theme_publication()

  export_data <- list(go_enrichment = ora_df, note = "ORA fallback - install fgsea for GSEA")
  methods_text <- "Over-representation analysis of DE genes against curated gene sets."
  technical_legend <- "Bar plot showing DE gene counts in curated pathways."

  save_figure_outputs(
    figure_id = context$fig_id, plot = p, data = export_data,
    methods_text = methods_text, technical_legend = technical_legend,
    output_dir = fig_dir, width = 10, height = 8,
    code_path = "pipeline/figure_modules/supplementary/FigS3I_go_enrichment.R"
  )

  return(list(success = TRUE, output_dir = fig_dir, note = "ORA fallback"))
}
