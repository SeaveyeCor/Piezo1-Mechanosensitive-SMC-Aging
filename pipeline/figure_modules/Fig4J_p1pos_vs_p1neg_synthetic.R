# ==============================================================================
# PAPER PANEL: Fig S4B   [internal module ID: Fig4J]
# PLOTS:  gene-programme scores in Piezo1+ vs Piezo1- synthetic SMCs
# CELLS:  CTRL synthetic-classified visceral SMCs only (Young + Aged CTRL)
# GROUPS COMPARED: Piezo1+ vs Piezo1- cells
# TEST:   Wilcoxon
# ==============================================================================
# Compares molecular profiles of Piezo1-expressing vs non-expressing synthetic
# SMCs across 5 gene panels: Synthetic, Cytoskeletal, NFAT/Ca2+, ECM, Contractile.
# Outputs: log2FC bar plot + module violin + stats XLSX
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

#' Generate Fig4J - P1+ vs P1- Synthetic SMC comparison
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

  gene_sets <- get_gene_sets()

  # Use fig_id from context (set by run_pipeline.R)
  fig_id <- context$fig_id

  # Get data
  seurat_obj <- context$seurat_obj
  smc_cells <- context$selection$smc
  phenotype_df <- context$phenotype
  metadata <- seurat_obj@meta.data

  # Expression matrix (SMC cells only)
  expr <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[, smc_cells])

  # --------------------------------------------------------------------------
  # Define populations: Ctrl bona fide SMCs, P1+ vs P1-
  # --------------------------------------------------------------------------
  smc_meta <- metadata[smc_cells, ]
  smc_meta$genotype <- ifelse(grepl("KO", smc_meta$sample), "KO", "Ctrl")
  ctrl_cells <- smc_cells[smc_meta$genotype == "Ctrl"]

  # P1+ vs P1-
  piezo1_raw <- as.numeric(expr["Piezo1", ctrl_cells])
  names(piezo1_raw) <- ctrl_cells
  p1_pos <- ctrl_cells[piezo1_raw > 0]
  p1_neg <- ctrl_cells[piezo1_raw == 0]

  # Synthetic cells from phenotype classification
  pheno_ctrl <- phenotype_df[phenotype_df$cell_barcode %in% ctrl_cells, ]
  synthetic_cells <- pheno_ctrl$cell_barcode[pheno_ctrl$phenotype == "Synthetic"]
  synthetic_p1pos <- intersect(synthetic_cells, p1_pos)
  synthetic_p1neg <- intersect(synthetic_cells, p1_neg)

  cat(sprintf("  Ctrl SMCs: %d (P1+: %d, P1-: %d)\n",
              length(ctrl_cells), length(p1_pos), length(p1_neg)))
  cat(sprintf("  Synthetic: %d (P1+: %d, P1-: %d)\n",
              length(synthetic_cells), length(synthetic_p1pos), length(synthetic_p1neg)))

  if (length(synthetic_p1pos) < 5 || length(synthetic_p1neg) < 5) {
    cat("  WARNING: Too few cells for meaningful comparison\n")
    return(list(success = FALSE, error = "Too few synthetic P1+ or P1- cells"))
  }

  # --------------------------------------------------------------------------
  # Gene panels
  # --------------------------------------------------------------------------
  panels <- list(
    "Synthetic Markers"  = list(genes = gene_sets$p1_comparison_synthetic, color = "#E41A1C"),
    "Cyto Remodeling"    = list(genes = gene_sets$p1_comparison_cyto, color = "#FF7F00"),
    "NFAT / Ca2+"        = list(genes = gene_sets$p1_comparison_nfat_ca, color = "#9467BD"),
    "ECM / Remodeling"   = list(genes = gene_sets$p1_comparison_ecm, color = "#17BECF"),
    "Contractile"        = list(genes = gene_sets$contractile_core, color = "#2CA02C")
  )

  # --------------------------------------------------------------------------
  # Per-gene statistics
  # --------------------------------------------------------------------------
  all_results <- data.frame()
  for (panel_name in names(panels)) {
    genes <- panels[[panel_name]]$genes
    genes <- genes[genes %in% rownames(expr)]

    for (g in genes) {
      vals_pos <- as.numeric(expr[g, synthetic_p1pos])
      vals_neg <- as.numeric(expr[g, synthetic_p1neg])

      mean_pos <- mean(vals_pos)
      mean_neg <- mean(vals_neg)
      det_pos <- 100 * mean(vals_pos > 0)
      det_neg <- 100 * mean(vals_neg > 0)
      log2fc <- log2(max(mean_pos, 0.001) / max(mean_neg, 0.001))

      pval <- NA
      if (length(vals_pos) >= 3 && length(vals_neg) >= 3) {
        pval <- wilcox.test(vals_pos, vals_neg)$p.value
      }
      sig <- ifelse(is.na(pval), "",
             ifelse(pval < 0.001, "***",
             ifelse(pval < 0.01, "**",
             ifelse(pval < 0.05, "*", "ns"))))

      all_results <- rbind(all_results, data.frame(
        panel = panel_name, gene = g,
        mean_p1pos = mean_pos, mean_p1neg = mean_neg,
        log2fc = log2fc, det_p1pos = det_pos, det_p1neg = det_neg,
        pval = pval, sig = sig,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --------------------------------------------------------------------------
  # Module-level statistics
  # --------------------------------------------------------------------------
  module_scores <- data.frame()
  for (panel_name in names(panels)) {
    genes <- panels[[panel_name]]$genes
    genes <- genes[genes %in% rownames(expr)]
    if (length(genes) == 0) next

    for (grp in c("P1+", "P1-")) {
      cells <- if (grp == "P1+") synthetic_p1pos else synthetic_p1neg
      scores <- colMeans(expr[genes, cells, drop = FALSE])
      module_scores <- rbind(module_scores, data.frame(
        panel = panel_name, group = grp, score = scores,
        stringsAsFactors = FALSE
      ))
    }
  }
  module_scores$panel <- factor(module_scores$panel, levels = names(panels))
  module_scores$group <- factor(module_scores$group, levels = c("P1+", "P1-"))

  # Bootstrap log2FC SE for each panel — resample cells with replacement
  # within each population, recompute the module-level log2FC, take the SD
  # across replicates as the standard error of the log2FC.
  set.seed(42)
  n_boot <- 1000
  boot_log2fc <- function(scores_pos, scores_neg, n) {
    out <- numeric(n)
    for (i in seq_len(n)) {
      sp <- sample(scores_pos, length(scores_pos), replace = TRUE)
      sn <- sample(scores_neg, length(scores_neg), replace = TRUE)
      out[i] <- log2(mean(sp) / max(mean(sn), 0.001))
    }
    out
  }

  mod_stats <- module_scores %>%
    group_by(panel) %>%
    summarize(
      p = wilcox.test(score[group == "P1+"], score[group == "P1-"])$p.value,
      mean_pos = mean(score[group == "P1+"]),
      mean_neg = mean(score[group == "P1-"]),
      .groups = "drop"
    ) %>%
    mutate(
      sig = ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
            ifelse(p < 0.05, "*", "ns"))),
      log2fc = log2(mean_pos / pmax(mean_neg, 0.001))
    )

  # Add bootstrap SE column
  mod_stats$log2fc_se <- NA_real_
  for (i in seq_len(nrow(mod_stats))) {
    pn <- as.character(mod_stats$panel[i])
    sp <- module_scores$score[module_scores$panel == pn & module_scores$group == "P1+"]
    sn <- module_scores$score[module_scores$panel == pn & module_scores$group == "P1-"]
    if (length(sp) >= 3 && length(sn) >= 3) {
      mod_stats$log2fc_se[i] <- sd(boot_log2fc(sp, sn, n_boot))
    }
  }

  # --------------------------------------------------------------------------
  # PRIMARY FIGURE (Fig4J main): Pathway-level log2FC summary (5 bars)
  # --------------------------------------------------------------------------
  # Compact 5-bar summary — one bar per gene-set panel — instead of the
  # ~50-bar per-gene plot. Per-gene detail moves to the supplementary figure
  # (Fig4J_supp_per_gene) and the supplementary table (Fig4J_data.xlsx ::
  # gene_level_stats sheet).
  pathway_bar_data <- mod_stats %>%
    mutate(panel = factor(panel, levels = names(panels)),
           panel_color = sapply(panel, function(x) panels[[as.character(x)]]$color),
           is_sig = sig %in% c("*", "**", "***"))

  # Errorbar tip x-coords (signed): bar end ± SE in the direction of the bar
  pathway_bar_data$err_lo <- pathway_bar_data$log2fc - pathway_bar_data$log2fc_se
  pathway_bar_data$err_hi <- pathway_bar_data$log2fc + pathway_bar_data$log2fc_se

  # x-axis padding accounting for error bar tips so annotations stay inside
  x_max <- max(pathway_bar_data$err_hi, na.rm = TRUE)
  x_min <- min(pathway_bar_data$err_lo, na.rm = TRUE)
  x_range <- max(abs(x_max), abs(x_min))
  x_pad <- x_range * 0.25
  # ns labels suppressed (only significant pathways get an asterisk).
  pathway_bar_data$display_label <- ifelse(pathway_bar_data$is_sig, pathway_bar_data$sig, "")
  # Position asterisk just beyond the error bar tip
  pathway_bar_data$label_x <- ifelse(
    pathway_bar_data$log2fc > 0,
    pathway_bar_data$err_hi + x_range * 0.03,
    pathway_bar_data$err_lo - x_range * 0.03
  )

  p_pathway <- ggplot(pathway_bar_data,
                       aes(x = log2fc, y = factor(panel, levels = rev(names(panels))),
                           fill = panel)) +
    geom_col(aes(alpha = is_sig), width = 0.6, color = "black", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = err_lo, xmax = err_hi),
                    height = 0.18, linewidth = 0.6, color = "black") +
    scale_fill_manual(values = sapply(panels, function(x) x$color), guide = "none") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.6) +
    geom_text(aes(label = display_label, x = label_x),
              size = 6, hjust = ifelse(pathway_bar_data$log2fc > 0, 0, 1),
              fontface = "bold") +
    scale_x_continuous(limits = c(x_min - x_pad, x_max + x_pad),
                        expand = expansion(mult = 0)) +
    labs(title = "P1+ vs P1- Synthetic SMCs: pathway-level log2 fold change",
         subtitle = sprintf("Module scores. P1+ (n=%d) vs P1- (n=%d). Wilcoxon test, faded = ns. See Fig4J_data.xlsx for per-gene values.",
                            length(synthetic_p1pos), length(synthetic_p1neg)),
         x = "log2FC (P1+ / P1-)") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "black"),
      axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  p_bar <- p_pathway  # used by save_figure_outputs below

  # --------------------------------------------------------------------------
  # SUPPLEMENTARY FIGURE (Fig4J_supp_per_gene): full per-gene log2FC bar plot
  # --------------------------------------------------------------------------
  gene_order <- all_results %>%
    arrange(factor(panel, levels = names(panels)), -abs(log2fc)) %>%
    pull(gene) %>% unique()

  per_gene_data <- all_results %>%
    mutate(gene = factor(gene, levels = rev(gene_order)),
           panel = factor(panel, levels = names(panels)),
           is_sig = sig %in% c("*", "**", "***"))

  p_per_gene <- ggplot(per_gene_data, aes(x = log2fc, y = gene, fill = panel)) +
    geom_col(aes(alpha = is_sig), width = 0.7) +
    scale_fill_manual(values = sapply(panels, function(x) x$color), name = "Panel") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    geom_text(aes(label = sig, x = ifelse(log2fc > 0, log2fc + 0.05, log2fc - 0.05)),
              size = 3.5, hjust = ifelse(per_gene_data$log2fc > 0, 0, 1)) +
    labs(title = "P1+ vs P1- Synthetic SMCs: per-gene log2 Fold Change (Supplement)",
         subtitle = sprintf("P1+ (n=%d) vs P1- (n=%d) | Wilcoxon test | Faded = ns",
                            length(synthetic_p1pos), length(synthetic_p1neg)),
         x = "log2FC (P1+ / P1-)") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "black"),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # --------------------------------------------------------------------------
  # MODULE-LEVEL VIOLIN (kept as a secondary panel)
  # --------------------------------------------------------------------------
  p_violin <- ggplot(module_scores, aes(x = group, y = score, fill = group)) +
    geom_violin(alpha = 0.6, linewidth = 0.3) +
    geom_boxplot(width = 0.2, outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("P1+" = "#D62728", "P1-" = "#7F7F7F"), guide = "none") +
    facet_wrap(~ panel, scales = "free_y", ncol = 5) +
    geom_text(data = mod_stats %>% mutate(panel = factor(panel, levels = names(panels))),
              aes(x = 1.5, y = Inf,
                  label = sprintf("%s\nlog2FC=%+.2f", sig, log2fc)),
              vjust = 1.2, size = 3.5, lineheight = 0.85, inherit.aes = FALSE) +
    labs(title = "Module Scores: P1+ vs P1- Synthetic SMCs",
         subtitle = sprintf("P1+ (n=%d, red) vs P1- (n=%d, grey) | Ctrl bona fide synthetic",
                            length(synthetic_p1pos), length(synthetic_p1neg)),
         y = "Mean Module Expression") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11, color = "black"),
      strip.text = element_text(face = "bold", size = 11),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # --------------------------------------------------------------------------
  # Save outputs
  # --------------------------------------------------------------------------
  fig_dir <- file.path(context$output_dir, fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  export_data <- list(
    gene_level_stats = all_results,
    module_level_stats = as.data.frame(mod_stats),
    population_summary = data.frame(
      population = c("Ctrl_bona_fide", "P1_pos", "P1_neg",
                     "Synthetic", "Synthetic_P1pos", "Synthetic_P1neg"),
      n_cells = c(length(ctrl_cells), length(p1_pos), length(p1_neg),
                  length(synthetic_cells), length(synthetic_p1pos), length(synthetic_p1neg))
    )
  )

  methods_text <- paste(
    "Piezo1-positive (raw count > 0) and Piezo1-negative synthetic SMCs were compared",
    "across 5 gene panels (Synthetic Markers, Cytoskeletal Remodeling, NFAT/Ca2+,",
    "ECM/Remodeling, Contractile). Synthetic cells were defined as contractile/synthetic",
    "expression ratio < 1.5. The primary panel (Fig S4B) shows pathway-level log2",
    "fold change (one bar per gene panel) computed from the mean module score in",
    "each population. Per-gene log2FC values are reported in the supplementary",
    "panel (Fig4J_supp_per_gene) and in the gene_level_stats sheet of the source",
    "data XLSX. Per-gene and per-module differences were tested by two-sided",
    "Wilcoxon rank-sum tests; analysis was restricted to Ctrl bona fide SMCs",
    "(clusters 1, 5, 15, 16)."
  )

  technical_legend <- paste(
    "Pathway-level log2 fold change (P1+/P1-) for the 5 gene panels in synthetic",
    "Ctrl bona fide SMCs. Bars show mean module expression log2FC; faded bars are",
    "non-significant (p >= 0.05). Per-gene values are in the supplementary",
    "Fig4J_supp_per_gene panel and in the gene_level_stats sheet of the source",
    "data file."
  )

  # Save pathway-level summary as the PRIMARY Fig4J output (compact 5-bar)
  save_figure_outputs(
    figure_id = fig_id,
    plot = p_bar,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 8, height = 4,
    code_path = "pipeline/figure_modules/Fig4J_p1pos_vs_p1neg_synthetic.R"
  )

  # Save the full per-gene plot as a supplementary companion.
  ggsave(file.path(fig_dir, paste0(fig_id, "_supp_per_gene.png")), p_per_gene,
         width = 10, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(fig_id, "_supp_per_gene.pdf")), p_per_gene,
         width = 10, height = 12)
  cat(sprintf("  Saved: %s_supp_per_gene.png/pdf  (per-gene supplement)\n", fig_id))

  # Save the module violin as an additional panel (kept for backward compat).
  ggsave(file.path(fig_dir, paste0(fig_id, "_modules.png")), p_violin,
         width = 16, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(fig_id, "_modules.pdf")), p_violin,
         width = 16, height = 6)
  cat(sprintf("  Saved: %s_modules.png\n", fig_id))

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    gene_stats = all_results,
    module_stats = mod_stats
  ))
}
