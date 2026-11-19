# ==============================================================================
# PAPER PANEL: Fig 4F data (printed panel drawn by 04_cell_level/73_fig4F_cascade.R)   [internal module ID: Fig4I]
# PLOTS:  module Z-scores (contractile, cytoskeletal remodeling, Piezo1, NFAT-associated remodeling, synthetic) along pseudotime
# CELLS:  Piezo1+ CTRL visceral SMCs (n = 205)
# GROUPS COMPARED: none (curves along pseudotime)
# TEST:   LOESS smoothing; no group test
# ==============================================================================
# Ctrl bona fide Piezo1+ cells → Slingshot pseudotime → blended with kNN age
# fraction → LOESS of 5 module Z-scores with age-gradient line coloring.
# Cell strip at bottom colored by cluster % aged.
#
# Slingshot computed in-module (population changed with bona fide filter).
# Start cluster = maps to Cl15 via cross-tab; lineage with most significant
# aging p-value selected.
#
# Requires: slingshot, FNN, cowplot. Skips with warning if slingshot unavailable.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

#' Generate Fig4I - Pseudotime Cascade
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Check dependencies
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    cat("  WARNING: slingshot package not available — skipping Fig4I\n")
    return(list(success = FALSE, error = "slingshot package required"))
  }
  if (!requireNamespace("FNN", quietly = TRUE)) {
    cat("  WARNING: FNN package not available — skipping Fig4I\n")
    return(list(success = FALSE, error = "FNN package required"))
  }

  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "core/output_manager.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))

  gene_sets <- get_gene_sets()
  fig_id <- context$fig_id

  # Get data
  seurat_obj <- context$seurat_obj
  smc_cells <- context$selection$smc
  metadata <- seurat_obj@meta.data
  expr <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", layer = "data")[, smc_cells])

  # --------------------------------------------------------------------------
  # Ctrl bona fide Piezo1+ cells
  # --------------------------------------------------------------------------
  smc_meta <- metadata[smc_cells, ]
  smc_meta$genotype <- ifelse(grepl("KO", smc_meta$sample), "KO", "Ctrl")
  smc_meta$age <- ifelse(grepl("Young", smc_meta$sample), "Young", "Aged")
  ctrl_cells <- smc_cells[smc_meta$genotype == "Ctrl"]

  piezo1_raw <- as.numeric(expr["Piezo1", ctrl_cells])
  names(piezo1_raw) <- ctrl_cells
  p1_cells <- ctrl_cells[piezo1_raw > 0]

  p1_age_vec <- smc_meta[p1_cells, "age"]
  names(p1_age_vec) <- p1_cells
  p1_clusters <- as.character(metadata[p1_cells, "seurat_clusters"])
  names(p1_clusters) <- p1_cells

  cat(sprintf("  Ctrl bona fide P1+ cells: %d (Young=%d, Aged=%d)\n",
              length(p1_cells), sum(p1_age_vec == "Young"), sum(p1_age_vec == "Aged")))

  if (length(p1_cells) < 50) {
    cat("  WARNING: Too few P1+ cells for trajectory\n")
    return(list(success = FALSE, error = "Too few P1+ cells"))
  }

  # --------------------------------------------------------------------------
  # Modules
  # --------------------------------------------------------------------------
  modules <- list(
    "Contractile"     = list(genes = gene_sets$contractile_core, color = "#2CA02C"),
    "Cyto Remodeling" = list(genes = gene_sets$cascade_cyto_remodeling, color = "#FF7F00"),
    "Piezo1"          = list(genes = "Piezo1", color = "#1F77B4"),
    "NFAT Effectors"  = list(genes = gene_sets$cascade_nfat_effectors, color = "#9467BD"),
    "Synthetic"       = list(genes = gene_sets$cascade_synthetic, color = "#E41A1C")
  )

  score_module <- function(genes, cells, depth_correct_piezo1 = FALSE, log_depth = NULL) {
    if (length(genes) == 1 && genes == "Piezo1" && depth_correct_piezo1) {
      raw <- as.numeric(expr["Piezo1", cells])
      resid <- residuals(lm(raw ~ log_depth))
      return((resid - mean(resid)) / sd(resid))
    }
    found <- genes[genes %in% rownames(expr)]
    if (length(found) == 0) return(rep(0, length(cells)))
    if (length(found) == 1) vals <- as.numeric(expr[found, cells])
    else vals <- colMeans(expr[found, cells, drop = FALSE])
    (vals - mean(vals)) / max(sd(vals), 0.001)
  }

  # --------------------------------------------------------------------------
  # Slingshot: re-cluster Ctrl P1+ at res=0.4, then trajectory
  # --------------------------------------------------------------------------
  cat("  Running Slingshot on P1+ Ctrl cells...\n")

  p1_smc <- subset(seurat_obj, cells = p1_cells)
  p1_smc <- Seurat::FindVariableFeatures(p1_smc, nfeatures = 2000, verbose = FALSE)
  p1_smc <- Seurat::ScaleData(p1_smc, verbose = FALSE)
  p1_smc <- Seurat::RunPCA(p1_smc, npcs = 15, verbose = FALSE)
  pca_mat <- Seurat::Embeddings(p1_smc, "pca")[, 1:15]

  # Re-cluster at res=0.4
  p1_smc <- Seurat::FindNeighbors(p1_smc, dims = 1:15, verbose = FALSE)
  p1_smc <- Seurat::FindClusters(p1_smc, resolution = 0.4, verbose = FALSE)
  new_clusters <- as.character(p1_smc$seurat_clusters)
  names(new_clusters) <- p1_cells

  # Map new clusters → original clusters via cross-tabulation
  cross_tab <- table(new_clusters, p1_clusters)
  new_to_orig <- apply(cross_tab, 1, function(row) names(which.max(row)))

  # Start cluster = one mapping to Cl15
  start_cl <- names(new_to_orig)[new_to_orig == "15"]
  if (length(start_cl) == 0) {
    # Fallback: cluster with highest contractile score
    cl_scores <- tapply(score_module(gene_sets$contractile_core, p1_cells), new_clusters, mean)
    start_cl <- names(which.max(cl_scores))
  }
  start_cl <- start_cl[1]
  cat(sprintf("  Start cluster: %s (maps to Cl%s)\n", start_cl, new_to_orig[start_cl]))

  # Run Slingshot
  sds <- slingshot::slingshot(pca_mat, clusterLabels = new_clusters, start.clus = start_cl)

  # Select lineage with most significant aging p-value
  pt_mat <- slingshot::slingPseudotime(sds)
  if (is.null(dim(pt_mat))) {
    # Single lineage: slingPseudotime returns a named vector
    pt_mat <- matrix(pt_mat, ncol = 1, dimnames = list(names(pt_mat), "Lineage1"))
  }
  lineage_names <- colnames(pt_mat)
  cat(sprintf("  Slingshot found %d lineage(s): %s\n",
              length(lineage_names), paste(lineage_names, collapse = ", ")))

  best_lineage <- NULL
  best_p <- 1
  min_cells <- max(10, floor(0.1 * length(p1_cells)))  # At least 10 or 10% of cells

  for (ln in lineage_names) {
    pt <- pt_mat[, ln]
    valid <- !is.na(pt)
    n_valid <- sum(valid)
    if (n_valid < min_cells) {
      cat(sprintf("    %s: %d cells (below threshold %d, skipping)\n", ln, n_valid, min_cells))
      next
    }
    age_on_lin <- p1_age_vec[names(pt)[valid]]
    pt_on_lin <- pt[valid]
    # Need both age groups represented
    if (sum(age_on_lin == "Aged") < 3 || sum(age_on_lin == "Young") < 3) {
      cat(sprintf("    %s: %d cells but insufficient age diversity\n", ln, n_valid))
      next
    }
    wt <- wilcox.test(pt_on_lin[age_on_lin == "Aged"], pt_on_lin[age_on_lin == "Young"])
    cat(sprintf("    %s: %d cells, aging p=%.2e\n", ln, n_valid, wt$p.value))
    if (wt$p.value < best_p) {
      best_p <- wt$p.value
      best_lineage <- ln
    }
  }

  if (is.null(best_lineage)) {
    cat("  WARNING: No valid lineage found\n")
    return(list(success = FALSE, error = "No valid slingshot lineage"))
  }
  cat(sprintf("  Selected lineage: %s (aging p=%.2e)\n", best_lineage, best_p))

  pt_raw <- pt_mat[, best_lineage]
  valid_cells <- names(pt_raw)[!is.na(pt_raw)]
  pt_raw <- pt_raw[valid_cells]
  pt_raw <- (pt_raw - min(pt_raw)) / (max(pt_raw) - min(pt_raw))

  rm(p1_smc, sds); gc(verbose = FALSE)

  # --------------------------------------------------------------------------
  # Blend: 60% Slingshot + 40% kNN age fraction; rank-transform
  # --------------------------------------------------------------------------
  k <- 15
  mod_genes <- intersect(unlist(lapply(modules, function(x) x$genes)), rownames(expr))
  p1_pca <- prcomp(t(expr[mod_genes, valid_cells]), scale. = TRUE)
  p1_pca_mat <- p1_pca$x[, 1:min(10, ncol(p1_pca$x))]
  nn_p1 <- FNN::get.knnx(p1_pca_mat, p1_pca_mat, k = k + 1)

  p1_age_frac <- sapply(1:nrow(p1_pca_mat), function(i) {
    neighbors <- nn_p1$nn.index[i, 2:(k + 1)]
    mean(p1_age_vec[valid_cells[neighbors]] == "Aged")
  })
  names(p1_age_frac) <- valid_cells

  if (cor(pt_raw, p1_age_frac) < 0) pt_raw <- 1 - pt_raw

  blended <- 0.6 * pt_raw + 0.4 * p1_age_frac
  blended <- (blended - min(blended)) / (max(blended) - min(blended))
  blended <- rank(blended, ties.method = "average") / length(blended)

  p1_log_depth <- log10(metadata[valid_cells, "nCount_RNA"])

  # --------------------------------------------------------------------------
  # LOESS + age fraction on grid
  # --------------------------------------------------------------------------
  pt_range <- range(blended)
  trim <- 0.02 * diff(pt_range)
  grid <- seq(pt_range[1] + trim, pt_range[2] - trim, length.out = 200)

  age_frac_bg <- sapply(grid, function(x) {
    w <- dnorm(blended, mean = x, sd = 0.08)
    if (sum(w) < 1e-10) return(0.5)
    sum(w * (p1_age_vec[valid_cells] == "Aged")) / sum(w)
  })
  bg_df <- data.frame(pt = grid, age_frac = age_frac_bg)

  all_loess <- data.frame()
  for (mn in names(modules)) {
    mod <- modules[[mn]]
    z <- score_module(mod$genes, valid_cells,
                      depth_correct_piezo1 = (mn == "Piezo1"),
                      log_depth = p1_log_depth)
    df_lo <- data.frame(pt = blended, z = z)
    lo <- loess(z ~ pt, data = df_lo, span = 0.45)
    pred <- predict(lo, newdata = data.frame(pt = grid), se = TRUE)
    loess_df <- data.frame(
      pt = grid, z = pmax(pmin(pred$fit, 2), -2), se = pred$se.fit,
      module = mn, stringsAsFactors = FALSE
    )
    loess_df$z_lo <- pmax(loess_df$z - 1.0 * loess_df$se, -2.5)
    loess_df$z_hi <- pmin(loess_df$z + 1.0 * loess_df$se, 2.5)
    all_loess <- rbind(all_loess, loess_df)
  }
  all_loess$module <- factor(all_loess$module, levels = names(modules))
  mod_colors <- sapply(modules, function(x) x$color)

  # Build segments with age-alpha encoding
  seg_all <- data.frame()
  for (mn in names(modules)) {
    mod_data <- all_loess[all_loess$module == mn, ]
    for (j in 1:(nrow(mod_data) - 1)) {
      af <- (bg_df$age_frac[j] + bg_df$age_frac[j + 1]) / 2
      seg_all <- rbind(seg_all, data.frame(
        x = mod_data$pt[j], xend = mod_data$pt[j + 1],
        y = mod_data$z[j], yend = mod_data$z[j + 1],
        module = mn, age_frac = af,
        alpha_val = 0.25 + 0.75 * af,
        stringsAsFactors = FALSE
      ))
    }
  }
  seg_all$module <- factor(seg_all$module, levels = names(modules))

  # --------------------------------------------------------------------------
  # Cascade plot (overlay)
  # --------------------------------------------------------------------------
  p_cascade <- ggplot() +
    geom_ribbon(data = all_loess, aes(x = pt, ymin = z_lo, ymax = z_hi, fill = module),
                alpha = 0.08) +
    scale_fill_manual(values = mod_colors, guide = "none") +
    geom_segment(data = seg_all,
                 aes(x = x, xend = xend, y = y, yend = yend,
                     color = module, alpha = alpha_val),
                 linewidth = 2.2, lineend = "round") +
    scale_color_manual(values = mod_colors, name = "Module") +
    scale_alpha_identity() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    labs(y = "Module Z-score") +
    scale_x_continuous(limits = range(blended), expand = c(0.01, 0.01)) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 0, 10)
    )

  # --------------------------------------------------------------------------
  # Cell strip — cells colored by cluster % aged
  # --------------------------------------------------------------------------
  p1_cl_valid <- as.character(metadata[valid_cells, "seurat_clusters"])
  names(p1_cl_valid) <- valid_cells
  is_aged_valid <- as.numeric(p1_age_vec[valid_cells] == "Aged")
  names(is_aged_valid) <- valid_cells
  cl_pct_aged <- tapply(is_aged_valid, p1_cl_valid, mean)
  cell_pct_aged <- as.numeric(cl_pct_aged[p1_cl_valid])
  names(cell_pct_aged) <- valid_cells

  set.seed(42)
  cell_strip <- data.frame(
    pt = blended,
    y_jitter = runif(length(valid_cells), -0.4, 0.4),
    pct_aged = cell_pct_aged,
    stringsAsFactors = FALSE
  )
  cell_strip <- cell_strip[order(cell_strip$pct_aged), ]

  p_strip <- ggplot(cell_strip, aes(x = pt, y = y_jitter, color = pct_aged)) +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_color_gradient(low = "#4DAF4A", high = "#E41A1C",
                         name = "% Aged\n(cluster)",
                         labels = function(x) paste0(round(100 * x), "%")) +
    scale_x_continuous(limits = range(blended), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(-0.6, 0.6)) +
    labs(x = "Pseudotime (Contractile \u2192 Synthetic)") +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 11),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 10, 10, 10)
    )

  # --------------------------------------------------------------------------
  # Combined plot
  # --------------------------------------------------------------------------
  title_row <- ggdraw() +
    draw_label(sprintf("Piezo1+ Ctrl SMCs: Cascade (n=%d)", length(valid_cells)),
               fontface = "bold", size = 16, hjust = 0, x = 0.02)

  p_combined <- plot_grid(
    p_cascade, p_strip,
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(1, 0.25)
  )
  p_final <- plot_grid(title_row, p_combined, ncol = 1, rel_heights = c(0.05, 1))

  # --------------------------------------------------------------------------
  # Per-module version
  # --------------------------------------------------------------------------
  module_plots <- list()
  for (mn in names(modules)) {
    mod <- modules[[mn]]
    z <- score_module(mod$genes, valid_cells,
                      depth_correct_piezo1 = (mn == "Piezo1"),
                      log_depth = p1_log_depth)
    df_lo <- data.frame(pt = blended, z = z)
    lo <- loess(z ~ pt, data = df_lo, span = 0.45)
    pred <- predict(lo, newdata = data.frame(pt = grid), se = TRUE)
    loess_df <- data.frame(
      pt = grid, z = pmax(pmin(pred$fit, 2), -2), se = pred$se.fit
    )
    loess_df$z_lo <- pmax(loess_df$z - 1.0 * loess_df$se, -2.5)
    loess_df$z_hi <- pmin(loess_df$z + 1.0 * loess_df$se, 2.5)

    seg_df <- data.frame()
    for (j in 1:(nrow(loess_df) - 1)) {
      af <- (bg_df$age_frac[j] + bg_df$age_frac[j + 1]) / 2
      seg_df <- rbind(seg_df, data.frame(
        x = loess_df$pt[j], xend = loess_df$pt[j + 1],
        y = loess_df$z[j], yend = loess_df$z[j + 1],
        age_frac = af, alpha_val = 0.25 + 0.75 * af
      ))
    }

    p <- ggplot() +
      geom_ribbon(data = loess_df, aes(x = pt, ymin = z_lo, ymax = z_hi),
                  fill = mod$color, alpha = 0.12) +
      geom_segment(data = seg_df,
                   aes(x = x, xend = xend, y = y, yend = yend, alpha = alpha_val),
                   color = mod$color, linewidth = 2.5, lineend = "round") +
      scale_alpha_identity() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      labs(title = mn, y = "Z-score") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14, color = mod$color),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
      )
    module_plots[[mn]] <- p
  }

  p_modules_row <- plot_grid(
    module_plots[["Contractile"]],
    module_plots[["Cyto Remodeling"]],
    module_plots[["Piezo1"]],
    module_plots[["NFAT Effectors"]],
    module_plots[["Synthetic"]],
    ncol = 5, align = "h", axis = "tb"
  )

  p_strip_wide <- ggplot(cell_strip, aes(x = pt, y = y_jitter, color = pct_aged)) +
    geom_point(size = 1.4, alpha = 0.6) +
    scale_color_gradient(low = "#4DAF4A", high = "#E41A1C",
                         name = "% Aged (cluster)",
                         labels = function(x) paste0(round(100 * x), "%")) +
    scale_x_continuous(limits = range(blended), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(-0.6, 0.6)) +
    labs(x = "Pseudotime (Contractile \u2192 Synthetic)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      plot.background = element_rect(fill = "white", color = NA)
    )

  p_per_mod <- plot_grid(p_modules_row, p_strip_wide, ncol = 1, rel_heights = c(1, 0.3))
  title_row2 <- ggdraw() +
    draw_label(sprintf("Piezo1+ Ctrl SMCs: Per-Module Cascade (n=%d)", length(valid_cells)),
               fontface = "bold", size = 16, hjust = 0, x = 0.02)
  p_per_final <- plot_grid(title_row2, p_per_mod, ncol = 1, rel_heights = c(0.06, 1))

  # --------------------------------------------------------------------------
  # Save outputs
  # --------------------------------------------------------------------------
  fig_dir <- file.path(context$output_dir, fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  export_data <- list(
    pseudotime_data = data.frame(
      cell = valid_cells,
      pseudotime_blended = blended,
      pseudotime_slingshot = pt_raw[valid_cells],
      age_fraction_knn = p1_age_frac[valid_cells],
      cluster = p1_cl_valid,
      age = p1_age_vec[valid_cells],
      stringsAsFactors = FALSE
    ),
    lineage_info = data.frame(
      lineage = best_lineage,
      aging_pvalue = best_p,
      n_cells = length(valid_cells),
      n_young = sum(p1_age_vec[valid_cells] == "Young"),
      n_aged = sum(p1_age_vec[valid_cells] == "Aged")
    )
  )

  methods_text <- paste(
    "Pseudotime trajectory was computed on Ctrl bona fide Piezo1+ SMCs using",
    "Slingshot. Cells were re-clustered at resolution 0.4 on 15 PCs, and the",
    "start cluster was mapped to Cl15 (Core Contractile) via cross-tabulation.",
    "The lineage with the most significant aging p-value was selected.",
    "Blended pseudotime = 60% Slingshot + 40% kNN (k=15) age fraction, then",
    "rank-transformed. Five module Z-scores were smoothed by LOESS (span=0.45).",
    "Line opacity encodes local age fraction (bold = aged-dominated).",
    "Cell strip shows individual cells colored by cluster % aged."
  )

  technical_legend <- paste(
    "Top: Overlay of 5 module LOESS curves along blended pseudotime.",
    "Line opacity increases with local age fraction.",
    "Bottom: Individual cells arranged by pseudotime, colored by cluster % aged",
    "(green = young-enriched, red = aged-enriched)."
  )

  # Save overlay as primary
  save_figure_outputs(
    figure_id = fig_id,
    plot = p_final,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 12, height = 8,
    code_path = "pipeline/figure_modules/Fig4I_pseudotime_cascade.R"
  )

  # Save per-module version
  ggsave(file.path(fig_dir, paste0(fig_id, "_per_module.png")), p_per_final,
         width = 20, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(fig_id, "_per_module.pdf")), p_per_final,
         width = 20, height = 7)
  cat(sprintf("  Saved: %s_per_module.png\n", fig_id))

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    n_cells = length(valid_cells),
    lineage = best_lineage,
    aging_p = best_p
  ))
}
