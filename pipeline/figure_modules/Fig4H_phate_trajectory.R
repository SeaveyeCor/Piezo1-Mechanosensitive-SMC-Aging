# ==============================================================================
# PAPER PANEL: Fig 4G   [internal module ID: Fig4H]
# PLOTS:  PHATE embedding of CTRL visceral SMCs with the trajectory through cluster centroids (Cl15 -> Cl16 -> Cl1 -> Cl5); coloured by cluster and by % aged
# CELLS:  CTRL visceral SMCs only (Young CTRL n = 469, Aged CTRL n = 477)
# GROUPS COMPARED: none (descriptive)
# TEST:   none
# ==============================================================================
# Custom PHATE embedding of Ctrl bona fide SMCs with cubic spline trajectory
# through cluster centroids + boundary waypoints (order: Cl15→Cl16→Cl1→Cl5).
# Multi-panel: cluster-colored view + age-heatmap view.
# Requires: FNN package for kNN graph construction.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

#' Generate Fig4H - PHATE Trajectory
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Check for FNN
  if (!requireNamespace("FNN", quietly = TRUE)) {
    cat("  WARNING: FNN package not available — skipping Fig4H\n")
    return(list(success = FALSE, error = "FNN package required for PHATE"))
  }

  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "core/output_manager.R"))

  fig_id <- context$fig_id

  # Get data
  seurat_obj <- context$seurat_obj
  smc_cells <- context$selection$smc
  metadata <- seurat_obj@meta.data

  # --------------------------------------------------------------------------
  # Ctrl bona fide SMCs only
  # --------------------------------------------------------------------------
  smc_meta <- metadata[smc_cells, ]
  smc_meta$genotype <- ifelse(grepl("KO", smc_meta$sample), "KO", "Ctrl")
  ctrl_cells <- smc_cells[smc_meta$genotype == "Ctrl"]
  orig_clusters <- as.character(metadata[ctrl_cells, "seurat_clusters"])
  names(orig_clusters) <- ctrl_cells
  age_vec <- ifelse(grepl("Young", metadata[ctrl_cells, "sample"]), "Young", "Aged")
  names(age_vec) <- ctrl_cells

  cat(sprintf("  Ctrl bona fide SMCs: %d\n", length(ctrl_cells)))

  # --------------------------------------------------------------------------
  # PCA on Ctrl subset
  # --------------------------------------------------------------------------
  ctrl_smc <- subset(seurat_obj, cells = ctrl_cells)
  ctrl_smc <- Seurat::FindVariableFeatures(ctrl_smc, nfeatures = 2000, verbose = FALSE)
  ctrl_smc <- Seurat::ScaleData(ctrl_smc, verbose = FALSE)
  ctrl_smc <- Seurat::RunPCA(ctrl_smc, npcs = 20, verbose = FALSE)
  pca_mat <- Seurat::Embeddings(ctrl_smc, "pca")[, 1:15]
  rm(ctrl_smc); gc(verbose = FALSE)

  # --------------------------------------------------------------------------
  # Custom PHATE: kNN → Gaussian kernel → diffuse t=5 → log → MDS
  # --------------------------------------------------------------------------
  cat("  Computing PHATE...\n")
  run_phate_r <- function(data_mat, ndim = 2, knn = 15, t = 5) {
    nn <- FNN::get.knnx(data_mat, data_mat, k = knn + 1)
    sigma <- nn$nn.dist[, knn + 1]
    sigma[sigma < 1e-10] <- 1e-10
    n <- nrow(data_mat)
    W <- matrix(0, n, n)
    for (i in 1:n) {
      idx <- nn$nn.index[i, 2:(knn + 1)]
      dists <- nn$nn.dist[i, 2:(knn + 1)]
      W[i, idx] <- exp(-dists^2 / (sigma[i] * sigma[idx]))
    }
    W <- (W + t(W)) / 2
    rs <- rowSums(W); rs[rs < 1e-10] <- 1e-10
    P <- W / rs
    Pt <- P
    for (step in 2:t) Pt <- Pt %*% P
    Pt[Pt < 1e-10] <- 1e-10
    Pt_log <- log(Pt)
    D <- as.matrix(dist(Pt_log))
    fit <- cmdscale(D, k = ndim)
    colnames(fit) <- paste0("PHATE", 1:ndim)
    rownames(fit) <- rownames(data_mat)
    return(fit)
  }
  phate_coords <- run_phate_r(pca_mat, ndim = 2, knn = 15, t = 5)

  # --------------------------------------------------------------------------
  # Age data
  # --------------------------------------------------------------------------
  is_aged <- as.numeric(age_vec == "Aged")
  names(is_aged) <- ctrl_cells
  cluster_pct_aged <- tapply(is_aged, orig_clusters, mean)
  age_cluster_avg <- cluster_pct_aged[orig_clusters]
  names(age_cluster_avg) <- ctrl_cells

  df <- data.frame(
    PHATE1 = phate_coords[, 1], PHATE2 = phate_coords[, 2],
    cluster = orig_clusters[rownames(phate_coords)],
    pct_aged = age_cluster_avg[rownames(phate_coords)],
    stringsAsFactors = FALSE
  )

  centroids <- df %>%
    group_by(cluster) %>%
    summarise(cx = median(PHATE1), cy = median(PHATE2), n = n(), .groups = "drop")

  cl_short <- c("1" = "Cl1", "5" = "Cl5", "15" = "Cl15", "16" = "Cl16")

  # --------------------------------------------------------------------------
  # Cubic spline trajectory: Cl15 → Cl16 → Cl1 → Cl5
  # --------------------------------------------------------------------------
  cat("  Building trajectory spline...\n")
  traj_order <- c("15", "16", "1", "5")
  ctrl_pts <- data.frame(x = numeric(0), y = numeric(0))

  for (i in seq_along(traj_order)) {
    cl <- traj_order[i]
    cx <- as.numeric(centroids$cx[centroids$cluster == cl])
    cy <- as.numeric(centroids$cy[centroids$cluster == cl])

    if (i < length(traj_order)) {
      next_cl <- traj_order[i + 1]
      cl_cells <- which(df$cluster == cl)
      next_cx <- as.numeric(centroids$cx[centroids$cluster == next_cl])
      next_cy <- as.numeric(centroids$cy[centroids$cluster == next_cl])

      dists <- sqrt((df$PHATE1[cl_cells] - next_cx)^2 +
                    (df$PHATE2[cl_cells] - next_cy)^2)
      boundary_cl <- cl_cells[order(dists)[1:max(1, floor(0.2 * length(cl_cells)))]]

      next_cells <- which(df$cluster == next_cl)
      dists2 <- sqrt((df$PHATE1[next_cells] - cx)^2 +
                     (df$PHATE2[next_cells] - cy)^2)
      boundary_next <- next_cells[order(dists2)[1:max(1, floor(0.2 * length(next_cells)))]]

      boundary_all <- c(boundary_cl, boundary_next)
      bx <- median(df$PHATE1[boundary_all])
      by <- median(df$PHATE2[boundary_all])

      ctrl_pts <- rbind(ctrl_pts, data.frame(x = cx, y = cy))
      ctrl_pts <- rbind(ctrl_pts, data.frame(x = bx, y = by))
    } else {
      ctrl_pts <- rbind(ctrl_pts, data.frame(x = cx, y = cy))
    }
  }

  t_param <- seq(0, 1, length.out = nrow(ctrl_pts))
  spline_x <- splinefun(t_param, ctrl_pts$x, method = "natural")
  spline_y <- splinefun(t_param, ctrl_pts$y, method = "natural")

  t_fine <- seq(0, 1, length.out = 300)
  smooth_curve <- data.frame(x = spline_x(t_fine), y = spline_y(t_fine))
  n_sc <- nrow(smooth_curve)
  arrow_start <- max(1, n_sc - 10)

  # Trajectory overlay helper
  add_trajectory <- function(p, lw = 1.8, al = 0.55) {
    p +
      geom_path(data = smooth_curve[1:(n_sc - 6), ],
                aes(x = x, y = y), inherit.aes = FALSE,
                linewidth = lw, color = "grey15", alpha = al) +
      geom_segment(data = data.frame(
        x = smooth_curve$x[arrow_start], y = smooth_curve$y[arrow_start],
        xend = smooth_curve$x[n_sc], yend = smooth_curve$y[n_sc]),
        aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
        arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
        linewidth = lw, color = "grey15", alpha = al)
  }

  # --------------------------------------------------------------------------
  # PLOT 1: Cluster-colored
  # --------------------------------------------------------------------------
  cluster_colors <- c("15" = "#1F77B4", "1" = "#2CA02C", "16" = "#9467BD", "5" = "#D62728")
  cl_names <- c("1" = "Trans. Hub", "5" = "Synthetic",
                "15" = "Core Contr.", "16" = "Adaptive")

  centroids_lab <- centroids %>%
    mutate(label = paste0(cl_short[cluster], "\n", cl_names[cluster]))

  p_cl <- ggplot(df, aes(x = PHATE1, y = PHATE2, color = cluster)) +
    geom_point(size = 1.4, alpha = 0.6) +
    scale_color_manual(values = cluster_colors,
                       labels = c("15" = "Cl15: Core Contractile",
                                  "16" = "Cl16: Adaptive",
                                  "1" = "Cl1: Transitional Hub",
                                  "5" = "Cl5: Synthetic"),
                       name = "Cluster") +
    geom_label(data = centroids_lab, aes(x = cx, y = cy, label = label),
               color = "black", size = 3.8, fontface = "bold",
               fill = "white", alpha = 0.8,
               lineheight = 0.85, inherit.aes = FALSE) +
    labs(subtitle = sprintf("n=%d Ctrl bona fide SMCs", nrow(df)),
         x = "PHATE 1", y = "PHATE 2") +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(size = 10.5, color = "black"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  p_cl <- add_trajectory(p_cl, lw = 2, al = 0.6)

  # --------------------------------------------------------------------------
  # PLOT 2: Age heatmap
  # --------------------------------------------------------------------------
  df_plot <- df[order(df$pct_aged), ]  # Plot low-aged first

  centroids_age <- centroids %>%
    mutate(
      label_name = cl_short[cluster],
      label_stats = sprintf("%.0f%% Aged\n(n=%d)", 100 * cluster_pct_aged[cluster], n)
    )

  p_age <- ggplot(df_plot, aes(x = PHATE1, y = PHATE2, color = pct_aged)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_gradient(low = "#4DAF4A", high = "#E41A1C",
                         name = "% Aged\n(cluster avg)",
                         labels = function(x) paste0(round(100 * x), "%")) +
    geom_text(data = centroids_age, aes(x = cx, y = cy, label = label_name),
              color = "black", size = 5, fontface = "bold",
              vjust = -1.2, inherit.aes = FALSE) +
    geom_text(data = centroids_age, aes(x = cx, y = cy, label = label_stats),
              color = "black", size = 3.5, vjust = 0.6, lineheight = 0.9,
              inherit.aes = FALSE) +
    labs(subtitle = "Green = young-enriched, Red = aged-enriched",
         x = "PHATE 1", y = "PHATE 2") +
    theme_minimal(base_size = 14) +
    theme(
      plot.subtitle = element_text(size = 10.5, color = "black"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
  p_age <- add_trajectory(p_age)

  # --------------------------------------------------------------------------
  # Save outputs
  # --------------------------------------------------------------------------
  fig_dir <- file.path(context$output_dir, fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  export_data <- list(
    phate_coordinates = data.frame(
      cell = rownames(phate_coords),
      PHATE1 = phate_coords[, 1],
      PHATE2 = phate_coords[, 2],
      cluster = orig_clusters[rownames(phate_coords)],
      age = age_vec[rownames(phate_coords)],
      pct_aged_cluster = age_cluster_avg[rownames(phate_coords)]
    ),
    cluster_summary = as.data.frame(centroids) %>%
      mutate(pct_aged = as.numeric(cluster_pct_aged[cluster])),
    spline_control_points = ctrl_pts
  )

  methods_text <- paste(
    "PHATE dimensionality reduction was performed on Ctrl bona fide SMCs",
    "(clusters 1, 5, 15, 16) using a custom implementation: PCA (2000 variable",
    "features, 15 PCs), kNN graph (k=15) with Gaussian kernel, diffusion (t=5",
    "steps), log-transform, then classical MDS to 2D. A cubic spline trajectory",
    "was fitted through cluster centroids and boundary waypoints in the order",
    "Cl15 (Core Contractile) -> Cl16 (Adaptive) -> Cl1 (Transitional Hub) ->",
    "Cl5 (Synthetic), representing the contractile-to-synthetic transition."
  )

  technical_legend <- paste(
    "PHATE embedding of Ctrl bona fide SMCs with trajectory spline.",
    "Left: cells colored by cluster identity. Right: cells colored by",
    "cluster-average % aged (green = young-enriched, red = aged-enriched).",
    "Arrow indicates contractile-to-synthetic direction."
  )

  # Save cluster view as primary
  save_multi_panel_figure(
    figure_id = fig_id,
    plots = list("_clusters" = p_cl, "_age" = p_age),
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 11, height = 8,
    code_path = "pipeline/figure_modules/Fig4H_phate_trajectory.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    n_cells = nrow(df),
    cluster_summary = as.data.frame(centroids)
  ))
}
