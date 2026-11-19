# ==============================================================================
# SMC ANALYSIS PIPELINE - MAIN RUNNER
# ==============================================================================
# Generates all figures for a given SMC selection criteria
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(FNN)
})

# ==============================================================================
# CONFIGURATION
# ==============================================================================

#' Get pipeline configuration
#'
#' Paths are resolved relative to the project root (working directory).
#' Set working directory to the project root before running the pipeline.
#'
#' @return List with pipeline settings
get_config <- function() {
  # Auto-detect project root from working directory
  project_root <- normalizePath(getwd(), winslash = "/")

  list(
    # Pipeline root directory (auto-detected)
    pipeline_root = file.path(project_root, "pipeline"),

    # Data path - deposited Seurat object (six mice)
    seurat_path = file.path(project_root, "data",
                            "seurat_paper_6samples.rds"),

    # Output settings
    base_output_dir = file.path(project_root, "pipeline_results"),

    # Figure dimensions (in inches)
    figure_sizes = list(
      standard = list(width = 8, height = 6),
      square = list(width = 6, height = 6),
      umap = list(width = 7, height = 6),
      heatmap = list(width = 10, height = 8),
      barplot = list(width = 6, height = 5)
    ),

    # Analysis parameters
    phenotype_ratio_threshold = 1.5,

    # Groups
    group_order = c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"),

    # Optional sample filter (none by default)
    default_exclude_samples = character(0)
  )
}

# ==============================================================================
# UMAP COMPUTATION
# ==============================================================================

#' Ensure UMAP coordinates exist in Seurat object
#'
#' Checks if UMAP reduction exists, and if not, computes it.
#' Will also run PCA if needed (prerequisite for UMAP).
#'
#' @param seurat_obj Seurat object
#' @param dims Number of PCA dimensions to use for UMAP (default: 30)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Seurat object with UMAP coordinates
ensure_umap <- function(seurat_obj, dims = 30, verbose = TRUE) {

  # Check if UMAP already exists
  if ("umap" %in% names(seurat_obj@reductions)) {
    if (verbose) cat("  UMAP coordinates found in Seurat object\n")
    return(seurat_obj)
  }

  if (verbose) cat("  UMAP not found - computing UMAP coordinates...\n")

  # Check if PCA exists and covers all cells
  pca_ok <- "pca" %in% names(seurat_obj@reductions) &&
    nrow(seurat_obj@reductions$pca@cell.embeddings) == ncol(seurat_obj)
  if (!pca_ok) {
    if (verbose) {
      if ("pca" %in% names(seurat_obj@reductions)) {
        n_pca <- nrow(seurat_obj@reductions$pca@cell.embeddings)
        cat(sprintf("    PCA exists for %d/%d cells - recomputing...\n", n_pca, ncol(seurat_obj)))
      } else {
        cat("    Running PCA (prerequisite for UMAP)...\n")
      }
    }

    # Check if variable features exist
    if (length(VariableFeatures(seurat_obj)) == 0) {
      if (verbose) cat("    Finding variable features...\n")
      seurat_obj <- FindVariableFeatures(seurat_obj,
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = FALSE)
    }

    # Check if data is scaled for ALL cells (not just a subset)
    scaled_data <- GetAssayData(seurat_obj, layer = "scale.data")
    n_scaled <- ncol(scaled_data)
    n_cells <- ncol(seurat_obj)
    if (length(scaled_data) == 0 || all(dim(scaled_data) == c(0, 0)) || n_scaled < n_cells) {
      if (verbose) {
        if (n_scaled > 0 && n_scaled < n_cells) {
          cat(sprintf("    Scale data exists for %d/%d cells - rescaling all...\n", n_scaled, n_cells))
        } else {
          cat("    Scaling data...\n")
        }
      }
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    }

    # Run PCA
    seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
    if (verbose) cat("    PCA complete (50 PCs computed)\n")
  } else {
    if (verbose) cat("    PCA found in Seurat object\n")
  }

  # Determine optimal number of dims (use available or requested, whichever is smaller)
  n_pcs <- ncol(seurat_obj@reductions$pca@cell.embeddings)
  dims_to_use <- min(dims, n_pcs)

  if (verbose) cat(sprintf("    Running UMAP with %d PCs...\n", dims_to_use))

  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj,
                         dims = 1:dims_to_use,
                         verbose = FALSE,
                         seed.use = 42)  # For reproducibility

  if (verbose) cat("    UMAP complete!\n")

  return(seurat_obj)
}

# ==============================================================================
# MAIN PIPELINE FUNCTION
# ==============================================================================

#' Run the full SMC analysis pipeline
#'
#' @param selection_module Path to selection criteria R file
#' @param figures Character vector of figure IDs to generate (default: all)
#' @param seurat_obj Optional: pre-loaded Seurat object (skips loading)
#' @param output_dir Optional: custom output directory
#' @param remove_epithelial Logical: remove epithelial-contaminated cells (default: FALSE)
#' @param exclude_samples Character vector of sample IDs to exclude (default: NULL)
#' @param depth_match Depth-matching filter. TRUE uses default (range_factor=2.0),
#'   or pass list(range_factor=N, target_samples=c(...)) to customize.
#'   Computes median-of-medians UMI across samples, then filters ONLY outlier
#'   samples (those whose median falls outside the range) to the matched window.
#'   Non-outlier samples are kept intact. Applied after sample exclusion +
#'   epithelial removal, before cell selection.
#'
#' @return List with run information and generated figure paths
run_pipeline <- function(
    selection_module,
    figures = "all",
    seurat_obj = NULL,
    output_dir = NULL,
    remove_epithelial = FALSE,
    exclude_samples = NULL,
    depth_match = NULL
) {
  config <- get_config()
  start_time <- Sys.time()

  cat("\n")
  cat("================================================================\n")
  cat("       SMC ANALYSIS PIPELINE\n")
  cat("================================================================\n")
  cat(sprintf("Started: %s\n", start_time))
  cat("\n")

  # --------------------------------------------------------------------------
  # 1. Load selection module
  # --------------------------------------------------------------------------
  cat("Loading selection criteria module...\n")
  # If relative path, make it relative to pipeline_root
  if (!file.exists(selection_module)) {
    selection_module <- file.path(config$pipeline_root, selection_module)
  }
  source(selection_module)

  method_info <- get_method_info()
  cat(sprintf("  Method: %s\n", method_info$method_name))
  cat(sprintf("  Description: %s\n", method_info$description))
  cat("\n")

  # --------------------------------------------------------------------------
  # 2. Load data
  # --------------------------------------------------------------------------
  if (is.null(seurat_obj)) {
    cat("Loading Seurat object...\n")
    seurat_obj <- readRDS(config$seurat_path)

    # Handle Seurat v5 assay format
    if (inherits(seurat_obj[["RNA"]], "Assay5")) {
      cat("  Joining Seurat v5 layers...\n")
      seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
    }

    cat(sprintf("  Loaded %d cells\n", ncol(seurat_obj)))
  } else {
    cat("Using provided Seurat object...\n")
    cat(sprintf("  %d cells\n", ncol(seurat_obj)))
  }

  # Cells with sample=NA take their sample from orig.ident (none in the deposited object)
  na_sample <- is.na(seurat_obj$sample)
  if (any(na_sample)) {
    seurat_obj$sample <- as.character(seurat_obj$sample)
    seurat_obj$sample[na_sample] <- as.character(seurat_obj$orig.ident[na_sample])
    cat(sprintf("  sample=NA for %d cells: taken from orig.ident\n", sum(na_sample)))
  }
  cat("\n")

  # Per-sample cell count and median UMI before any filtering (Table S0)
  raw_sample_stats <- data.frame(
    sample = sort(unique(as.character(seurat_obj$sample))),
    stringsAsFactors = FALSE)
  raw_sample_stats$n_cells <- sapply(raw_sample_stats$sample, function(s) sum(seurat_obj$sample == s))
  raw_sample_stats$median_umi <- sapply(raw_sample_stats$sample, function(s) median(seurat_obj$nCount_RNA[seurat_obj$sample == s]))

  # --------------------------------------------------------------------------
  # 2a. Ensure UMAP exists (compute if missing)
  # --------------------------------------------------------------------------
  seurat_obj <- ensure_umap(seurat_obj)

  # --------------------------------------------------------------------------
  # 2c. Exclude specific samples (uses config default if not specified)
  # --------------------------------------------------------------------------
  if (is.null(exclude_samples)) {
    exclude_samples <- config$default_exclude_samples
  }
  if (!is.null(exclude_samples) && length(exclude_samples) > 0) {
    cat("Excluding samples:", paste(exclude_samples, collapse = ", "), "\n")
    cells_before <- ncol(seurat_obj)
    cells_to_keep <- !seurat_obj@meta.data$orig.ident %in% exclude_samples
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    cells_excluded <- cells_before - ncol(seurat_obj)
    cat(sprintf("  Excluded %d cells from %d sample(s)\n", cells_excluded, length(exclude_samples)))
    cat(sprintf("  Remaining cells: %d\n", ncol(seurat_obj)))
    cat("\n")
  }

  # --------------------------------------------------------------------------
  # 2d. Remove epithelial-contaminated cells (optional)
  # --------------------------------------------------------------------------
  epithelial_removed <- 0
  goblet_removed <- 0
  celltype_removed <- 0
  if (remove_epithelial) {
    cat("Removing epithelial/goblet contaminated cells...\n")

    # Step 1: Remove cells annotated as non-muscularis cell types
    # Goblet/Epithelial = intestinal, Endothelial = vascular - not muscularis
    if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
      exclude_celltypes <- c("Goblet", "Epithelial", "Endothelial", "Endothelial_Lymphatic")
      cells_to_remove_by_type <- seurat_obj@meta.data$cell_type %in% exclude_celltypes
      celltype_removed <- sum(cells_to_remove_by_type)

      if (celltype_removed > 0) {
        # Count by type for reporting
        goblet_removed <- sum(seurat_obj@meta.data$cell_type == "Goblet")
        epithelial_type_removed <- sum(seurat_obj@meta.data$cell_type == "Epithelial")
        endothelial_removed <- sum(seurat_obj@meta.data$cell_type == "Endothelial")
        lymphatic_removed <- sum(seurat_obj@meta.data$cell_type == "Endothelial_Lymphatic")

        cat(sprintf("  Goblet cells removed: %d\n", goblet_removed))
        cat(sprintf("  Epithelial cells removed: %d\n", epithelial_type_removed))
        cat(sprintf("  Endothelial cells removed: %d\n", endothelial_removed))
        cat(sprintf("  Lymphatic Endothelial removed: %d\n", lymphatic_removed))

        cells_to_keep <- colnames(seurat_obj)[!cells_to_remove_by_type]
        seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
        cat(sprintf("  Cells after cell_type removal: %d\n", ncol(seurat_obj)))
      }
    }

    # Non-muscularis cells are removed by cell_type annotation.

    cat(sprintf("  TOTAL removed by cell_type: %d\n", celltype_removed))
    cat("  (Marker-based filtering disabled - ambient RNA contamination issue)\n")
    cat("\n")
  }

  # --------------------------------------------------------------------------
  # 2e. Depth-matched filtering (optional)
  # --------------------------------------------------------------------------
  depth_match_info <- NULL
  if (!is.null(depth_match)) {
    # Parse depth_match argument
    if (isTRUE(depth_match)) {
      range_factor <- 2.0
    } else if (is.list(depth_match)) {
      range_factor <- depth_match$range_factor %||% 2.0
    } else {
      stop("depth_match must be TRUE or a list with range_factor")
    }

    cat("Applying depth-matched filtering...\n")

    # Compute per-sample median UMI
    sample_col <- if ("sample" %in% colnames(seurat_obj@meta.data)) "sample" else "orig.ident"
    samples <- unique(seurat_obj@meta.data[[sample_col]])
    sample_medians <- sapply(samples, function(s) {
      median(seurat_obj$nCount_RNA[seurat_obj@meta.data[[sample_col]] == s])
    })
    names(sample_medians) <- samples

    # Target = median of per-sample medians (robust to outliers)
    target_depth <- median(sample_medians)
    depth_lo <- target_depth / range_factor
    depth_hi <- target_depth * range_factor

    # Identify outlier samples (median outside the depth window)
    outlier_samples <- names(sample_medians)[sample_medians < depth_lo | sample_medians > depth_hi]

    cat(sprintf("  Per-sample median UMIs:\n"))
    for (s in names(sort(sample_medians))) {
      flag <- if (s %in% outlier_samples) " <-- OUTLIER (will be filtered)" else ""
      cat(sprintf("    %-15s %6.0f%s\n", s, sample_medians[s], flag))
    }
    cat(sprintf("  Target depth (median of medians): %.0f\n", target_depth))
    cat(sprintf("  Depth window: [%.0f, %.0f] (range_factor=%.1f)\n",
                depth_lo, depth_hi, range_factor))
    cat(sprintf("  Outlier samples: %s\n",
                if (length(outlier_samples) > 0) paste(outlier_samples, collapse = ", ") else "none"))

    # Filter ONLY outlier samples to the depth window; keep all others intact
    cells_before <- ncol(seurat_obj)
    is_outlier_cell <- seurat_obj@meta.data[[sample_col]] %in% outlier_samples
    is_in_range <- seurat_obj$nCount_RNA >= depth_lo & seurat_obj$nCount_RNA <= depth_hi
    keep_mask <- !is_outlier_cell | (is_outlier_cell & is_in_range)

    # Report per-sample filtering
    cat(sprintf("  Per-sample retention:\n"))
    for (s in names(sort(sample_medians, decreasing = TRUE))) {
      s_mask <- seurat_obj@meta.data[[sample_col]] == s
      n_before <- sum(s_mask)
      n_after <- sum(s_mask & keep_mask)
      filtered <- if (s %in% outlier_samples) " [FILTERED]" else ""
      cat(sprintf("    %-15s %5d → %5d (%4.1f%%)%s\n",
                  s, n_before, n_after, 100 * n_after / max(n_before, 1), filtered))
    }

    seurat_obj_before_depth <- seurat_obj@meta.data[, sample_col, drop = FALSE]; colnames(seurat_obj_before_depth) <- "sample"
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep_mask])
    cells_removed <- cells_before - ncol(seurat_obj)
    cat(sprintf("  Total: %d → %d cells (removed %d from outlier samples)\n",
                cells_before, ncol(seurat_obj), cells_removed))
    cat("\n")

    depth_match_info <- list(
      target_depth = target_depth,
      range = c(depth_lo, depth_hi),
      range_factor = range_factor,
      cells_removed = cells_removed,
      sample_medians = sample_medians,
      outlier_samples = outlier_samples,
      n_before = sapply(outlier_samples, function(s) sum(is_outlier_cell & seurat_obj_before_depth$sample == s)),
      n_after = sapply(outlier_samples, function(s) sum(seurat_obj@meta.data[[sample_col]] == s))
    )
  }

  # --------------------------------------------------------------------------
  # 2f. KNN cluster assignment for cells with seurat_clusters=NA
  # --------------------------------------------------------------------------
  # Cells with seurat_clusters=NA (none in the deposited object) are assigned by kNN.
  # Use KNN (k=30) in PCA space to assign clusters from reference cells.
  na_cluster_mask <- is.na(seurat_obj$seurat_clusters)
  if (any(na_cluster_mask)) {
    n_na <- sum(na_cluster_mask)
    cat(sprintf("KNN cluster assignment: %d cells have seurat_clusters=NA\n", n_na))

    # Ensure PCA covers all cells
    pca_ok <- "pca" %in% names(seurat_obj@reductions) &&
      nrow(seurat_obj@reductions$pca@cell.embeddings) == ncol(seurat_obj)
    if (!pca_ok) {
      cat("  Computing PCA for KNN cluster assignment...\n")
      if (length(VariableFeatures(seurat_obj)) == 0) {
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                            nfeatures = 2000, verbose = FALSE)
      }
      scaled_data <- GetAssayData(seurat_obj, layer = "scale.data")
      if (length(scaled_data) == 0 || all(dim(scaled_data) == c(0, 0)) ||
          ncol(scaled_data) < ncol(seurat_obj)) {
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
      }
      seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
      cat("  PCA complete (50 PCs)\n")
    }

    pca_embeddings <- seurat_obj@reductions$pca@cell.embeddings
    ref_idx <- which(!na_cluster_mask)
    query_idx <- which(na_cluster_mask)

    ref_pca <- pca_embeddings[ref_idx, , drop = FALSE]
    query_pca <- pca_embeddings[query_idx, , drop = FALSE]
    ref_clusters <- seurat_obj$seurat_clusters[ref_idx]

    # KNN with k=30
    k <- 30
    knn_result <- FNN::get.knnx(data = ref_pca, query = query_pca, k = k)

    # Majority vote for each query cell
    assigned_clusters <- apply(knn_result$nn.index, 1, function(nn_ids) {
      neighbor_clusters <- ref_clusters[nn_ids]
      tbl <- table(neighbor_clusters)
      names(tbl)[which.max(tbl)]
    })

    # Assign clusters (preserve factor levels from existing clusters)
    existing_levels <- levels(seurat_obj$seurat_clusters)
    seurat_obj$seurat_clusters[query_idx] <- assigned_clusters
    if (!is.null(existing_levels)) {
      seurat_obj$seurat_clusters <- factor(seurat_obj$seurat_clusters, levels = existing_levels)
    }

    # Log cluster distribution for assigned cells
    assigned_tbl <- table(assigned_clusters)
    cat(sprintf("  Assigned %d cells to clusters via KNN (k=%d):\n", n_na, k))
    for (cl in names(sort(assigned_tbl, decreasing = TRUE))) {
      cat(sprintf("    Cluster %s: %d cells\n", cl, assigned_tbl[cl]))
    }

    # Report which samples were affected
    affected_samples <- unique(seurat_obj$orig.ident[query_idx])
    cat(sprintf("  Affected samples: %s\n", paste(affected_samples, collapse = ", ")))
    cat("\n")
  }

  # --------------------------------------------------------------------------
  # 3. Select cells
  # --------------------------------------------------------------------------
  cat("Selecting cells...\n")
  selection <- select_cells(seurat_obj)

  # --------------------------------------------------------------------------
  # 4. Create output directory
  # --------------------------------------------------------------------------
  if (is.null(output_dir)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    # Add suffix if epithelial removal was used
    method_suffix <- if (remove_epithelial) "_noEpithelial" else ""
    output_dir <- file.path(
      config$base_output_dir,
      paste0(method_info$method_name, method_suffix, "_", timestamp)
    )
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("Output directory: %s\n\n", output_dir))

  # --------------------------------------------------------------------------
  # 5. Classify SMC phenotypes
  # --------------------------------------------------------------------------
  cat("Classifying SMC phenotypes...\n")
  # Source gene_sets first so phenotype_classifier has access
  source(file.path(config$pipeline_root, "core/gene_sets.R"))
  source(file.path(config$pipeline_root, "core/phenotype_classifier.R"))
  phenotype_df <- classify_smc_phenotype(
    seurat_obj,
    selection$smc,
    ratio_threshold = config$phenotype_ratio_threshold
  )
  cat("\n")

  # --------------------------------------------------------------------------
  # 6. Prepare analysis context
  # --------------------------------------------------------------------------
  context <- list(
    seurat_obj = seurat_obj,
    raw_sample_stats = raw_sample_stats,
    depth_match_info = depth_match_info,
    selection = selection,
    phenotype = phenotype_df,
    method_info = method_info,
    config = config,
    output_dir = output_dir
  )

  # --------------------------------------------------------------------------
  # 7. Determine which figures to generate
  # --------------------------------------------------------------------------
  all_figures <- c(
    # Main figures (3E onwards - 3A-D are physiology figures in manuscript)
    "Fig3E", "Fig3F", "Fig3G", "Fig3H", "Fig3I", "Fig3J", "Fig3K",
    # Main figures (4E onwards - 4A-D are physiology figures in manuscript)
    "Fig4E", "Fig4F", "Fig4G", "Fig4H", "Fig4I", "Fig4J",
    # Supplementary panels (S3A, S3C/D, S3F, S3J, S4A in the manuscript) and score panels 3M, 3N
    "FigS3D", "FigS3F", "FigS3G", "FigS3H", "FigS3I", "FigS3K", "FigS4A", "FigS4H",
    # Supplemental Tables
    "SupplementalTables"
  )

  if (identical(figures, "all")) {
    figures_to_run <- all_figures
  } else {
    figures_to_run <- figures[figures %in% all_figures]
  }

  cat(sprintf("Generating %d figures...\n\n", length(figures_to_run)))

  # --------------------------------------------------------------------------
  # 8. Generate figures
  # --------------------------------------------------------------------------
  figure_results <- list()

  for (fig_id in figures_to_run) {
    cat(sprintf("--- %s ---\n", fig_id))

    result <- tryCatch({
      generate_figure(fig_id, context)
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      list(success = FALSE, error = e$message)
    })

    figure_results[[fig_id]] <- result
    cat("\n")
  }

  # --------------------------------------------------------------------------
  # 9. Collect all PNGs into single folder
  # --------------------------------------------------------------------------
  collect_pngs(output_dir)

  # --------------------------------------------------------------------------
  # 10. Generate run summary
  # --------------------------------------------------------------------------
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")

  # Create README
  generate_run_summary(
    output_dir = output_dir,
    method_info = method_info,
    selection = selection,
    phenotype = phenotype_df,
    figure_results = figure_results,
    duration = duration,
    remove_epithelial = remove_epithelial,
    epithelial_removed = epithelial_removed,
    goblet_removed = goblet_removed
  )

  cat("================================================================\n")
  cat("       PIPELINE COMPLETE\n")
  cat("================================================================\n")
  cat(sprintf("Duration: %.1f minutes\n", as.numeric(duration)))
  cat(sprintf("Output: %s\n", output_dir))
  cat(sprintf("Figures generated: %d/%d\n",
              sum(sapply(figure_results, function(x) isTRUE(x$success))),
              length(figures_to_run)))
  cat("================================================================\n\n")

  return(list(
    output_dir = output_dir,
    method_info = method_info,
    selection = selection,
    phenotype = phenotype_df,
    figure_results = figure_results,
    duration = duration
  ))
}

# ==============================================================================
# FIGURE GENERATION DISPATCHER
# ==============================================================================

#' Generate a specific figure
#'
#' @param fig_id Figure identifier
#' @param context Analysis context from run_pipeline
#'
#' @return List with success status and output paths
generate_figure <- function(fig_id, context) {
  # Get pipeline root directory - use config or derive from working directory
  pipeline_root <- context$config$pipeline_root

  # Map figure IDs to modules
  # NOTE: Fig3A-D and Fig4A-D are physiology figures in manuscript
  #       scRNA-seq figures start at Fig3E
  #       S3A-C are organotypic culture supplementary (not pipeline-generated)
  #       S3D-K are primary scRNA-seq panels (this pipeline)
  #       S4A-G are additional scRNA-seq supplementary (this pipeline)
  #       S5 = human data (not pipeline-generated)
  module_map <- list(
    # Main scRNA-seq figures (3E-K)
    "Fig3E" = "figure_modules/Fig3E_umap_phenotype.R",
    "Fig3F" = "figure_modules/Fig3F_piezo1_contractile_correlation.R",
    "Fig3G" = "figure_modules/Fig3G_piezo1_synthetic_correlation.R",
    "Fig3H" = "figure_modules/Fig3H_contractile_capacity_score.R",
    "Fig3I" = "figure_modules/Fig3I_contractile_rescue.R",
    "Fig3J" = "figure_modules/Fig3J_synthetic_rescue.R",
    "Fig3K" = "figure_modules/Fig3K_phenotype_proportions.R",
    # Main NFAT figures (4E-G)
    "Fig4E" = "figure_modules/Fig4E_nfatc3_remodeling.R",
    "Fig4F" = "figure_modules/Fig4F_nfat_piezo1_correlation.R",
    "Fig4G" = "figure_modules/Fig4G_nfat_synthetic_correlation.R",
    "Fig4H" = "figure_modules/Fig4H_phate_trajectory.R",
    "Fig4I" = "figure_modules/Fig4I_pseudotime_cascade.R",
    "Fig4J" = "figure_modules/Fig4J_p1pos_vs_p1neg_synthetic.R",
    # Supplementary modules (internal IDs; manuscript panels in the README)
    "FigS3D" = "figure_modules/supplementary/FigS3D_umap_celltypes.R",
    "FigS3F" = "figure_modules/supplementary/FigS3F_piezo1_log2fc.R",
    "FigS3G" = "figure_modules/supplementary/FigS3G_pathway_correlation.R",
    "FigS3H" = "figure_modules/supplementary/FigS3H_gene_correlations.R",
    "FigS3I" = "figure_modules/supplementary/FigS3I_go_enrichment.R",
    "FigS3K" = "figure_modules/supplementary/FigS3K_ecm_organization.R",
    "FigS4A" = "figure_modules/supplementary/FigS3L_ecm_heatmap.R",
    "FigS4H" = "figure_modules/supplementary/FigS4H_cyto_remodeling.R",
    # Supplemental Tables
    "SupplementalTables" = "figure_modules/supplementary/SupplementalTables.R"
  )

  if (!fig_id %in% names(module_map)) {
    return(list(success = FALSE, error = "Unknown figure ID"))
  }

  module_path <- file.path(pipeline_root, module_map[[fig_id]])

  if (!file.exists(module_path)) {
    cat(sprintf("  Module not yet implemented: %s\n", basename(module_path)))
    return(list(success = FALSE, error = "Module not implemented"))
  }

  # Pass fig_id through context so modules use it for output naming
  context$fig_id <- fig_id

  # Source and run the module
  source(module_path)
  result <- generate(context)

  return(result)
}

# ==============================================================================
# SUMMARY GENERATION
# ==============================================================================

#' Generate run summary README
#'
#' @param output_dir Output directory
#' @param method_info Method information
#' @param selection Cell selection results
#' @param phenotype Phenotype classification results
#' @param figure_results Results from figure generation
#' @param duration Run duration
#' @param remove_epithelial Whether epithelial removal was applied
#' @param epithelial_removed Number of epithelial cells removed (by marker)
#' @param goblet_removed Number of goblet cells removed (by cell_type)
generate_run_summary <- function(
    output_dir,
    method_info,
    selection,
    phenotype,
    figure_results,
    duration,
    remove_epithelial = FALSE,
    epithelial_removed = 0,
    goblet_removed = 0
) {
  lines <- c(
    "# SMC Analysis Pipeline Run Summary",
    "",
    sprintf("**Generated:** %s", Sys.time()),
    sprintf("**Duration:** %.1f minutes", as.numeric(duration)),
    ""
  )

  # Add epithelial removal info if applicable
  if (remove_epithelial) {
    lines <- c(lines,
      "## Non-Muscularis Cell Removal",
      "**Status:** ENABLED",
      sprintf("**Goblet cells removed:** %d", goblet_removed),
      sprintf("**Epithelial cells removed (by markers):** %d", epithelial_removed),
      "**Endothelial cells removed:** (see log)",
      "**Cell types excluded:** Goblet, Epithelial, Endothelial, Endothelial_Lymphatic",
      "**Markers used:** Reg3b, Reg3g, Lypd8, Epcam, Krt8, Krt19",
      "**Rationale:** Focus on muscularis-resident cells only",
      ""
    )
  }

  lines <- c(lines,
    "## Selection Method",
    sprintf("**Name:** %s", method_info$method_name),
    sprintf("**Description:** %s", method_info$description),
    "",
    "## Cell Selection Summary",
    sprintf("- **SMC cells:** %d", length(selection$smc)),
    sprintf("- **PDGFRa cells:** %d", length(selection$pdgfra)),
    sprintf("- **ICC cells:** %d", length(selection$icc)),
    "",
    "## Phenotype Classification",
    sprintf("- **Contractile SMCs:** %d (%.1f%%)",
            sum(phenotype$phenotype == "Contractile"),
            100 * mean(phenotype$phenotype == "Contractile")),
    sprintf("- **Synthetic SMCs:** %d (%.1f%%)",
            sum(phenotype$phenotype == "Synthetic"),
            100 * mean(phenotype$phenotype == "Synthetic")),
    "",
    "## Figures Generated",
    ""
  )

  # Add figure results
  for (fig_id in names(figure_results)) {
    result <- figure_results[[fig_id]]
    if (isTRUE(result$success)) {
      lines <- c(lines, sprintf("- [x] %s", fig_id))
    } else {
      lines <- c(lines, sprintf("- [ ] %s (failed: %s)", fig_id,
                                result$error %||% "unknown error"))
    }
  }

  lines <- c(lines,
    "",
    "## Directory Contents",
    "Each figure directory contains:",
    "- `FigXX.png` - Publication-quality figure (300 DPI)",
    "- `FigXX.pdf` - Vector format",
    "- `FigXX_data.xlsx` - Source data",
    "- `FigXX_methods.txt` - Methods description and technical legend",
    ""
  )

  writeLines(lines, file.path(output_dir, "README.md"))
}

# ==============================================================================
# PNG COLLECTION
# ==============================================================================

#' Collect all PNGs into a single folder for easy access
#'
#' @param output_dir The pipeline output directory
#' @return Invisible NULL
collect_pngs <- function(output_dir) {
  cat("Collecting PNGs into single folder...\n")

  # Create pngs directory
  pngs_dir <- file.path(output_dir, "pngs")
  dir.create(pngs_dir, showWarnings = FALSE, recursive = TRUE)

  # Find all PNG files in figure subdirectories
  png_files <- list.files(output_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)

  # Exclude any already in pngs folder
  png_files <- png_files[!grepl("/pngs/", png_files)]

  if (length(png_files) == 0) {
    cat("  No PNG files found to collect\n")
    return(invisible(NULL))
  }

  copied <- 0
  for (png_path in png_files) {
    # Get the filename
    filename <- basename(png_path)

    # Get the parent directory name (figure ID)
    parent_dir <- basename(dirname(png_path))

    # If filename doesn't start with Fig, prefix with parent directory
    # This handles cases like "trajectory.png" -> "FigS5_trajectory.png"
    if (!grepl("^Fig", filename)) {
      new_filename <- paste0(parent_dir, "_", filename)
    } else {
      new_filename <- filename
    }

    # Copy to pngs folder
    dest_path <- file.path(pngs_dir, new_filename)
    file.copy(png_path, dest_path, overwrite = TRUE)
    copied <- copied + 1
  }

  cat(sprintf("  Copied %d PNG files to %s\n", copied, pngs_dir))
  return(invisible(NULL))
}

# ==============================================================================
# QUICK RUN FUNCTIONS
# ==============================================================================

#' Run pipeline with cluster-based selection
#'
#' @param figures Figures to generate (default: all)
#' @return Pipeline results
run_cluster_based <- function(figures = "all") {
  pipeline_root <- dirname(sys.frame(1)$ofile)
  selection_module <- file.path(
    dirname(pipeline_root),
    "selection_criteria/cluster_based_smc.R"
  )
  run_pipeline(selection_module, figures = figures)
}

#' Run pipeline with Myh11+ selection
#'
#' @param figures Figures to generate (default: all)
#' @return Pipeline results
run_myh11_positive <- function(figures = "all") {
  pipeline_root <- dirname(sys.frame(1)$ofile)
  selection_module <- file.path(
    dirname(pipeline_root),
    "selection_criteria/myh11_positive.R"
  )
  run_pipeline(selection_module, figures = figures)
}
