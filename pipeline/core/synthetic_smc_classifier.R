# ==============================================================================
# SYNTHETIC SMC CLUSTER CLASSIFIER
# ==============================================================================
# Purpose: Identify high-confidence synthetic SMC clusters in addition to
#          high-confidence contractile SMCs, combining both for analysis
#
# Based on literature:
#   - Yap et al. 2021 "Six Shades of vSMCs" - 6 vSMC phenotypes
#   - Schneider et al. 2023 - Human colon snRNA-seq SIP syncytium
#   - Wolfson et al. 2025 - Mechanical stress in visceral SMCs
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ==============================================================================
# SYNTHETIC SMC MARKER GENE SETS
# ==============================================================================

#' Get synthetic SMC marker gene sets from literature
#' @return Named list of gene vectors for synthetic SMC subtypes
get_synthetic_smc_markers <- function() {
  list(
    # =========================================================================
    # CORE CONTRACTILE (for reference/comparison)
    # =========================================================================
    contractile_core = c("Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn"),

    # =========================================================================
    # SYNTHETIC/DEDIFFERENTIATION MASTER REGULATORS
    # =========================================================================
    # KLF4 is THE master TF driving SMC dedifferentiation (Yap et al. 2021)
    dedifferentiation_tfs = c(
      "Klf4",     # Master regulator of SMC dedifferentiation
      "Klf5",     # Proliferation/synthetic phenotype
      "Srf"       # Serum response factor (context-dependent)
    ),

    # =========================================================================
    # MESENCHYMAL-LIKE SMC MARKERS (Yap et al. 2021)
    # =========================================================================
    # KLF4-driven transition to mesenchymal phenotype
    # CD34+/CD44+ marks dedifferentiated/stem-like SMCs
    mesenchymal_like = c(
      "Klf4",     # Master TF
      "Cd34",     # Stem/progenitor marker
      "Cd44",     # Hyaluronan receptor, mesenchymal marker
      "Ly6a",     # Sca1 - stem cell antigen
      "Vim",      # Vimentin - mesenchymal marker
      "S100a4"    # FSP1 - fibroblast-specific protein 1
    ),

    # =========================================================================
    # FIBROBLAST-LIKE SMC MARKERS (Yap et al. 2021)
    # =========================================================================
    # ECM-producing, fibrotic phenotype
    fibroblast_like = c(
      "Fn1",      # Fibronectin 1
      "Bgn",      # Biglycan
      "Dcn",      # Decorin
      "Col1a1",   # Collagen type I alpha 1
      "Col1a2",   # Collagen type I alpha 2
      "Col3a1",   # Collagen type III alpha 1
      "Lum",      # Lumican
      "Postn",    # Periostin
      "Tnc"       # Tenascin C
    ),

    # =========================================================================
    # MECHANICAL STRESS-INDUCED MARKERS (Wolfson et al. 2025)
    # =========================================================================
    # Rapid gene expression changes in visceral SMCs under stress
    mechanical_stress = c(
      "Ereg",     # Epiregulin - EGFR ligand, immediate-early
      "Areg",     # Amphiregulin - EGFR ligand
      "Klf4",     # Also induced by mechanical stress
      "Pdgfa",    # PDGF-A - growth factor
      "Thbs1",    # Thrombospondin 1
      "Hbegf",    # Heparin-binding EGF-like growth factor
      "Fos",      # AP-1 component
      "Jun"       # AP-1 component
    ),

    # =========================================================================
    # INFLAMMATORY/MACROPHAGE-LIKE (Yap et al. 2021)
    # =========================================================================
    macrophage_like = c(
      "Cd68",     # Macrophage marker (low in SMC-derived)
      "Lgals3",   # Galectin-3
      "Aif1",     # Iba1
      "Ctss"      # Cathepsin S
    ),

    # =========================================================================
    # OSTEOGENIC SMC MARKERS (Yap et al. 2021)
    # =========================================================================
    osteogenic = c(
      "Runx2",    # Osteogenic master TF
      "Alpl",     # Alkaline phosphatase
      "Bglap",    # Osteocalcin
      "Spp1",     # Osteopontin
      "Ibsp",     # Bone sialoprotein
      "Mgp"       # Matrix Gla protein (can inhibit or promote)
    ),

    # =========================================================================
    # ADIPOCYTE-LIKE SMC MARKERS (Yap et al. 2021)
    # =========================================================================
    adipocyte_like = c(
      "Pparg",    # Adipogenic master TF
      "Fabp4",    # Fatty acid binding protein 4
      "Cebpa",    # C/EBP alpha
      "Adipoq",   # Adiponectin
      "Lpl"       # Lipoprotein lipase
    ),

    # =========================================================================
    # PROLIFERATION MARKERS
    # =========================================================================
    proliferation = c(
      "Mki67",    # Ki67
      "Pcna",     # Proliferating cell nuclear antigen
      "Top2a",    # Topoisomerase II alpha
      "Cdk1",     # Cyclin-dependent kinase 1
      "Ccnb1",    # Cyclin B1
      "Ccna2"     # Cyclin A2
    ),

    # =========================================================================
    # COMPOSITE SYNTHETIC SCORE (for cluster classification)
    # =========================================================================
    # Weighted combination of key synthetic markers
    synthetic_composite = c(
      # Dedifferentiation (high weight)
      "Klf4", "Klf5",
      # Mesenchymal (medium weight)
      "Vim", "Cd44", "S100a4",
      # Fibroblast-like (medium weight)
      "Fn1", "Col1a1", "Col3a1", "Bgn", "Dcn",
      # Stress-induced (lower weight)
      "Ereg", "Areg", "Thbs1"
    )
  )
}

# ==============================================================================
# CLUSTER-LEVEL CLASSIFICATION
# ==============================================================================

#' Calculate marker scores for each cluster
#' @param seurat_obj Seurat object with seurat_clusters
#' @param marker_sets Named list of marker gene vectors
#' @return Data frame with cluster scores for each marker set
calculate_cluster_marker_scores <- function(seurat_obj, marker_sets = NULL) {
  if (is.null(marker_sets)) {
    marker_sets <- get_synthetic_smc_markers()
  }

  expr_matrix <- GetAssayData(seurat_obj, layer = "data")
  clusters <- seurat_obj$seurat_clusters

  # Get unique cluster IDs (handle both factor and character)
  if (is.factor(clusters)) {
    cluster_ids <- levels(clusters)
  } else {
    cluster_ids <- sort(unique(as.character(clusters)))
  }

  # Initialize results
  results <- data.frame(cluster = cluster_ids, stringsAsFactors = FALSE)

  for (set_name in names(marker_sets)) {
    markers <- marker_sets[[set_name]]
    available <- markers[markers %in% rownames(expr_matrix)]

    if (length(available) == 0) {
      results[[set_name]] <- 0
      next
    }

    scores <- sapply(cluster_ids, function(clust) {
      cells <- which(clusters == clust)
      if (length(cells) == 0) return(0)
      mean(colMeans(as.matrix(expr_matrix[available, cells, drop = FALSE])))
    })

    results[[set_name]] <- scores
  }

  return(results)
}

#' Classify SMC clusters into contractile vs synthetic subtypes
#' @param cluster_scores Data frame from calculate_cluster_marker_scores
#' @param contractile_col Column name for contractile score
#' @param synthetic_col Column name for synthetic score
#' @param ratio_threshold Ratio threshold for classification
#' @return Data frame with SMC subtype classification
classify_smc_clusters <- function(
    cluster_scores,
    contractile_col = "contractile_core",
    synthetic_col = "synthetic_composite",
    ratio_threshold = 1.5,
    smc_threshold = 0.8  # Minimum contractile score to be considered SMC-origin
) {

  # Calculate ratio
  pseudo <- 0.001
  cluster_scores$cs_ratio <- (cluster_scores[[contractile_col]] + pseudo) /
    (cluster_scores[[synthetic_col]] + pseudo)

  # Classify
  cluster_scores$smc_subtype <- case_when(
    # High contractile score = SMC origin
    cluster_scores[[contractile_col]] >= smc_threshold ~ case_when(
      cluster_scores$cs_ratio > ratio_threshold ~ "SMC_Contractile",
      cluster_scores$cs_ratio <= 1/ratio_threshold ~ "SMC_Synthetic",
      TRUE ~ "SMC_Transitional"
    ),
    # Lower contractile but high synthetic = possible dedifferentiated SMC
    cluster_scores[[synthetic_col]] > 0.5 &
      cluster_scores$mesenchymal_like > 0.3 ~ "SMC_Synthetic_Dediff",
    # Default = not SMC
    TRUE ~ "Non_SMC"
  )

  return(cluster_scores)
}

#' Classify cells based on cluster assignments
#' @param seurat_obj Seurat object
#' @param cluster_classifications Data frame with cluster->subtype mapping
#' @return Seurat object with smc_subtype metadata
assign_smc_subtypes_to_cells <- function(seurat_obj, cluster_classifications) {
  # Create mapping
  cluster_to_subtype <- setNames(
    cluster_classifications$smc_subtype,
    as.character(cluster_classifications$cluster)
  )

  # Get cell clusters as character
  cell_clusters <- as.character(seurat_obj$seurat_clusters)

  # Apply mapping directly to metadata slot
  seurat_obj@meta.data$smc_subtype <- cluster_to_subtype[cell_clusters]

  # Report
  cat(sprintf("Assigned subtypes to %d cells\n", sum(!is.na(seurat_obj$smc_subtype))))

  return(seurat_obj)
}

# ==============================================================================
# HIGH-CONFIDENCE SMC SELECTION
# ==============================================================================

#' Get high-confidence SMC cells (contractile + synthetic combined)
#' @param seurat_obj Seurat object with smc_subtype metadata
#' @param include_transitional Whether to include transitional SMCs
#' @return Character vector of cell barcodes
get_high_confidence_smc_cells <- function(seurat_obj, include_transitional = FALSE) {
  subtypes <- c("SMC_Contractile", "SMC_Synthetic")

  if (include_transitional) {
    subtypes <- c(subtypes, "SMC_Transitional", "SMC_Synthetic_Dediff")
  }

  cells <- colnames(seurat_obj)[seurat_obj$smc_subtype %in% subtypes]
  return(cells)
}

#' Get SMC cells by specific subtype
#' @param seurat_obj Seurat object with smc_subtype metadata
#' @param subtype One of: "SMC_Contractile", "SMC_Synthetic", "SMC_Transitional"
#' @return Character vector of cell barcodes
get_smc_cells_by_subtype <- function(seurat_obj, subtype) {
  colnames(seurat_obj)[seurat_obj$smc_subtype == subtype]
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

#' Create SMC subtype UMAP
#' @param seurat_obj Seurat object with smc_subtype
#' @return ggplot object or NULL if UMAP not available
plot_smc_subtypes <- function(seurat_obj) {
  # Check if UMAP is available
  if (!"umap" %in% names(seurat_obj@reductions)) {
    cat("Note: UMAP not available, skipping visualization\n")
    return(NULL)
  }

  # Define colors
  subtype_colors <- c(
    "SMC_Contractile" = "#D32F2F",      # Red
    "SMC_Synthetic" = "#1976D2",         # Blue
    "SMC_Transitional" = "#7B1FA2",      # Purple
    "SMC_Synthetic_Dediff" = "#00796B",  # Teal
    "Non_SMC" = "#BDBDBD"                # Grey
  )

  p <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "smc_subtype",
    cols = subtype_colors,
    pt.size = 0.5,
    shuffle = TRUE
  ) +
    labs(
      title = "SMC Subtype Classification",
      subtitle = sprintf("Contractile: %d | Synthetic/Dediff: %d | Transitional: %d",
                         sum(seurat_obj$smc_subtype == "SMC_Contractile", na.rm = TRUE),
                         sum(grepl("Synthetic", seurat_obj$smc_subtype), na.rm = TRUE),
                         sum(seurat_obj$smc_subtype == "SMC_Transitional", na.rm = TRUE))
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(p)
}

# ==============================================================================
# MAIN CLASSIFICATION PIPELINE
# ==============================================================================

#' Run full SMC subtype classification pipeline
#' @param seurat_obj Seurat object (must have clustering)
#' @param output_dir Directory for output files
#' @param ratio_threshold Contractile/Synthetic ratio threshold
#' @return Seurat object with smc_subtype metadata
run_smc_subtype_classification <- function(
    seurat_obj,
    output_dir = "smc_subtype_classification",
    ratio_threshold = 1.5
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== SMC SUBTYPE CLASSIFICATION PIPELINE ===\n\n")

  # Step 1: Calculate cluster marker scores
  cat("Step 1: Calculating cluster marker scores...\n")
  marker_sets <- get_synthetic_smc_markers()
  cluster_scores <- calculate_cluster_marker_scores(seurat_obj, marker_sets)

  # Step 2: Classify SMC clusters
  cat("Step 2: Classifying SMC clusters...\n")
  cluster_classes <- classify_smc_clusters(
    cluster_scores,
    ratio_threshold = ratio_threshold
  )

  # Save cluster classifications
  write.csv(cluster_classes,
            file.path(output_dir, "cluster_smc_classifications.csv"),
            row.names = FALSE)
  cat("  Saved: cluster_smc_classifications.csv\n")

  # Step 3: Assign to cells
  cat("Step 3: Assigning subtypes to cells...\n")
  seurat_obj <- assign_smc_subtypes_to_cells(seurat_obj, cluster_classes)

  # Step 4: Summary statistics
  cat("\n=== CLASSIFICATION SUMMARY ===\n")
  subtype_counts <- table(seurat_obj$smc_subtype)
  for (subtype in names(subtype_counts)) {
    cat(sprintf("  %s: %d cells (%.1f%%)\n",
                subtype,
                subtype_counts[subtype],
                100 * subtype_counts[subtype] / ncol(seurat_obj)))
  }

  # High-confidence SMC count
  hc_smc <- get_high_confidence_smc_cells(seurat_obj, include_transitional = FALSE)
  cat(sprintf("\nHigh-confidence SMC (Contractile + Synthetic): %d cells\n", length(hc_smc)))

  # Step 5: Create visualization
  cat("\nStep 4: Creating visualizations...\n")
  p <- plot_smc_subtypes(seurat_obj)
  if (!is.null(p)) {
    ggsave(file.path(output_dir, "smc_subtype_umap.png"),
           p, width = 10, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "smc_subtype_umap.pdf"),
           p, width = 10, height = 8, bg = "white")
    cat("  Saved: smc_subtype_umap (PNG + PDF)\n")
  }

  # Save summary
  summary_df <- as.data.frame(subtype_counts)
  colnames(summary_df) <- c("subtype", "n_cells")
  summary_df$percentage <- 100 * summary_df$n_cells / ncol(seurat_obj)
  write.csv(summary_df,
            file.path(output_dir, "smc_subtype_summary.csv"),
            row.names = FALSE)

  cat("\n=== CLASSIFICATION COMPLETE ===\n")

  return(seurat_obj)
}

# ==============================================================================
# UTILITY: Compare with existing methods
# ==============================================================================

#' Compare new cluster-based subtype with existing Myh11+ method
#' @param seurat_obj Seurat object with both classifications
#' @return Comparison summary
compare_smc_methods <- function(seurat_obj) {
  if (!"is_smc" %in% colnames(seurat_obj@meta.data)) {
    cat("Note: is_smc (Myh11+) not found for comparison\n")
    return(NULL)
  }

  # Create comparison table
  comparison <- table(
    Myh11_positive = seurat_obj$is_smc,
    Cluster_subtype = seurat_obj$smc_subtype
  )

  cat("\n=== METHOD COMPARISON ===\n")
  print(comparison)

  # Calculate overlap
  myh11_cells <- colnames(seurat_obj)[seurat_obj$is_smc == TRUE]
  hc_smc_cells <- get_high_confidence_smc_cells(seurat_obj)

  overlap <- length(intersect(myh11_cells, hc_smc_cells))

  cat(sprintf("\nMyh11+ cells: %d\n", length(myh11_cells)))
  cat(sprintf("High-confidence SMC (cluster): %d\n", length(hc_smc_cells)))
  cat(sprintf("Overlap: %d (%.1f%% of Myh11+, %.1f%% of cluster)\n",
              overlap,
              100 * overlap / length(myh11_cells),
              100 * overlap / length(hc_smc_cells)))

  return(comparison)
}
