# ==============================================================================
# CLUSTER-BASED SMC SELECTION
# ==============================================================================
# Selects SMCs based on cell_type == "SMC_Contractile" from clustering
# This approach identifies 6,378 SMC cells from 8 pure contractile clusters
# ==============================================================================

#' Get method metadata
#'
#' @return List with method information
get_method_info <- function() {
  list(
    method_name = "cluster_based",
    description = paste(
      "Cluster-based bona fide SMC selection: SMCs identified from cell_type ==",
      "'SMC_Contractile', then filtered to 4 bona fide visceral SMC clusters",
      "(Cl1: Transitional Hub, Cl5: Rescued Synthetic, Cl15: Core Contractile,",
      "Cl16: Rescue Destination). Excludes Cl6 (lymphatic vessel SMC — Prox1/Flt4/Lyve1",
      "in 60-80%, NFAT dissociation) and Cl13 (NMJ SMC — Piezo1 at 5% det, no aging",
      "effect p=0.21). Doublet cells removed if doublet_cell_ids.txt is present."
    ),
    version = "2.0.0",
    n_expected_cells = NA,
    selection_column = "cell_type",
    selection_value = "SMC_Contractile",
    bona_fide_clusters = c("1", "5", "15", "16")
  )
}

#' Select cells for analysis
#'
#' @param seurat_obj Seurat object with all cells
#' @return List with smc, pdgfra, icc cell barcodes and metadata
select_cells <- function(seurat_obj) {
  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Select SMC cells based on cell_type
  if (!"cell_type" %in% colnames(metadata)) {
    stop("cell_type column not found in metadata. Required for cluster-based selection.")
  }

  smc_cells <- metadata$cell_barcode[metadata$cell_type == "SMC_Contractile"]

  # --------------------------------------------------------------------------
  # Bona fide cluster filter: keep only Cl1, Cl5, Cl15, Cl16
  # Excludes Cl6 (lymphatic vessel SMC) and Cl13 (NMJ SMC)
  # --------------------------------------------------------------------------
  bf_clusters <- c("1", "5", "15", "16")
  if ("seurat_clusters" %in% colnames(metadata)) {
    smc_cluster <- as.character(metadata[smc_cells, "seurat_clusters"])
    bf_mask <- smc_cluster %in% bf_clusters
    n_removed_bf <- sum(!bf_mask)
    smc_cells <- smc_cells[bf_mask]
    cat(sprintf("  Bona fide filter: kept %d SMCs in clusters %s (removed %d from Cl6/Cl13/other)\n",
                length(smc_cells), paste(bf_clusters, collapse = ","), n_removed_bf))
  } else {
    warning("seurat_clusters column not found — skipping bona fide cluster filter")
  }

  # --------------------------------------------------------------------------
  # Doublet removal (graceful skip if file not found)
  # --------------------------------------------------------------------------
  doublet_path <- file.path(getwd(), "data/doublet_cell_ids.txt")
  if (file.exists(doublet_path)) {
    doublet_ids <- readLines(doublet_path)
    n_doublets <- sum(smc_cells %in% doublet_ids)
    smc_cells <- setdiff(smc_cells, doublet_ids)
    cat(sprintf("  Doublet removal: removed %d cells (%d IDs in file)\n",
                n_doublets, length(doublet_ids)))
  } else {
    cat("  Doublet file not found — skipping doublet removal\n")
  }

  # Select PDGFRa+ cells
  pdgfra_cells <- if ("PDGFRa" %in% unique(metadata$cell_type)) {
    metadata$cell_barcode[metadata$cell_type == "PDGFRa"]
  } else if ("Fibroblast" %in% unique(metadata$cell_type)) {
    metadata$cell_barcode[metadata$cell_type == "Fibroblast"]
  } else {
    character(0)
  }

  # Select ICC cells
  icc_cells <- if ("ICC" %in% unique(metadata$cell_type)) {
    metadata$cell_barcode[metadata$cell_type == "ICC"]
  } else {
    character(0)
  }

  # Create classification metadata
  classification <- data.frame(
    cell_barcode = metadata$cell_barcode,
    is_smc = metadata$cell_barcode %in% smc_cells,
    is_pdgfra = metadata$cell_barcode %in% pdgfra_cells,
    is_icc = metadata$cell_barcode %in% icc_cells,
    cell_type = metadata$cell_type,
    stringsAsFactors = FALSE
  )

  # Add group info if available
  if ("group" %in% colnames(metadata)) {
    classification$group <- metadata$group
  }
  if ("age" %in% colnames(metadata)) {
    classification$age <- metadata$age
  }
  if ("genotype" %in% colnames(metadata)) {
    classification$genotype <- metadata$genotype
  }

  # Report selection
  cat("\n========================================\n")
  cat("CLUSTER-BASED SMC SELECTION\n")
  cat("========================================\n")
  cat(sprintf("Total cells in dataset: %d\n", nrow(metadata)))
  cat(sprintf("SMC cells selected: %d\n", length(smc_cells)))
  cat(sprintf("PDGFRa+ cells: %d\n", length(pdgfra_cells)))
  cat(sprintf("ICC cells: %d\n", length(icc_cells)))

  # Group breakdown if available
  if ("group" %in% colnames(classification)) {
    cat("\nSMC cells by group:\n")
    smc_by_group <- table(classification$group[classification$is_smc])
    for (grp in names(smc_by_group)) {
      cat(sprintf("  %s: %d\n", grp, smc_by_group[grp]))
    }
  }
  cat("========================================\n\n")

  return(list(
    smc = smc_cells,
    pdgfra = pdgfra_cells,
    icc = icc_cells,
    metadata = classification
  ))
}

#' Get detailed cluster information
#'
#' @return Data frame with cluster details
get_cluster_info <- function() {
  data.frame(
    cluster = c(1, 5, 15, 16, 6, 13),
    description = c(
      "Transitional Hub (bona fide)",
      "Rescued Synthetic (bona fide)",
      "Core Contractile (bona fide)",
      "Rescue Destination (bona fide)",
      "Lymphatic Vessel SMC (EXCLUDED — Prox1/Flt4/Lyve1 60-80%, NFAT dissociation)",
      "NMJ SMC (EXCLUDED — Piezo1 5% det, no aging effect p=0.21)"
    ),
    included = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}
