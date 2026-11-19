# ==============================================================================
# PHENOTYPE CLASSIFIER FOR SMC ANALYSIS PIPELINE
# ==============================================================================
# Classifies SMCs as Contractile or Synthetic based on marker expression
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# Source gene sets if not already loaded
if (!exists("get_gene_sets")) {
  # Try to find gene_sets.R relative to this file's location
  tryCatch({
    source(file.path(dirname(sys.frame(1)$ofile), "gene_sets.R"))
  }, error = function(e) {
    # Fallback: assume working directory is pipeline root
    source("core/gene_sets.R")
  })
}

#' Classify SMC phenotype based on contractile vs synthetic scores
#'
#' @param seurat_obj Seurat object with SMC cells
#' @param smc_cells Character vector of SMC cell barcodes
#' @param contractile_genes Character vector of contractile marker genes
#' @param synthetic_genes Character vector of synthetic marker genes
#' @param ratio_threshold Contractile/Synthetic ratio threshold for classification (default 1.5)
#' @param assay Assay to use for expression data (default "RNA")
#' @param slot Slot to use for expression data (default "data")
#'
#' @return Data frame with cell barcodes, scores, and phenotype classification
classify_smc_phenotype <- function(
    seurat_obj,
    smc_cells,
    contractile_genes = NULL,
    synthetic_genes = NULL,
    ratio_threshold = 1.5,
    assay = "RNA",
    slot = "data"
) {
  # Get default gene sets if not provided
  gene_sets <- get_gene_sets()

  if (is.null(contractile_genes)) {
    contractile_genes <- gene_sets$contractile_core
  }
  if (is.null(synthetic_genes)) {
    synthetic_genes <- gene_sets$synthetic_core
  }

  # Get expression matrix
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)

  # Filter to SMC cells
  smc_cells <- smc_cells[smc_cells %in% colnames(expr_matrix)]
  expr_smc <- expr_matrix[, smc_cells, drop = FALSE]

  # Filter genes to those available in the matrix
  available_genes <- rownames(expr_smc)
  contractile_avail <- contractile_genes[contractile_genes %in% available_genes]
  synthetic_avail <- synthetic_genes[synthetic_genes %in% available_genes]

  # Report gene availability
  cat(sprintf("Contractile genes: %d/%d available\n",
              length(contractile_avail), length(contractile_genes)))
  cat(sprintf("Synthetic genes: %d/%d available\n",
              length(synthetic_avail), length(synthetic_genes)))

  if (length(contractile_avail) == 0 || length(synthetic_avail) == 0) {
    stop("Insufficient genes available for phenotype classification")
  }

  # Calculate scores for each cell
  contractile_scores <- colMeans(as.matrix(expr_smc[contractile_avail, , drop = FALSE]))
  synthetic_scores <- colMeans(as.matrix(expr_smc[synthetic_avail, , drop = FALSE]))

  # Add pseudocount to avoid division by zero
  pseudo <- 0.001
  ratio <- (contractile_scores + pseudo) / (synthetic_scores + pseudo)

  # Classify phenotype
  phenotype <- ifelse(ratio > ratio_threshold, "Contractile", "Synthetic")

  # Create result data frame
  result <- data.frame(
    cell_barcode = smc_cells,
    contractile_score = contractile_scores,
    synthetic_score = synthetic_scores,
    score_ratio = ratio,
    phenotype = phenotype,
    stringsAsFactors = FALSE
  )

  # Report classification summary
  cat(sprintf("\nPhenotype classification (ratio threshold = %.2f):\n", ratio_threshold))
  cat(sprintf("  Contractile: %d (%.1f%%)\n",
              sum(phenotype == "Contractile"),
              100 * mean(phenotype == "Contractile")))
  cat(sprintf("  Synthetic: %d (%.1f%%)\n",
              sum(phenotype == "Synthetic"),
              100 * mean(phenotype == "Synthetic")))

  return(result)
}

#' Calculate module scores for gene sets
#'
#' @param seurat_obj Seurat object
#' @param cells Character vector of cell barcodes
#' @param gene_set Character vector of genes
#' @param assay Assay to use
#' @param slot Slot to use
#'
#' @return Named numeric vector of scores
calculate_module_score <- function(
    seurat_obj,
    cells,
    gene_set,
    assay = "RNA",
    slot = "data"
) {
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)

  cells <- cells[cells %in% colnames(expr_matrix)]
  expr_subset <- expr_matrix[, cells, drop = FALSE]

  available_genes <- gene_set[gene_set %in% rownames(expr_subset)]

  if (length(available_genes) == 0) {
    warning("No genes from gene set found in expression matrix")
    return(setNames(rep(NA, length(cells)), cells))
  }

  scores <- colMeans(as.matrix(expr_subset[available_genes, , drop = FALSE]))
  return(scores)
}

#' Get phenotype summary statistics by group
#'
#' @param phenotype_df Data frame from classify_smc_phenotype
#' @param metadata Data frame with cell metadata including group info
#' @param group_col Column name for group assignment
#'
#' @return Data frame with phenotype counts and percentages by group
get_phenotype_summary <- function(phenotype_df, metadata, group_col = "group") {
  # Merge phenotype with metadata
  merged <- merge(
    phenotype_df,
    metadata[, c("cell_barcode", group_col), drop = FALSE],
    by = "cell_barcode"
  )

  # Calculate summary by group
  summary_df <- merged %>%
    group_by(!!sym(group_col), phenotype) %>%
    summarize(
      n_cells = n(),
      mean_contractile_score = mean(contractile_score),
      mean_synthetic_score = mean(synthetic_score),
      .groups = "drop"
    ) %>%
    group_by(!!sym(group_col)) %>%
    mutate(
      total_cells = sum(n_cells),
      pct = 100 * n_cells / total_cells
    ) %>%
    ungroup()

  return(summary_df)
}

#' Classify cells and add to Seurat metadata
#'
#' @param seurat_obj Seurat object
#' @param smc_cells Character vector of SMC cell barcodes
#' @param ... Additional arguments passed to classify_smc_phenotype
#'
#' @return Seurat object with phenotype metadata added
add_phenotype_to_seurat <- function(seurat_obj, smc_cells, ...) {
  # Classify phenotypes
  phenotype_df <- classify_smc_phenotype(seurat_obj, smc_cells, ...)

  # Initialize columns with NA
  seurat_obj$smc_phenotype <- NA
  seurat_obj$contractile_score <- NA
  seurat_obj$synthetic_score <- NA
  seurat_obj$phenotype_ratio <- NA

  # Add to metadata for SMC cells
  idx <- match(phenotype_df$cell_barcode, colnames(seurat_obj))
  seurat_obj$smc_phenotype[idx] <- phenotype_df$phenotype
  seurat_obj$contractile_score[idx] <- phenotype_df$contractile_score
  seurat_obj$synthetic_score[idx] <- phenotype_df$synthetic_score
  seurat_obj$phenotype_ratio[idx] <- phenotype_df$score_ratio

  return(seurat_obj)
}
