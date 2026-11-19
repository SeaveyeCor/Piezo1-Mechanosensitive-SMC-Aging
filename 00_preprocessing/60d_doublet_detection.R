#!/usr/bin/env Rscript
# =============================================================================
# SYSTEMATIC DOUBLET DETECTION ON SMC POPULATION
# Uses scDblFinder to identify doublets across ALL clusters, not just Cl5.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
})

setwd(".")

# Load data
obj <- readRDS("data/seurat_paper_6samples.rds")
obj$sample <- as.character(obj$sample)
obj$sample[is.na(obj$sample)] <- as.character(obj$orig.ident[is.na(obj$sample)])

smc <- subset(obj, subset = cell_type == "SMC_Contractile")
meta_smc <- smc@meta.data
meta_smc$genotype <- ifelse(grepl("KO", meta_smc$sample), "KO", "Ctrl")
meta_smc$age <- ifelse(grepl("Young", meta_smc$sample), "Young", "Aged")

cat(sprintf("Total SMCs: %d\n", ncol(smc)))

# Convert to SCE for scDblFinder
sce <- as.SingleCellExperiment(smc)

# Run scDblFinder per sample (recommended approach)
set.seed(42)
sce <- scDblFinder(sce, samples = "sample")

# Extract results
dbl_class <- sce$scDblFinder.class
dbl_score <- sce$scDblFinder.score
names(dbl_class) <- colnames(sce)
names(dbl_score) <- colnames(sce)

cat(sprintf("\n=== DOUBLET DETECTION RESULTS ===\n"))
cat(sprintf("Singlets: %d\n", sum(dbl_class == "singlet")))
cat(sprintf("Doublets: %d (%.1f%%)\n", sum(dbl_class == "doublet"),
    100 * mean(dbl_class == "doublet")))

# Breakdown by cluster
bf_clusters <- c("1", "5", "15", "16")
cat("\n--- By Cluster (all cells) ---\n")
cat(sprintf("%-8s %8s %8s %8s\n", "Cluster", "Total", "Doublets", "%Dbl"))
for (cl in sort(unique(as.character(meta_smc$seurat_clusters)))) {
  cells_cl <- rownames(meta_smc)[as.character(meta_smc$seurat_clusters) == cl]
  n_dbl <- sum(dbl_class[cells_cl] == "doublet")
  bf_tag <- ifelse(cl %in% bf_clusters, " *", "")
  cat(sprintf("Cl%-6s %8d %8d %7.1f%%%s\n", cl, length(cells_cl), n_dbl,
      100 * n_dbl / length(cells_cl), bf_tag))
}

# Breakdown by cluster x genotype for bona fide
cat("\n--- Bona Fide Clusters by Group ---\n")
cat(sprintf("%-20s %8s %8s %8s\n", "Cluster_Group", "Total", "Doublets", "%Dbl"))
for (cl in bf_clusters) {
  for (geno in c("Ctrl", "KO")) {
    cells_cl <- rownames(meta_smc)[as.character(meta_smc$seurat_clusters) == cl &
                                    meta_smc$genotype == geno]
    if (length(cells_cl) == 0) next
    n_dbl <- sum(dbl_class[cells_cl] == "doublet")
    cat(sprintf("Cl%-3s %-14s %8d %8d %7.1f%%\n", cl, geno, length(cells_cl), n_dbl,
        100 * n_dbl / length(cells_cl)))
  }
}

# Check: do the scDblFinder doublets overlap with our k-means G2?
cat("\n--- Overlap with K-means ICC doublets ---\n")
umap_raw <- read.csv("R01_Piezo1_Aging_Gut/outputs/data/PITFALL_UMAP_coordinates.csv", row.names = 1)
cl5_ctrl <- rownames(meta_smc)[as.character(meta_smc$seurat_clusters) == "5" & meta_smc$genotype == "Ctrl"]
cl5_ctrl <- cl5_ctrl[cl5_ctrl %in% rownames(umap_raw)]
u1 <- umap_raw[cl5_ctrl, "UMAP_1"]
u2 <- umap_raw[cl5_ctrl, "UMAP_2"]
km <- kmeans(cbind(u1, u2), centers = 2, nstart = 10)
expr_all <- as.matrix(GetAssayData(smc, layer = "data"))
g1_kit <- mean(as.numeric(expr_all["Kit", cl5_ctrl[km$cluster == 1]]))
g2_kit <- mean(as.numeric(expr_all["Kit", cl5_ctrl[km$cluster == 2]]))
doublet_group <- ifelse(g1_kit > g2_kit, 1, 2)
kmeans_doublets <- cl5_ctrl[km$cluster == doublet_group]

scdblfinder_dbl_in_cl5 <- cl5_ctrl[dbl_class[cl5_ctrl] == "doublet"]
overlap <- sum(kmeans_doublets %in% scdblfinder_dbl_in_cl5)
cat(sprintf("K-means ICC doublets: %d\n", length(kmeans_doublets)))
cat(sprintf("scDblFinder doublets in Cl5 Ctrl: %d\n", length(scdblfinder_dbl_in_cl5)))
cat(sprintf("Overlap: %d (%.0f%% of k-means set)\n", overlap, 100 * overlap / length(kmeans_doublets)))

# ICC marker check on scDblFinder doublets vs singlets in Cl5
cat("\n--- ICC markers: scDblFinder doublets vs singlets in Cl5 Ctrl ---\n")
cl5_dbl <- cl5_ctrl[dbl_class[cl5_ctrl] == "doublet"]
cl5_sng <- cl5_ctrl[dbl_class[cl5_ctrl] == "singlet"]
cat(sprintf("%-10s %10s %10s\n", "Gene", "Doublets", "Singlets"))
for (g in c("Kit", "Ano1", "Etv1", "Prkcq", "Myh11", "Acta2", "Pdgfra", "Ptprc")) {
  vd <- as.numeric(expr_all[g, cl5_dbl])
  vs <- as.numeric(expr_all[g, cl5_sng])
  cat(sprintf("%-10s %10.3f %10.3f\n", g, mean(vd), mean(vs)))
}

# Save results
results <- data.frame(
  cell = colnames(sce),
  sample = meta_smc[colnames(sce), "sample"],
  cluster = as.character(meta_smc[colnames(sce), "seurat_clusters"]),
  dbl_class = as.character(dbl_class),
  dbl_score = dbl_score,
  stringsAsFactors = FALSE
)
write.csv(results, "output/doublets/scDblFinder_results.csv",
          row.names = FALSE)

# Save just the doublet cell IDs for easy filtering
doublet_ids <- results$cell[results$dbl_class == "doublet"]
writeLines(doublet_ids, "data/doublet_cell_ids.txt")

cat(sprintf("\nSaved %d doublet cell IDs\n", length(doublet_ids)))
cat("Done!\n")
