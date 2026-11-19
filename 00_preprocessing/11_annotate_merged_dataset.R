#!/usr/bin/env Rscript
# Annotate the merged dataset with all 3 Aged_KO samples

library(Seurat)
library(dplyr)

cat("===============================================\n")
cat("CELL TYPE ANNOTATION - MERGED DATASET (n=3)\n")
cat("===============================================\n\n")

# Load merged dataset
cat("Loading merged dataset with all 3 Aged_KO samples...\n")
seurat <- readRDS("Fulldataset-SMC_Piezo1-KO_Young-Aged_WITH_MAY15.rds")

cat("Dataset summary:\n")
cat("  Total cells:", ncol(seurat), "\n")
cat("  Total genes:", nrow(seurat), "\n\n")

cat("Sample distribution:\n")
print(table(seurat$orig.ident))
cat("\n")

cat("Group distribution:\n")
print(table(seurat$group))
cat("\n")

cat("Aged_KO sample breakdown:\n")
aged_ko <- subset(seurat, group == "Aged_KO")
print(table(aged_ko$orig.ident))
cat("Total Aged_KO:", ncol(aged_ko), "\n\n")

# Load existing cell type annotations
cat("Loading existing cell type annotations...\n")
annotations <- read.csv("FULL_CELL_TYPE_ANNOTATION/results/annotation_template_REVIEWED.csv")
cat("Annotations loaded for", nrow(annotations), "clusters\n\n")

# Apply annotations to cells based on seurat_clusters
cat("Applying cell type annotations...\n")
cluster_to_celltype <- setNames(
  annotations$manual_annotation,
  as.character(annotations$cluster)
)

# Direct assignment to avoid cell ID mismatch
seurat@meta.data$cell_type <- cluster_to_celltype[as.character(seurat@meta.data$seurat_clusters)]

cat("Cell type distribution:\n")
print(sort(table(seurat$cell_type), decreasing = TRUE))
cat("\n")

# Check cell type distribution by group
cat("Cell type distribution by experimental group:\n")
celltype_by_group <- seurat@meta.data %>%
  group_by(cell_type, group) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group, values_from = n_cells, values_fill = 0)

print(celltype_by_group)
cat("\n")

# Save annotated dataset
cat("Saving annotated merged dataset...\n")
saveRDS(seurat, "FULL_CELL_TYPE_ANNOTATION/seurat_merged_annotated_n3.rds")
cat("Saved to: FULL_CELL_TYPE_ANNOTATION/seurat_merged_annotated_n3.rds\n\n")

# Save cell type counts
write.csv(
  celltype_by_group,
  "FULL_CELL_TYPE_ANNOTATION/results/celltype_by_group_n3.csv",
  row.names = FALSE
)

cat("Cell type annotation complete!\n")
cat("Ready for phenotype classification.\n")
