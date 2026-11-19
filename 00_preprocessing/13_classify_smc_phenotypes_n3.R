#!/usr/bin/env Rscript
# SMC Phenotype Classification with all 3 Aged_KO samples (n=3)

library(Seurat)
library(dplyr)
library(ggplot2)

cat("==========================================================\n")
cat("SMC PHENOTYPE CLASSIFICATION - ALL 3 AGED_KO SAMPLES\n")
cat("==========================================================\n\n")

# Load fully annotated dataset (re-normalized for all cells)
cat("Loading fully annotated and normalized dataset (n=3)...\n")
seurat <- readRDS("FULL_CELL_TYPE_ANNOTATION/seurat_merged_annotated_complete_n3_normalized.rds")

cat("Total cells:", ncol(seurat), "\n")
cat("Total Aged_KO:", sum(seurat$group == "Aged_KO"), "\n\n")

# Define marker sets (same as before)
contractile_markers <- c("Myh11", "Acta2", "Tagln", "Cnn1", "Smtn", "Myocd")
synthetic_markers <- c(
  "Klf4", "Nfkb1", "Nfkb2", "Rela",
  "Areg", "Ereg",
  "Il1b", "Il6", "Il11", "Ptgs2", "Cxcl3", "Ccl20",
  "Pdgfa", "Thbs1", "Icam1"
)

# Check marker presence
contractile_present <- contractile_markers[contractile_markers %in% rownames(seurat)]
synthetic_present <- synthetic_markers[synthetic_markers %in% rownames(seurat)]

cat("Contractile markers available:", length(contractile_present), "/", length(contractile_markers), "\n")
cat("Synthetic markers available:", length(synthetic_present), "/", length(synthetic_markers), "\n\n")

# Calculate module scores
cat("Calculating module scores...\n")
seurat <- AddModuleScore(
  seurat,
  features = list(contractile_present),
  name = "Contractile_Score"
)
seurat <- AddModuleScore(
  seurat,
  features = list(synthetic_present),
  name = "Synthetic_Score"
)

# Calculate ratio
seurat@meta.data$contractile_synthetic_ratio <- seurat@meta.data$Contractile_Score1 / 
  (seurat@meta.data$Synthetic_Score1 + 0.1)

# Identify SMCs
# Extract Myh11 expression for all cells
myh11_expr <- FetchData(seurat, vars = "Myh11", layer = "data")[,1]
seurat@meta.data$myh11_expr <- myh11_expr
seurat@meta.data$is_smc <- myh11_expr > 0

# Classify SMC phenotype
cat("Classifying SMC phenotypes...\n")
seurat@meta.data$smc_phenotype <- "Non_SMC"
smc_mask <- seurat@meta.data$is_smc

seurat@meta.data$smc_phenotype[smc_mask] <- dplyr::case_when(
  seurat@meta.data$Contractile_Score1[smc_mask] > 0.5 &
    seurat@meta.data$contractile_synthetic_ratio[smc_mask] > 3 ~ "SMC_Contractile",
  
  seurat@meta.data$Synthetic_Score1[smc_mask] > 0.3 &
    seurat@meta.data$contractile_synthetic_ratio[smc_mask] < 1 &
    seurat@meta.data$myh11_expr[smc_mask] > 0.5 ~ "SMC_Synthetic",
  
  seurat@meta.data$contractile_synthetic_ratio[smc_mask] >= 1 &
    seurat@meta.data$contractile_synthetic_ratio[smc_mask] <= 3 ~ "SMC_Intermediate",
  
  TRUE ~ "SMC_Undifferentiated"
)

cat("\n=== OVERALL PHENOTYPE DISTRIBUTION ===\n")
print(table(seurat@meta.data$smc_phenotype))
cat("\n")

# By group
cat("=== PHENOTYPE DISTRIBUTION BY GROUP ===\n")
phenotype_by_group <- seurat@meta.data %>%
  filter(is_smc) %>%
  group_by(smc_phenotype, group) %>%
  summarise(n_cells = n(), .groups = "drop")

print(phenotype_by_group)
cat("\n")

# Calculate percentages within each group
cat("=== PHENOTYPE PERCENTAGES WITHIN EACH GROUP ===\n")
phenotype_pct <- seurat@meta.data %>%
  filter(is_smc) %>%
  group_by(group) %>%
  mutate(group_total = n()) %>%
  group_by(group, smc_phenotype) %>%
  summarise(
    n_cells = n(),
    pct = 100 * n() / first(group_total),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = smc_phenotype,
    names_from = group,
    values_from = c(n_cells, pct),
    values_fill = 0
  )

print(phenotype_pct)
cat("\n")

# Save results
cat("Saving results...\n")
saveRDS(seurat, "FULL_CELL_TYPE_ANNOTATION/seurat_smc_phenotypes_n3.rds")

write.csv(
  phenotype_by_group,
  "FULL_CELL_TYPE_ANNOTATION/results/phenotype_by_group_n3.csv",
  row.names = FALSE
)

write.csv(
  phenotype_pct,
  "FULL_CELL_TYPE_ANNOTATION/results/phenotype_percentages_n3.csv",
  row.names = FALSE
)

cat("\n==========================================================\n")
cat("SMC PHENOTYPE CLASSIFICATION COMPLETE (n=3)!\n")
cat("==========================================================\n")
