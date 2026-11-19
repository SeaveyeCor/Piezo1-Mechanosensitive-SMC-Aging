setwd(".")
suppressPackageStartupMessages({ library(Seurat); library(openxlsx); library(dplyr) })
source("pipeline/selection_criteria/cluster_based_smc.R")

cat("Loading Seurat...\n")
obj <- readRDS("data/seurat_paper_6samples.rds")
sel <- select_cells(obj)
smc <- sel$smc
meta <- obj@meta.data[smc, ]
ctrl <- smc[meta$genotype == "Ctrl"]
meta_ctrl <- meta[ctrl, ]
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
p1 <- as.numeric(expr["Piezo1", ctrl])
young <- meta_ctrl$age == "Young"
aged <- meta_ctrl$age == "Aged"

sd_dir <- "output/source_data"

# ==================================================================
# 1. Fig S4A source data (gene-level Piezo1 correlations)
# ==================================================================
cat("\n1. Creating Fig S4A source data...\n")

# Read Table S7 for the gene list
t7 <- read.csv("output/tables/Table_S7_Piezo1_Correlations.csv")
cat("  Table S7:", nrow(t7), "genes\n")

# Recompute per-gene correlations split by age (Table S7 is overall, we need by age)
results <- data.frame()
for (i in 1:nrow(t7)) {
  g <- t7$Gene[i]
  if (!(g %in% rownames(expr))) next
  gx <- as.numeric(expr[g, ctrl])

  ct_y <- cor.test(p1[young], gx[young], method = "spearman")
  ct_a <- cor.test(p1[aged], gx[aged], method = "spearman")
  ct_all <- cor.test(p1, gx, method = "spearman")

  pct_det <- round(100 * sum(gx > 0) / length(gx), 1)

  results <- rbind(results, data.frame(
    gene = g,
    category = t7$Category[i],
    pct_detecting = pct_det,
    rho_overall = round(ct_all$estimate, 4),
    p_overall = ct_all$p.value,
    rho_young = round(ct_y$estimate, 4),
    p_young = ct_y$p.value,
    rho_aged = round(ct_a$estimate, 4),
    p_aged = ct_a$p.value,
    mean_expr_young = round(mean(gx[young]), 4),
    sd_expr_young = round(sd(gx[young]), 4),
    mean_expr_aged = round(mean(gx[aged]), 4),
    sd_expr_aged = round(sd(gx[aged]), 4),
    stringsAsFactors = FALSE
  ))
}

wb <- createWorkbook()
addWorksheet(wb, "title")
writeData(wb, "title", data.frame(
  Field = c("Figure", "Title", "Manuscript"),
  Value = c("Fig S4A", "Gene-level Piezo1 Correlations",
            "Maladaptive Piezo1 Mechanotransduction Drives Smooth Muscle Aging in the Gut")
))
addWorksheet(wb, "gene_correlations")
writeData(wb, "gene_correlations", results)
addWorksheet(wb, "population_summary")
writeData(wb, "population_summary", data.frame(
  population = c("Ctrl_SMCs", "Young_Ctrl", "Aged_Ctrl"),
  n = c(length(ctrl), sum(young), sum(aged))
))
saveWorkbook(wb, file.path(sd_dir, "FigS4A_gene_piezo1_correlations_data.xlsx"), overwrite = TRUE)
cat("  Saved FigS4A_gene_piezo1_correlations_data.xlsx (", nrow(results), " genes)\n")

cat("\nDone.\n")
