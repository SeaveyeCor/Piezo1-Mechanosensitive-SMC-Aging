# ==============================================================================
# FIX SOURCE DATA: Add SD/SEM/N where missing + regenerate Fig3F
# ==============================================================================

setwd(".")

suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
  library(dplyr)
})

sd_dir <- "output/source_data"

# Load Seurat for cell-level SD calculations
cat("Loading Seurat object...\n")
obj <- readRDS("data/seurat_paper_6samples.rds")
source("pipeline/selection_criteria/cluster_based_smc.R")
source("pipeline/core/gene_sets.R")
# Use the analysed SMCs (from the pipeline's phenotype_by_cell sheet).
smc <- read.xlsx(file.path(sd_dir, "Fig3A_data.xlsx"), sheet = "phenotype_by_cell")$cell_barcode
stopifnot(length(smc) > 0, all(smc %in% colnames(obj)))
meta <- obj@meta.data[smc, ]
meta$group <- paste0(
  ifelse(meta$age == "Young", "Young", "Aged"), "_",
  ifelse(meta$genotype == "Ctrl", "Ctrl", "KO")
)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()

# ==============================================================================
# 1. REGENERATE Fig3F_data.xlsx (synthetic score bar plot with SD)
# ==============================================================================
cat("\n1. Regenerating Fig3F_data.xlsx (synthetic score bar)...\n")

synth_genes <- gs$synthetic_phenotype
synth_avail <- synth_genes[synth_genes %in% rownames(expr)]
score <- colMeans(as.matrix(expr[synth_avail, smc, drop = FALSE]))

df <- data.frame(score = score, group = meta$group)
stats <- df %>%
  group_by(group) %>%
  summarize(
    n = n(),
    mean = mean(score),
    sd = sd(score),
    sem = sd(score) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(genotype = ifelse(grepl("Ctrl", group), "Ctrl", "KO"))

# Pairwise tests
comps <- list(
  c("Young_Ctrl", "Aged_Ctrl", "Aging (WT)"),
  c("Young_KO", "Aged_KO", "Aging (KO)"),
  c("Young_Ctrl", "Young_KO", "KO effect (young)"),
  c("Aged_Ctrl", "Aged_KO", "KO rescue (aged)")
)

test_results <- data.frame()
for (comp in comps) {
  v1 <- df$score[df$group == comp[1]]
  v2 <- df$score[df$group == comp[2]]
  wt <- wilcox.test(v1, v2)
  p_adj <- min(wt$p.value * 4, 1)
  sig <- ifelse(p_adj < 0.001, "***", ifelse(p_adj < 0.01, "**", ifelse(p_adj < 0.05, "*", "ns")))
  test_results <- rbind(test_results, data.frame(
    comparison = comp[3], a = comp[1], b = comp[2],
    n_a = length(v1), n_b = length(v2),
    mean_a = mean(v1), sd_a = sd(v1),
    mean_b = mean(v2), sd_b = sd(v2),
    statistic = unname(wt$statistic),
    p_value = wt$p.value, p_adjusted = p_adj, significance = sig
  ))
}

wb <- createWorkbook()
addWorksheet(wb, "summary_statistics")
writeData(wb, "summary_statistics", stats)
addWorksheet(wb, "statistical_tests")
writeData(wb, "statistical_tests", test_results)
addWorksheet(wb, "cell_level_data")
writeData(wb, "cell_level_data", data.frame(cell_barcode = names(score), group = meta$group, synthetic_score = score))
addWorksheet(wb, "genes_used")
writeData(wb, "genes_used", data.frame(gene = synth_avail))
saveWorkbook(wb, file.path(sd_dir, "Fig3F_data.xlsx"), overwrite = TRUE)
cat("  Saved Fig3F_data.xlsx\n")

# ==============================================================================
# 2. FIX Fig3G_data.xlsx and Fig3H_data.xlsx: Add per-gene SD/N
# ==============================================================================
for (fig in c("Fig3G", "Fig3H")) {
  cat(sprintf("\n2. Fixing %s_data.xlsx (adding per-gene SD/N)...\n", fig))
  fpath <- file.path(sd_dir, paste0(fig, "_data.xlsx"))
  wb <- loadWorkbook(fpath)
  summ <- read.xlsx(wb, sheet = "summary")

  # Compute per-gene SD from cell-level expression
  for (i in 1:nrow(summ)) {
    gene <- summ$gene[i]
    if (!(gene %in% rownames(expr))) next
    gexpr <- as.numeric(expr[gene, smc])
    for (grp in c("Young_Ctrl", "Aged_Ctrl", "Aged_KO")) {
      vals <- gexpr[meta$group == grp]
      col_sd <- paste0("sd_", gsub("_", "", grp))
      col_n <- paste0("n_", gsub("_", "", grp))
      summ[i, col_sd] <- sd(vals)
      summ[i, col_n] <- length(vals)
    }
  }

  removeWorksheet(wb, "summary")
  addWorksheet(wb, "summary", gridLines = FALSE)
  writeData(wb, "summary", summ)
  worksheetOrder(wb) <- c(which(names(wb) == "summary"), which(names(wb) != "summary"))
  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat(sprintf("  Added sd_*/n_* columns to summary sheet\n"))
}

# ==============================================================================
# 3. FIX Fig3N_data.xlsx: Add SD/N to per_gene_summary
# ==============================================================================
cat("\n3. Fixing Fig3N_data.xlsx (adding per-gene SD/N)...\n")
fpath <- file.path(sd_dir, "Fig3N_data.xlsx")
wb <- loadWorkbook(fpath)
summ <- read.xlsx(wb, sheet = "per_gene_summary")

# Cyto remodeling genes
cyto_genes <- gs$cyto_remodeling
for (i in 1:nrow(summ)) {
  gene <- summ$gene[i]
  if (!(gene %in% rownames(expr))) next
  gexpr <- as.numeric(expr[gene, smc])
  for (grp in c("Young_Ctrl", "Aged_Ctrl", "Aged_KO")) {
    vals <- gexpr[meta$group == grp]
    col_sd <- paste0("sd_", gsub("_", "", grp))
    col_n <- paste0("n_", gsub("_", "", grp))
    summ[i, col_sd] <- sd(vals)
    summ[i, col_n] <- length(vals)
  }
}

removeWorksheet(wb, "per_gene_summary")
addWorksheet(wb, "per_gene_summary", gridLines = FALSE)
writeData(wb, "per_gene_summary", summ)
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("  Added sd_*/n_* columns\n")

# ==============================================================================
# 4. FIX Fig3B_data.xlsx: Add SD to proportions
# ==============================================================================
cat("\n4. Fixing Fig3B_data.xlsx (adding SD)...\n")
fpath <- file.path(sd_dir, "Fig3B_data.xlsx")
wb <- loadWorkbook(fpath)
prop <- read.xlsx(wb, sheet = "proportions")
# SD = SEM * sqrt(n)
prop$sd <- prop$sem * sqrt(prop$total)
removeWorksheet(wb, "proportions")
addWorksheet(wb, "proportions", gridLines = FALSE)
writeData(wb, "proportions", prop)
worksheetOrder(wb) <- c(which(names(wb) == "proportions"), which(names(wb) != "proportions"))
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("  Added SD column\n")

cat("\n=== ALL FIXES COMPLETE ===\n")
