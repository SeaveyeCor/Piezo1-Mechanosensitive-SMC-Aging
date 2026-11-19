setwd(".")
suppressPackageStartupMessages({ library(Seurat); library(openxlsx); library(dplyr) })
source("pipeline/core/gene_sets.R")
source("pipeline/core/run_pipeline.R")
source("pipeline/selection_criteria/cluster_based_smc.R")

config <- get_config()
cat("Loading Seurat...\n")
obj <- readRDS(config$seurat_path)
sel <- select_cells(obj)
smc <- sel$smc

meta <- obj@meta.data[smc, ]
keep <- rep(TRUE, nrow(meta))
cat("Before depth filter:", length(smc), "| After:", sum(keep), "\n")

# Use the analysed SMCs (4 groups)
meta$group <- paste0(
  ifelse(meta$age == "Young", "Young", "Aged"), "_",
  ifelse(meta$genotype == "Ctrl", "Ctrl", "KO")
)
cat("Group counts (all selected):\n")
print(table(meta$group))
cat("Total:", nrow(meta), "\n")

# The pipeline's analysed SMC set
# For Fig3F (bar plot), use same cells as other bar plots
# Check what Fig3E uses
sd_dir <- "output/source_data"
fig3e <- read.xlsx(file.path(sd_dir, "Fig3E_data.xlsx"), sheet="summary_statistics")
cat("\nFig3E cell counts:", paste0(fig3e$group, "=", fig3e$n, collapse=", "), "\n")
target_n <- sum(fig3e$n)
cat("Target total:", target_n, "\n")

# The analysed cells are the ones from the pipeline phenotype classification
# Read Fig3A to get the exact cell barcodes
fig3a <- read.xlsx(file.path(sd_dir, "Fig3A_data.xlsx"), sheet="phenotype_by_cell")
matched_cells <- fig3a$cell_barcode
cat("Fig3A cells:", length(matched_cells), "\n")

# Regenerate Fig3F with matched cells
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()
synth_genes <- gs$synthetic_phenotype
synth_avail <- synth_genes[synth_genes %in% rownames(expr)]

score <- colMeans(as.matrix(expr[synth_avail, matched_cells, drop = FALSE]))
meta_matched <- obj@meta.data[matched_cells, ]
meta_matched$group <- paste0(
  ifelse(meta_matched$age == "Young", "Young", "Aged"), "_",
  ifelse(meta_matched$genotype == "Ctrl", "Ctrl", "KO")
)

df <- data.frame(score = score, group = meta_matched$group)
stats <- df %>%
  group_by(group) %>%
  summarize(n=n(), mean=mean(score), sd=sd(score), sem=sd(score)/sqrt(n()), .groups="drop") %>%
  mutate(genotype = ifelse(grepl("Ctrl", group), "Ctrl", "KO"))

comps <- list(
  c("Young_Ctrl","Aged_Ctrl","Aging (WT)"),
  c("Young_KO","Aged_KO","Aging (KO)"),
  c("Young_Ctrl","Young_KO","KO effect (young)"),
  c("Aged_Ctrl","Aged_KO","KO rescue (aged)")
)
tests <- data.frame()
for (comp in comps) {
  v1 <- df$score[df$group == comp[1]]
  v2 <- df$score[df$group == comp[2]]
  wt <- wilcox.test(v1, v2)
  p_adj <- min(wt$p.value * 4, 1)
  sig <- ifelse(p_adj<0.001,"***",ifelse(p_adj<0.01,"**",ifelse(p_adj<0.05,"*","ns")))
  tests <- rbind(tests, data.frame(
    comparison=comp[3], a=comp[1], b=comp[2],
    n_a=length(v1), n_b=length(v2),
    mean_a=mean(v1), sd_a=sd(v1), mean_b=mean(v2), sd_b=sd(v2),
    statistic=unname(wt$statistic), p_value=wt$p.value, p_adjusted=p_adj, significance=sig
  ))
}

wb <- createWorkbook()
addWorksheet(wb, "title")
writeData(wb, "title", data.frame(
  Field=c("Figure","Title","Manuscript"),
  Value=c("Fig3F","Synthetic Score","Maladaptive Piezo1 Mechanotransduction Drives Smooth Muscle Aging in the Gut")
))
addWorksheet(wb, "summary_statistics")
writeData(wb, "summary_statistics", stats)
addWorksheet(wb, "statistical_tests")
writeData(wb, "statistical_tests", tests)
addWorksheet(wb, "cell_level_data")
writeData(wb, "cell_level_data", data.frame(cell_barcode=names(score), group=meta_matched$group, synthetic_score=score))
addWorksheet(wb, "genes_used")
writeData(wb, "genes_used", data.frame(gene=synth_avail))
saveWorkbook(wb, file.path(sd_dir, "Fig3F_data.xlsx"), overwrite=TRUE)
cat("\nFig3F_data.xlsx regenerated with n=", length(matched_cells), "\n")
cat("Group counts:", paste0(stats$group, "=", stats$n, collapse=", "), "\n")

# Fix Fig4A: rename "ECM Remodeling" -> "ECM Score" in correlation_summary
cat("\nFixing Fig4A pathway names...\n")
fpath <- file.path(sd_dir, "Fig4A_data.xlsx")
wb <- loadWorkbook(fpath)
for (sheet in c("correlation_summary", "plot_data")) {
  df <- read.xlsx(wb, sheet=sheet)
  df$pathway <- gsub("ECM Remodeling", "ECM Score", df$pathway)
  removeWorksheet(wb, sheet)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df)
}
saveWorkbook(wb, fpath, overwrite=TRUE)
cat("  Renamed ECM Remodeling -> ECM Score in Fig4A_data.xlsx\n")

cat("\nDone.\n")
