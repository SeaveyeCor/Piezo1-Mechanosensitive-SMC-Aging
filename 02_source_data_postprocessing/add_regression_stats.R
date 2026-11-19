setwd(".")
suppressPackageStartupMessages({ library(Seurat); library(openxlsx); library(dplyr) })
source("pipeline/selection_criteria/cluster_based_smc.R")
source("pipeline/core/gene_sets.R")

obj <- readRDS("data/seurat_paper_6samples.rds")
sel <- select_cells(obj)
smc <- sel$smc
meta <- obj@meta.data[smc, ]
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()

sd_dir <- "output/source_data"

# Get the analysed SMCs
fig3a <- read.xlsx(file.path(sd_dir, "Fig3A_phenotype_umap_data.xlsx"), sheet="phenotype_by_cell")
matched <- fig3a$cell_barcode
meta_m <- obj@meta.data[matched, ]
meta_m$group <- paste0(ifelse(meta_m$age=="Young","Young","Aged"), "_",
                        ifelse(meta_m$genotype=="Ctrl","Ctrl","KO"))
ctrl <- matched[meta_m$genotype == "Ctrl"]
meta_ctrl <- meta_m[ctrl, ]

# =====================================================================
# 1. FIX CORRELATION FILES: Add slope, intercept, R²
# =====================================================================
cat("=== FIXING CORRELATION FILES ===\n\n")

corr_files <- list(
  list(fig="Fig3C", file="Fig3C_piezo1_vs_contractile_data.xlsx",
       x_col="contractile_score", y_col="piezo1",
       x_label="Contractile Score", y_label="Piezo1 Expression"),
  list(fig="Fig3D", file="Fig3D_piezo1_vs_synthetic_data.xlsx",
       x_col="synthetic_score", y_col="piezo1",
       x_label="Synthetic Score", y_label="Piezo1 Expression"),
  list(fig="Fig4B", file="Fig4B_piezo1_vs_nfat_data.xlsx",
       x_col="nfat_score", y_col="piezo1",
       x_label="NFAT Pathway Score", y_label="Piezo1 Expression"),
  list(fig="Fig4C", file="Fig4C_synthetic_vs_nfat_data.xlsx",
       x_col="synthetic_score", y_col="nfat_score",
       x_label="Synthetic Score", y_label="NFAT Pathway Score")
)

for (info in corr_files) {
  fpath <- file.path(sd_dir, info$file)
  wb <- loadWorkbook(fpath)

  # Read cell-level data
  cl <- read.xlsx(wb, sheet="cell_level_data")

  # Compute slope, intercept, R² per age group
  cs <- read.xlsx(wb, sheet="correlation_summary")

  cs$slope <- NA
  cs$intercept <- NA
  cs$r_squared <- NA

  for (i in 1:nrow(cs)) {
    age_val <- cs$age[i]
    sub <- cl[cl$age == age_val, ]
    x <- sub[[info$x_col]]
    y <- sub[[info$y_col]]
    lm_fit <- lm(y ~ x)
    cs$slope[i] <- round(coef(lm_fit)[2], 6)
    cs$intercept[i] <- round(coef(lm_fit)[1], 6)
    cs$r_squared[i] <- round(summary(lm_fit)$r.squared, 6)
  }

  removeWorksheet(wb, "correlation_summary")
  addWorksheet(wb, "correlation_summary", gridLines=FALSE)
  writeData(wb, "correlation_summary", cs)

  # Reorder: title first, then correlation_summary
  ws <- names(wb)
  ti <- which(ws == "title")
  ci <- which(ws == "correlation_summary")
  others <- setdiff(seq_along(ws), c(ti, ci))
  worksheetOrder(wb) <- c(ti, ci, others)

  saveWorkbook(wb, fpath, overwrite=TRUE)
  cat(info$fig, ": added slope, intercept, R² to correlation_summary\n")
  for (i in 1:nrow(cs)) {
    cat("  ", cs$age[i], ": rho=", cs$rho[i], " slope=", cs$slope[i],
        " intercept=", cs$intercept[i], " R²=", cs$r_squared[i], "\n")
  }
  cat("\n")
}

# =====================================================================
# 2. FIX FIG3E: Contractile score as raw mean (not percentage)
# =====================================================================
cat("=== FIXING FIG3E CONTRACTILE SCORE ===\n\n")

contr_genes <- gs$contractile_core
contr_avail <- contr_genes[contr_genes %in% rownames(expr)]
contr_score <- colMeans(as.matrix(expr[contr_avail, matched, drop=FALSE]))

df <- data.frame(score = contr_score, group = meta_m$group)

stats <- df %>%
  group_by(group) %>%
  summarize(n=n(), mean=mean(score), sd=sd(score), sem=sd(score)/sqrt(n()), .groups="drop") %>%
  mutate(genotype = ifelse(grepl("Ctrl", group), "Ctrl", "KO"))

# Pairwise tests
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

fpath <- file.path(sd_dir, "Fig3E_contractile_score_data.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "title")
writeData(wb, "title", data.frame(
  Field=c("Figure","Title","Manuscript"),
  Value=c("Fig3E","Contractile Score","Maladaptive Piezo1 Mechanotransduction Drives Smooth Muscle Aging in the Gut")
))
addWorksheet(wb, "summary_statistics")
writeData(wb, "summary_statistics", stats)
addWorksheet(wb, "statistical_tests")
writeData(wb, "statistical_tests", tests)
addWorksheet(wb, "gene_panel")
writeData(wb, "gene_panel", data.frame(gene=contr_avail, gene_set="contractile_core"))
addWorksheet(wb, "cell_level_data")
writeData(wb, "cell_level_data", data.frame(
  cell=matched, score=contr_score, group=meta_m$group,
  genotype=ifelse(meta_m$genotype=="Ctrl","Ctrl","KO")
))
saveWorkbook(wb, fpath, overwrite=TRUE)

cat("Fig3E regenerated with raw mean expression (not percentage)\n")
cat("  Values:\n")
for (i in 1:nrow(stats)) {
  cat("  ", stats$group[i], ": mean=", round(stats$mean[i], 4),
      " SD=", round(stats$sd[i], 4), " SEM=", round(stats$sem[i], 4),
      " n=", stats$n[i], "\n")
}

cat("\nDone.\n")
