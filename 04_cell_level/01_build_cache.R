# 01_build_cache.R -- cache the bona fide visceral SMC set of the pipeline run.
# Cell set + group labels are taken from the pipeline output (Fig3H module), so no filter chain is re-derived here.
suppressMessages({library(Seurat); library(readxl)})
setwd(".")
OUT <- "output/cell_level"; dir.create(file.path(OUT, "data"), showWarnings=FALSE, recursive=TRUE)

ref <- as.data.frame(read_excel(
  "output/pipeline/Fig3H/Fig3H_data.xlsx",
  sheet = "cell_level_data"))
stopifnot(nrow(ref) > 0)
cat("target cells:", nrow(ref), "\n"); print(table(ref$group))

cat("loading full object ...\n")
obj <- readRDS("data/seurat_paper_6samples.rds")
cat("full object:", ncol(obj), "cells x", nrow(obj), "features\n")

missing <- setdiff(ref$cell, colnames(obj))
cat("target barcodes absent from object:", length(missing), "\n")
stopifnot(length(missing) == 0)

sub <- subset(obj, cells = ref$cell)
sub <- sub[, ref$cell]                        # enforce the reference ordering
stopifnot(identical(colnames(sub), ref$cell))

# group label comes from the pipeline output
sub$group <- factor(ref$group, levels = c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO"))
sub$genotype <- ref$genotype
sub$published_contractile_pct <- ref$score

cat("\ngroups in cached subset:\n"); print(table(sub$group))
cat("assay:", DefaultAssay(sub), " layers:", paste(Layers(sub), collapse=","), "\n")
saveRDS(sub, file.path(OUT, "data", "smc_cache.rds"))
cat("saved", file.path(OUT,"data","smc_cache.rds"), "\n")

# --- how many cells would the non-bona-fide (Cl6 + Cl13 retained) CTRL set have? ---
# supplement quotes Young CTRL n=773 / Aged CTRL n=846; neither pipeline run gives that.
md <- obj@meta.data
md$sample <- as.character(md$sample)
md$sample[is.na(md$sample)] <- as.character(md$orig.ident[is.na(md$sample)])
cat("\n--- provenance check for supplement's 773 / 846 ---\n")
if ("seurat_clusters" %in% names(md)) {
  ctrl <- md[grepl("WT|Ctrl|CTRL", md$orig.ident) | md$condition %in% c("Ctrl","CTRL","WT"), ]
  cat("clusters present:", paste(sort(unique(as.character(md$seurat_clusters))), collapse=","), "\n")
  for (cl in list(c("1","5","15","16"), c("1","5","6","13","15","16"))) {
    s <- ctrl[as.character(ctrl$seurat_clusters) %in% cl, ]
    cat("clusters {", paste(cl, collapse=","), "} CTRL cells by age:\n"); print(table(s$age))
  }
}
