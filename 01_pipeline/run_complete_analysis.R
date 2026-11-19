#!/usr/bin/env Rscript
# ==============================================================================
# COMPLETE ANALYSIS RUNNER
# ==============================================================================
# Single entry point for the single-cell figure modules and Tables S0-S12 of the manuscript.
# Module IDs (Fig3E, FigS3D, ...) are internal pipeline names; the README maps them to manuscript panels.
#
# Usage (from the repository root):
#   Rscript 01_pipeline/run_complete_analysis.R
#
# Output structure:
#   pipeline_results/<run_name>/
#     ├── <module ID>/          (one folder per figure module: PNG/PDF, *_data.xlsx, methods)
#     ├── Supplemental_Tables/  (S0-S7)
#     └── Tables_S8_S12/        (genome-wide tables)
#
# Requirements:
#   R 4.4+, Seurat, ggplot2, dplyr, tidyr, openxlsx, fgsea, msigdbr (see environment/R_sessionInfo.txt)
#
# Runtime: ~5-10 min
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("       COMPLETE ANALYSIS RUNNER\n")
cat("       Piezo1 Aging Manuscript — All Figures + Tables\n")
cat("================================================================\n")
overall_start_time <- Sys.time()
cat(sprintf("Started: %s\n", overall_start_time))
cat("\n")

# --- Validate working directory ---
if (!file.exists("pipeline/core/run_pipeline.R")) {
  stop("Working directory must be the smc-analysis project root.\n",
       "  Expected: pipeline/core/run_pipeline.R\n",
       "  Got: ", getwd())
}

if (!file.exists("data/seurat_paper_6samples.rds")) {
  stop("Seurat object not found at data/seurat_paper_6samples.rds")
}

# ==============================================================================
# PHASE 1: MAIN PIPELINE (figure modules + Tables S0-S7; module IDs are internal, see README)
# ==============================================================================
cat("================================================================\n")
cat("PHASE 1: Main Pipeline\n")
cat("================================================================\n\n")

source("pipeline/core/run_pipeline.R")

pipeline_result <- run_pipeline(
  selection_module = "selection_criteria/cluster_based_smc.R",
  figures = "all",
  remove_epithelial = TRUE,
  exclude_samples = character(0),
  depth_match = list(range_factor = 2.0)
)

pipeline_dir <- pipeline_result$output_dir
cat(sprintf("\nPhase 1 output: %s\n\n", pipeline_dir))

# ==============================================================================
# PHASE 2: GENOME-WIDE TABLES (Tables S8-S12)
# ==============================================================================
cat("================================================================\n")
cat("PHASE 2: Genome-wide tables (S8-S12)\n")
cat("================================================================\n\n")

tables_ext_dir <- file.path(pipeline_dir, "Tables_S8_S12")
dir.create(tables_ext_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(openxlsx)
})

# --- The visceral SMCs and their phenotype, exactly as analysed in Phase 1 (Fig 3A) ---
# Taken from the Phase 1 output (the cells, groups and phenotype calls analysed there).
sel <- read.xlsx(file.path(pipeline_dir, "Fig3E", "Fig3E_data.xlsx"), sheet = "phenotype_by_cell")
stopifnot(nrow(sel) > 0, !anyDuplicated(sel$cell_barcode)); cat("Phase 2 cells:", nrow(sel), "\n")

obj <- readRDS("data/seurat_paper_6samples.rds")
source("00_preprocessing/check_sample_identity.R"); check_sample_identity(obj)
stopifnot(all(sel$cell_barcode %in% colnames(obj)))
smc <- subset(obj, cells = sel$cell_barcode)
rm(obj)
smc$sample <- as.character(smc$sample)
smc$sample[is.na(smc$sample)] <- as.character(smc$orig.ident[is.na(smc$sample)])
smc$group <- paste0(
  ifelse(grepl("Young", smc$sample), "Young", "Aged"), "_",
  ifelse(grepl("KO", smc$sample), "KO", "Ctrl")
)
smc$phenotype <- sel$phenotype[match(colnames(smc), sel$cell_barcode)]

source("pipeline/core/gene_sets.R")
gene_sets <- get_gene_sets()
expr <- GetAssayData(smc, layer = "data")
smc$piezo1_status <- ifelse(as.numeric(expr["Piezo1", ]) > 0, "Piezo1+", "Piezo1-")

cat(sprintf("  Visceral SMCs: %d cells\n", ncol(smc)))
cat(sprintf("  Groups: %s\n", paste(names(table(smc$group)), table(smc$group),
                                      sep = "=", collapse = ", ")))
cat(sprintf("  Phenotype: %d contractile, %d synthetic\n",
            sum(smc$phenotype == "Contractile"), sum(smc$phenotype == "Synthetic")))
print(table(smc$group))
cat("\n")

# ==============================================================================
# Table S8: Contractile vs synthetic SMCs (genome-wide DE)
# ==============================================================================
cat("--- Table S8: Contractile vs Synthetic DE ---\n")

Idents(smc) <- "phenotype"
de_contr_synth <- FindMarkers(smc, ident.1 = "Contractile", ident.2 = "Synthetic",
                               test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0,
                               verbose = FALSE)
de_contr_synth$gene <- rownames(de_contr_synth)
write.csv(de_contr_synth %>% arrange(p_val_adj),
          file.path(tables_ext_dir, "Table_S8_Contractile_vs_Synthetic_DE.csv"),
          row.names = FALSE)

# ==============================================================================
# Tables S9-S10: Piezo1+ vs Piezo1- CTRL synthetic SMCs (DE + GSEA)
# ==============================================================================
cat("--- Tables S9-S10: Piezo1+ vs Piezo1- synthetic SMCs ---\n")

synth_ctrl <- subset(smc, cells = colnames(smc)[smc$phenotype == "Synthetic" &
                                                   smc$group %in% c("Young_Ctrl", "Aged_Ctrl")])
Idents(synth_ctrl) <- "piezo1_status"

n_p1pos <- sum(synth_ctrl$piezo1_status == "Piezo1+")
n_p1neg <- sum(synth_ctrl$piezo1_status == "Piezo1-")

if (n_p1pos >= 10 && n_p1neg >= 10) {
  de_p1 <- FindMarkers(synth_ctrl, ident.1 = "Piezo1+", ident.2 = "Piezo1-",
                        test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0,
                        verbose = FALSE)
  de_p1$gene <- rownames(de_p1)

  # Highlight NFAT and ECM genes (columns kept in Table S9)
  nfat_genes <- gene_sets$nfat_pathway
  ecm_genes <- c("Mmp2", "Mmp9", "Timp1", "Timp2", "Col1a1", "Fn1")
  de_p1$highlight <- case_when(
    de_p1$gene %in% nfat_genes ~ "NFAT",
    de_p1$gene %in% ecm_genes ~ "ECM",
    TRUE ~ "Other"
  )
  de_p1$label <- ifelse(de_p1$gene %in% c(nfat_genes, ecm_genes) & de_p1$p_val_adj < 0.1,
                          de_p1$gene, "")

  # Table S9
  write.csv(de_p1 %>% arrange(p_val_adj),
            file.path(tables_ext_dir, "Table_S9_P1pos_vs_P1neg_Synthetic_DE.csv"),
            row.names = FALSE)

  # Table S10: GSEA on the S9 ranking
  if (requireNamespace("fgsea", quietly = TRUE) &&
      requireNamespace("msigdbr", quietly = TRUE)) {
    suppressPackageStartupMessages({
      library(fgsea)
      library(msigdbr)
    })

    ranks <- setNames(de_p1$avg_log2FC, de_p1$gene)
    ranks <- sort(ranks[!is.na(ranks) & !duplicated(names(ranks))], decreasing = TRUE)

    gobp_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
    gobp_sets <- split(gobp_df$gene_symbol, gobp_df$gs_name)

    hallmark_df <- msigdbr(species = "Mus musculus", category = "H")
    hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

    gsea_res <- fgsea(pathways = c(gobp_sets, hallmark_sets), stats = ranks,
                       minSize = 15, maxSize = 300, nPermSimple = 10000)
    gsea_res <- as.data.frame(gsea_res[order(gsea_res$pval), ])
    gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, function(x) paste(x, collapse = ";"))

    write.csv(gsea_res, file.path(tables_ext_dir, "Table_S10_GSEA_P1pos_vs_P1neg.csv"),
              row.names = FALSE)
  }
} else {
  cat("  Skipped: insufficient Piezo1+ or Piezo1- synthetic cells\n")
}

# ==============================================================================
# Tables S11-S12: Genome-wide aging and rescue DE
# ==============================================================================
cat("--- Tables S11-S12: Genome-wide aging and rescue DE ---\n")

Idents(smc) <- "group"

# Table S11: Aging DE (Aged CTRL vs Young CTRL)
aging_de <- FindMarkers(smc, ident.1 = "Aged_Ctrl", ident.2 = "Young_Ctrl",
                         test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0,
                         verbose = FALSE)
aging_de$gene <- rownames(aging_de)
write.csv(aging_de %>% arrange(p_val_adj),
          file.path(tables_ext_dir, "Table_S11_Aging_DE_genome_wide.csv"), row.names = FALSE)
cat(sprintf("  S11: %d genes tested, %d significant\n",
            nrow(aging_de), sum(aging_de$p_val_adj < 0.05)))

# Table S12: Rescue DE (Aged P1KO vs Aged CTRL)
rescue_de <- FindMarkers(smc, ident.1 = "Aged_KO", ident.2 = "Aged_Ctrl",
                          test.use = "wilcox", min.pct = 0.05, logfc.threshold = 0,
                          verbose = FALSE)
rescue_de$gene <- rownames(rescue_de)
write.csv(rescue_de %>% arrange(p_val_adj),
          file.path(tables_ext_dir, "Table_S12_Rescue_DE_genome_wide.csv"), row.names = FALSE)
cat(sprintf("  S12: %d genes tested, %d significant\n",
            nrow(rescue_de), sum(rescue_de$p_val_adj < 0.05)))

cat("\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
end_time <- Sys.time()
total_duration <- as.numeric(difftime(end_time, overall_start_time, units = "mins"))

cat("\n")
cat("================================================================\n")
cat("       COMPLETE ANALYSIS — FINISHED\n")
cat("================================================================\n")
cat(sprintf("Total duration: %.1f minutes\n", total_duration))
cat(sprintf("Output: %s\n", pipeline_dir))
cat("\n")
cat("Phase 1 (Pipeline): figure modules, Tables S0-S7\n")
cat("Phase 2: Tables S8-S12\n")
cat("\n")
cat("Contents:\n")
cat(sprintf("  Main figures:          %d\n", length(list.files(pipeline_dir, pattern = "^Fig[34]", recursive = TRUE))))
cat(sprintf("  Supplementary:         %d\n", length(list.files(pipeline_dir, pattern = "^FigS[34]", recursive = TRUE))))
cat("================================================================\n")
