setwd(".")
suppressPackageStartupMessages({ library(Seurat); library(openxlsx); library(dplyr) })

sd <- "output/source_data"

cat("Loading Seurat...\n")
obj <- readRDS("data/seurat_paper_6samples.rds")
source("pipeline/selection_criteria/cluster_based_smc.R")
source("pipeline/core/gene_sets.R")
sel <- select_cells(obj)

# Get the analysed SMCs from Fig3A
fig3a <- read.xlsx(file.path(sd, "Fig3A_phenotype_umap_data.xlsx"), sheet = "phenotype_by_cell")
matched_cells <- fig3a$cell_barcode
meta <- obj@meta.data[matched_cells, ]
meta$group <- paste0(
  ifelse(meta$age == "Young", "Young", "Aged"), "_",
  ifelse(meta$genotype == "Ctrl", "Ctrl", "KO")
)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()

# =====================================================================
# 1. Fig 3M (ECM score): add cell_level_data sheet
# =====================================================================
cat("\n1. Fig 3M: adding cell_level_data...\n")
ecm_genes <- gs$ecm_structural
ecm_avail <- ecm_genes[ecm_genes %in% rownames(expr)]
ecm_score <- colMeans(as.matrix(expr[ecm_avail, matched_cells, drop = FALSE]))

fpath <- file.path(sd, "Fig3M_ecm_score_data.xlsx")
wb <- loadWorkbook(fpath)
if ("cell_level_data" %in% names(wb)) removeWorksheet(wb, "cell_level_data")
addWorksheet(wb, "cell_level_data")
writeData(wb, "cell_level_data", data.frame(
  cell_barcode = matched_cells,
  ecm_score = ecm_score,
  group = meta$group
))
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("  Added cell_level_data (", length(matched_cells), " cells)\n")

# =====================================================================
# 2. FIX FigS3A: Add cell_type, group, sample to UMAP coordinates
# =====================================================================
cat("\n2. Fixing FigS3A (adding cell_type, group, sample)...\n")
fpath <- file.path(sd, "FigS3A_celltype_umap_data.xlsx")
wb <- loadWorkbook(fpath)
umap <- read.xlsx(wb, sheet = "umap_coordinates")

# Get metadata for these cells
all_meta <- obj@meta.data
umap_cells <- umap$cell_barcode

# Add columns from metadata
if ("cell_type" %in% colnames(all_meta)) {
  umap$cell_type <- all_meta[umap_cells, "cell_type"]
} else if ("cell_type_label" %in% colnames(all_meta)) {
  umap$cell_type <- all_meta[umap_cells, "cell_type_label"]
}

umap$group <- paste0(
  ifelse(all_meta[umap_cells, "age"] == "Young", "Young", "Aged"), "_",
  ifelse(all_meta[umap_cells, "genotype"] == "Ctrl", "Ctrl", "KO")
)

# Try to get sample ID
for (col in c("orig.ident", "sample", "Sample")) {
  if (col %in% colnames(all_meta)) {
    umap$sample <- all_meta[umap_cells, col]
    break
  }
}

removeWorksheet(wb, "umap_coordinates")
addWorksheet(wb, "umap_coordinates")
writeData(wb, "umap_coordinates", umap)
# Keep as first data sheet
ws_order <- seq_along(names(wb))
ui <- which(names(wb) == "umap_coordinates")
ti <- which(names(wb) == "title")
others <- setdiff(ws_order, c(ti, ui))
worksheetOrder(wb) <- c(ti, ui, others)
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("  Added cell_type, group, sample columns (", nrow(umap), " cells)\n")
cat("  Cell types found: ", paste(unique(umap$cell_type), collapse = ", "), "\n")

# =====================================================================
# 3. FIX Fig3B: Fix the n column (garbage number)
# =====================================================================
cat("\n3. Fixing Fig3B (n column)...\n")
fpath <- file.path(sd, "Fig3B_phenotype_proportions_data.xlsx")
wb <- loadWorkbook(fpath)
prop <- read.xlsx(wb, sheet = "proportions")
cat("  Before fix - n values: ", paste(prop$n, collapse = ", "), "\n")

# The n should be the count of cells with that phenotype in that group
# Recalculate from total and pct
prop$n <- round(prop$pct / 100 * prop$total)

cat("  After fix  - n values: ", paste(prop$n, collapse = ", "), "\n")

removeWorksheet(wb, "proportions")
addWorksheet(wb, "proportions")
writeData(wb, "proportions", prop)
ws_order <- seq_along(names(wb))
pi <- which(names(wb) == "proportions")
ti <- which(names(wb) == "title")
others <- setdiff(ws_order, c(ti, pi))
worksheetOrder(wb) <- c(ti, pi, others)
saveWorkbook(wb, fpath, overwrite = TRUE)
cat("  Fixed n column\n")

cat("\nAll fixes done.\n")
