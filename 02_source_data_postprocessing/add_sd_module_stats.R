setwd(".")
suppressPackageStartupMessages({ library(Seurat); library(openxlsx); library(dplyr) })
source("pipeline/core/gene_sets.R")
source("pipeline/selection_criteria/cluster_based_smc.R")

obj <- readRDS("data/seurat_paper_6samples.rds")
sel <- select_cells(obj)
smc <- sel$smc
meta <- obj@meta.data[smc, ]
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()

# Classify
contr_sc <- colMeans(as.matrix(expr[gs$contractile_core[gs$contractile_core %in% rownames(expr)], smc, drop=FALSE]))
synth_sc <- colMeans(as.matrix(expr[gs$synthetic_core[gs$synthetic_core %in% rownames(expr)], smc, drop=FALSE]))
ratio <- (contr_sc + 0.001) / (synth_sc + 0.001)
ctrl_synth <- smc[meta$genotype == "Ctrl" & ratio <= 1.5]
p1x <- as.numeric(expr["Piezo1", ctrl_synth])
p1_pos <- ctrl_synth[p1x > 0]
p1_neg <- ctrl_synth[p1x == 0]

# Fix original Fig4F: add SD to gene_level_stats and module_level_stats
sd_dir <- "output/source_data"
fpath <- file.path(sd_dir, "FigS4B_p1pos_vs_p1neg_data.xlsx")
wb <- loadWorkbook(fpath)

# Fix gene_level_stats
gl <- read.xlsx(wb, sheet="gene_level_stats")
for (i in 1:nrow(gl)) {
  g <- gl$gene[i]
  if (!(g %in% rownames(expr))) next
  gp <- as.numeric(expr[g, p1_pos])
  gn <- as.numeric(expr[g, p1_neg])
  gl$sd_p1pos[i] <- round(sd(gp), 4)
  gl$sd_p1neg[i] <- round(sd(gn), 4)
}
removeWorksheet(wb, "gene_level_stats")
addWorksheet(wb, "gene_level_stats")
writeData(wb, "gene_level_stats", gl)

# Fix module_level_stats
ml <- read.xlsx(wb, sheet="module_level_stats")
ml$sd_pos <- NA
ml$sd_neg <- NA
# Recompute from genes
panels <- list(
  "Synthetic Markers" = c("Klf4", "Klf5", "Vim", "Fn1", "Col1a1", "Spp1", "Lgals3", "Cdkn1a", "Icam1", "Vcan", "S100a4"),
  "Cyto Remodeling" = gs$cyto_remodeling,
  "NFAT / Ca2+" = c("Nfatc1","Nfatc2","Nfatc4","Ppp3ca","Calm1","Itpr2","Itpr3","Trpm7","Rcan1"),
  "ECM / Remodeling" = c("Bgn","Col6a1","Dcn","Sparc","Timp2","Mmp2","Mmp9","Timp1"),
  "Contractile" = gs$contractile_core
)

for (j in 1:nrow(ml)) {
  panel_name <- ml$panel[j]
  if (panel_name %in% names(panels)) {
    genes <- panels[[panel_name]]
    avail <- genes[genes %in% rownames(expr)]
    if (length(avail) >= 2) {
      sp <- colMeans(as.matrix(expr[avail, p1_pos, drop=FALSE]))
      sn <- colMeans(as.matrix(expr[avail, p1_neg, drop=FALSE]))
      ml$sd_pos[j] <- round(sd(sp), 4)
      ml$sd_neg[j] <- round(sd(sn), 4)
    }
  }
}

removeWorksheet(wb, "module_level_stats")
addWorksheet(wb, "module_level_stats")
writeData(wb, "module_level_stats", ml)

# Reorder
ws <- names(wb)
ti <- which(ws == "title")
others <- setdiff(seq_along(ws), ti)
worksheetOrder(wb) <- c(ti, others)

saveWorkbook(wb, fpath, overwrite=TRUE)
cat("Fig S4B: added SD to gene_level_stats and module_level_stats\n")
