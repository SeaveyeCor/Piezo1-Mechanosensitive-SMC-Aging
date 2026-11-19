# 37: cluster label per cell of the SMC cache (output/cell_level/data/clusters.rds), read by 52, 52b and 53
# (Fig S3E). Every cell in the deposited object has a stored cluster label.
# Usage (from the repository root, after 01_build_cache.R): Rscript 04_cell_level/37_cluster_labels.R
suppressMessages(library(Seurat))
o <- readRDS("output/cell_level/data/smc_cache.rds")
cl <- as.character(o$seurat_clusters); stopifnot(!any(is.na(cl)))
saveRDS(cl, "output/cell_level/data/clusters.rds"); cat(length(cl), "cluster labels saved\n")
