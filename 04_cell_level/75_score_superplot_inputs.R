# PAPER PANEL: Fig 3E, 3F, 3M, 3N, 4D: per-cell scores and statistics
# PLOTS:  contractile, synthetic, ECM, cytoskeletal and NFAT-associated remodeling scores
# CELLS:  visceral SMCs (n = 1,604: Young CTRL 469, Aged CTRL 477, Young P1KO 222, Aged P1KO 436)
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Mann-Whitney U across cells, Bonferroni x4
#
# 75: superplot inputs for 3E (contractile), 3F (synthetic), 3M (ECM), 3N (cytoskeletal) from per_cell_scores.csv, and
# 4D from 71's output; same statistics as the figure legends (two-sided Mann-Whitney U across cells, Bonferroni m = 4).
# Usage (from the repository root, after 01_build_cache.R, 02_scores_singlecell.R and 71_score_panels.R):
#   Rscript 04_cell_level/75_score_superplot_inputs.R <out dir>
suppressMessages(library(Seurat))
OUT <- commandArgs(TRUE)[1]; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sc <- read.csv("output/cell_level/data/per_cell_scores.csv")
o <- readRDS("output/cell_level/data/smc_cache.rds")
smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)]); names(smp) <- colnames(o)
sc$sample <- smp[sc$cell]; stopifnot(!any(is.na(sc$sample)))
P <- c("3E" = "contractile_capacity", "3F" = "synthetic_supplement", "3M" = "ecm_structural", "3N" = "cytoskeletal_remod")
stopifnot(all(P %in% colnames(sc)))
pc <- do.call(rbind, lapply(names(P), function(p) data.frame(panel = p, metric = P[[p]], cell = sc$cell, sample = sc$sample, group = sc$group, value = sc[[P[[p]]]])))
d4 <- read.csv("output/cell_level/data/score_panels/fig4D_superplot_inputs/per_cell_values.csv")
pc <- rbind(pc, d4[, c("panel", "metric", "cell", "sample", "group", "value")])
write.csv(pc, file.path(OUT, "per_cell_values.csv"), row.names = FALSE)
pa <- aggregate(value ~ panel + metric + group + sample, pc, function(v) c(m = mean(v), n = length(v)))
pa <- data.frame(pa[, 1:4], animal_mean = pa$value[, "m"], n_cells = pa$value[, "n"])
write.csv(pa, file.path(OUT, "per_animal_values.csv"), row.names = FALSE)
CMP <- list(c("Young_Ctrl","Aged_Ctrl"), c("Young_KO","Aged_KO"), c("Young_Ctrl","Young_KO"), c("Aged_Ctrl","Aged_KO"))
st <- do.call(rbind, lapply(unique(pc$panel), function(p) do.call(rbind, lapply(CMP, function(cc) {
  x <- pc$value[pc$panel == p & pc$group == cc[1]]; y <- pc$value[pc$panel == p & pc$group == cc[2]]
  pv <- wilcox.test(x, y)$p.value; pb <- min(1, 4 * pv)
  data.frame(panel = p, group1 = cc[1], group2 = cc[2], n1 = length(x), n2 = length(y), p = pv, p_bonf = pb,
             sig = ifelse(pb < 1e-4, "****", ifelse(pb < 1e-3, "***", ifelse(pb < .01, "**", ifelse(pb < .05, "*", "ns"))))) }))))
write.csv(st, file.path(OUT, "stats_cell_level.csv"), row.names = FALSE)
print(aggregate(value ~ panel + group, pc, mean)); print(st[st$group1 == "Aged_Ctrl", ]); print(pa[grepl("Aged", pa$group), ])
