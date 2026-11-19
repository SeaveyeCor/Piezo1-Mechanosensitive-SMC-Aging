# PAPER PANEL: Fig S3B
# PLOTS:  aging change in channel expression within expressing cells, vs detection-matched background genes
# CELLS:  CTRL visceral SMCs only (Young CTRL 469, Aged CTRL 477)
# GROUPS COMPARED: Aged CTRL vs Young CTRL
# TEST:   empirical p vs 120 detection-matched background genes (depth-matched bootstrap)
#
# 48: single-cell violins for the conditional result:
# violin + every cell as a quasirandom dot + crossbar, house
# palette, theme_publication(), house brackets. Expressing control cells only.
# Crossbar = MEAN (the quantity the background test is computed on), stated in caption.
suppressMessages({library(Seurat); library(ggplot2); library(dplyr); library(ggbeeswarm)})
setwd("."); set.seed(1)
source("pipeline/utils/plotting_theme.R"); source("pipeline/utils/statistics.R")
OUT <- "output/cell_level"; DIR <- "output/cell_level/channels_aging"
o <- readRDS(file.path(OUT,"data","smc_cache.rds")); DefaultAssay(o) <- "RNA"
cnt <- GetAssayData(o,layer="counts"); grp <- as.character(o$group)
smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)]); tot <- Matrix::colSums(cnt)
res <- read.csv(file.path(OUT,"data","conditional_positive_channels.csv"))
CH <- res$gene[order(-res$log2FC_pos)]
colors <- get_colors(); dot_cols <- c(Ctrl="#3F3F3F", KO="#c2306b")
stars <- function(p) ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*","ns")))

cells <- do.call(rbind, lapply(CH, function(g){ x <- as.numeric(cnt[g,])
  keep <- grp %in% c("Young_Ctrl","Aged_Ctrl") & x>0
  data.frame(gene=g, group=factor(grp[keep], levels=c("Young_Ctrl","Aged_Ctrl")),
             genotype="Ctrl", counts=x[keep], animal=smp[keep],
             cell_barcode=colnames(cnt)[keep], total_UMI=tot[keep]) }))
cells$gene <- factor(cells$gene, levels=CH)
cat("expressing control cells per channel:\n"); print(table(cells$gene, cells$group))

violin <- function(d, title, ylab, base=18, dot=1.5) {
  ggplot(d, aes(group, counts, fill=genotype)) +
    geom_violin(colour="black", linewidth=0.35, alpha=0.32, scale="width", width=0.80, adjust=1.6) +
    ggbeeswarm::geom_quasirandom(aes(colour=genotype), width=0.33, size=dot, alpha=0.75,
                                 stroke=0, show.legend=FALSE) +
    stat_summary(fun=mean, geom="crossbar", width=0.5, linewidth=0.4, colour="black", fatten=1.4) +
    scale_fill_manual(values=colors$genotype, guide="none") +
    scale_colour_manual(values=dot_cols) +
    scale_x_discrete(labels=c(Young_Ctrl="3", Aged_Ctrl="24")) +
    scale_y_continuous(trans="log2", breaks=c(1,2,4,8,16,32)) +
    labs(title=title, x="Age (months)", y=ylab) +
    theme_publication(base_size=base) + theme(legend.position="none")
}

## Piezo1 standalone, 7x6 like Fig3H
d1 <- cells[cells$gene=="Piezo1",]; ytop <- max(d1$counts)
tests <- data.frame(a="Young_Ctrl", b="Aged_Ctrl",
                    p_adjusted=res$p_emp[res$gene=="Piezo1"],
                    significance=stars(res$p_emp[res$gene=="Piezo1"]))
p1 <- violin(d1, "Piezo1 in Expressing SMCs", "Piezo1 UMI per cell\n(expressing cells, log2 axis)") +
  expand_limits(y=ytop*2.6)
p1 <- add_pairwise_brackets(p1, tests, c("Young_Ctrl","Aged_Ctrl"), ytop*1.15, step=ytop*0.9)
ggsave(file.path(DIR,"figures","Fig_S3B_Piezo1_expressing_cells_violin.png"), p1, width=7, height=6, dpi=300, bg="white")
ggsave(file.path(DIR,"figures","Fig_S3B_Piezo1_expressing_cells_violin.pdf"), p1, width=7, height=6, bg="white")

## all nine, faceted; star per facet from the background test (ns omitted, house rule)
ann <- data.frame(gene=factor(res$gene, levels=CH), lab=stars(res$p_emp))
ann <- ann[ann$lab!="ns",]
ytops <- cells %>% group_by(gene) %>% summarise(y=max(counts)*1.6, .groups="drop")
ann <- merge(ann, ytops, by="gene")
pall <- violin(cells, "Channel Expression in Expressing SMCs", "UMI per cell (expressing cells, log2 axis)",
               base=13, dot=1.0) +
  geom_text(data=ann, aes(x=1.5, y=y, label=lab), inherit.aes=FALSE, size=6, fontface="bold") +
  facet_wrap(~gene, ncol=3, scales="free_y") +
  theme(strip.background=element_blank(), strip.text=element_text(face="bold", size=13))
ggsave(file.path(DIR,"figures","Fig_S3B_all_channels_expressing_cells_violins.png"), pall, width=11, height=10, dpi=300, bg="white")
ggsave(file.path(DIR,"figures","Fig_S3B_all_channels_expressing_cells_violins.pdf"), pall, width=11, height=10, bg="white")
write.csv(cells, file.path(DIR,"data","expressing_cells_counts_long.csv"), row.names=FALSE)
cat("wrote violins\n")
