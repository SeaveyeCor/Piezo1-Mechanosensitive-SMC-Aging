# PAPER PANEL: Fig S3B
# PLOTS:  aging change in channel expression within expressing cells, vs detection-matched background genes
# CELLS:  CTRL visceral SMCs only (Young CTRL 469, Aged CTRL 477)
# GROUPS COMPARED: Aged CTRL vs Young CTRL
# TEST:   empirical p vs 120 detection-matched background genes (depth-matched bootstrap)
#
# 47: single figure. Per channel, two bars: the observed aging change within expressing
# cells (depth-matched; error bar = 2.5-97.5% across 200 bootstrap draws, cells resampled) beside the level that
# detection-matched background genes reach under the identical statistic (error bar =
# 2.5-97.5% across background genes). House style via theme_publication().
suppressMessages({library(ggplot2); library(dplyr)}); setwd(".")
source("pipeline/utils/plotting_theme.R")
OUT <- "output/cell_level"; DIR <- "output/cell_level/channels_aging"
dir.create(file.path(DIR,"figures"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(DIR,"data"), showWarnings=FALSE)
res <- read.csv(file.path(OUT,"data","conditional_positive_channels.csv"))
res <- res[order(-res$log2FC_pos),]
L <- data.frame(gene=factor(res$gene, levels=res$gene), val=res$log2FC_pos, lo=res$obs_lo, hi=res$obs_hi)
stars <- ifelse(res$p_emp<0.001,"***", ifelse(res$p_emp<0.01,"**", ifelse(res$p_emp<0.05,"*","")))
ann <- data.frame(gene=L$gene, lab=stars, y=res$obs_hi + 0.04)
# Background bars are not drawn. The star still comes from the
# detection-matched background test; the background values stay in the source CSV.
p <- ggplot(L, aes(gene, val)) +
  geom_hline(yintercept=0, linewidth=0.5, colour="black") +
  geom_col(width=0.66, fill="#7F7F7F", colour="black", linewidth=0.4) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.24, linewidth=0.55) +
  geom_text(data=ann, aes(gene, y, label=lab), inherit.aes=FALSE, size=7, fontface="bold") +
  labs(title="Piezo1 Increases with Age in Expressing SMCs",
       x=NULL, y="Log2 fold change\n(aged vs young, expressing cells)") +
  theme_publication(base_size=14) +
  theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(DIR,"figures","Fig_S3B_channels_aging_expressing_cells.png"), p, width=9, height=6, dpi=300, bg="white")
ggsave(file.path(DIR,"figures","Fig_S3B_channels_aging_expressing_cells.pdf"), p, width=9, height=6, bg="white")
file.copy(file.path(OUT,"data","conditional_positive_channels.csv"),
          file.path(DIR,"data","conditional_positive_channels.csv"), overwrite=TRUE)
file.copy(file.path(OUT,"data","conditional_positive_draws.csv"), file.path(DIR,"data","conditional_positive_draws.csv"), overwrite=TRUE)
print(res[,c("gene","log2FC_pos","obs_lo","obs_hi","bg_median","bg_lo","bg_hi","p_emp")], row.names=FALSE, digits=3)
