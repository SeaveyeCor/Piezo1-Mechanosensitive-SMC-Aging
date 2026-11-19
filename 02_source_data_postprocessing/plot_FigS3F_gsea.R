# PAPER PANEL: Fig S3F
# PLOTS:  dot plot of the 12 GO pathways in output/source_data/FigS3F_data.xlsx (sheet go_enrichment):
#         x = -log10(adjusted p), dot size = gene count, red = higher in P1KO (rescue), blue = higher in CTRL (aging)
# CELLS:  Aged CTRL and Aged P1KO visceral SMCs
# GROUPS COMPARED: Aged P1KO vs Aged CTRL
# TEST:   GSEA (fgsea); adjusted p
suppressPackageStartupMessages({ library(ggplot2); library(openxlsx) })
d <- read.xlsx("output/source_data/FigS3F_data.xlsx", sheet = "go_enrichment")
d$mlog10 <- -log10(d$padj)
d$pathway_clean <- factor(d$pathway_clean, levels = rev(d$pathway_clean))
d$direction <- factor(d$direction, levels = c("Higher in KO (rescue)", "Higher in Ctrl (aging)"))
p <- ggplot(d, aes(mlog10, pathway_clean, size = size, colour = direction)) +
  geom_point() +
  scale_size_continuous(name = "Gene Count", breaks = c(50, 100, 150, 200), range = c(2, 9)) +
  scale_colour_manual(values = c("Higher in KO (rescue)" = "#E8212B", "Higher in Ctrl (aging)" = "#2B7BBA"), name = NULL) +
  labs(x = expression(-log[10]~"(adjusted p value)"), y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(colour = "black", size = 13), legend.text = element_text(size = 13),
        legend.title = element_text(size = 14), plot.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_legend(override.aes = list(size = 5), order = 2), size = guide_legend(order = 1))
out <- "output/figures/FigS3F"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out, "FigS3F_gsea.png"), p, width = 10, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out, "FigS3F_gsea.pdf"), p, width = 10, height = 5.5, bg = "white", device = cairo_pdf)
cat("Fig S3F saved to", out, "\n")
