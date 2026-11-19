# PAPER PANEL: Fig 3E, 3F, 3M, 3N, 4D (violins) and 4H-4J (velocity bars): drawing
# PLOTS:  per-cell values with group mean and SEM
# CELLS:  visceral SMCs (n = 1,604: Young CTRL 469, Aged CTRL 477, Young P1KO 222, Aged P1KO 436); velocity panels: 1,368 SMCs
# GROUPS COMPARED: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Mann-Whitney U across cells, Bonferroni x4 (inputs from 75_ and the velocity source data)
#
# One consistent "superplot" format for Fig 3M, 3N, 4D, 4H, 4I, 4J.
# Bar = group mean over cells (the quantity the legend reports); error bar = SEM over cells;
# small dots = every cell; large outlined markers = per-mouse means (shape = mouse 1 / mouse 2).
# 4I is a per-cell yes/no (positive Piezo1 velocity), so per-cell dots carry no information there:
# that panel shows the pooled % with SEM across cells and per-mouse % markers only.
# Stars = the cell-level two-sided Mann-Whitney U, Bonferroni m=4, as in the figure legends
# (stats_cell_level.csv, written by 01_compute_per_cell.py). Nothing typed.
suppressPackageStartupMessages({library(ggplot2); library(ggbeeswarm); library(dplyr); library(patchwork)})
D <- commandArgs(TRUE)[1]  # data dir as argument
source("pipeline/utils/plotting_theme.R")  # run from the smc-analysis root
cols <- get_colors()$genotype; glab <- get_group_labels()
GROUPS <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
pc <- read.csv(file.path(D,"per_cell_values.csv")); pa <- read.csv(file.path(D,"per_animal_values.csv"))
st <- read.csv(file.path(D,"stats_cell_level.csv"))
prep <- function(d) { d$group <- factor(d$group, levels=GROUPS); d$genotype <- ifelse(grepl("KO", d$group), "KO", "Ctrl"); d }
pc <- prep(pc); pa <- prep(pa)
pa$mouse <- paste("Mouse", ave(as.character(pa$sample), pa$group, FUN = function(s) match(s, sort(unique(s)))))
META <- list(
  "3E"=list(t="Contractile score", y="Contractile capacity score\n(mean expression)"),
  "3F"=list(t="Synthetic score", y="Synthetic phenotype score\n(mean expression)"),
  "3M"=list(t="ECM score", y="ECM score\n(mean expression)"),
  "3N"=list(t="Cytoskeletal remodeling score", y="Cytoskeletal score\n(mean expression)"),
  "4D"=list(t="NFAT activity score", y="NFAT activity score\n(mean expression)"),
  "4H"=list(t="Contractile velocity index", y="Contractile genes with\npositive velocity (%)"),
  "4I"=list(t="Piezo1 upregulating", y="SMCs with positive\nPiezo1 velocity (%)"),
  "4J"=list(t="Contractile to synthetic", y="Transition probability (%)"),
  "S4C"=list(t="Myocd upregulating", y="SMCs with positive\nMyocd velocity (%)"),
  # 4D, 5-gene NFAT-associated remodeling score; files Fig4D_rec_points
  "4D_rec"=list(lab="4D", t="NFAT-associated remodeling", y="NFAT-associated remodeling\nscore (mean expression)"))
BIN <- c("4I","S4C")  # per-cell yes/no panels: no per-cell dots, SEM = sqrt(p(1-p)/n)
POS <- setNames(seq_along(GROUPS), GROUPS)
brackets <- function(s, ytop) {
  s <- s[s$sig != "ns", ]; if (!nrow(s)) return(NULL)
  s$x1 <- POS[s$group1]; s$x2 <- POS[s$group2]; s$w <- s$x2 - s$x1
  s <- s[order(s$w, s$x1), ]; s$y <- ytop * (1.04 + 0.085 * (seq_len(nrow(s)) - 1)); s
}
panel <- function(p) {
  m <- META[[p]]; d <- pc[pc$panel == p, ]; a <- pa[pa$panel == p, ]
  g <- d %>% group_by(group, genotype) %>% summarise(n=n(), mean=mean(value), sem=sd(value)/sqrt(n()), .groups="drop")
  # every bar is mean +/- SEM across cells, including 4I. For a yes/no
  # value the SEM is sqrt(p(1-p)/n).
  g$lo <- g$mean - g$sem; g$hi <- g$mean + g$sem
  ytop <- if (p %in% BIN) max(c(g$hi, a$animal_mean)) else max(d$value)
  b <- brackets(st[st$panel == p, ], ytop)
  pl <- ggplot(g, aes(group, mean)) +
    geom_col(aes(fill=genotype), width=0.72, alpha=0.35, colour="black", linewidth=0.4)
  if (!(p %in% BIN)) pl <- pl + geom_quasirandom(data=d, aes(group, value, colour=genotype), size=0.45, alpha=0.45, width=0.3, show.legend=FALSE)
  pl <- pl + geom_errorbar(aes(ymin=lo, ymax=hi), width=0.22, linewidth=0.6) +
    geom_point(data=a, aes(group, animal_mean, shape=mouse, fill=genotype), size=4.2, stroke=0.9, colour="black",
               position=position_dodge(width=0.35)) +
    scale_fill_manual(values=cols, breaks=c("Ctrl","KO"), labels=c("CTRL","P1KO"), name=NULL) +
    scale_colour_manual(values=c(Ctrl="#4D4D4D", KO="#C2185B"), guide="none") +
    scale_shape_manual(values=c("Mouse 1"=21, "Mouse 2"=24), name=NULL) +
    scale_x_discrete(labels=glab) + labs(title=paste0(if (is.null(m$lab)) p else m$lab, "  ", m$t), x="Age (months)", y=m$y) +
    theme_publication(base_size=13) + theme(legend.position="bottom", legend.key.size=unit(0.5,"cm"),
      plot.title=element_text(hjust=0, size=14))
  if (!is.null(b)) pl <- pl +
    geom_segment(data=b, aes(x=x1, xend=x2, y=y, yend=y), inherit.aes=FALSE, linewidth=0.4) +
    geom_text(data=b, aes(x=(x1+x2)/2, y=y, label=sig), inherit.aes=FALSE, vjust=-0.1, size=4)
  pl <- pl + guides(fill=guide_legend(override.aes=list(shape=22, size=5, alpha=1)))
  if (p %in% c("4H","4J",BIN)) pl + scale_y_continuous(breaks=seq(0,100,25), expand=expansion(mult=c(0.02, 0.08)))
  else pl + scale_y_continuous(expand=expansion(mult=c(0.02, 0.08)))
}
PAN <- c("3M","3N","4D","4H","4I","4J")
# optional panel list as 2nd arg (e.g. S4C, Myocd velocity, from 65_S4C_myocd_velocity.py); overview only for the default six
if (length(commandArgs(TRUE)) > 1) PAN <- strsplit(commandArgs(TRUE)[2], ",")[[1]]
for (p in PAN) { pl <- panel(p)
  ggsave(file.path(D, paste0("Fig", p, "_points.png")), pl, width=5.2, height=5.4, dpi=300, bg="white")
  ggsave(file.path(D, paste0("Fig", p, "_points.pdf")), pl, width=5.2, height=5.4, bg="white") }
if (identical(PAN, c("3M","3N","4D","4H","4I","4J"))) {
all <- wrap_plots(lapply(PAN, panel), ncol=3) + plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave(file.path(D, "P1KO_Fig3M3N4D4HIJ_points_overview.png"), all, width=15.5, height=11.5, dpi=200, bg="white")
}
cat("wrote", length(PAN), "panels + overview\n")
