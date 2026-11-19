# PAPER PANEL: Fig 3D, 4A, 4B, 4C (correlations) and the Fig 4D score
# PLOTS:  score-score and Piezo1-score correlations; NFAT-associated remodeling score
# CELLS:  correlations: CTRL visceral SMCs only (469 young, 477 aged); 4D: visceral SMCs (n = 1,604: Young CTRL 469, Aged CTRL 477, Young P1KO 222, Aged P1KO 436)
# GROUPS COMPARED: correlations: Young CTRL and Aged CTRL separately; 4D: Aged vs Young CTRL; Aged vs Young P1KO; P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   Spearman correlation; 4D: Mann-Whitney U, Bonferroni x4
#
# 71: score panels 3D, 4A, 4B, 4C, 4D.
#  - 4D: NFAT-associated remodeling score (5 genes). Superplot inputs for 02_superplots.R (panel "4D_rec").
#  - 4B, 4C: NFAT pathway score without Nfat5 (calcineurin-independent), 12 genes.
#  - 3D, 4C: synthetic phenotype score = the SI Methods set (Vim, Col1a2, Dcn, Igfbp6, Efemp1), as in 3F.
#  - 4A: Spearman rho of Piezo1 vs each pathway score.
# The alternative gene sets (13-gene NFAT pathway; Klf4, Vim, Mmp2, Mmp9, Timp1 synthetic set) are computed alongside
# for comparison. Scores = mean log-normalized RNA data (zeros included). CTRL cells only for correlations (469 / 477).
# Usage (from .): Rscript 04_cell_level/71_score_panels.R
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(patchwork)})
source("pipeline/utils/plotting_theme.R"); source("pipeline/core/gene_sets.R"); GS <- get_gene_sets()
OUT <- "output/cell_level/data/score_panels"; FIG <- "output/cell_level/panels_publication/score_panels"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
o <- readRDS("output/cell_level/data/smc_cache.rds"); x <- GetAssayData(o, assay="RNA", layer="data")
grp <- as.character(o$group); smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)])
sc <- function(g, k=TRUE) colMeans(as.matrix(x[g, k, drop=FALSE]))
SYN_ALT <- c("Klf4","Vim","Mmp2","Mmp9","Timp1"); SYN_NEW <- GS$synthetic_phenotype
NFATP13 <- GS$nfat_pathway; NFATP12 <- setdiff(NFATP13, "Nfat5")
ACT6 <- GS$nfatc3_remodeling; ACT5 <- setdiff(ACT6, "Fas")
stopifnot(identical(SYN_NEW, c("Vim","Col1a2","Dcn","Igfbp6","Efemp1")), length(ACT6) == 6)
star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "ns")))
ctrl <- grp %in% c("Young_Ctrl","Aged_Ctrl"); age <- ifelse(grp == "Young_Ctrl", "Young", "Aged")
sp <- function(a, b, k) { t <- cor.test(a[k], b[k], method="spearman", exact=FALSE); c(rho=unname(t$estimate), p=t$p.value, n=sum(k)) }

# ---------- correlations: alternative vs final gene sets ----------
p1 <- x["Piezo1", ]
COR <- list(
  "3D" = list(xl="Synthetic phenotype score", yl="Piezo1 expression", t="Piezo1 vs synthetic score",
              alt=list(x=sc(SYN_ALT), y=p1), paper=list(x=sc(SYN_NEW), y=p1),
              alt_genes=paste("synthetic:", paste(SYN_ALT, collapse=", ")), paper_genes=paste("synthetic:", paste(SYN_NEW, collapse=", "))),
  "4B" = list(xl="Piezo1 expression", yl="NFAT pathway score", t="Piezo1 vs NFAT pathway score",
              alt=list(x=p1, y=sc(NFATP13)), paper=list(x=p1, y=sc(NFATP12)),
              alt_genes=paste("NFAT pathway:", paste(NFATP13, collapse=", ")), paper_genes=paste("NFAT pathway:", paste(NFATP12, collapse=", "))),
  "4C" = list(xl="Synthetic phenotype score", yl="NFAT pathway score", t="NFAT pathway vs synthetic score",
              alt=list(x=sc(SYN_ALT), y=sc(NFATP13)), paper=list(x=sc(SYN_NEW), y=sc(NFATP12)),
              alt_genes=paste("NFAT pathway 13 genes; synthetic:", paste(SYN_ALT, collapse=", ")),
              paper_genes=paste("NFAT pathway 12 genes (no Nfat5); synthetic:", paste(SYN_NEW, collapse=", "))))
# reference values for the alternative gene sets (pipeline source data)
REF_ALT <- list("3D"=c(0.1430, 0.1631), "4B"=c(0.0908, 0.2078), "4C"=c(0.0286, 0.1936))
cs <- list(); cells <- list()
for (pn in names(COR)) {
  for (v in c("alt","paper")) for (a in c("Young","Aged")) {
    k <- ctrl & age == a; r <- sp(COR[[pn]][[v]]$x, COR[[pn]][[v]]$y, k)
    cs[[length(cs)+1]] <- data.frame(panel=pn, version=ifelse(v=="alt","alternative","paper"), age=a, n=r["n"], rho=r["rho"], p=r["p"],
                                     star=star(r["p"]), genes=COR[[pn]][[paste0(v,"_genes")]])
    if (v == "alt") stopifnot(abs(r["rho"] - REF_ALT[[pn]][match(a, c("Young","Aged"))]) < 5e-4)  # matches the pipeline source data
  }
  cells[[pn]] <- data.frame(panel=pn, cell=colnames(x)[ctrl], sample=smp[ctrl], age=age[ctrl],
                            x_paper=COR[[pn]]$paper$x[ctrl], y_paper=COR[[pn]]$paper$y[ctrl],
                            x_alternative=COR[[pn]]$alt$x[ctrl], y_alternative=COR[[pn]]$alt$y[ctrl])
}
cs <- do.call(rbind, cs); rownames(cs) <- NULL
write.csv(cs, file.path(OUT, "correlations_paper_vs_alternative.csv"), row.names=FALSE)
write.csv(do.call(rbind, cells), file.path(OUT, "correlation_panels_per_cell.csv"), row.names=FALSE)

corplot <- function(pn) {
  m <- COR[[pn]]; d <- cells[[pn]]; d$age <- factor(d$age, levels=c("Young","Aged"))
  r <- cs[cs$panel == pn & cs$version == "paper", ]
  lab <- sprintf("%s:  rho = %+.3f  %s", c("Young (3mo)","Aged (24mo)"), r$rho[match(c("Young","Aged"), r$age)], r$star[match(c("Young","Aged"), r$age)])
  ggplot(d, aes(x_paper, y_paper, colour=age)) +
    geom_smooth(method="lm", formula=y ~ x, se=TRUE, linewidth=1.5, fill="grey70") +
    scale_colour_manual(values=c(Young="#4DAF4A", Aged="#E41A1C"), guide="none") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.05, vjust=c(1.6, 3.4), label=lab, colour=c("#4DAF4A","#E41A1C"), fontface="bold", size=4) +
    labs(title=paste0(pn, "  ", m$t), x=m$xl, y=m$yl) + theme_publication(base_size=13) + theme(plot.title=element_text(hjust=0, size=14))
}
for (pn in names(COR)) { pl <- corplot(pn)
  ggsave(file.path(FIG, paste0("Fig", pn, ".png")), pl, width=5.6, height=5, dpi=300, bg="white")
  ggsave(file.path(FIG, paste0("Fig", pn, ".pdf")), pl, width=5.6, height=5, bg="white") }

# ---------- 4A: Piezo1 vs pathway scores ----------
P4A_alt <- list("NFAT Pathway"=NFATP13, "TGF-beta/Stress"=GS$tgfbeta_pathway, "Calcium SOCE/IP3"=GS$calcium_pathway,
                "Core Contractile"=GS$contractile_core, "Dedifferentiation"=GS$dedifferentiation, "NFAT-Responsive Remodeling"=ACT6,
                "ECM Remodeling"=GS$ecm_remodeling, "Synthetic Phenotype"=SYN_ALT)
P4A_paper <- P4A_alt; P4A_paper[["NFAT Pathway"]] <- NFATP12; P4A_paper[["Synthetic Phenotype"]] <- SYN_NEW
P4A_paper[["NFAT-Responsive Remodeling"]] <- ACT5; names(P4A_paper)[names(P4A_paper) == "NFAT-Responsive Remodeling"] <- "NFAT-Associated Remodeling"
SUB4A <- c("NFAT Pathway"=0.0908, "TGF-beta/Stress"=0.1132, "Calcium SOCE/IP3"=-0.1053, "Core Contractile"=-0.2402,
           "Dedifferentiation"=0.1282, "NFAT-Responsive Remodeling"=0.1463, "ECM Remodeling"=0.0501, "Synthetic Phenotype"=0.1430)  # young rho, FigS3G_data.xlsx
a4 <- list()
for (v in c("alt","paper")) { L <- if (v == "alt") P4A_alt else P4A_paper
  for (nm in names(L)) for (a in c("Young","Aged")) { k <- ctrl & age == a; r <- sp(p1, sc(L[[nm]]), k)
    if (v == "alt" && a == "Young") stopifnot(abs(r["rho"] - SUB4A[[nm]]) < 5e-4)
    a4[[length(a4)+1]] <- data.frame(version=ifelse(v=="alt","alternative","paper"), pathway=nm, age=a, n=r["n"], rho=r["rho"], p=r["p"],
                                     star=star(r["p"]), differs_from_alternative=(v=="paper" && !identical(L[[nm]], P4A_alt[[nm]]) ) | nm == "NFAT-Associated Remodeling",
                                     genes=paste(L[[nm]], collapse=", ")) } }
a4 <- do.call(rbind, a4); rownames(a4) <- NULL
write.csv(a4, file.path(OUT, "fig4A_paper_vs_alternative.csv"), row.names=FALSE)
d4 <- a4[a4$version == "paper", ]
ord <- c("Core Contractile","Calcium SOCE/IP3","Synthetic Phenotype","ECM Remodeling","Dedifferentiation","NFAT Pathway","TGF-beta/Stress","NFAT-Associated Remodeling")
d4$pathway <- factor(d4$pathway, levels=rev(ord)); d4$age <- factor(ifelse(d4$age=="Young","3mo","24mo"), levels=c("3mo","24mo"))
pA <- ggplot(d4, aes(rho, pathway, colour=age)) + geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_point(size=3.2) + geom_text(aes(label=star, vjust=ifelse(age=="3mo", -0.8, 1.8)), size=3.2, show.legend=FALSE) +  # 3mo star above, 24mo below
  scale_colour_manual(values=c("3mo"="#4DAF4A","24mo"="#E41A1C"), name="Age") +
  labs(title="4A  Piezo1 vs pathway scores", x="Spearman rho with Piezo1", y=NULL) + theme_publication(base_size=12) +
  theme(plot.title=element_text(hjust=0, size=14), legend.position="bottom")
ggsave(file.path(FIG, "Fig4A.png"), pA, width=6, height=5.2, dpi=300, bg="white")
ggsave(file.path(FIG, "Fig4A.pdf"), pA, width=6, height=5.2, bg="white")

# ---------- 4D: NFAT-associated remodeling score, 5 genes (superplot inputs) ----------
G4 <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO"); s5 <- sc(ACT5); s6 <- sc(ACT6)
stopifnot(abs(mean(s6[grp=="Aged_Ctrl"]) - 0.16713) < 1e-4)  # 6-gene score, aged CTRL mean
star4 <- function(p) ifelse(p < 1e-4, "****", star(p))
D4 <- file.path(OUT, "fig4D_superplot_inputs"); dir.create(D4, showWarnings=FALSE)
write.csv(data.frame(panel="4D_rec", metric="nfat_associated_remodeling_5genes", cell=colnames(x), sample=smp, group=grp, value=s5),
          file.path(D4, "per_cell_values.csv"), row.names=FALSE)
pa <- aggregate(value ~ sample + group, data.frame(sample=smp, group=grp, value=s5), mean); pa$n_cells <- as.integer(table(smp)[pa$sample])
write.csv(data.frame(panel="4D_rec", metric="nfat_associated_remodeling_5genes", group=pa$group, sample=pa$sample, animal_mean=pa$value, n_cells=pa$n_cells),
          file.path(D4, "per_animal_values.csv"), row.names=FALSE)
CMP <- list(c("Young_Ctrl","Aged_Ctrl"), c("Young_KO","Aged_KO"), c("Young_Ctrl","Young_KO"), c("Aged_Ctrl","Aged_KO"))
st <- do.call(rbind, lapply(CMP, function(c) { p <- wilcox.test(s5[grp==c[1]], s5[grp==c[2]])$p.value; p0 <- wilcox.test(s6[grp==c[1]], s6[grp==c[2]])$p.value
  data.frame(panel="4D_rec", group1=c[1], group2=c[2], n1=sum(grp==c[1]), n2=sum(grp==c[2]), p=p, p_bonf=min(1,4*p), sig=star4(min(1,4*p)),
             p_bonf_6genes=min(1,4*p0), sig_6genes=star4(min(1,4*p0))) }))
write.csv(st[,1:8], file.path(D4, "stats_cell_level.csv"), row.names=FALSE)
write.csv(st, file.path(OUT, "fig4D_stats_paper_vs_alternative.csv"), row.names=FALSE)
bars <- do.call(rbind, lapply(G4, function(g) data.frame(group=g, n=sum(grp==g), mean_5genes=mean(s5[grp==g]), sem_5genes=sd(s5[grp==g])/sqrt(sum(grp==g)),
                                                         mean_6genes=mean(s6[grp==g]), sem_6genes=sd(s6[grp==g])/sqrt(sum(grp==g)))))
write.csv(bars, file.path(OUT, "fig4D_bars_paper_vs_alternative.csv"), row.names=FALSE)
options(width=200); print(cs, digits=3); print(a4[,1:7], digits=3); print(st, digits=3); print(bars, digits=4)
