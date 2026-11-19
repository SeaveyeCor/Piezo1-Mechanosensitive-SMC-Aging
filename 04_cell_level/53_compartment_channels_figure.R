# PAPER PANEL: Fig S3E
# PLOTS:  % of contractile SMCs (outside the synthetic cluster Cl5) detecting each mechanosensitive channel
# CELLS:  contractile visceral SMCs (n = 1,302)
# GROUPS COMPARED: P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   logistic regression adjusted for log(genes detected), Bonferroni x8 channels
#
# 53: figure for 52 (channel detection in contractile SMCs, CTRL vs P1KO) for all 9 channels, plus a forest
# of the nFeature-matched aged P1KO - aged CTRL difference. Also exports per-cell raw channel
# counts.
suppressMessages({library(Seurat); library(ggplot2); library(patchwork); library(Matrix)})
setwd("."); source("pipeline/utils/plotting_theme.R")
OUT <- "output/cell_level"; FIG <- file.path(OUT,"panels_publication"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
o <- readRDS(file.path(OUT,"data","smc_cache.rds")); DefaultAssay(o) <- "RNA"
cnt <- GetAssayData(o,layer="counts"); grp <- as.character(o$group)
smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)])
cl <- as.character(readRDS(file.path(OUT,"data","clusters.rds")))
res <- read.csv(file.path(OUT,"data","channels_contractile_compartment_agedgeno.csv"))
res <- res[res$gene != "Piezo1",]  # Piezo1 is the deleted gene, not a compensation candidate
CH <- res$gene; inC <- cl != "5"
GR <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
GLAB <- c(Young_Ctrl="CTRL, 3 mo", Aged_Ctrl="CTRL, 24 mo", Young_KO="P1KO, 3 mo", Aged_KO="P1KO, 24 mo")
GCOL <- c(Young_Ctrl="#CCCCCC", Aged_Ctrl="#7F7F7F", Young_KO="#ffc2da", Aged_KO="#ff7fb2")

# A: forest of matched differences
res$gene <- factor(res$gene, levels=rev(CH))
pA <- ggplot(res, aes(matched_diff_pp, gene)) +
  geom_vline(xintercept=0, linewidth=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, linewidth=0.6) +
  geom_point(size=3.2, shape=21, fill="#ff7fb2", colour="#fe0267", stroke=0.9) +
  labs(title="Aged P1KO vs Aged CTRL",
       subtitle="Contractile SMCs, matched on genes detected per cell (300 bootstrap draws)",
       x="Difference in % SMCs detecting\n(percentage points, median and 95% interval)", y=NULL) +
  theme_publication(base_size=13) + theme(plot.subtitle=element_text(size=9.5, colour="grey25"))

# B: per-group bars with per-animal circles, contractile compartment
gd <- do.call(rbind, lapply(CH, function(g){ x <- as.numeric(cnt[g,])>0
  do.call(rbind, lapply(GR, function(k){ m <- grp==k & inC
    data.frame(gene=g, grp=k, pct=100*mean(x[m]), n=sum(m)) })) }))
ad <- do.call(rbind, lapply(CH, function(g){ x <- as.numeric(cnt[g,])>0
  do.call(rbind, lapply(sort(unique(smp)), function(s){ m <- smp==s & inC
    data.frame(gene=g, grp=grp[smp==s][1], sample=s, pct=100*mean(x[m]), n=sum(m)) })) }))
for (d in c("gd","ad")) { z <- get(d); z$grp <- factor(z$grp, levels=GR); z$gene <- factor(z$gene, levels=CH); assign(d, z) }
lab <- data.frame(gene=factor(CH, levels=CH), y=sapply(CH, function(g) max(ad$pct[ad$gene==g]))*1.12 + 1,
                  txt=ifelse(res$logit_p[match(CH, as.character(res$gene))] * length(CH) < 0.05, "*", "ns"))
pB <- ggplot(gd, aes(grp, pct, fill=grp)) +
  geom_col(width=0.7, colour="black", linewidth=0.3) +
  geom_point(data=ad, aes(grp, pct), inherit.aes=FALSE, shape=21, fill="white", size=1.8, stroke=0.6) +
  geom_segment(data=lab, aes(x=2, xend=4, y=y, yend=y), inherit.aes=FALSE, linewidth=0.4) +
  geom_text(data=lab, aes(x=3, y=y, label=txt), inherit.aes=FALSE, vjust=-0.4, size=3.6, fontface="bold") +
  facet_wrap(~gene, ncol=3, scales="free_y") +
  scale_fill_manual(values=GCOL, labels=GLAB, name=NULL) +
  scale_x_discrete(labels=c(Young_Ctrl="3", Aged_Ctrl="24", Young_KO="3", Aged_KO="24")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Channel Detection in Contractile SMCs",
       subtitle="Bars: raw rate, circles: animals. Bracket: aged P1KO vs aged CTRL after matching on genes detected\n(stars: depth-adjusted logistic p, Bonferroni-corrected for 8 channels). Unmatched gaps reflect library complexity.", x="Age (months)", y="SMCs with detected transcript (%)") +
  theme_publication(base_size=12) +
  theme(legend.position="bottom", strip.background=element_blank(), strip.text=element_text(face="bold"), plot.subtitle=element_text(size=9, colour="grey25"))

# Paper-style variant: bars = pooled % with SEM across cells (sqrt(p(1-p)/n)),
# no per-animal circles; genotype fill as in the manuscript (grey CTRL, pink P1KO)
gd$se <- 100*sqrt((gd$pct/100)*(1-gd$pct/100)/gd$n)
gd$geno <- factor(ifelse(grepl("KO", gd$grp), "P1KO", "CTRL"), levels=c("CTRL","P1KO"))
gd2 <- gd; gd2$grp <- factor(as.character(gd2$grp), levels=c("Young_Ctrl","Young_KO","Aged_Ctrl","Aged_KO"))  # 3 CTRL, 3 P1KO, 24 CTRL, 24 P1KO
# both within-age genotype comparisons, young from 52b (same method)
ry <- read.csv(file.path(OUT,"data","channels_contractile_compartment_younggeno.csv")); ry <- ry[match(CH, ry$gene),]
# one p per comparison = depth-adjusted logistic p, Bonferroni-corrected for the 8 channels
# tested at each age
star <- function(r) { q <- pmin(1, r$logit_p * length(CH)); ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", "ns"))) }
top <- sapply(CH, function(g) max((gd$pct+gd$se)[gd$gene==g]))
lab2 <- rbind(data.frame(gene=factor(CH, levels=CH), x0=1, x1=2, y=top*1.10 + 0.6, txt=star(ry)),
              data.frame(gene=factor(CH, levels=CH), x0=3, x1=4, y=top*1.10 + 0.6, txt=star(res[match(CH, as.character(res$gene)),])))
fr <- rbind(data.frame(ry[,c("gene","matched_diff_pp","lo","hi")], age="3 mo"),
            data.frame(res[,c("gene","matched_diff_pp","lo","hi")], age="24 mo"))
fr$gene <- factor(as.character(fr$gene), levels=rev(CH)); fr$age <- factor(fr$age, levels=c("24 mo","3 mo"))
pA2 <- ggplot(fr, aes(matched_diff_pp, gene, fill=age)) +
  geom_vline(xintercept=0, linewidth=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25, linewidth=0.6, position=position_dodge(width=0.6)) +
  geom_point(size=3, shape=21, colour="#fe0267", stroke=0.9, position=position_dodge(width=0.6)) +
  scale_fill_manual(values=c("3 mo"="white", "24 mo"="#ff7fb2"), breaks=c("3 mo","24 mo"), name=NULL) +
  labs(title="P1KO vs CTRL, Within Age",
       subtitle="Contractile SMCs, matched on genes detected per cell (300 bootstrap draws)",
       x="Difference in % SMCs detecting, P1KO \u2212 CTRL\n(percentage points, median and 95% interval)", y=NULL) +
  theme_publication(base_size=13) + theme(legend.position="bottom", plot.subtitle=element_text(size=9.5, colour="grey25"))
pB2 <- ggplot(gd2, aes(grp, pct, fill=geno, colour=geno)) +
  geom_col(width=0.7, linewidth=0.4) +
  geom_errorbar(aes(ymin=pmax(pct-se,0), ymax=pct+se), width=0.25, linewidth=0.5, colour="black") +
  geom_segment(data=lab2, aes(x=x0, xend=x1, y=y, yend=y), inherit.aes=FALSE, linewidth=0.4) +
  geom_text(data=lab2, aes(x=(x0+x1)/2, y=y, label=txt), inherit.aes=FALSE, vjust=-0.4, size=3.6, fontface="bold") +
  facet_wrap(~gene, ncol=3, scales="free_y") +
  scale_fill_manual(values=c(CTRL="#7F7F7F", P1KO="#ff7fb2"), name=NULL) +
  scale_colour_manual(values=c(CTRL="black", P1KO="#fe0267"), name=NULL) +
  scale_x_discrete(labels=c(Young_Ctrl="3", Aged_Ctrl="24", Young_KO="3", Aged_KO="24")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Channel Detection in Contractile SMCs",
       subtitle="Mean \u00b1 SEM across SMCs (unmatched). Brackets: P1KO vs CTRL within age, after matching on genes detected\n(stars: depth-adjusted logistic p, Bonferroni-corrected for 8 channels). Unmatched gaps reflect library complexity.", x="Age (months)", y="SMCs with detected transcript (%)") +
  theme_publication(base_size=12) +
  theme(legend.position="bottom", strip.background=element_blank(), strip.text=element_text(face="bold"), plot.subtitle=element_text(size=9, colour="grey25"))
p2 <- pA2 + pB2 + plot_layout(widths=c(1,1.5))
ggsave(file.path(FIG,"Fig_S3E_no_compensation_contractile_sem.png"), p2, width=15, height=8.5, dpi=300, bg="white")
ggsave(file.path(FIG,"Fig_S3E_no_compensation_contractile_sem.pdf"), p2, width=15, height=8.5, bg="white")
p <- pA + pB + plot_layout(widths=c(1,1.5))
ggsave(file.path(FIG,"Fig_S3E_no_compensation_contractile.png"), p, width=15, height=8.5, dpi=300, bg="white")
ggsave(file.path(FIG,"Fig_S3E_no_compensation_contractile.pdf"), p, width=15, height=8.5, bg="white")
write.csv(merge(gd, aggregate(cbind(pct,n)~gene+grp, ad, function(v) paste(round(v,2), collapse="; ")),
                by=c("gene","grp"), suffixes=c("_pooled","_per_animal")),
          file.path(OUT,"data","channels_contractile_compartment_by_group.csv"), row.names=FALSE)

# per-cell raw export: every bona fide SMC x 20 channels (UMI counts) + metadata
ALL <- c("Piezo1","Piezo2","Tmem63a","Tmem63b","Tmem120a","Trpv1","Trpv2","Trpv4","Trpc1","Trpc3",
         "Trpc6","Trpm4","Trpm7","Pkd1","Pkd2","Kcnk2","Kcnk4","Kcnk10","Asic1","Asic3")
ALL <- ALL[ALL %in% rownames(cnt)]
raw <- data.frame(cell_barcode=colnames(cnt), group=grp, animal=smp, smc_subcluster=cl,
                  compartment=ifelse(inC,"contractile (non-Cl5)","synthetic (Cl5)"),
                  total_UMI=Matrix::colSums(cnt), genes_detected=o$nFeature_RNA,
                  as.matrix(t(cnt[ALL,])), check.names=FALSE)
write.csv(raw, file.path(OUT,"data","per_cell_channel_counts_SMC.csv"), row.names=FALSE)
cat(sprintf("raw export: %d cells x %d channels\n", nrow(raw), length(ALL)))
