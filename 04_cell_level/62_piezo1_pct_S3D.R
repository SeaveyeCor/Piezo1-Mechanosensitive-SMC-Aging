# PAPER PANEL: Fig S3D
# PLOTS:  % of visceral SMCs with Piezo1 detected (>= 1 UMI)
# CELLS:  visceral SMCs (n = 1,604: Young CTRL 469, Aged CTRL 477, Young P1KO 222, Aged P1KO 436)
# GROUPS COMPARED: P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   logistic regression adjusted for log(genes detected); depth-matched bootstrap interval
#
# 62: Fig S3D analysed like Fig S3E. Bars = % SMCs with
# Piezo1 detected (>=1 UMI), pooled; error bars = SEM across SMCs, sqrt(p(1-p)/n). P1KO vs CTRL
# within age: (i) cells resampled with replacement and matched 1:1 on genes detected (within 10%),
# 300 bootstrap draws, 95% interval; (ii) logistic regression detected ~ genotype + log(genes
# detected). Stars from the logistic p (* <0.05, ** <0.01, *** <0.001); the
# interval does not gate the star. Population: all bona fide SMCs (deletion validation), whereas
# S3E uses contractile SMCs (see 52).
suppressMessages({library(Seurat); library(Matrix); library(ggplot2)}); setwd(".")
source("pipeline/utils/plotting_theme.R"); set.seed(1)
OUT <- "output/cell_level"; FIG <- file.path(OUT,"panels_publication"); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)
o <- readRDS(file.path(OUT,"data","smc_cache.rds")); DefaultAssay(o) <- "RNA"
cnt <- GetAssayData(o,layer="counts"); grp <- as.character(o$group); nF <- o$nFeature_RNA
smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)])
x <- as.numeric(cnt["Piezo1",] > 0)
test <- function(age){ rf <- which(grp==paste0(age,"_Ctrl")); ts <- which(grp==paste0(age,"_KO"))
  v <- replicate(300, { rb <- sample(rf, replace=TRUE); tb <- sample(ts, replace=TRUE)
    ur <- integer(); ut <- integer(); used <- rep(FALSE, length(tb))
    for (k in rb) { ok <- which(!used & abs(nF[tb]-nF[k])/nF[k] <= 0.10)
      if (!length(ok)) next; j <- ok[sample.int(length(ok),1)]; used[j] <- TRUE; ur <- c(ur,k); ut <- c(ut,tb[j]) }
    100*(mean(x[ut]) - mean(x[ur])) })
  keep <- c(rf,ts); G <- as.numeric(keep %in% ts)
  f <- summary(glm(x[keep] ~ G + log(nF[keep]), family=binomial))$coefficients[2,]
  data.frame(gene="Piezo1", age=age, pct_ctrl=100*mean(x[rf]), pct_ko=100*mean(x[ts]),
             matched_diff_pp=median(v), lo=quantile(v,.025), hi=quantile(v,.975),
             logit_beta=f[1], logit_p=f[4]) }
st <- rbind(test("Young"), test("Aged")); rownames(st) <- NULL
# one p per comparison; the star comes from the logistic p alone (interval = effect size only).
# Labels unchanged (both intervals exclude 0), so the figure and CSV did not need re-running.
st$label <- ifelse(st$logit_p < .001, "***", ifelse(st$logit_p < .01, "**", ifelse(st$logit_p < .05, "*", "ns")))
print(st, digits=3); write.csv(st, file.path(OUT,"data","piezo1_pct_S3D_tests.csv"), row.names=FALSE)
GR <- c("Young_Ctrl","Young_KO","Aged_Ctrl","Aged_KO")
gd <- do.call(rbind, lapply(GR, function(g){ m <- grp==g; p <- mean(x[m])
  data.frame(grp=g, n=sum(m), k=sum(x[m]), pct=100*p, se=100*sqrt(p*(1-p)/sum(m))) }))
gd$grp <- factor(gd$grp, levels=GR); gd$geno <- factor(ifelse(grepl("KO",gd$grp),"P1KO","CTRL"), levels=c("CTRL","P1KO"))
ad <- do.call(rbind, lapply(sort(unique(smp)), function(s){ m <- smp==s
  data.frame(grp=grp[m][1], animal=s, pct=100*mean(x[m])) }))
ad$grp <- factor(ad$grp, levels=GR)
ad <- do.call(rbind, lapply(split(ad, ad$grp), function(d){ d$dx <- if (nrow(d)==1) 0 else c(-.13,.13); d }))
mk <- function(points){
  yb <- max(c(gd$pct+gd$se, if (points) ad$pct))*1.10
  br <- data.frame(x0=c(1,3), x1=c(2,4), y=yb, lab=st$label)
  p <- ggplot(gd, aes(grp, pct, fill=geno, colour=geno)) +
    geom_col(width=0.7, linewidth=0.5) +
    geom_errorbar(aes(ymin=pmax(pct-se,0), ymax=pct+se), width=0.25, linewidth=0.6, colour="black")
  if (points) p <- p + geom_point(data=ad, aes(as.numeric(grp)+dx, pct), inherit.aes=FALSE, shape=21, size=2.6, fill="white", stroke=.8)
  p + geom_segment(data=br, aes(x=x0, xend=x1, y=y, yend=y), inherit.aes=FALSE, linewidth=0.5) +
    geom_text(data=br, aes(x=(x0+x1)/2, y=y, label=lab), inherit.aes=FALSE, vjust=-0.3, size=6, fontface="bold") +
    scale_fill_manual(values=c(CTRL="#7F7F7F", P1KO="#ff7fb2"), name=NULL) +
    scale_colour_manual(values=c(CTRL="black", P1KO="#fe0267"), name=NULL) +
    scale_x_discrete(labels=c("3","3","24","24")) +
    scale_y_continuous(limits=c(0, yb*1.15), expand=c(0,0)) +
    labs(x="Age (months)", y="% SMCs with Piezo1 detected") +
    theme_publication(base_size=16) + theme(legend.position="top") }
for (ext in c("png","pdf")) {
  ggsave(file.path(FIG, paste0("Fig_S3D_piezo1_pct_sem.",ext)), mk(FALSE), width=4.8, height=5.6, dpi=300, bg="white")
  ggsave(file.path(FIG, paste0("Fig_S3D_piezo1_pct_sem_permouse.",ext)), mk(TRUE), width=4.8, height=5.6, dpi=300, bg="white") }
write.csv(merge(gd, aggregate(pct~grp, ad, function(v) paste(round(v,2), collapse="; ")), by="grp", suffixes=c("","_per_animal")),
          file.path(OUT,"data","piezo1_pct_S3D_by_group.csv"), row.names=FALSE)
