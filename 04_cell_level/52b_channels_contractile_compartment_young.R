# PAPER PANEL: Fig S3E
# PLOTS:  % of contractile SMCs (outside the synthetic cluster Cl5) detecting each mechanosensitive channel
# CELLS:  contractile visceral SMCs (n = 1,302)
# GROUPS COMPARED: P1KO vs CTRL at 3 mo; P1KO vs CTRL at 24 mo
# TEST:   logistic regression adjusted for log(genes detected), Bonferroni x8 channels
#
# 52b: script 52 run on the young contrast (young CTRL vs young P1KO), same method, for the within-age genotype comparison at 3 months.
# Piezo2: contractile compartment only (non-Cl5),
# 1:1 nFeature matching at 10% tolerance over 300 bootstrap draws (cells resampled with replacement), plus logistic with nFeature covariate.
suppressMessages({library(Seurat); library(Matrix)}); setwd("."); set.seed(1)
OUT <- "output/cell_level"; dir.create(file.path(OUT, "data"), showWarnings=FALSE, recursive=TRUE)
o <- readRDS(file.path(OUT,"data","smc_cache.rds")); DefaultAssay(o) <- "RNA"
cnt <- GetAssayData(o,layer="counts"); grp <- as.character(o$group); nF <- o$nFeature_RNA
smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)])
cl <- readRDS(file.path(OUT,"data","clusters.rds"))
CH <- read.csv(file.path(OUT,"data","conditional_positive_channels.csv"))$gene
inC <- cl != "5"
rf <- which(grp=="Young_Ctrl" & inC); ts <- which(grp=="Young_KO" & inC)
cat(sprintf("contractile compartment: young CTRL %d, young KO %d; median nFeature %g vs %g\n",
    length(rf), length(ts), median(nF[rf]), median(nF[ts])))
# Bootstrap. Each draw resamples cells of both groups with replacement before
# matching, so the 95% interval carries cell-sampling noise.
DR <- replicate(300, { rb <- sample(rf, replace=TRUE); tb <- sample(ts, replace=TRUE)
  ur <- integer(); ut <- integer(); used <- rep(FALSE, length(tb))
  for (k in rb) { ok <- which(!used & abs(nF[tb]-nF[k])/nF[k] <= 0.10)
    if (!length(ok)) next; j <- ok[sample.int(length(ok),1)]; used[j] <- TRUE; ur <- c(ur,k); ut <- c(ut,tb[j]) }
  list(r=ur,t=ut) }, simplify=FALSE)
cat(sprintf("median matched pairs %g\n", median(sapply(DR, function(D) length(D$r)))))
res <- do.call(rbind, lapply(CH, function(g){ x <- as.numeric(cnt[g,]>0)
  v <- sapply(DR, function(D) 100*(mean(x[D$t])-mean(x[D$r])))
  vr <- sapply(DR, function(D) 100*mean(x[D$r])); vt <- sapply(DR, function(D) 100*mean(x[D$t]))   # matched % per group
  keep <- c(rf,ts); y <- x[keep]; G <- factor(grp[keep], levels=c("Young_Ctrl","Young_KO"))
  f <- summary(glm(y ~ G + log(nF[keep]), family=binomial))$coefficients[2,]
  pa <- sapply(sort(unique(smp[keep])), function(s) round(100*mean(x[keep][smp[keep]==s]),1))
  data.frame(gene=g, pct_ctrl=round(100*mean(x[rf]),1), pct_ko=round(100*mean(x[ts]),1),
    matched_diff_pp=median(v), lo=quantile(v,.025), hi=quantile(v,.975),
    matched_pct_ctrl=median(vr), matched_pct_ctrl_lo=quantile(vr,.025), matched_pct_ctrl_hi=quantile(vr,.975),
    matched_pct_ko=median(vt), matched_pct_ko_lo=quantile(vt,.025), matched_pct_ko_hi=quantile(vt,.975),
    logit_beta=f[1], logit_p=f[4], per_animal=paste(names(pa), pa, sep="=", collapse="; ")) }))
rownames(res) <- NULL
print(res[,1:8], digits=3, row.names=FALSE); cat("\n"); print(res[,c("gene","per_animal")], row.names=FALSE)
write.csv(res, file.path(OUT,"data","channels_contractile_compartment_younggeno.csv"), row.names=FALSE)
