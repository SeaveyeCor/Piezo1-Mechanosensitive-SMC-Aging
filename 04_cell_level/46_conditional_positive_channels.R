# PAPER PANEL: Fig S3B
# PLOTS:  aging change in channel expression within expressing cells, vs detection-matched background genes
# CELLS:  CTRL visceral SMCs only (Young CTRL 469, Aged CTRL 477)
# GROUPS COMPARED: Aged CTRL vs Young CTRL
# TEST:   empirical p vs 120 detection-matched background genes (depth-matched bootstrap)
#
# 46: within cells that EXPRESS a channel, does expression
# rise with age? This conditions on detection, which is the quantity RNAscope measures
# (transcripts per expressing cell) and which the whole-population mean destroys via the
# ~80% dropout. Depth-matched within positives; each channel tested against its own
# detection-matched background using the identical statistic.
suppressMessages({library(Seurat); library(Matrix)}); setwd("."); set.seed(1)
OUT <- "output/cell_level"; dir.create(file.path(OUT, "data"), showWarnings=FALSE, recursive=TRUE)
o <- readRDS(file.path(OUT,"data","smc_cache.rds")); DefaultAssay(o) <- "RNA"
ex <- GetAssayData(o,layer="data"); cnt <- GetAssayData(o,layer="counts")
grp <- as.character(o$group); tot <- Matrix::colSums(cnt)
yc <- which(grp=="Young_Ctrl"); ac <- which(grp=="Aged_Ctrl"); keep <- c(yc,ac)
det <- Matrix::rowMeans(ex[,keep]>0)*100
# mechanically activated ion channels surveyed
CH <- c("Trpm7","Pkd1","Pkd2","Trpc1","Piezo1","Tmem63a","Trpc3","Piezo2","Tmem63b","Trpv2","Trpm4","Tmem120a",
        "Trpv4","Kcnk2","Trpc6","Kcnk4","Kcnk10","Trpv1","Asic1","Asic3")
CH <- CH[CH %in% rownames(cnt) & det[CH] >= 5]            # need enough positive cells
cat(sprintf("channels tested (>=5%% detection): %d\n", length(CH)))

# Bootstrap. Each draw resamples
# young and aged cells with replacement before matching, so the 95% interval carries cell-sampling
# noise. The background test compares medians.
mkdraw <- function(){ yb <- sample(yc, replace=TRUE); ab <- sample(ac, replace=TRUE)
  up <- integer(); un <- integer(); used <- rep(FALSE, length(ab))
  for (k in yb) { ok <- which(!used & abs(tot[ab]-tot[k])/tot[k] <= 0.10)
    if (!length(ok)) next; j <- ok[sample.int(length(ok),1)]; used[j] <- TRUE; up <- c(up,k); un <- c(un,ab[j]) }
  list(y=up,a=un) }
DR <- replicate(200, mkdraw(), simplify=FALSE)
stat <- function(g){ x <- as.numeric(cnt[g,])
  v <- sapply(DR, function(D){ py<-D$y[x[D$y]>0]; pa<-D$a[x[D$a]>0]
    if (length(py)<20 || length(pa)<20) return(NA)
    log2(mean(x[pa])/mean(x[py])) })
  c(median(v,na.rm=TRUE), sum(!is.na(v)), quantile(v,.025,na.rm=TRUE), quantile(v,.975,na.rm=TRUE)) }

res <- do.call(rbind, lapply(CH, function(g){
  s <- stat(g); if (is.na(s[1])) return(NULL)
  pool <- names(det)[abs(det-det[g]) <= 2 & !grepl(paste0("^",g,"$"), names(det))]
  if (length(pool) < 40) return(NULL)
  set.seed(2); bgv <- sapply(sample(pool, min(120,length(pool))), function(b) stat(b)[1])
  bgv <- bgv[is.finite(bgv)]
  data.frame(gene=g, det=round(det[g],1), log2FC_pos=s[1], obs_lo=s[3], obs_hi=s[4],
             bg_median=median(bgv), bg_lo=quantile(bgv,.025), bg_hi=quantile(bgv,.975),
             pctile=100*mean(bgv<s[1]),
             p_emp=2*min(mean(bgv>=s[1]),mean(bgv<=s[1])), n_bg=length(bgv)) }))
res <- res[order(-res$log2FC_pos),]
res$verdict <- ifelse(res$p_emp<0.05 & res$log2FC_pos>0, "UP beyond background",
               ifelse(res$p_emp<0.05 & res$log2FC_pos<0, "DOWN beyond background", "within background"))
print(res[,c("gene","det","log2FC_pos","bg_lo","bg_hi","pctile","p_emp","verdict")],
      row.names=FALSE, digits=3)
write.csv(res, file.path(OUT,"data","conditional_positive_channels.csv"), row.names=FALSE)
# per-draw values behind each channel's error bar (no RNG used here; DR is fixed above)
dr <- do.call(rbind, lapply(as.character(res$gene), function(g){ x <- as.numeric(cnt[g,])
  do.call(rbind, lapply(seq_along(DR), function(i){ D <- DR[[i]]
    py <- D$y[x[D$y]>0]; pa <- D$a[x[D$a]>0]
    data.frame(gene=g, draw=i, n_young_expressing=length(py), n_aged_expressing=length(pa),
               mean_young=mean(x[py]), mean_aged=mean(x[pa]),
               log2FC=if (length(py)<20 || length(pa)<20) NA else log2(mean(x[pa])/mean(x[py]))) })) }))
write.csv(dr, file.path(OUT,"data","conditional_positive_draws.csv"), row.names=FALSE)
