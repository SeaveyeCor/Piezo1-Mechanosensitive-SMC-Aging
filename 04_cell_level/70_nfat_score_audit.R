# 70: robustness of the Fig 4D NFAT-associated score to removing individual members (Fas, Wnt5a, Timp1),
# and of the NFAT pathway score (Fig 4B, C) to removing Nfat5.
# Score = colMeans of log-normalized RNA data (zeros included), exactly as 02_scores_singlecell.R / Fig4E module.
# Tests as Fig 4D: two-sided Mann-Whitney U across cells, Bonferroni x4. Usage (from .):
# Rscript 04_cell_level/70_nfat_score_audit.R
suppressMessages({library(Seurat); library(Matrix)})
OUT <- "output/cell_level/data/nfat_audit"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
o <- readRDS("output/cell_level/data/smc_cache.rds"); DefaultAssay(o) <- "RNA"
x <- GetAssayData(o, layer="data"); cnt <- GetAssayData(o, layer="counts")
grp <- as.character(o$group); smp <- as.character(o$sample); smp[is.na(smp)] <- as.character(o$orig.ident[is.na(smp)])
G4 <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
CMP <- list(c("Young_Ctrl","Aged_Ctrl"), c("Young_KO","Aged_KO"), c("Young_Ctrl","Young_KO"), c("Aged_Ctrl","Aged_KO"))
full <- c("Wnt5a","Mmp2","Mmp9","Vcam1","Fas","Timp1")
sets <- list("All 6"=full,
             "Without Fas"=setdiff(full,"Fas"), "Without Wnt5a"=setdiff(full,"Wnt5a"), "Without Timp1"=setdiff(full,"Timp1"),
             "Without Mmp2"=setdiff(full,"Mmp2"), "Without Mmp9"=setdiff(full,"Mmp9"), "Without Vcam1"=setdiff(full,"Vcam1"),
             "Without Fas, Wnt5a, Timp1 (Mmp2, Mmp9, Vcam1)"=c("Mmp2","Mmp9","Vcam1"),
             "Mmp2 + Mmp9 only"=c("Mmp2","Mmp9"))
star <- function(p) ifelse(p<1e-4,"****",ifelse(p<1e-3,"***",ifelse(p<.01,"**",ifelse(p<.05,"*","ns"))))
res <- list(); pm <- list()
for (nm in names(sets)) {
  s <- colMeans(as.matrix(x[sets[[nm]], , drop=FALSE]))
  m <- sapply(G4, function(g) mean(s[grp==g]))
  for (c in CMP) {
    p <- wilcox.test(s[grp==c[1]], s[grp==c[2]])$p.value; pb <- min(1, 4*p)
    res[[length(res)+1]] <- data.frame(set=nm, genes=paste(sets[[nm]], collapse=", "), comparison=paste(c[2],"vs",c[1]),
      mean1=m[c[1]], mean2=m[c[2]], fold=m[c[2]]/m[c[1]], p_bonf=pb, sig=star(pb))
  }
  for (a in sort(unique(smp))) pm[[length(pm)+1]] <- data.frame(set=nm, sample=a, group=grp[smp==a][1], n=sum(smp==a), mean=mean(s[smp==a]))
}
res <- do.call(rbind, res); pm <- do.call(rbind, pm)
write.csv(res, file.path(OUT,"nfat_activity_gene_removal_tests.csv"), row.names=FALSE)
write.csv(pm, file.path(OUT,"nfat_activity_gene_removal_per_mouse.csv"), row.names=FALSE)
# per-gene behaviour: detection and mean log-normalized expression by group
pg <- do.call(rbind, lapply(c(full,"Fasl","Nfat5"), function(g) { if (!g %in% rownames(x)) return(data.frame(gene=g, group=G4, pct_detect=NA, mean_lognorm=NA))
  data.frame(gene=g, group=G4, pct_detect=sapply(G4, function(k) 100*mean(cnt[g,grp==k]>0)), mean_lognorm=sapply(G4, function(k) mean(x[g,grp==k]))) }))
write.csv(pg, file.path(OUT,"nfat_genes_per_group.csv"), row.names=FALSE)
# NFAT pathway score (Fig 4B: vs Piezo1; Fig 4C: vs synthetic score), with and without Nfat5, CTRL cells by age
path <- c("Nfatc1","Nfatc2","Nfatc3","Nfatc4","Nfat5","Ppp3ca","Ppp3cb","Ppp3cc","Rcan1","Rcan2","Rcan3","Calm1","Calm3")
# Synthetic score: the SI Methods set (Vim, Col1a2, Dcn, Igfbp6, Efemp1) and, for comparison, an alternative
# set (Klf4, Vim, Mmp2, Mmp9, Timp1). Both reported.
SYN <- list(alternative_set=c("Klf4","Vim","Mmp2","Mmp9","Timp1"), methods_set=c("Vim","Col1a2","Dcn","Igfbp6","Efemp1"))
p1 <- x["Piezo1",]
pc <- list()
for (nm in c("13 genes","12 genes, without Nfat5")) {
  gs <- if (grepl("without", nm)) setdiff(path,"Nfat5") else path; sc <- colMeans(as.matrix(x[gs,]))
  for (g in c("Young_Ctrl","Aged_Ctrl")) { k <- grp==g
    a <- cor.test(p1[k], sc[k], method="spearman", exact=FALSE)
    b <- cor.test(sc[k], colMeans(as.matrix(x[SYN[[1]],k])), method="spearman", exact=FALSE)
    d <- cor.test(sc[k], colMeans(as.matrix(x[SYN[[2]],k])), method="spearman", exact=FALSE)
    pc[[length(pc)+1]] <- data.frame(set=nm, group=g, n=sum(k), rho_vs_Piezo1_4B=a$estimate, p_4B=a$p.value,
      rho_4C_alternative_synthetic=b$estimate, p_4C_alternative=b$p.value, rho_4C_methods_synthetic=d$estimate, p_4C_methods=d$p.value) }
  for (c in CMP) { p <- wilcox.test(sc[grp==c[1]], sc[grp==c[2]])$p.value
    pc[[length(pc)+1]] <- data.frame(set=nm, group=paste(c[2],"vs",c[1],"(fold, Bonferroni p)"), n=NA, rho_vs_Piezo1_4B=mean(sc[grp==c[2]])/mean(sc[grp==c[1]]), p_4B=min(1,4*p),
      rho_4C_alternative_synthetic=NA, p_4C_alternative=NA, rho_4C_methods_synthetic=NA, p_4C_methods=NA) }
}
d3 <- do.call(rbind, lapply(names(SYN), function(nm) do.call(rbind, lapply(c("Young_Ctrl","Aged_Ctrl"), function(g) { k <- grp==g
  t <- cor.test(p1[k], colMeans(as.matrix(x[SYN[[nm]],k])), method="spearman", exact=FALSE); data.frame(synthetic_set=nm, group=g, n=sum(k), rho_3D=t$estimate, p=t$p.value) }))))
write.csv(d3, file.path(OUT,"fig3D_synthetic_definition_check.csv"), row.names=FALSE); print(d3, digits=3)
pc <- do.call(rbind, pc); write.csv(pc, file.path(OUT,"nfat_pathway_with_without_Nfat5.csv"), row.names=FALSE)
options(width=220); print(res[,c("set","comparison","fold","p_bonf","sig")], digits=3); print(pg, digits=3); print(pc, digits=3)
print(reshape(pm[,c("set","sample","mean")], idvar="set", timevar="sample", direction="wide"), digits=3)
