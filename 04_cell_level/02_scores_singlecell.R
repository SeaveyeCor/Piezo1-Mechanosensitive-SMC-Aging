# 02_scores_singlecell.R -- per-cell score values + single-cell distribution panels.
# Gene panels as listed in the SI Methods.
suppressMessages({library(Seurat); library(ggplot2); library(dplyr); library(ggbeeswarm)
                  library(openxlsx); library(readxl)})
setwd(".")
OUT <- "output/cell_level"; dir.create(file.path(OUT, "data"), showWarnings=FALSE, recursive=TRUE); dir.create(file.path(OUT, "panels"), showWarnings=FALSE, recursive=TRUE)
obj <- readRDS(file.path(OUT,"data","smc_cache.rds"))
DefaultAssay(obj) <- "RNA"
expr <- GetAssayData(obj, layer = "data")          # log-normalized
cat("expr matrix:", nrow(expr), "genes x", ncol(expr), "cells\n")

PANELS <- list(
  contractile_capacity = c("Myh11","Acta2","Tagln","Cnn1","Mylk","Myocd","Smtn"),
  contractile_corr19   = c("Acta2","Myh11","Myl9","Cnn1","Cald1","Tagln","Tpm1","Tpm2",
                           "Calm1","Calm2","Calm3","Mylk","Rock1","Rock2","Prkca","Prkcb",
                           "Prkcd","Dapk3","Ppp1r14a"),
  synthetic_supplement = c("Vim","Col1a2","Dcn","Igfbp6","Efemp1"),
  synthetic_pipeline   = c("Klf4","Vim","Mmp2","Mmp9","Timp1"),
  nfat_pathway         = c("Nfatc1","Nfatc2","Nfatc3","Nfatc4","Nfat5","Ppp3ca","Ppp3cb",
                           "Ppp3cc","Rcan1","Rcan2","Rcan3","Calm1","Calm3"),
  nfat_activity        = c("Wnt5a","Mmp2","Mmp9","Vcam1","Fas","Timp1"),
  ecm_structural       = c("Bgn","Col6a1","Dcn","Sparc","Timp2"),
  cytoskeletal_remod   = c("Ezr","Syne2","Lmna","Vim")
)

# ---- gene availability audit: absence in a lookup is not absence in the data ----
allg <- rownames(expr)
audit <- do.call(rbind, lapply(names(PANELS), function(p) {
  g <- PANELS[[p]]
  data.frame(panel = p, gene = g,
             in_object = g %in% allg,
             case_insensitive_match = sapply(g, function(x) {
               h <- allg[tolower(allg) == tolower(x)]; if (length(h)) h[1] else NA_character_ }),
             detected_pct = sapply(g, function(x)
               if (x %in% allg) round(100*mean(expr[x,] > 0),2) else NA_real_),
             stringsAsFactors = FALSE)
}))
cat("\n--- genes NOT found as written ---\n")
print(audit[!audit$in_object, c("panel","gene","case_insensitive_match")], row.names = FALSE)
cat("\n--- genes present but detected in <1% of cells ---\n")
print(audit[audit$in_object & audit$detected_pct < 1, c("panel","gene","detected_pct")], row.names = FALSE)

score_of <- function(genes) {
  g <- intersect(genes, allg)
  if (!length(g)) return(rep(NA_real_, ncol(expr)))
  Matrix::colMeans(expr[g, , drop = FALSE])
}
sc <- as.data.frame(lapply(PANELS, score_of))
sc$cell <- colnames(expr); sc$group <- obj$group

# ---- reproducibility check against the published Fig3H panel ----
yc <- mean(sc$contractile_capacity[sc$group == "Young_Ctrl"])
sc$contractile_capacity_pct <- 100 * sc$contractile_capacity / yc
chk <- data.frame(mine = sc$contractile_capacity_pct, published = obj$published_contractile_pct)
cat("\n--- recomputed vs pipeline contractile score (module Fig3H = Fig 3E) ---\n")
cat("max abs diff:", round(max(abs(chk$mine - chk$published)), 6),
    " | Pearson r:", round(cor(chk$mine, chk$published), 6), "\n")

# ---- stats: the four valid pairwise comparisons, Bonferroni m = 4 ----
CMP <- list(c("Young_Ctrl","Aged_Ctrl"), c("Young_KO","Aged_KO"),
            c("Young_Ctrl","Young_KO"),  c("Aged_Ctrl","Aged_KO"))
stars <- function(p) ifelse(p < 1e-4,"****", ifelse(p < 1e-3,"***",
                     ifelse(p < 1e-2,"**", ifelse(p < 0.05,"*","ns"))))
stat_tbl <- do.call(rbind, lapply(names(PANELS), function(p) {
  do.call(rbind, lapply(CMP, function(cc) {
    a <- sc[[p]][sc$group == cc[1]]; b <- sc[[p]][sc$group == cc[2]]
    w <- suppressWarnings(wilcox.test(a, b))
    data.frame(panel = p, group1 = cc[1], group2 = cc[2], n1 = length(a), n2 = length(b),
               median1 = median(a), median2 = median(b), W = unname(w$statistic),
               p_value = w$p.value, p_bonferroni = min(1, w$p.value*4),
               sig = stars(min(1, w$p.value*4)))
  }))
}))
write.csv(stat_tbl, file.path(OUT,"data","score_stats.csv"), row.names = FALSE)
write.csv(sc,       file.path(OUT,"data","per_cell_scores.csv"), row.names = FALSE)
write.csv(audit,    file.path(OUT,"data","gene_panel_audit.csv"), row.names = FALSE)

# ---- single-cell panels: every cell drawn, distribution + median on top ----
GCOL <- c(Young_Ctrl="#7F7F7F", Aged_Ctrl="#7F7F7F", Young_KO="#ff7fb2", Aged_KO="#ff7fb2")
LAB  <- c(Young_Ctrl="Young\nCTRL", Aged_Ctrl="Aged\nCTRL",
          Young_KO="Young\nP1KO",  Aged_KO="Aged\nP1KO")
TITLES <- c(contractile_capacity="Contractile capacity score (7 genes)",
            contractile_corr19  ="Contractile score, correlation set (19 genes)",
            synthetic_supplement="Synthetic phenotype score - supplement panel (5 genes)",
            synthetic_pipeline  ="Synthetic phenotype score - pipeline panel (5 genes)",
            nfat_pathway        ="NFAT pathway score (13 upstream genes)",
            nfat_activity       ="NFAT activity score (6 downstream effectors)",
            ecm_structural      ="ECM score (5 structural genes)",
            cytoskeletal_remod  ="Cytoskeletal remodeling score (4 genes)")

make_panel <- function(p) {
  d <- data.frame(group = sc$group, y = sc[[p]])
  nlab <- d %>% count(group) %>% mutate(lab = paste0(LAB[as.character(group)], "\nn=", n))
  labv <- setNames(nlab$lab, nlab$group)
  st <- stat_tbl[stat_tbl$panel == p, ]
  st$x <- match(st$group1, levels(d$group)); st$xend <- match(st$group2, levels(d$group))
  st <- st[order(abs(st$xend - st$x), pmin(st$x, st$xend)), ]   # short spans lowest
  rng <- range(d$y); pad <- diff(rng) * 0.11
  brs <- data.frame(x = st$x, xend = st$xend,
                    y = max(rng) + pad * seq_len(nrow(st)), lab = st$sig)
  ggplot(d, aes(group, y)) +
    geom_violin(aes(fill = group), colour = NA, alpha = 0.30, scale = "width", width = 0.85) +
    ggbeeswarm::geom_quasirandom(aes(colour = group), size = 0.45, alpha = 0.42,
                                 width = 0.33, stroke = 0) +
    stat_summary(fun = median, geom = "crossbar", width = 0.55,
                 linewidth = 0.5, colour = "black", fatten = 1.6) +
    geom_segment(data = brs, aes(x = x, xend = xend, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.3, colour = "grey35") +
    geom_text(data = brs, aes(x = (x + xend)/2, y = y, label = lab),
              inherit.aes = FALSE, vjust = -0.25, size = 3.1, colour = "grey20") +
    scale_fill_manual(values = GCOL, guide = "none") +
    scale_colour_manual(values = GCOL, guide = "none") +
    scale_x_discrete(labels = labv) +
    labs(title = TITLES[p], x = NULL,
         y = "mean log-normalized expression",
         caption = "each dot = one bona fide visceral SMC; bar = median; Mann-Whitney, Bonferroni m=4") +
    expand_limits(y = max(rng) + pad * (nrow(st) + 0.9)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11.5),
          plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0),
          axis.text.x = element_text(size = 9, lineheight = 0.95))
}

for (p in names(PANELS)) {
  g <- make_panel(p)
  ggsave(file.path(OUT,"panels",paste0("SC_", p, ".png")), g,
         width = 6.0, height = 4.6, dpi = 300, bg = "white")
  ggsave(file.path(OUT,"panels",paste0("SC_", p, ".pdf")), g, width = 6.0, height = 4.6)
}
cat("\nwrote", length(PANELS), "single-cell score panels\n")

wb <- createWorkbook()
addWorksheet(wb,"per_cell_scores");  writeData(wb,"per_cell_scores", sc)
addWorksheet(wb,"statistics");       writeData(wb,"statistics", stat_tbl)
addWorksheet(wb,"gene_panel_audit"); writeData(wb,"gene_panel_audit", audit)
saveWorkbook(wb, file.path(OUT,"data","SourceData_singlecell_scores.xlsx"), overwrite = TRUE)
cat("wrote source data workbook\n")
print(stat_tbl[stat_tbl$panel %in% c("ecm_structural","cytoskeletal_remod","nfat_activity"), ], row.names=FALSE)
