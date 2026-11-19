# PAPER PANEL: Fig 4F
# PLOTS:  module Z-scores along pseudotime (contractile, cytoskeletal remodeling, Piezo1, NFAT-associated remodeling, synthetic)
# CELLS:  Piezo1+ CTRL visceral SMCs (n = 205)
# GROUPS COMPARED: none (curves along pseudotime)
# TEST:   LOESS smoothing; no group test
#
# 73: Fig 4F pseudotime cascade.
#  Pseudotime from pipeline module Fig4I (pipeline/figure_modules/Fig4I_pseudotime_cascade.R):
#  205 Piezo1+ CTRL SMCs, Slingshot pseudotime blended with kNN age fraction, five module Z-scores smoothed by LOESS.
#  Steps: (1) module curves from the exported pseudotime (Fig4I_data.xlsx, sheet pseudotime_data) and the analysed-SMC
#  cache, (2) plot with the 5-gene NFAT-associated remodeling module, (3) export everything needed to replot:
#  per-cell pseudotime + module scores, LOESS grid, cascade order.
# Usage (from .): Rscript 04_cell_level/73_fig4F_cascade.R
suppressMessages({library(Seurat); library(Matrix); library(ggplot2); library(cowplot); library(openxlsx)})
source("pipeline/core/gene_sets.R"); GS <- get_gene_sets()
OUT <- "output/cell_level/data/score_panels/fig4F"; FIG <- "output/cell_level/panels_publication/score_panels"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE); dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

# ---------- pseudotime (exported by the pipeline) ----------
pt_df <- read.xlsx("output/pipeline/Fig4I/Fig4I_data.xlsx", sheet="pseudotime_data")
li <- read.xlsx("output/pipeline/Fig4I/Fig4I_data.xlsx", sheet="lineage_info")
stopifnot(nrow(pt_df) == 205, li$n_cells == 205)
o <- readRDS("output/cell_level/data/smc_cache.rds"); x <- GetAssayData(o, assay="RNA", layer="data")
grp <- as.character(o$group)
stopifnot(all(pt_df$cell %in% colnames(x)))
# the 205 cells are exactly the Piezo1+ CTRL SMCs of the analysed SMC set
p1pos_ctrl <- colnames(x)[grp %in% c("Young_Ctrl","Aged_Ctrl") & x["Piezo1", ] > 0]
stopifnot(setequal(p1pos_ctrl, pt_df$cell))
cells <- pt_df$cell; blended <- pt_df$pseudotime_blended; names(blended) <- cells
age <- pt_df$age; stopifnot(sum(age=="Young") == 107, sum(age=="Aged") == 98)
log_depth <- log10(o$nCount_RNA[cells])

ACT6 <- GS$cascade_nfat_effectors; ACT5 <- setdiff(ACT6, "Fas")
stopifnot(identical(ACT6, GS$nfatc3_remodeling), length(ACT5) == 5)
modules_alt <- list("Contractile"=list(genes=GS$contractile_core, color="#2CA02C"),
                    "Cyto Remodeling"=list(genes=GS$cascade_cyto_remodeling, color="#FF7F00"),
                    "Piezo1"=list(genes="Piezo1", color="#1F77B4"),
                    "NFAT Effectors"=list(genes=ACT6, color="#9467BD"),
                    "Synthetic"=list(genes=GS$cascade_synthetic, color="#E41A1C"))
modules_paper <- modules_alt; names(modules_paper)[4] <- "NFAT-assoc. remodeling"; modules_paper[[4]]$genes <- ACT5

# same scoring as the pipeline module: mean log-normalized expression -> Z across the 205 cells; Piezo1 depth-corrected
raw_score <- function(genes) { g <- intersect(genes, rownames(x)); stopifnot(length(g) == length(genes)); colMeans(as.matrix(x[g, cells, drop=FALSE])) }
z_score <- function(genes, piezo1=FALSE) {
  if (piezo1) { r <- residuals(lm(as.numeric(x["Piezo1", cells]) ~ log_depth)); return((r - mean(r)) / sd(r)) }
  v <- raw_score(genes); (v - mean(v)) / max(sd(v), 0.001)
}
# grid + age fraction exactly as the pipeline
pt_range <- range(blended); trim <- 0.02 * diff(pt_range); grid <- seq(pt_range[1] + trim, pt_range[2] - trim, length.out=200)
age_frac_bg <- sapply(grid, function(g) { w <- dnorm(blended, mean=g, sd=0.08); sum(w * (age == "Aged")) / sum(w) })
bg_df <- data.frame(pt=grid, age_frac=age_frac_bg)

fit_modules <- function(modules, version) {
  lo_all <- list(); pc <- data.frame(cell=cells, sample=as.character(o$sample[cells]), age=age, pseudotime_blended=blended,
                                     pseudotime_slingshot=pt_df$pseudotime_slingshot, age_fraction_knn=pt_df$age_fraction_knn, cluster=pt_df$cluster)
  pc$sample[is.na(pc$sample)] <- as.character(pc$orig.ident[is.na(pc$sample)])
  for (mn in names(modules)) {
    m <- modules[[mn]]; z <- z_score(m$genes, piezo1 = (mn == "Piezo1"))
    pc[[paste0(mn, " | mean expr")]] <- if (mn == "Piezo1") as.numeric(x["Piezo1", cells]) else raw_score(m$genes)
    pc[[paste0(mn, " | Z")]] <- z
    lo <- loess(z ~ pt, data=data.frame(pt=blended, z=z), span=0.45)
    pr <- predict(lo, newdata=data.frame(pt=grid), se=TRUE)
    d <- data.frame(version=version, module=mn, pt=grid, z=pmax(pmin(pr$fit, 2), -2), se=pr$se.fit, age_frac=age_frac_bg)
    d$z_lo <- pmax(d$z - d$se, -2.5); d$z_hi <- pmin(d$z + d$se, 2.5); d$z_unclipped <- pr$fit
    lo_all[[mn]] <- d
  }
  list(loess=do.call(rbind, lo_all), per_cell=pc)
}
ALT <- fit_modules(modules_alt, "alternative"); PAPER <- fit_modules(modules_paper, "paper")

# cascade order: pseudotime where each LOESS curve first crosses Z = 0 (down for contractile, up for the rest), and its peak
cross <- function(d, up=TRUE) { s <- if (up) d$z > 0 else d$z < 0; i <- which(s)[1]; if (is.na(i) || i == 1) return(d$pt[1]); approx(d$z[c(i-1, i)], d$pt[c(i-1, i)], xout=0)$y }
order_tab <- function(L) do.call(rbind, lapply(split(L$loess, L$loess$module), function(d) {
  d <- d[order(d$pt), ]; up <- d$module[1] != "Contractile"
  data.frame(version=d$version[1], module=d$module[1], genes=NA, pt_cross_zero=cross(d, up), pt_peak=d$pt[which.max(if (up) d$z else -d$z)],
             z_start=d$z[1], z_end=d$z[nrow(d)], spearman_rho_vs_pseudotime=NA, p=NA)
}))
ot <- rbind(order_tab(ALT), order_tab(PAPER)); rownames(ot) <- NULL
for (i in seq_len(nrow(ot))) {
  mods <- if (ot$version[i] == "alternative") modules_alt else modules_paper; mn <- ot$module[i]
  ot$genes[i] <- paste(mods[[mn]]$genes, collapse=", ")
  pcv <- if (ot$version[i] == "alternative") ALT$per_cell else PAPER$per_cell
  ct <- cor.test(pcv[[paste0(mn, " | Z")]], blended, method="spearman", exact=FALSE); ot$spearman_rho_vs_pseudotime[i] <- unname(ct$estimate); ot$p[i] <- ct$p.value
}
ot <- ot[order(ot$version, ot$pt_cross_zero), ]
# the NFAT curve with and without Fas: how different?
a6 <- ALT$per_cell[["NFAT Effectors | Z"]]; a5 <- PAPER$per_cell[["NFAT-assoc. remodeling | Z"]]
l6 <- ALT$loess[ALT$loess$module == "NFAT Effectors", ]; l5 <- PAPER$loess[PAPER$loess$module == "NFAT-assoc. remodeling", ]
chk <- data.frame(metric=c("Pearson r, per-cell Z (6 vs 5 genes)", "Spearman rho, per-cell Z (6 vs 5 genes)", "max |LOESS difference| on grid (Z units)",
                           "mean |LOESS difference| on grid (Z units)", "Fas detected in (of 205 Piezo1+ CTRL SMCs)", "Fas mean log-expr, 205 cells"),
                  value=c(cor(a6, a5), cor(a6, a5, method="spearman"), max(abs(l6$z - l5$z)), mean(abs(l6$z - l5$z)),
                          sum(x["Fas", cells] > 0), mean(x["Fas", cells])))
print(ot[, c("version","module","pt_cross_zero","pt_peak","spearman_rho_vs_pseudotime","p")], digits=3); print(chk, digits=3)

write.csv(ot, file.path(OUT, "fig4F_cascade_order_paper_vs_alternative.csv"), row.names=FALSE)
write.csv(rbind(ALT$loess, PAPER$loess), file.path(OUT, "fig4F_loess_grid.csv"), row.names=FALSE)
write.csv(merge(ALT$per_cell, PAPER$per_cell[, c("cell", "NFAT-assoc. remodeling | mean expr", "NFAT-assoc. remodeling | Z")], by="cell"),
          file.path(OUT, "fig4F_per_cell_pseudotime_and_modules.csv"), row.names=FALSE)
write.csv(chk, file.path(OUT, "fig4F_nfat_module_with_without_Fas.csv"), row.names=FALSE)
write.csv(bg_df, file.path(OUT, "fig4F_age_fraction_grid.csv"), row.names=FALSE)

# ---------- plot ----------
cascade_plot <- function(L, modules, title) {
  lo <- L$loess; lo$module <- factor(lo$module, levels=names(modules)); cols <- sapply(modules, function(m) m$color)
  seg <- do.call(rbind, lapply(names(modules), function(mn) {
    d <- lo[lo$module == mn, ]; j <- seq_len(nrow(d) - 1)
    data.frame(x=d$pt[j], xend=d$pt[j+1], y=d$z[j], yend=d$z[j+1], module=mn, alpha_val=0.25 + 0.75 * (bg_df$age_frac[j] + bg_df$age_frac[j+1]) / 2)
  })); seg$module <- factor(seg$module, levels=names(modules))
  p_c <- ggplot() +
    geom_ribbon(data=lo, aes(x=pt, ymin=z_lo, ymax=z_hi, fill=module), alpha=0.08) + scale_fill_manual(values=cols, guide="none") +
    geom_segment(data=seg, aes(x=x, xend=xend, y=y, yend=yend, color=module, alpha=alpha_val), linewidth=2.2, lineend="round") +
    scale_color_manual(values=cols, name="Module") + scale_alpha_identity() +
    guides(color=guide_legend(override.aes=list(linewidth=4, alpha=1), keywidth=unit(2.2, "lines"), keyheight=unit(1.6, "lines"))) +
    geom_hline(yintercept=0, linetype="dashed", color="grey60", linewidth=0.3) + labs(y="Module Z-score") +
    scale_x_continuous(limits=range(blended), expand=c(0.01, 0.01)) + theme_minimal(base_size=14) +
    theme(legend.position="right", legend.text=element_text(size=22), legend.title=element_text(size=24, face="bold"), panel.grid.minor=element_blank(), axis.title.x=element_blank(),
          axis.text.x=element_blank(), plot.background=element_rect(fill="white", color=NA), plot.margin=margin(10, 10, 0, 10))
  cl <- pt_df$cluster; is_aged <- as.numeric(age == "Aged"); pct <- tapply(is_aged, cl, mean)[cl]
  set.seed(42); strip <- data.frame(pt=blended, y=runif(length(cells), -0.4, 0.4), pct_aged=as.numeric(pct)); strip <- strip[order(strip$pct_aged), ]
  p_s <- ggplot(strip, aes(pt, y, color=pct_aged)) + geom_point(size=1.2, alpha=0.6) +
    scale_color_gradient(low="#4DAF4A", high="#E41A1C", name="% Aged (cluster)", labels=function(v) paste0(round(100 * v), "%"), breaks=c(0.2, 0.4, 0.6)) +
    scale_x_continuous(limits=range(blended), expand=c(0.01, 0.01)) + scale_y_continuous(limits=c(-0.6, 0.6)) +
    labs(x="Pseudotime (Contractile → Synthetic)") + theme_minimal(base_size=14) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid=element_blank(),
          legend.position="right", legend.text=element_text(size=15), legend.title=element_text(size=16, face="bold"), legend.key.height=unit(0.75, "lines"), legend.key.width=unit(1.4, "lines"),
          plot.background=element_rect(fill="white", color=NA), plot.margin=margin(0, 10, 10, 10))
  tr <- ggdraw() + draw_label(title, fontface="bold", size=16, hjust=0, x=0.02)
  out <- plot_grid(tr, plot_grid(p_c, p_s, ncol=1, align="v", axis="lr", rel_heights=c(1, 0.25)), ncol=1, rel_heights=c(0.05, 1))
  attr(out, "legends") <- list(module=cowplot::get_legend(p_c), aged=cowplot::get_legend(p_s))
  out
}
p_paper <- cascade_plot(PAPER, modules_paper, sprintf("Piezo1+ CTRL SMCs: cascade (n=%d)", length(cells)))
p_alt <- cascade_plot(ALT, modules_alt, sprintf("6-gene NFAT module: Piezo1+ CTRL SMCs: cascade (n=%d)", length(cells)))
ggsave(file.path(FIG, "Fig4F.png"), p_paper, width=12, height=8, dpi=300, bg="white")
ggsave(file.path(FIG, "Fig4F.pdf"), p_paper, width=12, height=8, bg="white", device=cairo_pdf)
# the two legends on their own (vector PDF and 600 dpi PNG), to place at any size in the assembled figure
for (nm in names(attr(p_paper, "legends"))) {
  lg <- attr(p_paper, "legends")[[nm]]
  w <- grid::convertWidth(sum(lg$widths), "in", valueOnly=TRUE) + 0.3; h <- grid::convertHeight(sum(lg$heights), "in", valueOnly=TRUE) + 0.3
  for (ext in c("png", "pdf")) ggsave(file.path(FIG, sprintf("Fig4F_legend_%s.%s", nm, ext)), ggdraw(lg), width=w, height=h, dpi=600, bg="white",
                                      device=if (ext == "pdf") cairo_pdf else "png")
}
ggsave(file.path(FIG, "Fig4F_6gene_NFAT_module.png"), p_alt, width=12, height=8, dpi=150, bg="white")
# NFAT module only, 5 and 6 genes overlaid
ov <- rbind(transform(l6, which="6 genes, with Fas"), transform(l5, which="5 genes, without Fas"))
p_ov <- ggplot(ov, aes(pt, z, color=which, fill=which)) + geom_ribbon(aes(ymin=z_lo, ymax=z_hi), alpha=0.12, color=NA) + geom_line(linewidth=1.6) +
  scale_color_manual(values=c("#9467BD", "#4B2A80"), name=NULL) + scale_fill_manual(values=c("#9467BD", "#4B2A80"), name=NULL) +
  geom_hline(yintercept=0, linetype="dashed", color="grey60") + labs(x="Pseudotime (Contractile → Synthetic)", y="Module Z-score",
  title="Fig 4F NFAT module with and without Fas", subtitle=sprintf("205 Piezo1+ CTRL SMCs. Per-cell scores correlate r = %.2f;\nlargest difference between the curves %.2f Z", chk$value[1], chk$value[3])) +
  theme_minimal(base_size=13) + theme(legend.position="bottom", plot.background=element_rect(fill="white", color=NA))
ggsave(file.path(FIG, "Fig4F_NFAT_module_with_without_Fas.png"), p_ov, width=8, height=5.5, dpi=200, bg="white")
cat("done\n")
