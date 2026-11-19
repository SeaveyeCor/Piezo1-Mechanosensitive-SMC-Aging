# ==============================================================================
# PAPER PANEL: Fig S4B
# PLOTS:  GO biological-process programme scores, Piezo1+ vs Piezo1- cells
# CELLS:  CTRL synthetic-classified visceral SMCs (Young + Aged CTRL)
# GROUPS COMPARED: Piezo1+ vs Piezo1-
# TEST:   Wilcoxon
# ==============================================================================
# GO term scores (unbiased programme definitions)
# ==============================================================================

setwd(".")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
})

source("pipeline/utils/plotting_theme.R")
source("pipeline/core/gene_sets.R")
source("pipeline/selection_criteria/cluster_based_smc.R")

cat("Loading Seurat...\n")
obj <- readRDS("data/seurat_paper_6samples.rds")
sel <- select_cells(obj)
smc <- sel$smc
meta <- obj@meta.data[smc, ]
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
gs <- get_gene_sets()

# Classify phenotype
contr_score <- colMeans(as.matrix(expr[gs$contractile_core[gs$contractile_core %in% rownames(expr)], smc, drop = FALSE]))
synth_score <- colMeans(as.matrix(expr[gs$synthetic_core[gs$synthetic_core %in% rownames(expr)], smc, drop = FALSE]))
ratio <- (contr_score + 0.001) / (synth_score + 0.001)
is_synthetic <- ratio <= 1.5

# Ctrl synthetic SMCs, P1+ vs P1-
ctrl_synth <- smc[meta$genotype == "Ctrl" & is_synthetic]
p1_expr <- as.numeric(expr["Piezo1", ctrl_synth])
p1_pos <- ctrl_synth[p1_expr > 0]
p1_neg <- ctrl_synth[p1_expr == 0]

cat("Ctrl Synthetic: ", length(ctrl_synth), " (P1+=", length(p1_pos), ", P1-=", length(p1_neg), ")\n\n")

# GO gene sets
go_sets <- list(
  "Muscle contraction\n(GO:0006936)" = c(
    "Myh11", "Acta2", "Actg2", "Myl9", "Myl6", "Mylk", "Cnn1", "Cald1",
    "Tagln", "Tpm1", "Tpm2", "Des", "Smtn", "Lmod1", "Tnnt2", "Tnnc1"
  ),
  "Ca2+ transport\n(GO:0006816)" = c(
    "Atp2a2", "Atp2b1", "Itpr1", "Itpr2", "Itpr3", "Ryr2",
    "Orai1", "Stim1", "Stim2", "Slc8a1", "Trpc1", "Trpc6", "Trpv4", "Trpm7"
  ),
  "Calcineurin-NFAT\n(GO:0033173)" = c(
    "Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4", "Nfat5",
    "Ppp3ca", "Ppp3cb", "Ppp3cc", "Rcan1", "Rcan2", "Rcan3",
    "Calm1", "Calm2", "Calm3"
  ),
  "ECM organization\n(GO:0030198)" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col6a1", "Fn1", "Bgn", "Dcn",
    "Sparc", "Lum", "Mmp2", "Mmp9", "Timp1", "Timp2", "Adamts1"
  ),
  "Cell migration\n(GO:0030334)" = c(
    "Vim", "Cd44", "S100a4", "Vcam1", "Icam1",
    "Mmp2", "Mmp9", "Fn1", "Spp1", "Lgals3"
  ),
  "Cell proliferation\n(GO:0008283)" = c(
    "Mki67", "Pcna", "Top2a", "Klf4", "Klf5",
    "Cdkn1a", "Ccnd1", "Ccnd2", "Myc"
  ),
  "Inflammatory response\n(GO:0006954)" = c(
    "C3", "C4b", "Il6", "Cxcl1", "Ccl2", "Ccl7",
    "Vcam1", "Icam1", "Fas", "Tnfrsf12a"
  ),
  "Cytoskeleton organization\n(GO:0007010)" = c(
    "Ezr", "Msn", "Rdx", "Vim", "Lmna", "Lmnb1",
    "Syne1", "Syne2", "Actn1", "Actn4", "Flna", "Flnb",
    "Myl9", "Myl6", "Cfl1", "Pfn1", "Tpm1", "Tpm2"
  ),
  "Cell adhesion\n(GO:0007155)" = c(
    "Vcam1", "Icam1", "Cd44", "Fn1", "Col1a1", "Col1a2",
    "Itga5", "Itgb1", "Itgav", "Spp1"
  )
)

# Clean names for data export (no linebreaks)
go_names_clean <- gsub("\n", " ", names(go_sets))

# Compute log2FC and stats
results <- data.frame()
gene_details <- data.frame()

for (i in seq_along(go_sets)) {
  go_name <- names(go_sets)[i]
  go_clean <- go_names_clean[i]
  genes <- go_sets[[i]]
  avail <- genes[genes %in% rownames(expr)]
  if (length(avail) < 3) next

  sc_pos <- colMeans(as.matrix(expr[avail, p1_pos, drop = FALSE]))
  sc_neg <- colMeans(as.matrix(expr[avail, p1_neg, drop = FALSE]))

  m_pos <- mean(sc_pos)
  m_neg <- mean(sc_neg)
  sd_pos <- sd(sc_pos)
  sd_neg <- sd(sc_neg)
  l2fc <- log2((m_pos + 0.01) / (m_neg + 0.01))
  wt <- wilcox.test(sc_pos, sc_neg)
  sig <- ifelse(wt$p.value < 0.001, "***", ifelse(wt$p.value < 0.01, "**", ifelse(wt$p.value < 0.05, "*", "ns")))

  # SE of log2FC via bootstrap
  set.seed(42)
  boot_l2fc <- replicate(1000, {
    bp <- sample(sc_pos, replace = TRUE)
    bn <- sample(sc_neg, replace = TRUE)
    log2((mean(bp) + 0.01) / (mean(bn) + 0.01))
  })
  l2fc_se <- sd(boot_l2fc)

  results <- rbind(results, data.frame(
    go_term = go_name,
    go_term_clean = go_clean,
    n_genes = length(avail),
    genes_used = paste(avail, collapse = ", "),
    mean_p1pos = round(m_pos, 4),
    sd_p1pos = round(sd_pos, 4),
    mean_p1neg = round(m_neg, 4),
    sd_p1neg = round(sd_neg, 4),
    log2fc = round(l2fc, 3),
    log2fc_se = round(l2fc_se, 3),
    p_value = wt$p.value,
    significance = sig,
    stringsAsFactors = FALSE
  ))

  # Per-gene detail
  for (g in avail) {
    gp <- as.numeric(expr[g, p1_pos])
    gn <- as.numeric(expr[g, p1_neg])
    gl2 <- log2((mean(gp) + 0.01) / (mean(gn) + 0.01))
    gw <- wilcox.test(gp, gn)
    gsig <- ifelse(gw$p.value < 0.001, "***", ifelse(gw$p.value < 0.01, "**", ifelse(gw$p.value < 0.05, "*", "ns")))
    gene_details <- rbind(gene_details, data.frame(
      go_term = go_clean,
      gene = g,
      mean_p1pos = round(mean(gp), 4),
      sd_p1pos = round(sd(gp), 4),
      mean_p1neg = round(mean(gn), 4),
      sd_p1neg = round(sd(gn), 4),
      log2fc = round(gl2, 3),
      p_value = gw$p.value,
      sig = gsig
    ))
  }
}

# =====================================================================
# PLOT
# =====================================================================

plot_df <- results
plot_df$alpha_val <- ifelse(plot_df$significance != "ns", 1, 0.4)

# Create rich text labels: bold name, plain GO ID
plot_df$label_rich <- gsub("^(.+)\n\\((.+)\\)$", "**\\1**<br><span style='font-size:9pt'>(\\2)</span>", plot_df$go_term)
plot_df$label_rich <- factor(plot_df$label_rich, levels = rev(plot_df$label_rich))
plot_df$go_term_fct <- factor(plot_df$go_term, levels = rev(results$go_term))

# Assign distinct colors per GO term
go_colors <- c(
  "Muscle contraction\n(GO:0006936)" = "#66C2A5",
  "Ca2+ transport\n(GO:0006816)" = "#FC8D62",
  "Calcineurin-NFAT\n(GO:0033173)" = "#8DA0CB",
  "ECM organization\n(GO:0030198)" = "#E78AC3",
  "Cell migration\n(GO:0030334)" = "#A6D854",
  "Cell proliferation\n(GO:0008283)" = "#FFD92F",
  "Inflammatory response\n(GO:0006954)" = "#E5C494",
  "Cytoskeleton organization\n(GO:0007010)" = "#B3B3B3",
  "Cell adhesion\n(GO:0007155)" = "#66C2A5"
)

p <- ggplot(plot_df, aes(x = log2fc, y = label_rich)) +
  geom_col(aes(fill = go_term_fct, alpha = alpha_val),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = log2fc - log2fc_se, xmax = log2fc + log2fc_se),
                 height = 0.2, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_text(aes(label = ifelse(significance != "ns", significance, ""),
                x = ifelse(log2fc > 0, log2fc + log2fc_se + 0.02,
                           log2fc - log2fc_se - 0.02)),
            hjust = ifelse(plot_df$log2fc > 0, 0, 1),
            size = 5, fontface = "bold") +
  scale_fill_manual(values = go_colors) +
  scale_alpha_identity() +
  guides(fill = "none") +
  labs(
    x = expression(log[2]*FC ~ "(P1+ / P1-)"),
    y = NULL
  ) +
  theme_publication(base_size = 14) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 11, lineheight = 1.2),
    plot.margin = margin(10, 20, 10, 10)
  )

# Save
out_dir <- "output/figures/FigS4B"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "FigS4B_go_programs.png"), p,
       width = 8, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "FigS4B_go_programs.pdf"), p,
       width = 8, height = 5, bg = "white")


# Save data xlsx
wb <- createWorkbook()
addWorksheet(wb, "title")
writeData(wb, "title", data.frame(
  Field = c("Figure", "Title", "Manuscript"),
  Value = c("Fig S4B", "GO programs, Piezo1+ vs Piezo1- CTRL synthetic SMCs",
            "Maladaptive Piezo1 Mechanotransduction Drives Smooth Muscle Aging in the Gut")
))
addWorksheet(wb, "module_level_stats")
writeData(wb, "module_level_stats", results[, c("go_term_clean", "n_genes", "genes_used",
  "mean_p1pos", "sd_p1pos", "mean_p1neg", "sd_p1neg",
  "log2fc", "log2fc_se", "p_value", "significance")])
addWorksheet(wb, "gene_level_stats")
writeData(wb, "gene_level_stats", gene_details)
addWorksheet(wb, "population_summary")
writeData(wb, "population_summary", data.frame(
  population = c("Ctrl_synthetic", "P1_pos", "P1_neg"),
  n_cells = c(length(ctrl_synth), length(p1_pos), length(p1_neg))
))

sd_out <- "output/source_data"
saveWorkbook(wb, file.path(sd_out, "FigS4B_go_programs_data.xlsx"), overwrite = TRUE)

cat("\nSaved:\n")
cat("  Figure: output/figures/FigS4B/FigS4B_go_programs.png/pdf\n")
cat("  Data:   FigS4B_go_programs_data.xlsx\n")
cat("\nDone.\n")
