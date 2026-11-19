library(Seurat)
cat("=== SMC METADATA EXPORT FOR RNA VELOCITY ===\n")

obj <- readRDS("./data/seurat_paper_6samples.rds")
obj$sample <- as.character(obj$sample)
obj$sample[is.na(obj$sample)] <- as.character(obj$orig.ident[is.na(obj$sample)])

# ── 1. Subset to the 6 samples ─────────────────────────────────────────────
used <- c("WT_Young_1", "WT_Young_2", "WT_Aged_1", "WT_Aged_2",
          "KO_Young_1", "KO_Aged")
obj <- subset(obj, sample %in% used)
cat(sprintf("After sample selection: %d cells\n", ncol(obj)))

# ── 2. Depth match ──────────────────────────────────────────────────────────
sample_medians <- sapply(used, function(s) median(obj$nCount_RNA[obj$sample == s]))
target <- median(sample_medians)
lo <- target / 2; hi <- target * 2
outliers <- names(sample_medians)[sample_medians < lo | sample_medians > hi]

cat(sprintf("Depth match: target=%.0f, range=[%.0f, %.0f]\n", target, lo, hi))
cat(sprintf("Outliers: %s\n", paste(outliers, collapse=", ")))

is_outlier <- obj$sample %in% outliers
in_range <- obj$nCount_RNA >= lo & obj$nCount_RNA <= hi
keep <- !is_outlier | (is_outlier & in_range)
obj <- subset(obj, cells = colnames(obj)[keep])
cat(sprintf("After depth match: %d cells\n", ncol(obj)))
for (s in used) {
  cat(sprintf("  %s: %d cells\n", s, sum(obj$sample == s)))
}

# ── 3. Normalize (log-normalized data for the phenotype scores below) ───────
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
obj <- NormalizeData(obj, verbose = FALSE)

# ── 4. Bona fide filter ─────────────────────────────────────────────────────
cat("\n=== BONA FIDE SMC SELECTION ===\n")
smc <- obj$cell_type == "SMC_Contractile"
cl <- as.character(obj$seurat_clusters)
bonafide <- smc & !is.na(cl) & cl %in% c("1", "5", "15", "16")

cat("Bona fide SMCs per sample:\n")
for (s in used) {
  n_bf <- sum(bonafide & obj$sample == s)
  n_smc <- sum(smc & obj$sample == s)
  cat(sprintf("  %s: %d bona fide / %d SMC\n", s, n_bf, n_smc))
}
cat(sprintf("Total bona fide: %d\n", sum(bonafide)))

# Doublet removal placeholder (we don't have the doublet file here, skip)
# The bonafide pipeline removes 63 doublets; we note this

# ── 5. Phenotype check ──────────────────────────────────────────────────────
cat("\n=== PHENOTYPE ===\n")
grp <- ifelse(grepl("^WT_Young", obj$sample), "Young_Ctrl",
       ifelse(grepl("^WT_Aged", obj$sample), "Aged_Ctrl",
       ifelse(grepl("^KO_Young", obj$sample), "Young_KO", "Aged_KO")))

for (g in c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")) {
  mask <- bonafide & grp == g
  n <- sum(mask)
  # Quick phenotype classification
  contr_genes <- c("Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn")
  synth_genes <- c("Col1a1", "Fn1", "Spp1", "Klf4", "Vim")
  avail_c <- contr_genes[contr_genes %in% rownames(obj)]
  avail_s <- synth_genes[synth_genes %in% rownames(obj)]

  barcodes <- colnames(obj)[mask]
  if (length(barcodes) > 0) {
    c_score <- colMeans(as.matrix(GetAssayData(obj, layer = "data")[avail_c, barcodes]))
    s_score <- colMeans(as.matrix(GetAssayData(obj, layer = "data")[avail_s, barcodes]))
    pct_contr <- 100 * sum(c_score / (s_score + 0.01) > 1.5) / length(c_score)
    cat(sprintf("  %s: %d cells, %.1f%% contractile\n", g, n, pct_contr))
  }
}

# ── 6. Export metadata ───────────────────────────────────────────────────────
cat("\n=== EXPORT ===\n")
bc_clean <- gsub("^.*_([ACGT]+)-[0-9]+$", "\\1", colnames(obj))

meta <- data.frame(
  seurat_barcode = colnames(obj),
  barcode = bc_clean,
  sample = obj$sample,
  group = grp,
  cell_type = obj$cell_type,
  cluster = as.character(obj$seurat_clusters),
  is_smc = smc,
  is_bonafide = bonafide,
  nCount_RNA = obj$nCount_RNA,
  nFeature_RNA = obj$nFeature_RNA,
  stringsAsFactors = FALSE
)

outfile <- "data/velocity/scvelo_ready/velocity_metadata_corrected.csv"
write.csv(meta, outfile, row.names = FALSE)

vel_smc <- meta[meta$is_bonafide, ]
cat(sprintf("Exported %d cells (%d bona fide SMCs) to %s\n",
            nrow(meta), nrow(vel_smc), outfile))

cat(sprintf("\nBona fide per group:\n"))
print(table(vel_smc$group))
cat(sprintf("\nBona fide per sample:\n"))
print(table(vel_smc$sample))
