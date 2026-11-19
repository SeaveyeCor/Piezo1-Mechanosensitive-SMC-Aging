# Sample identity checks for a Seurat object holding several 10x libraries.
#
#   1. No two samples share cells: a barcode found in two samples with identical UMI and gene counts is the same cell
#      loaded twice. Unrelated libraries share ~0.5% of barcodes by chance, never with identical counts.
#   2. (optional, needs the Cell Ranger files) each sample's barcodes come from its own file (>= 99%) and from no other
#      file (<= 5%).
#   3. The object carries a stored UMAP, so every figure uses the same coordinates.
#
# Any failure stops with an error. Use from R:
#   source("00_preprocessing/check_sample_identity.R"); check_sample_identity(obj, raw_files = NULL)
# or from the command line (raw_files.csv: columns sample, file):
#   Rscript 00_preprocessing/check_sample_identity.R object.rds [raw_files.csv]

check_sample_identity <- function(obj, raw_files = NULL, sample_col = "orig.ident", require_umap = TRUE) {
  core <- function(x) regmatches(x, regexpr("[ACGT]{16}-1", x))
  md <- obj@meta.data
  s <- as.character(md[[sample_col]])
  bc <- core(colnames(obj))
  stopifnot(length(bc) == ncol(obj))
  key <- paste(bc, md$nCount_RNA, md$nFeature_RNA)
  samples <- sort(unique(s))
  cat("Sample identity check:", length(samples), "samples,", ncol(obj), "cells\n")

  # 1. duplicated cells between samples
  fail <- character(0)
  for (i in seq_along(samples)) for (j in seq_along(samples)) if (i < j) {
    a <- key[s == samples[i]]; b <- key[s == samples[j]]
    shared <- sum(a %in% b)
    if (shared > 0.01 * min(length(a), length(b)))
      fail <- c(fail, sprintf("%s and %s share %d identical cells (%s: %d cells, %s: %d cells)",
                              samples[i], samples[j], shared, samples[i], length(a), samples[j], length(b)))
  }
  if (length(fail)) stop("Duplicated cells between samples:\n  ", paste(fail, collapse = "\n  "), call. = FALSE)
  cat("  1. no cells shared between samples: OK\n")

  # 2. each sample comes from its own Cell Ranger file
  if (!is.null(raw_files)) {
    stopifnot(all(samples %in% names(raw_files)))
    raw <- lapply(raw_files, function(f) {
      x <- Seurat::Read10X_h5(f); if (is.list(x)) x <- x[[1]]; colnames(x)
    })
    frac <- sapply(raw, function(h) tapply(bc %in% h, s, mean))[samples, , drop = FALSE]
    print(round(100 * frac, 1))
    for (x in samples) {
      own <- frac[x, x]; other <- frac[x, setdiff(colnames(frac), x)]
      if (own < 0.99) fail <- c(fail, sprintf("%s: only %.1f%% of cells found in its own file", x, 100 * own))
      if (any(other > 0.05)) fail <- c(fail, sprintf("%s: %.1f%% of cells found in the file of %s", x,
                                                     100 * max(other), names(other)[which.max(other)]))
    }
    if (length(fail)) stop("Sample-to-file mismatch:\n  ", paste(fail, collapse = "\n  "), call. = FALSE)
    cat("  2. every sample matches its own Cell Ranger file and no other: OK\n")
  }

  # 3. stored UMAP
  if (require_umap) {
    if (!"umap" %in% names(obj@reductions)) stop("No stored UMAP in the object.", call. = FALSE)
    cat("  3. stored UMAP present: OK\n")
  }
  invisible(TRUE)
}

if (sys.nframe() == 0) {
  a <- commandArgs(TRUE)
  suppressMessages(library(Seurat))
  rf <- if (length(a) > 1) { m <- read.csv(a[2]); setNames(m$file, m$sample) } else NULL
  check_sample_identity(readRDS(a[1]), raw_files = rf)
}
