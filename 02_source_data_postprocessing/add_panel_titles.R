library(openxlsx)

sd <- "output/source_data"

# Source-data file -> panel title (paper panel letters)
titles <- list(
  "Fig3A_data.xlsx" = "SMC Phenotype UMAP",
  "Fig3B_data.xlsx" = "Phenotype Proportions",
  "Fig3C_data.xlsx" = "Piezo1 vs Contractile Score",
  "Fig3D_data.xlsx" = "Piezo1 vs Synthetic Score",
  "Fig3E_data.xlsx" = "Contractile Score",
  "Fig3F_data.xlsx" = "Synthetic Score",
  "Fig3G_data.xlsx" = "Contractile Gene Changes",
  "Fig3H_data.xlsx" = "Synthetic Gene Changes",
  "Fig3M_data.xlsx" = "ECM Score",
  "Fig3N_data.xlsx" = "Cytoskeletal Remodeling Score",
  "Fig4A_data.xlsx" = "Pathway Correlations with Piezo1",
  "Fig4B_data.xlsx" = "Piezo1 vs NFAT Pathway Score",
  "Fig4C_data.xlsx" = "Synthetic vs NFAT Pathway Score",
  "Fig4D_data.xlsx" = "NFAT-associated Remodeling Score",
  "Fig4F_data.xlsx" = "Gene-module Dynamics on Pseudotime",
  "Fig4G_data.xlsx" = "PHATE Trajectory",
  "FigS3A_data.xlsx" = "Cell Type UMAP",
  "FigS3C_data.xlsx" = "Piezo1 Log2 Fold Change",
  "FigS3D_data.xlsx" = "Piezo1+ SMC Fraction",
  "FigS3F_data.xlsx" = "GSEA, Aged P1KO vs Aged CTRL",
  "FigS3J_data.xlsx" = "ECM Gene Heatmap",
  "FigS4B_data.xlsx" = "Piezo1+ vs Piezo1- Synthetic SMCs"
)

for (fn in names(titles)) {
  fpath <- file.path(sd, fn)
  if (!file.exists(fpath)) {
    cat("SKIP (not found):", fn, "\n")
    next
  }

  title <- titles[[fn]]
  fig_id <- gsub("_data\\.xlsx$", "", fn)

  wb <- tryCatch(loadWorkbook(fpath), error = function(e) NULL)
  if (is.null(wb)) {
    cat("SKIP (load error):", fn, "\n")
    next
  }

  # Add or replace title sheet
  if ("title" %in% names(wb)) removeWorksheet(wb, "title")
  addWorksheet(wb, "title", gridLines = FALSE)
  writeData(wb, "title", data.frame(
    Field = c("Figure", "Title", "Manuscript"),
    Value = c(fig_id, title,
              "Maladaptive Piezo1 Mechanotransduction Drives Smooth Muscle Aging in the Gut")
  ))

  # Move title to first position
  sheet_order <- seq_along(names(wb))
  title_idx <- which(names(wb) == "title")
  other_idx <- setdiff(sheet_order, title_idx)
  worksheetOrder(wb) <- c(title_idx, other_idx)

  saveWorkbook(wb, fpath, overwrite = TRUE)
  cat("Added title:", fn, "->", title, "\n")
}
