# Post-pipeline source-data steps. Run from the repository root after 01_pipeline/run_complete_analysis.R.
# Order matters: add_sd_sem_n.R must run before source_data_synthetic_score.R, so that 3F uses the analysed SMCs.
library(openxlsx)
run_dir <- tail(sort(Sys.glob("pipeline_results/cluster_based_noEpithelial_*")), 1)
stopifnot(length(run_dir) == 1)
sd <- "output/source_data"
dir.create(sd, recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures/FigS4B", recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(run_dir, "Supplemental_Tables/Table_S7_Piezo1_Correlations.csv"), "output/tables/", overwrite = TRUE)
# the cell-level scripts (04_cell_level/) read the pipeline workbooks from output/pipeline/<module>/
for (f in Sys.glob(file.path(run_dir, "*", "*_data.xlsx"))) {
  d <- file.path("output/pipeline", basename(dirname(f))); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(f, d, overwrite = TRUE)
}

# manuscript panel -> pipeline module folder, and the long workbook name
panels <- data.frame(
  panel  = c("Fig3A","Fig3B","Fig3C","Fig3D","Fig3E","Fig3G","Fig3H","Fig3M","Fig3N","Fig4A","Fig4B","Fig4C","Fig4D",
             "Fig4F","FigS4B","Fig4G","FigS3A","FigS3C","FigS3D","FigS3F","FigS3J"),
  module = c("Fig3E","Fig3K","Fig3F","Fig3G","Fig3H","Fig3I","Fig3J","FigS3K","FigS4H","FigS3G","Fig4F","Fig4G","Fig4E",
             "Fig4I","Fig4J","Fig4H","FigS3D","FigS3F","FigS3F","FigS3I","FigS4A"),
  long   = c("phenotype_umap","phenotype_proportions","piezo1_vs_contractile","piezo1_vs_synthetic","contractile_score",
             "contractile_gene_rescue","synthetic_gene_rescue","ecm_score","cyto_remodeling","pathway_correlations",
             "piezo1_vs_nfat","synthetic_vs_nfat","nfat_activity_score","pseudotime_cascade","p1pos_vs_p1neg",
             "phate_trajectory","celltype_umap","piezo1_ko_validation","piezo1_pct_expressing","gsea_ko_vs_ctrl",
             "ecm_gene_heatmap"))
for (i in seq_len(nrow(panels)))
  file.copy(file.path(run_dir, panels$module[i], paste0(panels$module[i], "_data.xlsx")),
            file.path(sd, paste0(panels$panel[i], "_data.xlsx")), overwrite = TRUE)

step <- function(f) { cat("\n====", f, "\n"); source(file.path("02_source_data_postprocessing", f), local = new.env()) }
# Stage 1 works on the short names (FigXX_data.xlsx)
step("source_data_FigS3F_gsea.R")
step("add_sd_sem_n.R")
step("source_data_synthetic_score.R")
step("add_panel_titles.R")
# Stage 2 works on the long names
long <- setNames(c(panels$long, "synthetic_score"), c(panels$panel, "Fig3F"))
for (k in names(long)) {
  s <- file.path(sd, paste0(k, "_data.xlsx"))
  if (file.exists(s)) file.copy(s, file.path(sd, sprintf("%s_%s_data.xlsx", k, long[[k]])), overwrite = TRUE)
}
for (f in c("add_cell_level_sheets.R", "source_data_supplementary.R", "add_sd_module_stats.R",
            "add_regression_stats.R", "plot_FigS4B_go_programs.R", "plot_FigS3F_gsea.R")) step(f)
cat("\nSource data written to", sd, "\n")
