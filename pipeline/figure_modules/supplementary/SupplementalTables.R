# ==============================================================================
# SUPPLEMENTAL TABLES GENERATOR
# ==============================================================================
# Generates the supplementary tables:
# - Table S1: Cell type markers used for identification
# - Table S2a: Phenotype classification genes (7 contractile + 5 synthetic)
# - Table S2b: Gene scoring panels (contractile capacity + synthetic phenotype)
# - Table S3: Pathway gene sets (NFAT, ECM scores)
# - Table S4: Gene expression by experimental group
# - Table S5: Contractile gene rescue statistics
# - Table S6: Synthetic gene rescue statistics
# - Table S7: Piezo1 gene correlations
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

#' Generate all supplemental tables
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/statistics.R"))
  source(file.path(pipeline_root, "core/gene_sets.R"))

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Create output directory
  table_dir <- file.path(context$output_dir, "Supplemental_Tables")
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

  cat("  Generating Supplemental Tables...\n")

  # ==========================================================================
  # TABLE S1: Cell Type Markers
  # ==========================================================================
  cat("  Creating Table S1: Cell type markers...\n")

  table_s1 <- create_table_s1_markers(expr, smc_cells, metadata)

  # ==========================================================================
  # TABLE S2a: Phenotype Classification Genes
  # ==========================================================================
  cat("  Creating Table S2a: Phenotype classification genes...\n")

  gene_sets <- get_gene_sets()
  table_s2a <- create_table_s2a_classification(gene_sets)

  # ==========================================================================
  # TABLE S2b: Gene Scoring Panels
  # ==========================================================================
  cat("  Creating Table S2b: Gene scoring panels...\n")

  table_s2b <- create_table_s2b_scoring(gene_sets)

  # ==========================================================================
  # TABLE S3: Pathway Gene Sets
  # ==========================================================================
  cat("  Creating Table S3: Pathway gene sets...\n")

  table_s3 <- create_table_s3_pathway_genes(expr, smc_cells, metadata, gene_sets)

  # ==========================================================================
  # TABLE S4: Gene Expression by Group
  # ==========================================================================
  cat("  Creating Table S4: Gene expression levels...\n")

  table_s4 <- create_table_s4_expression(expr, smc_cells, metadata, gene_sets)

  # ==========================================================================
  # TABLE S5: Contractile Gene Rescue Statistics
  # ==========================================================================
  cat("  Creating Table S5: Contractile gene rescue...\n")

  table_s5 <- create_table_s5_contractile(expr, smc_cells, metadata, gene_sets)

  # ==========================================================================
  # TABLE S6: Synthetic Gene Rescue Statistics
  # ==========================================================================
  cat("  Creating Table S6: Synthetic gene rescue...\n")

  table_s6 <- create_table_s6_synthetic(expr, smc_cells, metadata, gene_sets)

  # ==========================================================================
  # TABLE S7: Piezo1 Gene Correlations
  # ==========================================================================
  cat("  Creating Table S7: Piezo1 gene correlations...\n")

  table_s7 <- create_table_s7_piezo1_correlations(expr, smc_cells, metadata, gene_sets)

  # ==========================================================================
  # Tables: S0 Sample QC, S0b Piezo1 KO Validation,
  # S1c Cluster x Group, S4b Panel Scores, S4c Phenotype Proportions
  # ==========================================================================
  cat("  Creating Table S0: Sample inclusion + QC...\n")
  table_s0  <- create_table_s0_sample_inclusion(seurat_obj, context)
  cat("  Creating Table S0b: Piezo1 KO validation per sample...\n")
  table_s0b <- create_table_s0b_piezo1_validation(expr, smc_cells, metadata)
  cat("  Creating Table S1c: Cluster composition by group...\n")
  table_s1c <- create_table_s1c_cluster_by_group(metadata)
  cat("  Creating Table S4b: Panel scores by group...\n")
  table_s4b <- create_table_s4b_panel_scores(expr, smc_cells, metadata, gene_sets)
  cat("  Creating Table S4c: Phenotype proportions...\n")
  table_s4c <- create_table_s4c_phenotype_proportions(metadata, smc_cells, context)

  # ==========================================================================
  # Save all tables to single Excel workbook
  # ==========================================================================
  cat("  Saving tables to Excel...\n")

  wb <- createWorkbook()

  # Table S0 (Sample Inclusion + QC)
  addWorksheet(wb, "S0 Sample Inclusion")
  writeData(wb, "S0 Sample Inclusion", table_s0$data)

  # Table S0b (Piezo1 KO Validation per sample)
  addWorksheet(wb, "S0b Piezo1 KO Validation")
  writeData(wb, "S0b Piezo1 KO Validation", table_s0b$data)

  # Table S1
  addWorksheet(wb, "S1 CellType Markers")
  writeData(wb, "S1 CellType Markers", table_s1$data)

  # Table S1c (Cluster by Group)
  addWorksheet(wb, "S1c Cluster by Group")
  writeData(wb, "S1c Cluster by Group", table_s1c$data)

  # Table S2a (Classification Genes)
  addWorksheet(wb, "S2a Classification Genes")
  writeData(wb, "S2a Classification Genes", table_s2a$data)

  # Table S2b (Scoring Panels)
  addWorksheet(wb, "S2b Scoring Panels")
  writeData(wb, "S2b Scoring Panels", table_s2b$data)

  # Table S3 (Pathway Gene Sets)
  addWorksheet(wb, "S3 Pathway GeneSets")
  writeData(wb, "S3 Pathway GeneSets", table_s3$data)

  # Table S4 (Gene Expression)
  addWorksheet(wb, "S4 GeneExpression")
  writeData(wb, "S4 GeneExpression", table_s4$data)

  # Table S4b (Panel Scores by Group)
  addWorksheet(wb, "S4b Panel Scores")
  writeData(wb, "S4b Panel Scores", table_s4b$data)

  # Table S4c (Phenotype Proportions)
  addWorksheet(wb, "S4c Phenotype Proportions")
  writeData(wb, "S4c Phenotype Proportions", table_s4c$data)

  # Table S5 (Contractile Rescue)
  addWorksheet(wb, "S5 Contractile Rescue")
  writeData(wb, "S5 Contractile Rescue", table_s5$data)

  # Table S6 (Synthetic Rescue)
  addWorksheet(wb, "S6 Synthetic Rescue")
  writeData(wb, "S6 Synthetic Rescue", table_s6$data)

  # Table S7 (Piezo1 Correlations)
  addWorksheet(wb, "S7 Piezo1 Correlations")
  writeData(wb, "S7 Piezo1 Correlations", table_s7$data)

  # Table legends sheet
  addWorksheet(wb, "Table Legends")
  legends <- data.frame(
    Table = c("S0", "S0b", "S1", "S1c", "S2a", "S2b", "S3", "S4", "S4b", "S4c", "S5", "S6", "S7"),
    Title = c(
      "Sample Inclusion + QC",
      "Piezo1 KO Validation per Sample",
      "Cell Type Markers",
      "Cluster Composition by Group",
      "Phenotype Classification Genes",
      "Gene Scoring Panels",
      "Pathway Gene Sets (NFAT, ECM)",
      "Gene Expression by Experimental Group",
      "Panel Scores by Group",
      "Phenotype Proportions per Group + Cl5 Sub-cluster Split",
      "Contractile Gene Rescue in High-Confidence SMCs",
      "Synthetic Gene Rescue in High-Confidence SMCs",
      "Piezo1 Gene Correlations"
    ),
    Description = c(
      table_s0$legend,
      table_s0b$legend,
      table_s1$legend,
      table_s1c$legend,
      table_s2a$legend,
      table_s2b$legend,
      table_s3$legend,
      table_s4$legend,
      table_s4b$legend,
      table_s4c$legend,
      table_s5$legend,
      table_s6$legend,
      table_s7$legend
    )
  )
  writeData(wb, "Table Legends", legends)

  saveWorkbook(wb, file.path(table_dir, "Supplemental_Tables.xlsx"), overwrite = TRUE)

  # Tables S0-S7 in the manuscript layout
  source(file.path(pipeline_root, "figure_modules/supplementary/PaperTables_S0_S7.R"), local = TRUE)
  write_paper_tables_S0_S7(seurat_obj, context, expr, smc_cells, metadata, gene_sets, table_s0b, table_s1c, table_s4c, table_dir)

  # Also save individual CSV files
  write.csv(table_s0$data,  file.path(table_dir, "Table_S0_Sample_Inclusion.csv"), row.names = FALSE)
  write.csv(table_s0b$data, file.path(table_dir, "Table_S0b_Piezo1_KO_Validation.csv"), row.names = FALSE)
  write.csv(table_s1$data,  file.path(table_dir, "Table_S1_CellType_Markers.csv"), row.names = FALSE)
  write.csv(table_s1c$data, file.path(table_dir, "Table_S1c_Cluster_by_Group.csv"), row.names = FALSE)
  write.csv(table_s2a$data, file.path(table_dir, "Table_S2a_Classification_Genes.csv"), row.names = FALSE)
  write.csv(table_s2b$data, file.path(table_dir, "Table_S2b_Scoring_Panels.csv"), row.names = FALSE)
  write.csv(table_s3$data,  file.path(table_dir, "Table_S3_Pathway_GeneSets.csv"), row.names = FALSE)
  write.csv(table_s4$data,  file.path(table_dir, "Table_S4_GeneExpression.csv"), row.names = FALSE)
  write.csv(table_s4b$data, file.path(table_dir, "Table_S4b_Panel_Scores_by_Group.csv"), row.names = FALSE)
  write.csv(table_s4c$data, file.path(table_dir, "Table_S4c_Phenotype_Proportions.csv"), row.names = FALSE)
  write.csv(table_s5$data,  file.path(table_dir, "Table_S5_Contractile_Rescue.csv"), row.names = FALSE)
  write.csv(table_s6$data,  file.path(table_dir, "Table_S6_Synthetic_Rescue.csv"), row.names = FALSE)
  write.csv(table_s7$data,  file.path(table_dir, "Table_S7_Piezo1_Correlations.csv"), row.names = FALSE)

  # Save methods text
  methods_text <- paste(
    "SUPPLEMENTAL TABLES METHODS",
    "",
    "Table S1: Cell Type Markers",
    table_s1$methods,
    "",
    "Table S2a: Phenotype Classification Genes",
    table_s2a$methods,
    "",
    "Table S2b: Gene Scoring Panels",
    table_s2b$methods,
    "",
    "Table S3: Pathway Gene Sets",
    table_s3$methods,
    "",
    "Table S4: Gene Expression by Group",
    table_s4$methods,
    "",
    "Table S5: Contractile Gene Rescue",
    table_s5$methods,
    "",
    "Table S6: Synthetic Gene Rescue",
    table_s6$methods,
    "",
    "Table S7: Piezo1 Gene Correlations",
    table_s7$methods,
    sep = "\n"
  )

  writeLines(methods_text, file.path(table_dir, "Supplemental_Tables_methods.txt"))

  cat("  Saved: Supplemental_Tables.xlsx\n")
  cat("  Saved: Individual CSV files\n")
  cat("  Saved: Supplemental_Tables_methods.txt\n")
  cat("  Done!\n")

  return(list(
    success = TRUE,
    output_dir = table_dir,
    tables = list(
      S1 = table_s1,
      S2a = table_s2a,
      S2b = table_s2b,
      S3 = table_s3,
      S4 = table_s4,
      S5 = table_s5,
      S6 = table_s6,
      S7 = table_s7
    )
  ))
}

# ==============================================================================
# TABLE S1: Cell Type Markers
# ==============================================================================

create_table_s1_markers <- function(expr, smc_cells, metadata) {
  # Define marker genes used for cluster-based cell type classification
  # These match the marker sets used for the cluster-based cell-type annotation (00_preprocessing/11_)
  celltype_markers <- list(
    "SMC" = list(
      cell_type_col = "SMC_Contractile",
      genes = c("Myh11", "Acta2", "Tagln", "Cnn1", "Mylk")
    ),
    "ICC" = list(
      cell_type_col = "ICC",
      genes = c("Kit", "Ano1", "Slc8a1")
    ),
    "Fibroblast (PDGFRa+)" = list(
      cell_type_col = "Fibroblast",
      genes = c("Pdgfra", "Cd34", "Col1a1")
    ),
    "Macrophage" = list(
      cell_type_col = "Immune_Macrophage",
      genes = c("Cd68", "Adgre1", "Cx3cr1")
    ),
    "T Cell" = list(
      cell_type_col = "Immune_TCell",
      genes = c("Cd3e", "Cd3d", "Cd4", "Cd8a")
    ),
    "B Cell" = list(
      cell_type_col = "Immune_BCell",
      genes = c("Cd19", "Ms4a1", "Cd79a")
    ),
    "Enteric Glia" = list(
      cell_type_col = "Glial",
      genes = c("Sox10", "Plp1", "S100b")
    ),
    "Enteric Neuron" = list(
      cell_type_col = "Neuronal",
      genes = c("Tubb3", "Snap25", "Elavl4")
    ),
    "Endothelial" = list(
      cell_type_col = "Endothelial",
      genes = c("Pecam1", "Cdh5", "Vwf")
    ),
    "Endothelial Lymphatic" = list(
      cell_type_col = "Endothelial_Lymphatic",
      genes = c("Prox1", "Lyve1", "Flt4")
    ),
    "Pericyte" = list(
      cell_type_col = "Pericyte",
      genes = c("Pdgfrb", "Rgs5", "Cspg4")
    ),
    "Mesothelial" = list(
      cell_type_col = "Mesothelial",
      genes = c("Msln", "Upk3b", "Wt1")
    ),
    "Epithelial" = list(
      cell_type_col = "Epithelial",
      genes = c("Epcam", "Krt8", "Krt18")
    ),
    "Goblet" = list(
      cell_type_col = "Goblet",
      genes = c("Muc2", "Tff3", "Clca3")
    ),
    "Erythrocyte" = list(
      cell_type_col = "Erythrocyte",
      genes = c("Hba-a1", "Hbb-bt")
    )
  )

  # Calculate expression statistics for each marker
  results <- data.frame()
  all_cells <- colnames(expr)

  for (ct_name in names(celltype_markers)) {
    ct_info <- celltype_markers[[ct_name]]
    ct_col_val <- ct_info$cell_type_col
    genes <- ct_info$genes

    # Get cells of this type and all other cells
    type_cells <- rownames(metadata)[metadata$cell_type == ct_col_val]
    type_cells <- type_cells[type_cells %in% all_cells]
    other_cells <- setdiff(all_cells, type_cells)

    n_type <- length(type_cells)

    for (gene in genes) {
      if (!gene %in% rownames(expr)) next

      if (n_type > 0) {
        type_expr <- as.numeric(expr[gene, type_cells])
        mean_type <- mean(type_expr, na.rm = TRUE)
        pct_type <- 100 * mean(type_expr > 0)
      } else {
        mean_type <- NA
        pct_type <- NA
      }

      other_expr <- as.numeric(expr[gene, other_cells])
      mean_other <- mean(other_expr, na.rm = TRUE)
      pct_other <- 100 * mean(other_expr > 0)

      # Log2FC (add pseudocount)
      if (!is.na(mean_type)) {
        log2fc <- log2((mean_type + 0.01) / (mean_other + 0.01))
      } else {
        log2fc <- NA
      }

      results <- rbind(results, data.frame(
        `Cell Type` = ct_name,
        Gene = gene,
        `N Cells` = n_type,
        `Mean Expression` = round(mean_type, 4),
        `Pct Expressing` = round(pct_type, 1),
        `Mean Other` = round(mean_other, 4),
        `Pct Other` = round(pct_other, 1),
        `Log2FC vs Other` = round(log2fc, 2),
        check.names = FALSE
      ))
    }
  }

  legend <- paste(
    "Marker genes used for cluster-based cell type classification.",
    "Each cell type was identified by calculating mean expression of its marker set",
    "per Seurat cluster (resolution 0.5) and assigning clusters to the cell type",
    "with the highest marker score (threshold > 0.2).",
    "Mean Expression = normalized log counts in the identified cell type.",
    "Pct Expressing = percentage of cells in that type with expression > 0.",
    "Mean Other / Pct Other = same metrics for all other cells.",
    "Log2FC vs Other = log2 fold change of mean expression in the cell type vs all others."
  )

  methods <- paste(
    "Cell type classification used a cluster-based approach.",
    "Seurat clustering (resolution 0.5, 30 PCs) identified transcriptional clusters.",
    "Each cluster was assigned to the cell type with the highest mean marker gene score.",
    "SMC clusters were further validated by high expression of Myh11, Acta2, Tagln, Cnn1, Mylk",
    "and absence of fibroblast (Pdgfra), ICC (Kit, Ano1), and immune (Ptprc, Cd68) markers.",
    "Epithelial-contaminated cells expressing Epcam, Krt8, Krt19, Reg3b, Reg3g, or Lypd8",
    "were removed prior to SMC analysis."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S2a: Phenotype Classification Genes (7 contractile + 5 synthetic)
# ==============================================================================

create_table_s2a_classification <- function(gene_sets) {
  # 12 genes used for contractile:synthetic ratio classification
  gene_info <- list(
    list(gene="Myh11",  phenotype="Contractile", fn="Smooth muscle myosin heavy chain (major contractile protein)",
         refs="Owens et al. Physiol Rev 2004"),
    list(gene="Acta2",  phenotype="Contractile", fn="Smooth muscle alpha-actin (thin filament)",
         refs="Owens et al. Physiol Rev 2004"),
    list(gene="Tagln",  phenotype="Contractile", fn="Transgelin / SM22alpha (actin-binding, cytoskeletal)",
         refs="Owens et al. Physiol Rev 2004"),
    list(gene="Cnn1",   phenotype="Contractile", fn="Calponin 1 (thin filament regulatory)",
         refs="Owens et al. Physiol Rev 2004"),
    list(gene="Mylk",   phenotype="Contractile", fn="Myosin light chain kinase (calcium-dependent activation)",
         refs="Owens et al. Physiol Rev 2004"),
    list(gene="Myocd",  phenotype="Contractile", fn="Myocardin (master transcriptional regulator)",
         refs="Wang et al. PNAS 2003"),
    list(gene="Smtn",   phenotype="Contractile", fn="Smoothelin (late differentiation marker)",
         refs="van Eys et al. Trends Cardiovasc Med 2007"),
    list(gene="Col1a1", phenotype="Synthetic", fn="Collagen type I alpha 1 (ECM production)",
         refs="Yap 2021; Schneider 2023; Wolfson 2025"),
    list(gene="Fn1",    phenotype="Synthetic", fn="Fibronectin 1 (ECM glycoprotein, migration)",
         refs="Yap 2021; Wolfson 2025"),
    list(gene="Spp1",   phenotype="Synthetic", fn="Osteopontin (ECM, migration, calcification)",
         refs="Yap 2021"),
    list(gene="Klf4",   phenotype="Synthetic", fn="Kruppel-like factor 4 (dedifferentiation TF, represses Myocd)",
         refs="Yap 2021; Wolfson 2025"),
    list(gene="Vim",    phenotype="Synthetic", fn="Vimentin (mesenchymal intermediate filament)",
         refs="Schneider 2023; Wolfson 2025")
  )

  results <- do.call(rbind, lapply(gene_info, function(g) {
    data.frame(
      Phenotype = g$phenotype,
      Gene = g$gene,
      Function = g$fn,
      References = g$refs,
      stringsAsFactors = FALSE
    )
  }))

  legend <- paste(
    "Genes used for SMC phenotype classification via the contractile:synthetic ratio.",
    "7 contractile markers and 5 synthetic markers define the two phenotypic poles.",
    "Module scores are calculated as mean normalized expression across each panel;",
    "cells with ratio > 1.5 are classified as contractile (Figures 3E, 3K).",
    "All synthetic panel genes are supported by published SMC phenotype literature",
    "(Yap 2021, Schneider 2023, Wolfson 2025)."
  )

  methods <- paste(
    "SMC phenotype classification used 7 contractile markers (Myh11, Acta2, Tagln,",
    "Cnn1, Mylk, Myocd, Smtn) and 5 synthetic markers (Col1a1, Fn1, Spp1, Klf4, Vim).",
    "Module scores were calculated as the mean normalized expression across each panel.",
    "The contractile:synthetic ratio was used to classify individual cells (ratio > 1.5)."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S2b: Gene Scoring Panels (contractile capacity + synthetic phenotype)
# ==============================================================================

create_table_s2b_scoring <- function(gene_sets) {
  gene_info <- list(
    # Contractile Capacity Score (7 genes -> Fig 3E)
    list(panel="Contractile Capacity Score", gene="Myh11",
         fn="Smooth muscle myosin heavy chain (major contractile protein)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Acta2",
         fn="Smooth muscle alpha-actin (thin filament)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Tagln",
         fn="Transgelin / SM22alpha (actin-binding, cytoskeletal)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Cnn1",
         fn="Calponin 1 (thin filament regulatory)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Mylk",
         fn="Myosin light chain kinase (calcium-dependent activation)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Myocd",
         fn="Myocardin (master transcriptional regulator)",
         figure="Fig 3H", refs=""),
    list(panel="Contractile Capacity Score", gene="Smtn",
         fn="Smoothelin (late differentiation marker)",
         figure="Fig 3H", refs=""),
    # Synthetic Phenotype Score (5 genes -> Figs 3G, 4G)
    list(panel="Synthetic Phenotype Score", gene="Klf4",
         fn="Kruppel-like factor 4 (dedifferentiation TF, represses Myocd)",
         figure="Figs 3G, 4G",
         refs="Shankman et al. Nat Med 2015; Alencar et al. ATVB 2021"),
    list(panel="Synthetic Phenotype Score", gene="Vim",
         fn="Vimentin (mesenchymal intermediate filament)",
         figure="Figs 3G, 4G",
         refs="Liu et al. Cell Commun Signal 2022; Rensen et al. Cardiovasc Res 2007"),
    list(panel="Synthetic Phenotype Score", gene="Mmp2",
         fn="Matrix metalloproteinase 2 (ECM degradation)",
         figure="Figs 3G, 4G",
         refs="Yap et al. Circulation 2020; Liu et al. Cell Commun Signal 2022"),
    list(panel="Synthetic Phenotype Score", gene="Mmp9",
         fn="Matrix metalloproteinase 9 (ECM degradation)",
         figure="Figs 3G, 4G",
         refs="Liu et al. Cell Commun Signal 2022; Nagase et al. Matrix Biol 2006"),
    list(panel="Synthetic Phenotype Score", gene="Timp1",
         fn="Tissue inhibitor of metalloproteinases 1",
         figure="Figs 3G, 4G",
         refs="Baker et al. Circ Res 1996; Brew et al. BBA 2010"),
    # NFAT Activity Score (5 genes -> Fig 4D) — downstream NFAT effectors
    list(panel="NFAT Activity Score", gene="Wnt5a",
         fn="Wnt5a (non-canonical Wnt/Ca2+; upstream NFAT amplifier)",
         figure="Fig 4E",
         refs="Dejmek et al. Mol Cell Biol 2006; Dodge et al. JBC 2012"),
    list(panel="NFAT Activity Score", gene="Mmp2",
         fn="Matrix metalloproteinase 2 (direct NFATc3 transcriptional target)",
         figure="Fig 4E",
         refs="Zhang et al. 2025 (ChIP + luciferase)"),
    list(panel="NFAT Activity Score", gene="Mmp9",
         fn="Matrix metalloproteinase 9 (direct NFATc3 transcriptional target)",
         figure="Fig 4E",
         refs="Zhang et al. 2025 (ChIP + luciferase)"),
    list(panel="NFAT Activity Score", gene="Vcam1",
         fn="Vascular cell adhesion molecule 1 (NFAT-family responsive)",
         figure="Fig 4E",
         refs=""),
    list(panel="NFAT Activity Score", gene="Fas",
         fn="Apoptosis receptor (NFAT-associated, FasL co-regulation)",
         figure="Fig 4E",
         refs=""),
    list(panel="NFAT Activity Score", gene="Timp1",
         fn="Tissue inhibitor of metalloproteinases 1 (co-regulated with MMP2/9)",
         figure="Fig 4E",
         refs=""),
    # ECM Score (5 genes -> Fig 3M) — structural ECM panel
    list(panel="ECM Score", gene="Bgn",
         fn="Biglycan (small leucine-rich proteoglycan)",
         figure="Fig S3K", refs=""),
    list(panel="ECM Score", gene="Col6a1",
         fn="Collagen type VI alpha 1 (basement-membrane-associated)",
         figure="Fig S3K", refs=""),
    list(panel="ECM Score", gene="Dcn",
         fn="Decorin (small leucine-rich proteoglycan, collagen-binding)",
         figure="Fig S3K", refs=""),
    list(panel="ECM Score", gene="Sparc",
         fn="Secreted protein acidic and rich in cysteine (matricellular)",
         figure="Fig S3K", refs=""),
    list(panel="ECM Score", gene="Timp2",
         fn="Tissue inhibitor of metalloproteinases 2 (the only protease-axis gene without NFAT-Activity overlap)",
         figure="Fig S3K", refs=""),
    # Cytoskeletal Remodeling Score (4 genes -> Fig 3N).
    list(panel="Cytoskeletal Remodeling Score", gene="Ezr",
         fn="Ezrin (actin-membrane crosslinker, mesenchymal)",
         figure="Fig S4H", refs=""),
    list(panel="Cytoskeletal Remodeling Score", gene="Syne2",
         fn="Nesprin-2 (nuclear envelope linker, cytoskeletal coupling)",
         figure="Fig S4H", refs=""),
    list(panel="Cytoskeletal Remodeling Score", gene="Lmna",
         fn="Lamin A/C (nuclear lamina, mechanotransduction)",
         figure="Fig S4H", refs=""),
    list(panel="Cytoskeletal Remodeling Score", gene="Vim",
         fn="Vimentin (intermediate filament; also in Synthetic Phenotype Score)",
         figure="Fig S4H", refs="")
  )

  results <- do.call(rbind, lapply(gene_info, function(g) {
    data.frame(
      Score.Panel = g$panel,
      Gene = g$gene,
      Function = g$fn,
      Figure = g$figure,
      References = g$refs,
      stringsAsFactors = FALSE
    )
  }))

  legend <- paste(
    "Gene scoring panels used for quantitative analyses across the manuscript.",
    "Each score is the mean normalized expression of the listed genes within each",
    "cell, computed across the bona fide SMC population.",
    "(1) Contractile Capacity Score (7 genes; Fig 3E): canonical SMC contractile",
    "markers — same set as Table S2a.",
    "(2) Synthetic Phenotype Score (5 genes; Figs 3D, 3F, 4A, 4C): dedifferentiation and",
    "matrix remodeling, used for correlation analyses with Piezo1 and NFAT.",
    "(3) NFAT Activity Score (5 genes; Fig 4D): downstream NFAT-responsive",
    "effectors. Reads out NFAT being active rather than just present (upstream",
    "NFAT machinery is in NFAT Pathway, Table S3).",
    "(4) ECM Score (5 genes; Fig 3M): structural ECM panel; it shares no gene with",
    "NFAT Activity (Mmp2/Mmp9/Timp1 are counted only there).",
    "(5) Cytoskeletal Remodeling Score (4 genes; Fig 3N): cytoskeletal step of the",
    "contractile-to-synthetic dedifferentiation cascade. Vim is shared with the",
    "Synthetic Phenotype Score by design (mesenchymal driver of both phenotypes).",
    "All scores use Bonferroni-adjusted (m=4) Mann-Whitney U for pairwise group",
    "comparisons across {Young CTRL, Aged CTRL, Young KO, Aged KO}, excluding the",
    "two diagonal cross-group pairs that conflate aging + genotype."
  )

  methods <- paste(
    "Module scores are computed as the mean of log-normalized expression across",
    "the gene set within each bona fide SMC, then group-averaged across all cells",
    "in each experimental condition. Group comparisons use two-sided Mann-Whitney",
    "U tests with Bonferroni correction (m = 4 valid pairwise comparisons:",
    "Aging WT, Aging KO, KO@young, KO rescue@aged; the two diagonal cross-group",
    "comparisons that conflate age and genotype are excluded). The Synthetic",
    "Phenotype Score (5 genes) was selected from 10 candidate gene sets by data-",
    "driven validation requiring significant Spearman correlation with Piezo1 in",
    "aged SMCs (rho > 0.10, p < 0.01), significant correlation with NFAT pathway",
    "score in aged SMCs (rho > 0.25, p < 0.001), age-dependent score increase",
    "(p < 0.001), and Piezo1 KO rescue (p < 0.001)."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S3: Pathway Gene Sets (NFAT, ECM)
# ==============================================================================

create_table_s3_pathway_genes <- function(expr, smc_cells, metadata, gene_sets) {
  # NFAT Pathway genes
  nfat_genes <- gene_sets$nfat_pathway
  if (is.null(nfat_genes) || length(nfat_genes) == 0) {
    nfat_genes <- c("Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4", "Nfat5",
                    "Ppp3ca", "Ppp3cb", "Ppp3cc", "Rcan1", "Rcan2", "Rcan3",
                    "Calm1", "Calm3")
  }

  # ECM Score — UNIFIED ECM SCORE (5 genes; no overlap with NFAT scores)
  # Replaces ecm_remodeling_validated + ecm_organization_validated.
  ecm_structural <- gene_sets$ecm_structural
  if (is.null(ecm_structural)) ecm_structural <- c("Bgn", "Col6a1", "Dcn", "Sparc", "Timp2")

  # NFAT Activity genes (downstream effectors; was "NFAT-Responsive Remodeling")
  nfatc3_targets <- gene_sets$nfatc3_remodeling
  if (is.null(nfatc3_targets)) nfatc3_targets <- c("Wnt5a", "Mmp2", "Mmp9", "Vcam1", "Fas", "Timp1")

  # Cytoskeletal Remodeling, parallel to ECM Score and NFAT Activity
  cyto_remodeling <- gene_sets$cyto_remodeling
  if (is.null(cyto_remodeling)) cyto_remodeling <- c("Ezr", "Syne2", "Lmna", "Vim")

  # Build combined table
  all_pathway_genes <- data.frame(
    Pathway = c(
      rep("NFAT Pathway", length(nfat_genes)),
      rep("NFAT Activity", length(nfatc3_targets)),
      rep("ECM Score", length(ecm_structural)),
      rep("Cytoskeletal Remodeling", length(cyto_remodeling))
    ),
    Gene = c(nfat_genes, nfatc3_targets, ecm_structural, cyto_remodeling),
    stringsAsFactors = FALSE
  )

  # Remove duplicates within the same pathway
  all_pathway_genes <- all_pathway_genes[!duplicated(all_pathway_genes[, c("Pathway", "Gene")]), ]

  # Track which genes appear in multiple pathways for cross-referencing
  gene_pathway_counts <- table(all_pathway_genes$Gene)
  multi_pathway_genes <- names(gene_pathway_counts[gene_pathway_counts > 1])

  # Build cross-reference lookup: for each gene, list OTHER pathways it appears in
  cross_ref <- list()
  for (g in multi_pathway_genes) {
    pathways_for_gene <- all_pathway_genes$Pathway[all_pathway_genes$Gene == g]
    cross_ref[[g]] <- pathways_for_gene
  }

  # Get cells by group
  young_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Young_Ctrl"]]
  aged_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Aged_Ctrl"]]
  young_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Young_KO"]]
  aged_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Aged_KO"]]

  results <- data.frame()

  for (i in 1:nrow(all_pathway_genes)) {
    gene <- all_pathway_genes$Gene[i]
    pathway <- all_pathway_genes$Pathway[i]

    if (!gene %in% rownames(expr)) {
      results <- rbind(results, data.frame(
        Pathway = pathway,
        Gene = gene,
        `Mean Young Ctrl` = NA, `Mean Aged Ctrl` = NA,
        `Mean Young KO` = NA, `Mean Aged KO` = NA,
        `Aging Effect` = NA, `Rescue Effect` = NA,
        `Aging P-value` = NA, `Rescue P-value` = NA,
        `Also In Pathway` = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      ))
      next
    }

    yc_expr <- as.numeric(expr[gene, young_ctrl_cells])
    ac_expr <- as.numeric(expr[gene, aged_ctrl_cells])
    yko_expr <- as.numeric(expr[gene, young_ko_cells])
    ako_expr <- as.numeric(expr[gene, aged_ko_cells])

    mean_yc <- mean(yc_expr, na.rm = TRUE)
    mean_ac <- mean(ac_expr, na.rm = TRUE)
    mean_yko <- mean(yko_expr, na.rm = TRUE)
    mean_ako <- mean(ako_expr, na.rm = TRUE)

    test_aging <- tryCatch(wilcox.test(ac_expr, yc_expr), error = function(e) list(p.value = NA))
    test_rescue <- tryCatch(wilcox.test(ako_expr, ac_expr), error = function(e) list(p.value = NA))

    # Build cross-reference note for genes in multiple pathways
    also_in <- ""
    if (gene %in% multi_pathway_genes) {
      other_pathways <- setdiff(cross_ref[[gene]], pathway)
      if (length(other_pathways) > 0) {
        also_in <- paste("Also in:", paste(other_pathways, collapse = ", "))
      }
    }

    results <- rbind(results, data.frame(
      Pathway = pathway,
      Gene = gene,
      `Mean Young Ctrl` = round(mean_yc, 4),
      `Mean Aged Ctrl` = round(mean_ac, 4),
      `Mean Young KO` = round(mean_yko, 4),
      `Mean Aged KO` = round(mean_ako, 4),
      `Aging Effect` = round(mean_ac - mean_yc, 4),
      `Rescue Effect` = round(mean_ako - mean_ac, 4),
      `Aging P-value` = ifelse(is.na(test_aging$p.value), NA, sprintf("%.2e", test_aging$p.value)),
      `Rescue P-value` = ifelse(is.na(test_rescue$p.value), NA, sprintf("%.2e", test_rescue$p.value)),
      `Also In Pathway` = also_in,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  # Sort by pathway then gene
  results <- results %>% arrange(Pathway, Gene)

  legend <- paste(
    "Gene sets for pathway score calculations.",
    "NFAT Pathway (13 genes): Nfatc1-4, Nfat5, Ppp3ca/cb/cc, Rcan1-3, Calm1, Calm3 —",
    "upstream signaling machinery used for cell-level correlation analyses (Figures 4F, 4G);",
    "transcriptionally near-invariant across groups (NFAT activity is regulated post-",
    "translationally via calcineurin-dependent nuclear translocation).",
    "NFAT Activity (6 genes): Wnt5a, Mmp2, Mmp9, Vcam1, Fas, Timp1 — downstream effectors",
    "with >50% aging increase and >30% KO decrease (both p<0.05). Reads out NFAT being",
    "active (post-translational regulation), not just expressed.",
    "ECM Score (5 genes): Bgn, Col6a1, Dcn, Sparc, Timp2 — unified ECM score with",
    "no gene overlap with NFAT Activity.",
    "Cytoskeletal Remodeling (4 genes): Ezr, Syne2, Lmna, Vim — cytoskeletal step of",
    "the contractile-to-synthetic dedifferentiation cascade (Fig 3N). Vim overlaps",
    "with the Synthetic Phenotype Score by design (mesenchymal driver of both phenotypes)."
  )

  methods <- paste(
    "Pathway gene sets were curated from published literature and validated",
    "against our scRNA-seq data. NFAT pathway includes transcription factors",
    "(Nfatc1-4, Nfat5), calcineurin phosphatase subunits (Ppp3ca/cb/cc),",
    "RCAN regulators (Rcan1-3), and calmodulin (Calm1, Calm3).",
    "NFAT Activity genes (downstream effectors) were validated to show the expected pattern",
    "(UP with aging, DOWN with Piezo1 KO) in our data.",
    "ECM gene sets were validated with >20% change criteria in both directions.",
    "Expression values are normalized log counts.",
    "P-values from Wilcoxon rank-sum tests."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S4: Gene Expression by Group
# ==============================================================================

create_table_s4_expression <- function(expr, smc_cells, metadata, gene_sets) {
  # Combine all relevant genes
  all_genes <- unique(c(
    gene_sets$contractile_core,
    gene_sets$contractile_score,
    gene_sets$synthetic_core,
    gene_sets$synthetic_score,
    "Piezo1"
  ))

  available_genes <- all_genes[all_genes %in% rownames(expr)]

  # Get cells by group
  groups <- c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")

  results <- data.frame()

  for (gene in available_genes) {
    row_data <- data.frame(Gene = gene)

    for (grp in groups) {
      cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == grp]]
      gene_expr <- as.numeric(expr[gene, cells])

      grp_label <- gsub("_", " ", grp)
      row_data[[paste0(grp_label, " Mean")]] <- round(mean(gene_expr, na.rm = TRUE), 4)
      row_data[[paste0(grp_label, " SD")]] <- round(sd(gene_expr, na.rm = TRUE), 4)
      row_data[[paste0(grp_label, " Pct Expressing")]] <- round(100 * mean(gene_expr > 0), 1)
    }

    # Determine gene category with clear labels (T3: Classification vs Extended)
    if (gene %in% gene_sets$contractile_core) {
      row_data$Category <- "Contractile (Classification Panel)"
    } else if (gene %in% gene_sets$contractile_score) {
      row_data$Category <- "Contractile (Extended Panel)"
    } else if (gene %in% gene_sets$synthetic_core) {
      row_data$Category <- "Synthetic (Classification Panel)"
    } else if (gene %in% gene_sets$synthetic_score) {
      row_data$Category <- "Synthetic (Extended Panel)"
    } else {
      row_data$Category <- "Other"
    }

    results <- rbind(results, row_data)
  }

  # Reorder columns
  results <- results[, c("Gene", "Category",
                         "Young Ctrl Mean", "Young Ctrl SD", "Young Ctrl Pct Expressing",
                         "Aged Ctrl Mean", "Aged Ctrl SD", "Aged Ctrl Pct Expressing",
                         "Young KO Mean", "Young KO SD", "Young KO Pct Expressing",
                         "Aged KO Mean", "Aged KO SD", "Aged KO Pct Expressing")]

  # Compute dynamic counts from data
  n_total <- length(smc_cells)
  group_counts <- table(metadata[smc_cells[smc_cells %in% rownames(metadata)], "group"])

  legend <- paste(
    "Mean normalized expression levels for contractile, synthetic, and pathway",
    sprintf("genes across all four experimental groups in bona fide SMCs (n=%s).", format(n_total, big.mark=",")),
    sprintf("Cell counts: Young CTRL n=%s, Aged CTRL n=%s, Young KO n=%s, Aged KO n=%s.",
            format(group_counts["Young_Ctrl"], big.mark=","),
            format(group_counts["Aged_Ctrl"], big.mark=","),
            format(group_counts["Young_KO"], big.mark=","),
            format(group_counts["Aged_KO"], big.mark=",")),
    "Category: 'Classification Panel' = genes used in the contractile:synthetic ratio",
    "for phenotype classification (Figs 3E/3K; see Table S2a);",
    "'Extended Panel' = additional markers tracked for broader analysis.",
    "Pct Expressing = percentage of cells with expression > 0."
  )

  methods <- paste(
    "Expression values are log-normalized counts from the Seurat pipeline.",
    "Statistics were calculated for SMCs only (bona fide cluster-based selection).",
    "Groups: Young_Ctrl, Aged_Ctrl, Young_KO, Aged_KO."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S5: Contractile Gene Rescue Statistics
# ==============================================================================

create_table_s5_contractile <- function(expr, smc_cells, metadata, gene_sets) {
  # Get all contractile-related genes
  contractile_genes <- unique(c(
    gene_sets$contractile_core,
    gene_sets$contractile_score,
    gene_sets$muscle_contraction
  ))

  results <- calculate_rescue_stats(expr, smc_cells, metadata, contractile_genes, "Contractile")

  legend <- paste(
    "Mean normalized expression levels of contractile marker genes across",
    "experimental groups in high-confidence SMCs.",
    "Aging Effect = Aged Ctrl - Young Ctrl (negative = decline).",
    "Rescue Effect = Aged KO - Aged Ctrl (positive = restoration).",
    "Pct Rescued = (Rescue Effect / |Aging Effect|) * 100.",
    "Passes All Criteria = Shows expected aging decline AND significant KO rescue.",
    "P-values from Wilcoxon rank-sum tests."
  )

  methods <- paste(
    "Contractile gene rescue was assessed by comparing expression between groups.",
    "Genes showing significant rescue meet three criteria: (1) decline with aging",
    "(Aged_Ctrl < Young_Ctrl), (2) rescued by Piezo1 KO (Aged_KO > Aged_Ctrl),",
    "(3) rescue is statistically significant (Wilcoxon p < 0.05).",
    "The contractile score (Figure 3E) uses seven canonical SMC contractile markers",
    "(Myh11, Acta2, Tagln, Cnn1, Mylk, Myocd, Smtn; see Table S2, Used_For column)."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S6: Synthetic Gene Rescue Statistics
# ==============================================================================

create_table_s6_synthetic <- function(expr, smc_cells, metadata, gene_sets) {
  # Get all synthetic-related genes
  synthetic_genes <- unique(c(
    gene_sets$synthetic_core,
    gene_sets$synthetic_score,
    gene_sets$ecm_remodeling
  ))

  results <- calculate_rescue_stats(expr, smc_cells, metadata, synthetic_genes, "Synthetic")

  legend <- paste(
    "Mean normalized expression levels of synthetic/ECM marker genes across",
    "experimental groups in high-confidence SMCs.",
    "Aging Effect = Aged Ctrl - Young Ctrl (positive = increase with aging).",
    "Rescue Effect = Aged KO - Aged Ctrl (negative = rescue/reduction).",
    "Pct Rescued = (|Rescue Effect| / Aging Effect) * 100.",
    "P-values from Wilcoxon rank-sum tests."
  )

  methods <- paste(
    "Synthetic gene rescue was assessed by comparing expression between groups.",
    "For synthetic genes, rescue means reduction in expression (Aged_KO < Aged_Ctrl)",
    "following age-related increase (Aged_Ctrl > Young_Ctrl).",
    "Genes meeting all criteria are marked as passing."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# TABLE S7: Piezo1 Gene Correlations
# ==============================================================================

create_table_s7_piezo1_correlations <- function(expr, smc_cells, metadata, gene_sets) {
  # Check if Piezo1 is available
  if (!"Piezo1" %in% rownames(expr)) {
    return(list(
      data = data.frame(Note = "Piezo1 not found in expression matrix"),
      legend = "Piezo1 gene correlations could not be computed.",
      methods = "Piezo1 was not detected in the expression matrix."
    ))
  }

  # Get Piezo1 expression in SMC cells
  piezo1_expr <- as.numeric(expr["Piezo1", smc_cells])

  # Correlate with all contractile, synthetic, and pathway genes
  target_genes <- unique(c(
    gene_sets$contractile_core,
    gene_sets$contractile_extended,
    gene_sets$synthetic_core,
    gene_sets$synthetic_extended,
    gene_sets$nfat_pathway,
    gene_sets$nfatc3_remodeling,
    gene_sets$ecm_structural,
    "Tgfb1", "Klf4"  # Added for manuscript correlation claims
  ))

  available_targets <- target_genes[target_genes %in% rownames(expr)]

  results <- data.frame()

  for (gene in available_targets) {
    gene_expr <- as.numeric(expr[gene, smc_cells])

    # Spearman correlation
    cor_result <- tryCatch({
      cor.test(piezo1_expr, gene_expr, method = "spearman")
    }, error = function(e) list(estimate = NA, p.value = NA))

    # Determine gene category — priority: contractile > synthetic_core > pathway/target > synthetic_extended
    # NFAT-responsive remodeling and ECM sets take priority over synthetic_extended
    if (gene %in% gene_sets$contractile_core) {
      category <- "Contractile Core"
    } else if (gene %in% gene_sets$contractile_extended) {
      category <- "Contractile Extended"
    } else if (gene %in% gene_sets$synthetic_core) {
      category <- "Synthetic Core"
    } else if (gene %in% gene_sets$nfat_pathway) {
      category <- "NFAT Pathway"
    } else if (gene %in% gene_sets$nfatc3_remodeling) {
      category <- "NFAT Activity"
    } else if (gene %in% gene_sets$ecm_structural) {
      category <- "ECM Score"
    } else if (gene %in% gene_sets$synthetic_extended) {
      category <- "Synthetic Extended"
    } else {
      category <- "Other"
    }

    results <- rbind(results, data.frame(
      Gene = gene,
      Category = category,
      `Spearman Rho` = round(as.numeric(cor_result$estimate), 4),
      `P-value` = ifelse(is.na(cor_result$p.value), NA, sprintf("%.2e", cor_result$p.value)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  # Sort by signed rho (most positive first, most negative last)
  results <- results %>%
    arrange(desc(`Spearman Rho`))

  n_cells <- length(smc_cells)
  legend <- paste(
    "Spearman correlations between Piezo1 expression and contractile, synthetic,",
    "and pathway genes in high-confidence SMCs (n =", format(n_cells, big.mark = ","), "cells).",
    "Positive rho indicates co-expression with Piezo1;",
    "negative rho indicates inverse relationship.",
    "Sorted by signed rho (most positive first, most negative last)."
  )

  methods <- paste(
    "Spearman rank correlations were calculated between Piezo1 expression and",
    "each target gene across all SMC cells. Target genes include contractile markers,",
    "synthetic markers, NFAT pathway components, and ECM genes.",
    "Correlations were computed on log-normalized expression values."
  )

  return(list(data = results, legend = legend, methods = methods))
}

# ==============================================================================
# Helper function: Calculate rescue statistics
# ==============================================================================

calculate_rescue_stats <- function(expr, smc_cells, metadata, genes, gene_type) {
  # Get cells by group
  young_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Young_Ctrl"]]
  aged_ctrl_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Aged_Ctrl"]]
  aged_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Aged_KO"]]
  young_ko_cells <- smc_cells[smc_cells %in% rownames(metadata)[metadata$group == "Young_KO"]]

  results <- data.frame()

  for (gene in genes) {
    if (!gene %in% rownames(expr)) next

    yc_expr <- as.numeric(expr[gene, young_ctrl_cells])
    ac_expr <- as.numeric(expr[gene, aged_ctrl_cells])
    ako_expr <- as.numeric(expr[gene, aged_ko_cells])
    yko_expr <- as.numeric(expr[gene, young_ko_cells])

    mean_yc <- mean(yc_expr, na.rm = TRUE)
    mean_ac <- mean(ac_expr, na.rm = TRUE)
    mean_ako <- mean(ako_expr, na.rm = TRUE)
    mean_yko <- mean(yko_expr, na.rm = TRUE)

    aging_effect <- mean_ac - mean_yc
    rescue_effect <- mean_ako - mean_ac

    # Statistical tests
    test_aging <- tryCatch({
      wilcox.test(ac_expr, yc_expr)
    }, error = function(e) list(p.value = NA))

    test_rescue <- tryCatch({
      wilcox.test(ako_expr, ac_expr)
    }, error = function(e) list(p.value = NA))

    # Determine if passes criteria
    if (gene_type == "Contractile") {
      shows_decline <- aging_effect < 0
      shows_rescue <- rescue_effect > 0
    } else {
      # Synthetic: increase with aging, decrease with KO
      shows_decline <- aging_effect > 0  # "decline" = increase for synthetic
      shows_rescue <- rescue_effect < 0   # "rescue" = decrease for synthetic
    }

    rescue_significant <- !is.na(test_rescue$p.value) && test_rescue$p.value < 0.05
    passes <- shows_decline & shows_rescue & rescue_significant

    # Calculate percent rescued
    if (aging_effect != 0) {
      pct_rescued <- 100 * abs(rescue_effect) / abs(aging_effect)
    } else {
      pct_rescued <- NA
    }

    results <- rbind(results, data.frame(
      Gene = gene,
      `Mean Young Ctrl` = round(mean_yc, 4),
      `Mean Aged Ctrl` = round(mean_ac, 4),
      `Mean Young KO` = round(mean_yko, 4),
      `Mean Aged KO` = round(mean_ako, 4),
      `Aging Effect` = round(aging_effect, 4),
      `Rescue Effect` = round(rescue_effect, 4),
      `Pct Rescued` = round(pct_rescued, 1),
      `Aging P-value` = ifelse(is.na(test_aging$p.value), NA, sprintf("%.2e", test_aging$p.value)),
      `Rescue P-value` = ifelse(is.na(test_rescue$p.value), NA, sprintf("%.2e", test_rescue$p.value)),
      `Passes All Criteria` = passes,
      check.names = FALSE
    ))
  }

  # Sort: passing genes first (by Pct Rescued desc), then non-passing (by rescue p-value)
  results <- results %>%
    mutate(.sort_key = ifelse(
      `Passes All Criteria`,
      -`Pct Rescued`,
      1e6 + as.numeric(`Rescue P-value`)
    )) %>%
    arrange(desc(`Passes All Criteria`), .sort_key) %>%
    select(-.sort_key)

  return(results)
}

# ==============================================================================
# TABLE S0: Sample Inclusion + QC
# ==============================================================================
# Per-sample cell counts, median UMI and bona fide SMC contribution (working table).
create_table_s0_sample_inclusion <- function(seurat_obj, context) {
  md <- seurat_obj@meta.data
  smc_cells <- context$selection$smc

  all_samples <- c("WT_Young_1","WT_Young_2","KO_Young_1","WT_Aged_1","WT_Aged_2","KO_Aged")
  group_map <- c(WT_Young_1="Young_Ctrl", WT_Young_2="Young_Ctrl", KO_Young_1="Young_KO",
                 WT_Aged_1="Aged_Ctrl",   WT_Aged_2="Aged_Ctrl",   KO_Aged="Aged_KO")

  # Cell counts in current (post-filter) seurat
  current_counts <- table(md$sample)
  umi_med_current <- sapply(all_samples, function(s) {
    m <- seurat_obj$nCount_RNA[md$sample == s]
    if (length(m) == 0) NA_real_ else median(m)
  })

  # Bona fide SMC counts per sample (intersection with smc_cells)
  bf_per_sample <- table(md[smc_cells, "sample"])
  bf_count <- sapply(all_samples, function(s) {
    v <- bf_per_sample[s]; if (is.na(v) || length(v) == 0) 0L else as.integer(v)
  })

  status <- rep("INCLUDED", length(all_samples))
  reason <- rep("Pass all QC criteria", length(all_samples))

  results <- data.frame(
    Sample = all_samples,
    Group  = unname(group_map[all_samples]),
    `Cells (post-filter)` = as.integer(ifelse(is.na(current_counts[all_samples]), 0, current_counts[all_samples])),
    `Median UMI (post-filter)` = round(umi_med_current),
    `Bona fide SMCs` = bf_count,
    Status = status,
    Reason = reason,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  legend <- paste("Per-sample cell count, median UMI and bona fide SMC contribution (working table; the manuscript Table S0 is in Tables_S0_S7.xlsx).")
  return(list(data = results, legend = legend))
}

# ==============================================================================
# TABLE S0b: Piezo1 KO Validation per Sample
# ==============================================================================
create_table_s0b_piezo1_validation <- function(expr, smc_cells, metadata) {
  if (!"Piezo1" %in% rownames(expr)) {
    return(list(data = data.frame(Note = "Piezo1 not in expression matrix"),
                legend = "Piezo1 KO validation not computed (gene missing)."))
  }
  p1 <- as.numeric(expr["Piezo1", smc_cells])
  md_bf <- metadata[smc_cells, ]
  md_bf$piezo1 <- p1

  # Reference: Aged_Ctrl mean
  ac_mean <- mean(md_bf$piezo1[md_bf$group == "Aged_Ctrl"], na.rm = TRUE)
  ac_p1 <- md_bf$piezo1[md_bf$group == "Aged_Ctrl"]

  samples <- sort(unique(md_bf$sample))
  rows <- lapply(samples, function(s) {
    e <- md_bf$piezo1[md_bf$sample == s]
    grp <- unique(md_bf$group[md_bf$sample == s])[1]
    data.frame(
      Sample = s,
      Group  = grp,
      `n cells` = length(e),
      `Mean Piezo1` = round(mean(e, na.rm = TRUE), 4),
      `Pct Piezo1+` = round(100 * mean(e > 0, na.rm = TRUE), 1),
      `Log2FC vs A_Ctrl mean` = round(log2((mean(e, na.rm = TRUE) + 1e-6) / (ac_mean + 1e-6)), 3),
      `Wilcoxon p (vs A_Ctrl)` = if (grp == "Aged_Ctrl") "n/a (reference)" else
        signif(suppressWarnings(wilcox.test(e, ac_p1)$p.value), 3),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, rows)

  legend <- paste(
    "Per-sample Piezo1 expression statistics within the bona fide SMC pool.",
    "log2FC and Wilcoxon p computed against the Aged_Ctrl group average (reference baseline).",
    "Confirms genuine Piezo1 knockout in KO samples (significantly reduced mean and detection)."
  )
  return(list(data = results, legend = legend))
}

# ==============================================================================
# TABLE S1c: Cluster Composition by Group
# ==============================================================================
create_table_s1c_cluster_by_group <- function(metadata) {
  md <- metadata
  md$cluster_str <- paste0("Cl", as.character(md$seurat_clusters))
  md$cluster_str[md$cluster_str == "ClNA"] <- "Cl_NA"

  groups <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
  tab <- as.data.frame.matrix(table(md$cluster_str, md$group))
  for (g in groups) if (!g %in% colnames(tab)) tab[[g]] <- 0
  tab <- tab[, groups]
  tab$Cluster <- rownames(tab)
  tab$Total <- rowSums(tab[, groups])

  # Bona fide cluster flag
  tab$`Bona fide SMC pool` <- ifelse(tab$Cluster %in% c("Cl1","Cl5","Cl15","Cl16"), "Yes", "No")

  # Best-effort cell type label (from majority cell_type per cluster)
  cell_type_by_cluster <- sapply(rownames(tab), function(cl) {
    cells <- md$cell_type[md$cluster_str == cl]
    if (length(cells) == 0) return(NA_character_)
    tt <- sort(table(cells), decreasing = TRUE)
    if (length(tt) > 0) names(tt)[1] else NA_character_
  })
  tab$`Cell Type` <- cell_type_by_cluster

  # Order: bona fide first, then total desc
  tab <- tab[order(tab$`Bona fide SMC pool` == "No", -tab$Total), ]
  tab <- tab[, c("Cluster","Cell Type",groups,"Total","Bona fide SMC pool")]

  legend <- paste(
    "Per-cluster cell counts split by experimental group on the post-filter dataset.",
    "Bona fide SMC clusters (Cl1, Cl5, Cl15, Cl16) listed first.",
    "Use this table to see which clusters expand/contract with aging or Piezo1 KO."
  )
  return(list(data = tab, legend = legend))
}

# ==============================================================================
# TABLE S4b: Panel Scores by Group
# ==============================================================================
# Per-panel mean ± SD across the 9 active scoring panels, plus aging + rescue effects
# and Wilcoxon p-values. The headline numbers from Fig 3G/3H/4B-E/S3J/K.
create_table_s4b_panel_scores <- function(expr, smc_cells, metadata, gene_sets) {
  panels <- list(
    "Contractile Capacity"        = gene_sets$contractile_core,
    "Synthetic Phenotype" = gene_sets$synthetic_phenotype %||% c("Vim","Col1a2","Dcn","Igfbp6","Efemp1"),
    "Synthetic Classification"    = gene_sets$synthetic_core,
    "NFAT Pathway (upstream)"     = gene_sets$nfat_pathway,
    "NFAT Activity (downstream)"  = gene_sets$nfatc3_remodeling,
    "ECM Score"                   = gene_sets$ecm_structural %||% c("Bgn","Col6a1","Dcn","Sparc","Timp2"),
    "Cytoskeletal Remodeling"     = gene_sets$cascade_cyto_remodeling %||% c("Ezr","Syne2","Lmna","Vim"),
    "Mechanosensors (Piezo1)"     = "Piezo1",
    "Mechanosensors (Piezo2)"     = "Piezo2"
  )
  groups <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
  md_bf <- metadata[smc_cells, ]

  rows <- list()
  for (panel_name in names(panels)) {
    genes <- intersect(panels[[panel_name]], rownames(expr))
    if (length(genes) == 0) next
    score <- colMeans(as.matrix(expr[genes, smc_cells, drop = FALSE]))
    fmt <- function(g) {
      s <- score[md_bf$group == g]
      sprintf("%.4f ± %.4f (n=%d)", mean(s, na.rm=TRUE), sd(s, na.rm=TRUE), length(s))
    }
    yc <- score[md_bf$group == "Young_Ctrl"]
    ac <- score[md_bf$group == "Aged_Ctrl"]
    ak <- score[md_bf$group == "Aged_KO"]
    aging_p  <- tryCatch(suppressWarnings(wilcox.test(ac, yc)$p.value), error = function(e) NA)
    rescue_p <- tryCatch(suppressWarnings(wilcox.test(ak, ac)$p.value), error = function(e) NA)
    rows[[length(rows)+1]] <- data.frame(
      Panel = panel_name,
      Genes = paste(genes, collapse = ", "),
      `Young Ctrl` = fmt("Young_Ctrl"),
      `Aged Ctrl`  = fmt("Aged_Ctrl"),
      `Young KO`   = fmt("Young_KO"),
      `Aged KO`    = fmt("Aged_KO"),
      `Aging Effect (A-Y Ctrl)`   = round(mean(ac, na.rm=TRUE) - mean(yc, na.rm=TRUE), 4),
      `Rescue Effect (AKO-ACtrl)` = round(mean(ak, na.rm=TRUE) - mean(ac, na.rm=TRUE), 4),
      `Aging P-value`  = signif(aging_p, 3),
      `Rescue P-value` = signif(rescue_p, 3),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  results <- do.call(rbind, rows)

  legend <- paste(
    "Aggregated panel SCORES (mean of normalized expression across panel genes, computed",
    "per cell, then summarized per experimental group). One row per panel x 9 panels.",
    "Reports per-group mean +- SD with n cells, plus aging effect, rescue effect,",
    "and Wilcoxon p-values (scores shown in Fig 3E, 3F, 3M, 3N, 4D)."
  )
  return(list(data = results, legend = legend))
}

# ==============================================================================
# TABLE S4c: Phenotype Proportions per Group + Cl5 sub-cluster split
# ==============================================================================
create_table_s4c_phenotype_proportions <- function(metadata, smc_cells, context) {
  phenotype <- context$phenotype
  if (is.null(phenotype) || nrow(phenotype) == 0) {
    return(list(data = data.frame(Note = "Phenotype classification not in context"),
                legend = "S4c not generated."))
  }

  md_bf <- metadata[smc_cells, ]
  md_bf$cell_barcode <- rownames(md_bf)
  if (!"phenotype" %in% colnames(phenotype)) {
    return(list(data = data.frame(Note = "phenotype column missing"),
                legend = "S4c not generated."))
  }
  ph <- phenotype[, c("cell_barcode","phenotype")]
  md_bf <- merge(md_bf, ph, by = "cell_barcode", all.x = TRUE)

  groups <- c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO")
  # Whole-pool proportions
  prop_rows <- lapply(groups, function(g) {
    sub <- md_bf[md_bf$group == g, ]
    n_total <- nrow(sub)
    n_c <- sum(sub$phenotype == "Contractile", na.rm = TRUE)
    n_s <- sum(sub$phenotype == "Synthetic", na.rm = TRUE)
    data.frame(
      Group = g,
      `Total cells` = n_total,
      `Contractile (n)` = n_c,
      `Contractile (%)` = round(100 * n_c / max(n_total, 1), 1),
      `Synthetic (n)` = n_s,
      `Synthetic (%)` = round(100 * n_s / max(n_total, 1), 1),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  prop_df <- do.call(rbind, prop_rows)

  # Cl5 split
  cl5 <- md_bf[as.character(md_bf$seurat_clusters) == "5", ]
  cl5_rows <- lapply(groups, function(g) {
    sub <- cl5[cl5$group == g, ]
    n_c <- sum(sub$phenotype == "Contractile", na.rm = TRUE)
    n_s <- sum(sub$phenotype == "Synthetic", na.rm = TRUE)
    data.frame(
      Group = g,
      `Cl5 Contractile` = n_c,
      `Cl5 Synthetic` = n_s,
      `Cl5 total` = n_c + n_s,
      `Cl5 % Synthetic` = round(100 * n_s / max(n_c + n_s, 1), 1),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  cl5_df <- do.call(rbind, cl5_rows)

  # Stitch with section headers (matches build_supp_tables.py output)
  hdr1 <- as.data.frame(matrix("", nrow = 1, ncol = ncol(prop_df) + ncol(cl5_df) - 1),
                         stringsAsFactors = FALSE)
  colnames(hdr1) <- c(colnames(prop_df), colnames(cl5_df)[-1])
  hdr1[1, "Group"] <- "=== Whole bona fide pool: phenotype classification (per-cell ratio > 1.5 = Contractile) ==="

  prop_pad <- prop_df
  for (c in colnames(cl5_df)[-1]) prop_pad[[c]] <- ""

  blank_row <- hdr1; blank_row[1, ] <- ""

  hdr2 <- hdr1; hdr2[1, ] <- ""
  hdr2[1, "Group"] <- "=== Cl5 sub-cluster phenotype split (the rescued-synthetic cluster) ==="

  cl5_pad <- as.data.frame(matrix("", nrow = nrow(cl5_df), ncol = ncol(prop_df) - 1),
                            stringsAsFactors = FALSE)
  colnames(cl5_pad) <- colnames(prop_df)[-1]
  cl5_pad <- cbind(Group = cl5_df$Group, cl5_pad, cl5_df[, -1])

  combined <- rbind(hdr1, prop_pad, blank_row, hdr2, cl5_pad)

  legend <- paste(
    "Two related blocks. TOP: percentage of bona fide SMCs classified as Contractile vs",
    "Synthetic per group (= Fig 3B data). BOTTOM: Cl5-specific split — Cl5 is the",
    "'Rescued Synthetic' cluster contributing nearly all Synthetic-classified cells.",
    "Aging shifts Aged Ctrl to 23.9% Synthetic; Piezo1 KO restores Aged to 96.8% Contractile."
  )
  return(list(data = combined, legend = legend))
}

# null-coalesce helper used above
`%||%` <- function(a, b) if (!is.null(a)) a else b
