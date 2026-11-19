# ==============================================================================
# Tables S0-S7 in the manuscript layout (Tables_S0_S7.xlsx)
# ==============================================================================
# Called at the end of SupplementalTables.R generate(). Uses the analysed SMCs
# and the gene sets of the manuscript figures:
#   NFAT pathway = 12 genes (gene_sets$nfat_pathway without Nfat5), as in Fig 4A-4C
#   NFAT activity / NFAT-associated remodeling = 5 genes (gene_sets$nfatc3_remodeling
#   without Fas), as in Fig 4A and 4D
# Descriptive text (Table S1 cell-type references, Table S2 panel descriptions) and the
# two Table S0 rows for libraries not in the deposited object are read from data/.
# Each sheet: row 1 = table title, row 2 = column headers, data from row 3.
# ==============================================================================

paper_gene_sets <- function(gene_sets) {
  list(
    contractile          = gene_sets$contractile_core,
    synthetic_phenotype  = gene_sets$synthetic_phenotype,
    synthetic_class      = gene_sets$synthetic_core,
    nfat_pathway         = setdiff(gene_sets$nfat_pathway, "Nfat5"),
    nfat_activity        = setdiff(gene_sets$nfatc3_remodeling, "Fas"),
    ecm                  = gene_sets$ecm_structural,
    cyto                 = gene_sets$cascade_cyto_remodeling
  )
}

p3 <- function(p) ifelse(is.na(p), NA, signif(p, 3))

group_cells <- function(smc_cells, metadata) {
  g <- c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")
  setNames(lapply(g, function(k) smc_cells[metadata[smc_cells, "group"] == k]), g)
}

# per-gene means, effects and Wilcoxon p (aging: Aged_Ctrl vs Young_Ctrl; rescue: Aged_KO vs Aged_Ctrl)
gene_stats <- function(expr, gene, G) {
  v <- lapply(G, function(cc) as.numeric(expr[gene, cc]))
  m <- sapply(v, mean)
  pa <- tryCatch(wilcox.test(v$Aged_Ctrl, v$Young_Ctrl)$p.value, error = function(e) NA)
  pr <- tryCatch(wilcox.test(v$Aged_KO, v$Aged_Ctrl)$p.value, error = function(e) NA)
  list(m = m, aging = unname(m["Aged_Ctrl"] - m["Young_Ctrl"]), rescue = unname(m["Aged_KO"] - m["Aged_Ctrl"]),
       pa = if (is.nan(pa)) NA else pa, pr = if (is.nan(pr)) NA else pr)
}

paper_s0 <- function(context, metadata, smc_cells) {
  rs <- context$raw_sample_stats; dm <- context$depth_match_info
  order <- c("WT_Young_1", "WT_Young_2", "KO_Young_1", "WT_Aged_1", "WT_Aged_2", "KO_Aged")
  order <- order[order %in% rs$sample]
  grp <- sapply(order, function(s) as.character(metadata$group[metadata$sample == s][1]))
  bf <- table(metadata[smc_cells, "sample"])
  out <- data.frame(Sample = order, Group = grp,
                    `Raw cells` = rs$n_cells[match(order, rs$sample)],
                    `Median UMI` = round(rs$median_umi[match(order, rs$sample)]),
                    `Bona fide SMCs` = as.integer(bf[order]),
                    Status = "INCLUDED", Reason = "Pass all QC criteria",
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(dm)) for (s in dm$outlier_samples) {
    i <- match(s, out$Sample)
    out$Status[i] <- "INCLUDED (depth-matched)"
    out$Reason[i] <- sprintf(paste0("Median UMI %s after removal of non-muscularis cells, outside the [%s, %s] window ",
                                    "(%.1f-%.1fx the median of per-sample medians, %s); cells restricted to that window, ",
                                    "%s of %s retained"),
                             format(round(dm$sample_medians[[s]]), big.mark = ","), format(round(dm$range[1]), big.mark = ","),
                             format(round(dm$range[2]), big.mark = ","), 1 / dm$range_factor, dm$range_factor,
                             format(round(dm$target_depth), big.mark = ","),
                             format(dm$n_after[[s]], big.mark = ","), format(dm$n_before[[s]], big.mark = ","))
  }
  out
}

paper_s1b <- function(s1c) {
  lab <- read.csv("data/table_S1b_cluster_labels.csv", check.names = FALSE, stringsAsFactors = FALSE)
  d <- s1c
  hit <- match(d$Cluster, lab$Cluster)
  d$`Cell Type` <- ifelse(is.na(hit), d$`Cell Type`, lab$`Cell Type`[hit])
  d
}

paper_s3 <- function(expr, G, P) {
  panels <- list("Cytoskeletal Remodeling" = P$cyto, "ECM Score" = P$ecm, "NFAT Activity" = P$nfat_activity,
                 "NFAT Pathway" = P$nfat_pathway, "Synthetic Phenotype" = P$synthetic_phenotype)
  rows <- list()
  for (pw in names(panels)) for (g in sort(panels[[pw]])) {
    s <- gene_stats(expr, g, G)
    others <- setdiff(names(panels)[sapply(panels, function(x) g %in% x)], pw)
    rows[[length(rows) + 1]] <- data.frame(Pathway = pw, Gene = g,
      `Mean Young Ctrl` = round(s$m[["Young_Ctrl"]], 4), `Mean Aged Ctrl` = round(s$m[["Aged_Ctrl"]], 4),
      `Mean Young KO` = round(s$m[["Young_KO"]], 4), `Mean Aged KO` = round(s$m[["Aged_KO"]], 4),
      `Aging Effect` = round(s$aging, 4), `Rescue Effect` = round(s$rescue, 4),
      `Aging P-value` = p3(s$pa), `Rescue P-value` = p3(s$pr),
      `Also In Pathway` = if (length(others)) paste("Also in:", paste(others, collapse = ", ")) else NA,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

paper_s4 <- function(expr, G) {
  gl <- read.csv("data/table_S4_genes.csv", check.names = FALSE, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(gl)), function(i) {
    g <- gl$Gene[i]; r <- data.frame(Gene = g, Category = gl$Category[i], check.names = FALSE, stringsAsFactors = FALSE)
    for (k in names(G)) { v <- as.numeric(expr[g, G[[k]]]); lb <- gsub("_", " ", k)
      r[[paste(lb, "Mean")]] <- round(mean(v), 4); r[[paste(lb, "SD")]] <- round(sd(v), 4)
      r[[paste(lb, "Pct Expressing")]] <- round(100 * mean(v > 0), 1) }
    r })
  do.call(rbind, rows)
}

paper_s4b <- function(expr, G, P) {
  panels <- list("Contractile Capacity" = P$contractile, "Synthetic Phenotype" = P$synthetic_phenotype,
                 "Synthetic Classification" = P$synthetic_class, "NFAT Pathway (upstream)" = P$nfat_pathway,
                 "NFAT Activity (downstream)" = P$nfat_activity, "ECM Score" = P$ecm, "Cytoskeletal Remodeling" = P$cyto,
                 "Mechanosensors (Piezo1)" = "Piezo1", "Mechanosensors (Piezo2)" = "Piezo2")
  rows <- lapply(names(panels), function(nm) {
    g <- panels[[nm]]
    sc <- lapply(G, function(cc) colMeans(as.matrix(expr[g, cc, drop = FALSE])))
    f <- function(k) sprintf("%.4f ± %.4f (n=%d)", mean(sc[[k]]), sd(sc[[k]]), length(sc[[k]]))
    data.frame(Panel = nm, Genes = paste(g, collapse = ", "),
               `Young Ctrl` = f("Young_Ctrl"), `Aged Ctrl` = f("Aged_Ctrl"), `Young KO` = f("Young_KO"), `Aged KO` = f("Aged_KO"),
               `Aging Effect (A-Y Ctrl)` = round(mean(sc$Aged_Ctrl) - mean(sc$Young_Ctrl), 4),
               `Rescue Effect (AKO-ACtrl)` = round(mean(sc$Aged_KO) - mean(sc$Aged_Ctrl), 4),
               `Aging P-value` = p3(wilcox.test(sc$Aged_Ctrl, sc$Young_Ctrl)$p.value),
               `Rescue P-value` = p3(wilcox.test(sc$Aged_KO, sc$Aged_Ctrl)$p.value),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Tables S5/S6: calculate_rescue_stats() (SupplementalTables.R), Pct Rescued capped at 100 with the raw value kept
paper_rescue <- function(expr, smc_cells, metadata, genes, type) {
  d <- calculate_rescue_stats(expr, smc_cells, metadata, genes, type)
  raw <- d$`Pct Rescued`
  d$`Pct Rescued` <- pmin(raw, 100)
  for (k in c("Aging P-value", "Rescue P-value")) d[[k]] <- p3(suppressWarnings(as.numeric(d[[k]])))
  d$Capped <- !is.na(raw) & raw > 100
  d$`Pct Rescued (raw)` <- raw
  d
}

paper_s7 <- function(expr, smc_cells, gene_sets, P) {
  targets <- unique(c(gene_sets$contractile_core, gene_sets$contractile_extended, gene_sets$synthetic_core,
                      gene_sets$synthetic_extended, gene_sets$nfat_pathway, gene_sets$nfatc3_remodeling,
                      gene_sets$ecm_structural, "Tgfb1", "Klf4", P$cyto, P$synthetic_phenotype))
  targets <- targets[targets %in% rownames(expr)]
  p1 <- as.numeric(expr["Piezo1", smc_cells])
  cat_of <- function(g) {
    if (g %in% gene_sets$contractile_core) "Contractile Core"
    else if (g %in% gene_sets$contractile_extended) "Contractile Extended"
    else if (g %in% gene_sets$synthetic_core) "Synthetic Core"
    else if (g %in% gene_sets$nfat_pathway) "NFAT Pathway"
    else if (g %in% gene_sets$nfatc3_remodeling) "NFAT Activity"
    else if (g %in% gene_sets$ecm_structural) "ECM Score"
    else if (g %in% gene_sets$synthetic_extended) "Synthetic Extended"
    else if (g %in% P$cyto) "Cytoskeletal Remodeling"
    else if (g %in% P$synthetic_phenotype) "Synthetic Phenotype"
    else "Other"
  }
  d <- do.call(rbind, lapply(targets, function(g) {
    ct <- suppressWarnings(cor.test(p1, as.numeric(expr[g, smc_cells]), method = "spearman"))
    data.frame(Gene = g, Category = cat_of(g), `Spearman Rho` = round(unname(ct$estimate), 4), `P-value` = p3(ct$p.value),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  d[order(-d$`Spearman Rho`), ]
}

# Table S1: descriptive text from data/table_S1_cell_types.csv; every number is computed here.
#   {n:ClA,ClB}     cells in those clusters in the deposited object (all 7 libraries)
#   {pct:Gene:ClX}  % of cells in cluster X with >= 1 UMI of Gene
#   {contractile}, {synthetic}, {..._Cl5}, {Cl5_total}, {Cl5_median_ratio}, {Cl5_pct_...}: phenotype classification of this run
fill_table_s1 <- function(phenotype) {
  obj <- readRDS("data/seurat_paper_6samples.rds")
  cl <- paste0("Cl", as.character(obj$seurat_clusters))
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  fmt <- function(x) format(x, big.mark = ",", trim = TRUE)
  ph <- phenotype; ph$cluster <- cl[match(ph$cell_barcode, colnames(obj))]; stopifnot(!anyNA(ph$cluster))
  c5 <- ph[ph$cluster == "Cl5", ]
  v <- c(contractile = fmt(sum(ph$phenotype == "Contractile")), synthetic = fmt(sum(ph$phenotype == "Synthetic")),
         contractile_Cl5 = fmt(sum(c5$phenotype == "Contractile")), synthetic_Cl5 = fmt(sum(c5$phenotype == "Synthetic")),
         Cl5_total = fmt(nrow(c5)), Cl5_median_ratio = sprintf("%.2f", median(c5$score_ratio)),
         Cl5_pct_contractile = round(100 * mean(c5$phenotype == "Contractile")),
         Cl5_pct_synthetic = round(100 * mean(c5$phenotype == "Synthetic")))
  fill <- function(txt) {
    for (m in unique(regmatches(txt, gregexpr("\\{[^}]+\\}", txt))[[1]])) {
      key <- substr(m, 2, nchar(m) - 1); parts <- strsplit(key, ":")[[1]]
      val <- if (parts[1] == "n") fmt(sum(cl %in% strsplit(parts[2], ",")[[1]]))
             else if (parts[1] == "pct") round(100 * mean(counts[parts[2], cl == parts[3]] > 0))
             else v[[key]]
      txt <- gsub(m, val, txt, fixed = TRUE)
    }
    txt
  }
  t <- read.csv("data/table_S1_cell_types.csv", check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_along(t)) t[[j]] <- vapply(t[[j]], fill, character(1), USE.NAMES = FALSE)
  stopifnot(!any(grepl("\\{", unlist(t))))
  t
}

write_paper_tables_S0_S7 <- function(seurat_obj, context, expr, smc_cells, metadata, gene_sets, s0b, s1c, s4c, table_dir) {
  P <- paper_gene_sets(gene_sets)
  stopifnot(length(P$nfat_pathway) == 12, length(P$nfat_activity) == 5)
  G <- group_cells(smc_cells, metadata)
  s6_genes <- unique(c(gene_sets$synthetic_core, gene_sets$synthetic_score, gene_sets$ecm_remodeling, P$synthetic_phenotype))
  s5 <- paper_rescue(expr, smc_cells, metadata, unique(c(gene_sets$contractile_core, gene_sets$contractile_score,
                                                           gene_sets$muscle_contraction)), "Contractile")
  s5 <- s5[, c("Gene", "Mean Young Ctrl", "Mean Aged Ctrl", "Mean Young KO", "Mean Aged KO", "Aging Effect", "Rescue Effect",
               "Pct Rescued", "Aging P-value", "Rescue P-value", "Pct Rescued (raw)")]
  s6 <- paper_rescue(expr, smc_cells, metadata, s6_genes, "Synthetic")
  # S6 order: genes passing all criteria first, each block by Pct Rescued (raw), highest first
  s6 <- s6[order(!s6$`Passes All Criteria`, -ifelse(is.na(s6$`Pct Rescued (raw)`), -Inf, s6$`Pct Rescued (raw)`)), ]
  s6 <- s6[, c("Gene", "Mean Young Ctrl", "Mean Aged Ctrl", "Mean Young KO", "Mean Aged KO", "Aging Effect", "Rescue Effect",
               "Pct Rescued", "Aging P-value", "Rescue P-value", "Passes All Criteria", "Capped", "Pct Rescued (raw)")]
  s0b_d <- s0b$data; s0b_d$`Wilcoxon p (vs A_Ctrl)` <- as.character(s0b_d$`Wilcoxon p (vs A_Ctrl)`)
  tabs <- list(
    "S0 Sample Inclusion"         = list("S0: Sample Inclusion + QC", paper_s0(context, metadata, smc_cells)),
    "S0b Piezo1 KO Validation"    = list("S0b: Piezo1 KO Validation per Sample", s0b_d),
    "S1 Cell Type Markers"        = list("S1: Cell Type Markers (canonical reference)",
                                         fill_table_s1(context$phenotype)),
    "S1b Cluster by Group"        = list("S1b: Cluster Composition by Group", paper_s1b(s1c$data)),
    "S2 Gene Panels"              = list("S2: Gene Panels (consolidated classification + scoring panels)",
                                         read.csv("data/table_S2_gene_panels.csv", check.names = FALSE, stringsAsFactors = FALSE)),
    "S3 Pathway GeneSets"         = list("S3: Pathway Gene Sets ; per-gene aging+rescue stats", paper_s3(expr, G, P)),
    "S4 Gene Expression by Group" = list("S4: Gene Expression by Group", paper_s4(expr, G)),
    "S4b Panel Scores by Group"   = list("S4b: Panel Scores by Group", paper_s4b(expr, G, P)),
    "S4c Phenotype Proportions"   = list("S4c: Phenotype proportions", s4c$data),
    "S5 Contractile Rescue"       = list("S5: Contractile Rescue", s5),
    "S6 Synthetic Rescue"         = list("S6: Synthetic Rescue", s6),
    "S7 Piezo1 Correlations"      = list("S7: Piezo1 Correlations", paper_s7(expr, smc_cells, gene_sets, P)))
  # S2 gene lists must match the panels used
  s2 <- tabs[["S2 Gene Panels"]][[2]]
  chk <- c("NFAT Pathway" = paste(P$nfat_pathway, collapse = ", "), "NFAT Activity" = paste(P$nfat_activity, collapse = ", "),
           "Synthetic Phenotype Score" = paste(P$synthetic_phenotype, collapse = ", "), "ECM Score" = paste(P$ecm, collapse = ", "),
           "Cytoskeletal Remodeling" = paste(P$cyto, collapse = ", "), "Contractile" = paste(P$contractile, collapse = ", "))
  for (k in names(chk)) stopifnot(setequal(strsplit(s2$Genes[s2$Panel == k], ", ")[[1]], strsplit(chk[[k]], ", ")[[1]]))
  wb <- createWorkbook()
  for (nm in names(tabs)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, tabs[[nm]][[1]], startRow = 1, startCol = 1)
    d <- tabs[[nm]][[2]]
    if (nm == "S4c Phenotype Proportions") {
      for (j in seq_along(d)) { x <- as.character(d[[j]]); x[!is.na(x) & (x == "" | startsWith(x, "==="))] <- NA; d[[j]] <- x }
      for (j in 2:ncol(d)) d[[j]] <- suppressWarnings(as.numeric(d[[j]]))
    }
    writeData(wb, nm, d, startRow = 2, startCol = 1, keepNA = FALSE)
  }
  saveWorkbook(wb, file.path(table_dir, "Tables_S0_S7.xlsx"), overwrite = TRUE)
  cat("  Saved: Tables_S0_S7.xlsx (manuscript layout)\n")
  invisible(tabs)
}
