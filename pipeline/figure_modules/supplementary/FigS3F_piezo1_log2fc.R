# ==============================================================================
# PAPER PANEL: Fig S3C (Table S0b)   [internal module ID: FigS3F]
# PLOTS:  log2 fold change of Piezo1 expression, P1KO vs CTRL, at each age
# CELLS:  visceral SMCs (bona fide clusters 1/5/15/16, doublets removed; n = 1,604: 469 / 477 / 222 / 436)
# GROUPS COMPARED: Young P1KO vs Young CTRL; Aged P1KO vs Aged CTRL
# TEST:   Wilcoxon
# ==============================================================================
# Shows Piezo1 KO validation in Young and Aged SMCs
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

#' Generate FigS3F
#'
#' @param context Analysis context from run_pipeline
#' @return List with success status and output paths
generate <- function(context) {
  # Source utilities
  pipeline_root <- context$config$pipeline_root
  source(file.path(pipeline_root, "utils/plotting_theme.R"))
  source(file.path(pipeline_root, "utils/statistics.R"))
  source(file.path(pipeline_root, "core/output_manager.R"))

  colors <- get_colors()

  # Get data
  seurat_obj <- context$seurat_obj
  selection <- context$selection
  smc_cells <- selection$smc

  metadata <- seurat_obj@meta.data
  metadata$cell_barcode <- rownames(metadata)

  # Get expression data
  expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

  # Get Piezo1 expression in SMCs
  if (!"Piezo1" %in% rownames(expr)) {
    return(list(success = FALSE, error = "Piezo1 not found in expression matrix"))
  }

  piezo1_expr <- as.numeric(expr["Piezo1", smc_cells])

  # Create data frame
  plot_data <- data.frame(
    cell_barcode = smc_cells,
    piezo1 = piezo1_expr,
    age = metadata[smc_cells, "age"],
    genotype = metadata[smc_cells, "genotype"],
    group = metadata[smc_cells, "group"],
    stringsAsFactors = FALSE
  )

  # Calculate mean expression by age and genotype
  summary_by_age <- plot_data %>%
    group_by(age, genotype) %>%
    summarize(
      mean_expr = mean(piezo1, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  # Calculate log2FC for each age
  young_ctrl <- summary_by_age$mean_expr[summary_by_age$age == "Young" &
                                          summary_by_age$genotype == "Ctrl"]
  young_ko <- summary_by_age$mean_expr[summary_by_age$age == "Young" &
                                        summary_by_age$genotype == "KO"]
  aged_ctrl <- summary_by_age$mean_expr[summary_by_age$age == "Aged" &
                                         summary_by_age$genotype == "Ctrl"]
  aged_ko <- summary_by_age$mean_expr[summary_by_age$age == "Aged" &
                                       summary_by_age$genotype == "KO"]

  pseudo <- 0.01
  log2fc_young <- log2((young_ko + pseudo) / (young_ctrl + pseudo))
  log2fc_aged <- log2((aged_ko + pseudo) / (aged_ctrl + pseudo))

  # Statistical tests
  test_young <- wilcox_test(
    plot_data$piezo1, plot_data$group,
    "Young_Ctrl", "Young_KO"
  )

  test_aged <- wilcox_test(
    plot_data$piezo1, plot_data$group,
    "Aged_Ctrl", "Aged_KO"
  )

  # Bootstrap SE for the per-age log2FC (1000 resamples within each cell pool)
  set.seed(42)
  n_boot <- 1000
  boot_log2fc_se <- function(ctrl_vals, ko_vals) {
    out <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      ck <- sample(ctrl_vals, length(ctrl_vals), replace = TRUE)
      kk <- sample(ko_vals,   length(ko_vals),   replace = TRUE)
      out[i] <- log2((mean(kk) + pseudo) / (mean(ck) + pseudo))
    }
    sd(out)
  }
  young_ctrl_vals <- plot_data$piezo1[plot_data$group == "Young_Ctrl"]
  young_ko_vals   <- plot_data$piezo1[plot_data$group == "Young_KO"]
  aged_ctrl_vals  <- plot_data$piezo1[plot_data$group == "Aged_Ctrl"]
  aged_ko_vals    <- plot_data$piezo1[plot_data$group == "Aged_KO"]
  se_young <- boot_log2fc_se(young_ctrl_vals, young_ko_vals)
  se_aged  <- boot_log2fc_se(aged_ctrl_vals,  aged_ko_vals)

  # Create bar plot data — factor order forces Young (3) on left, Aged (24) on right
  bar_data <- data.frame(
    age = factor(c("Young", "Aged"), levels = c("Young", "Aged")),
    age_label = c("3", "24"),
    log2fc = c(log2fc_young, log2fc_aged),
    log2fc_se = c(se_young, se_aged),
    significance = c(test_young$significance, test_aged$significance),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------------
  # Panel A: Log2FC bar plot — KO vs Ctrl per age, with bootstrap SE
  # --------------------------------------------------------------------------
  p_log2fc <- ggplot(bar_data, aes(x = age, y = log2fc, fill = age)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = log2fc - log2fc_se, ymax = log2fc + log2fc_se),
                  width = 0.2, linewidth = 0.7, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(aes(label = significance,
                  y = ifelse(log2fc < 0,
                             log2fc - log2fc_se - 0.3,
                             log2fc + log2fc_se + 0.3)),
              size = 7, fontface = "bold") +
    scale_fill_manual(values = colors$age) +
    scale_x_discrete(labels = c("Young" = "3", "Aged" = "24")) +
    labs(
      title = "Piezo1 KO Validation",
      x = NULL,
      y = expression(Log[2]~Fold~Change~"(KO/Ctrl)")
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # --------------------------------------------------------------------------
  # Panel B: % Piezo1 expressing cells per group
  # Order: paired by age (Young Ctrl, Young KO, Aged Ctrl, Aged KO)
  # so the KO effect is visible side-by-side at each age.
  # X-axis labels: just "3" / "24" — bar color encodes Ctrl vs KO.
  # --------------------------------------------------------------------------
  group_order <- c("Young_Ctrl", "Young_KO", "Aged_Ctrl", "Aged_KO")
  age_labels <- c("Young_Ctrl" = "3", "Young_KO" = "3",
                  "Aged_Ctrl"  = "24", "Aged_KO"  = "24")

  pct_data <- plot_data %>%
    mutate(group = factor(group, levels = group_order)) %>%
    group_by(group) %>%
    summarize(
      n_total = n(),
      n_expressing = sum(piezo1 > 0),
      pct_expressing = 100 * n_expressing / n_total,
      .groups = "drop"
    )

  # Bar colors (Ctrl = grey, KO = pink) — same scheme as the rest of the pipeline
  bar_colors <- c("Young_Ctrl" = "#7F7F7F", "Aged_Ctrl" = "#7F7F7F",
                  "Young_KO"   = "#ff7fb2", "Aged_KO"   = "#ff7fb2")

  p_pct <- ggplot(pct_data, aes(x = group, y = pct_expressing, fill = group)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", pct_expressing)),
              vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    scale_x_discrete(labels = age_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "% Piezo1 Expressing",
      x = NULL,
      y = "% Cells Expressing Piezo1"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  # Combine panels
  p <- cowplot::plot_grid(p_log2fc, p_pct, ncol = 2, rel_widths = c(1, 1.2),
                          labels = c("", ""), align = "h", axis = "bt")

  # Create output directory
  fig_dir <- file.path(context$output_dir, context$fig_id)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data for export
  export_data <- list(
    log2fc_summary = bar_data,
    pct_expressing = pct_data,
    expression_summary = summary_by_age,
    statistical_tests = data.frame(
      comparison = c("Young KO vs Ctrl", "Aged KO vs Ctrl"),
      statistic = c(test_young$statistic, test_aged$statistic),
      p_value = c(test_young$p_value, test_aged$p_value),
      significance = c(test_young$significance, test_aged$significance)
    )
  )

  # Methods text
  methods_text <- paste(
    "Piezo1 knockout efficiency was validated by comparing Piezo1 expression",
    "between KO and Control SMCs at each age. Left panel: Log2 fold change",
    "calculated as log2(KO mean / Ctrl mean) of log-normalized expression.",
    "Right panel: percentage of SMCs with detectable Piezo1 expression (>0",
    "log-normalized counts) in each experimental group.",
    "Statistical significance assessed using Wilcoxon rank-sum test."
  )

  technical_legend <- paste(
    "Left: Bar plot showing log2 fold change of Piezo1 expression (KO vs Ctrl).",
    "Green = Young (3 mo), Red = Aged (24 mo).",
    "Negative values indicate successful Piezo1 knockdown.",
    "Right: Percentage of SMCs expressing Piezo1 (>0 log-normalized counts)",
    "per experimental group. Grey = Ctrl, Pink = KO.",
    "Significance: * p<0.05, ** p<0.01, *** p<0.001."
  )

  # Save outputs
  save_figure_outputs(
    figure_id = context$fig_id,
    plot = p,
    data = export_data,
    methods_text = methods_text,
    technical_legend = technical_legend,
    output_dir = fig_dir,
    width = 10,
    height = 5,
    code_path = "pipeline/figure_modules/supplementary/FigS3F_piezo1_log2fc.R"
  )

  return(list(
    success = TRUE,
    output_dir = fig_dir,
    log2fc = bar_data,
    pct_expressing = pct_data
  ))
}
