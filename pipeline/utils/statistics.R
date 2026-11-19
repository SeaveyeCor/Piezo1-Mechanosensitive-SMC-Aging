# ==============================================================================
# STATISTICAL FUNCTIONS FOR SMC ANALYSIS
# ==============================================================================

#' Perform Wilcoxon test between two groups
#' @param values Numeric vector of values
#' @param groups Factor/character vector of group labels
#' @param group1 Name of first group
#' @param group2 Name of second group
#' @return List with statistic, p-value, and formatted result
wilcox_test <- function(values, groups, group1, group2) {
  vals1 <- values[groups == group1]
  vals2 <- values[groups == group2]

  # Handle edge cases

  if (length(vals1) < 2 || length(vals2) < 2) {
    return(list(
      statistic = NA,
      p_value = NA,
      significance = "ns",
      n1 = length(vals1),
      n2 = length(vals2),
      mean1 = mean(vals1, na.rm = TRUE),
      mean2 = mean(vals2, na.rm = TRUE)
    ))
  }

  test <- tryCatch({
    wilcox.test(vals1, vals2)
  }, error = function(e) {
    list(statistic = NA, p.value = NA)
  })

  list(
    statistic = unname(test$statistic),
    p_value = test$p.value,
    significance = format_significance(test$p.value),
    n1 = length(vals1),
    n2 = length(vals2),
    mean1 = mean(vals1, na.rm = TRUE),
    mean2 = mean(vals2, na.rm = TRUE)
  )
}

#' Calculate Spearman correlation
#' @param x Numeric vector
#' @param y Numeric vector
#' @return List with rho, p-value, and formatted result
spearman_cor <- function(x, y) {
  # Remove NAs
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]

  if (length(x) < 3) {
    return(list(
      rho = NA,
      p_value = NA,
      significance = "ns",
      n = length(x)
    ))
  }

  test <- cor.test(x, y, method = "spearman")

  list(
    rho = unname(test$estimate),
    p_value = test$p.value,
    significance = format_significance(test$p.value),
    n = length(x)
  )
}

#' Calculate log2 fold change with pseudocount
#' @param mean1 Mean of numerator group
#' @param mean2 Mean of denominator group
#' @param pseudo Pseudocount to add (default 0.01)
#' @return Numeric log2FC value
calc_log2fc <- function(mean1, mean2, pseudo = 0.01) {
  log2((mean1 + pseudo) / (mean2 + pseudo))
}

#' Calculate mean and SEM by group
#' @param values Numeric vector
#' @param groups Group labels
#' @return Data frame with group statistics
group_stats <- function(values, groups) {
  df <- data.frame(value = values, group = groups)

  stats <- df %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      n = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      sem = sd(value, na.rm = TRUE) / sqrt(dplyr::n()),
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )

  return(stats)
}

#' Calculate percentage of cells expressing a gene (value > 0)
#' @param values Numeric vector of expression values
#' @param groups Optional group labels
#' @return Numeric percentage or data frame by group
pct_expressing <- function(values, groups = NULL) {
  if (is.null(groups)) {
    return(100 * sum(values > 0, na.rm = TRUE) / length(values))
  }

  df <- data.frame(value = values, group = groups)

  df %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      n_total = dplyr::n(),
      n_expressing = sum(value > 0, na.rm = TRUE),
      pct_expressing = 100 * n_expressing / n_total,
      .groups = "drop"
    )
}

#' Calculate module score (mean of gene set)
#' @param expr_matrix Expression matrix (genes x cells)
#' @param genes Character vector of genes in module
#' @return Numeric vector of scores per cell
calc_module_score <- function(expr_matrix, genes) {
  available <- genes[genes %in% rownames(expr_matrix)]

  if (length(available) == 0) {
    warning("No genes from module found in expression matrix")
    return(rep(NA, ncol(expr_matrix)))
  }

  if (length(available) == 1) {
    return(as.numeric(expr_matrix[available, ]))
  }

  colMeans(as.matrix(expr_matrix[available, ]), na.rm = TRUE)
}

#' Format significance stars
#' @param p_value Numeric p-value
#' @return Character string
format_significance <- function(p_value) {
  if (is.na(p_value)) return("ns")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  return("ns")
}

# ==============================================================================
# PAIRWISE GROUP TESTS — 4 valid comparisons across the standard 4 groups
# ==============================================================================
# Standard policy across the pipeline:
#   Tests 4 pairwise Mann-Whitney U comparisons:
#     1. Young_Ctrl vs Aged_Ctrl  — aging effect in WT
#     2. Young_KO   vs Aged_KO    — aging effect under KO
#     3. Young_Ctrl vs Young_KO   — KO effect at young
#     4. Aged_Ctrl  vs Aged_KO    — KO rescue at aged
#   Excludes the two diagonal cross-group pairs that conflate age + genotype:
#     - Young_Ctrl vs Aged_KO
#     - Young_KO   vs Aged_Ctrl
# These are documented in pipeline/utils/statistics.R::PAIRWISE_COMPARISONS.

#' Standard 4-group pairwise comparison list.
#' Each row: a, b (group names), label (display label).
get_pairwise_comparisons <- function() {
  data.frame(
    a = c("Young_Ctrl", "Young_KO", "Young_Ctrl", "Aged_Ctrl"),
    b = c("Aged_Ctrl",  "Aged_KO",  "Young_KO",   "Aged_KO"),
    label = c("Aging (WT)", "Aging (KO)", "KO effect (young)", "KO rescue (aged)"),
    stringsAsFactors = FALSE
  )
}

#' Run all 4 valid pairwise Mann-Whitney tests on a per-cell scalar metric.
#'
#' Returns both raw and Bonferroni-adjusted p-values (m = number of
#' comparisons). The `significance` column uses the adjusted p-value, so
#' display brackets are conservative by default.
#'
#' @param values Numeric vector (one value per cell)
#' @param groups Character/factor vector of group labels (same length)
#' @param comparisons Optional override for comparison list (default = standard 4)
#' @param adjust Multiple-testing adjustment method (default "bonferroni"; pass
#'               "none" to skip correction)
#' @return Data frame with columns: comparison, a, b, n_a, n_b, mean_a, mean_b,
#'         statistic, p_value, p_adjusted, significance
pairwise_group_tests <- function(values, groups, comparisons = NULL,
                                  adjust = "bonferroni") {
  if (is.null(comparisons)) comparisons <- get_pairwise_comparisons()
  rows <- list()
  for (i in seq_len(nrow(comparisons))) {
    a <- comparisons$a[i]
    b <- comparisons$b[i]
    res <- wilcox_test(values, groups, a, b)
    rows[[i]] <- data.frame(
      comparison = comparisons$label[i],
      a = a, b = b,
      n_a = res$n1, n_b = res$n2,
      mean_a = res$mean1, mean_b = res$mean2,
      statistic = res$statistic,
      p_value = res$p_value,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$p_adjusted <- p.adjust(out$p_value, method = adjust)
  out$significance <- vapply(out$p_adjusted, format_significance, character(1))
  out
}

#' Add pairwise Mann-Whitney brackets to a ggplot bar chart.
#'
#' Brackets are stacked from narrowest span outward so wider brackets sit
#' above narrower ones without overlap. The y-coordinate is computed
#' relative to `y_top` (the highest bar+SEM in the plot).
#'
#' @param p ggplot object
#' @param tests Data frame from pairwise_group_tests()
#' @param group_order Character vector giving the x-position of each group
#'                    (e.g. c("Young_Ctrl","Aged_Ctrl","Young_KO","Aged_KO"))
#' @param y_top Numeric — the top of the existing data (max bar + SEM)
#' @param step  Numeric — vertical spacing between brackets (default y_top * 0.10)
#' @param drop_ns Logical — if TRUE, omit non-significant brackets
#' @return Modified ggplot
add_pairwise_brackets <- function(p, tests, group_order, y_top,
                                   step = NULL, drop_ns = TRUE) {
  if (is.null(step)) step <- y_top * 0.10
  g2x <- setNames(seq_along(group_order), group_order)
  # Compute span and sort narrowest first so brackets nest
  tests$x1   <- pmin(g2x[tests$a], g2x[tests$b])
  tests$x2   <- pmax(g2x[tests$a], g2x[tests$b])
  tests$span <- tests$x2 - tests$x1
  tests <- tests[order(tests$span, tests$x1), ]
  if (drop_ns) tests <- tests[tests$significance != "ns", ]
  if (nrow(tests) == 0) return(p)

  y_cur <- y_top + step * 0.4
  for (i in seq_len(nrow(tests))) {
    x1 <- tests$x1[i]; x2 <- tests$x2[i]
    h  <- step * 0.18
    p <- p +
      annotate("segment", x = x1, xend = x1, y = y_cur, yend = y_cur + h,
               color = "black", linewidth = 0.4) +
      annotate("segment", x = x2, xend = x2, y = y_cur, yend = y_cur + h,
               color = "black", linewidth = 0.4) +
      annotate("segment", x = x1, xend = x2, y = y_cur + h, yend = y_cur + h,
               color = "black", linewidth = 0.4) +
      annotate("text", x = (x1 + x2) / 2, y = y_cur + h * 1.4,
               label = tests$significance[i],
               size = 3.6, fontface = "bold")
    y_cur <- y_cur + step
  }
  p
}
