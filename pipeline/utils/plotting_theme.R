# ==============================================================================
# STANDARDIZED PLOTTING THEME FOR SMC ANALYSIS
# ==============================================================================
# Consistent styling across all figures
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ==============================================================================
# COLOR PALETTES
# ==============================================================================

#' Get standardized color palette
#' @return Named list of color vectors
get_colors <- function() {
  list(
    # Genotype colors
    genotype = c(
      "Ctrl" = "#7F7F7F",
      "KO" = "#ff7fb2"
    ),

    # Group colors (all 4 groups)
    group = c(
      "Young_Ctrl" = "#7F7F7F",
      "Aged_Ctrl" = "#7F7F7F",
      "Young_KO" = "#ff7fb2",
      "Aged_KO" = "#ff7fb2"
    ),

    # Age colors (for correlation plots)
    age = c(
      "Young" = "#4DAF4A",
      "Aged" = "#E41A1C"
    ),

    # SMC phenotype colors
    phenotype = c(
      "Contractile" = "#E41A1C",
      "Synthetic" = "#377EB8"
    ),

    # Border/accent colors
    borders = c(
      "KO_border" = "#fe0267",
      "default" = "black"
    ),

    # Heatmap gradient
    heatmap = c(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B"
    ),

    # Cell type colors (for UMAP)
    cell_types = c(
      "SMC_Contractile" = "#E41A1C",
      "Fibroblast" = "#377EB8",
      "ICC" = "#4DAF4A",
      "Goblet" = "#984EA3",
      "Glial" = "#FF7F00",
      "Endothelial" = "#A65628",
      "Immune_Macrophage" = "#F781BF",
      "Immune_TCell" = "#999999",
      "Epithelial" = "#66C2A5",
      "Pericyte" = "#FC8D62",
      "Mesothelial" = "#8DA0CB",
      "Other" = "#CCCCCC"
    ),

    # Greyed out cells
    inactive = "#CCCCCC"
  )
}

#' Get group ordering for x-axis
#' @return Character vector of ordered group names
get_group_order <- function() {
  c("Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO")
}

#' Get x-axis labels for 4-group bar plots — age in months only.
#' Genotype is conveyed by bar color (Ctrl grey, P1KO pink) plus the legend,
#' so the x-axis labels only need the age. Changed from the
#' "CTRL 3" / "P1KO 24" form to keep bar charts uncluttered.
#'
#' Note: UMAP per-panel titles (Fig3E, FigS3E) and correlation panels
#' (Fig3F/G, Fig4F/G) intentionally use the fuller form with explicit
#' genotype prefix because those panels don't carry a genotype color legend.
#' @return Named character vector
get_group_labels <- function() {
  c(
    "Young_Ctrl" = "3",
    "Aged_Ctrl"  = "24",
    "Young_KO"   = "3",
    "Aged_KO"    = "24"
  )
}

# ==============================================================================
# GGPLOT THEMES
# ==============================================================================

#' Standard theme for publication figures
#' @param base_size Base font size (default 14)
#' @return ggplot theme object
theme_smc <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      # Title
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 2, hjust = 0.5, color = "black"),

      # Axes
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),

      # Legend
      legend.title = element_text(face = "bold", size = base_size - 2),
      legend.text = element_text(size = base_size - 3),
      legend.position = "bottom",

      # Panel
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.x = element_blank(),

      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}

#' PUBLICATION THEME - STANDARD FOR ALL MANUSCRIPT FIGURES
#' Sized for 6-8" wide exports at 300 DPI. Text renders at ~8-12pt when
#' printed at single-column width (3.3"). Meets Nature/Cell minimum (5-7pt).
#' @param base_size Base font size (default 18 for publication readability)
#' @return ggplot theme object
theme_publication <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      # Title: LARGE and BOLD, centered
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,  # 20pt
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = base_size - 4,  # 14pt
        hjust = 0.5,
        color = "black",
        margin = margin(b = 8)
      ),

      # Axes: Bold titles, large readable text
      axis.title = element_text(
        face = "bold",
        size = base_size      # 18pt
      ),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(
        size = base_size - 2,  # 16pt
        color = "black"
      ),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),

      # Legend: Bold title, readable text
      legend.title = element_text(
        face = "bold",
        size = base_size - 2  # 16pt
      ),
      legend.text = element_text(
        size = base_size - 4  # 14pt
      ),
      legend.position = "bottom",
      legend.key.size = unit(1, "cm"),

      # Panel
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),

      # Margins
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

#' Theme for UMAP plots
#' @param base_size Base font size
#' @return ggplot theme object
theme_umap <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "black"),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 2, color = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      aspect.ratio = 1
    )
}

#' Theme for heatmaps
#' @param base_size Base font size
#' @return ggplot theme object
theme_heatmap <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "black"),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "italic", size = base_size - 2),
      axis.text.x = element_text(face = "bold", size = base_size - 1),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Format significance stars
#' @param p_value Numeric p-value
#' @return Character string with stars or "ns"
format_significance <- function(p_value) {
  ifelse(p_value < 0.001, "***",
         ifelse(p_value < 0.01, "**",
                ifelse(p_value < 0.05, "*", "ns")))
}

#' Vectorized significance formatting
#' @param p_values Numeric vector of p-values
#' @return Character vector
format_significance_vec <- function(p_values) {
  sapply(p_values, format_significance)
}

#' Create color scale for genotype
#' @return ggplot scale
scale_fill_genotype <- function() {
  colors <- get_colors()
  scale_fill_manual(values = colors$genotype, name = "Genotype")
}

#' Create color scale for age
#' @return ggplot scale
scale_color_age <- function() {
  colors <- get_colors()
  scale_color_manual(values = colors$age, name = "Age")
}

#' Create color scale for phenotype
#' @return ggplot scale
scale_fill_phenotype <- function() {
  colors <- get_colors()
  scale_fill_manual(values = colors$phenotype, name = "Phenotype")
}
