# ==============================================================================
# OUTPUT MANAGER FOR SMC ANALYSIS PIPELINE
# ==============================================================================
# Handles standardized saving of PNG, PDF, XLSX, and methods.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(openxlsx)
})

#' Save figure outputs in standardized format
#'
#' @param figure_id Figure identifier (e.g., "Fig3A", "Fig3B1")
#' @param plot ggplot object
#' @param data List of data frames for XLSX sheets
#' @param methods_text Character string describing methods
#' @param technical_legend Character string for technical legend
#' @param output_dir Directory to save outputs
#' @param width Plot width in inches (default 8)
#' @param height Plot height in inches (default 6)
#' @param code_path Path to the generating code file
#'
#' @return Invisible NULL
save_figure_outputs <- function(
    figure_id,
    plot,
    data,
    methods_text,
    technical_legend = NULL,
    output_dir,
    width = 8,
    height = 6,
    code_path = NULL
) {
  # Create output directory if needed
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Base filename
  base_name <- figure_id

  # Save PNG
  png_path <- file.path(output_dir, paste0(base_name, ".png"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s\n", basename(png_path)))

 # Save PDF
  pdf_path <- file.path(output_dir, paste0(base_name, ".pdf"))
  ggsave(pdf_path, plot, width = width, height = height)
  cat(sprintf("  Saved: %s\n", basename(pdf_path)))

  # Save XLSX
  xlsx_path <- file.path(output_dir, paste0(base_name, "_data.xlsx"))
  save_xlsx(data, xlsx_path)
  cat(sprintf("  Saved: %s\n", basename(xlsx_path)))

  # Save methods.txt
  methods_path <- file.path(output_dir, paste0(base_name, "_methods.txt"))
  save_methods_txt(
    figure_id = figure_id,
    methods_text = methods_text,
    technical_legend = technical_legend,
    code_path = code_path,
    output_path = methods_path
  )
  cat(sprintf("  Saved: %s\n", basename(methods_path)))

  invisible(NULL)
}

#' Save multiple plots for a single figure (e.g., Fig3B1, Fig3B2)
#'
#' @param figure_id Base figure identifier (e.g., "Fig3B")
#' @param plots Named list of ggplot objects (names become suffixes)
#' @param data List of data frames for XLSX
#' @param methods_text Methods description
#' @param technical_legend Technical legend
#' @param output_dir Output directory
#' @param width Plot width
#' @param height Plot height
#' @param code_path Code file path
#'
#' @return Invisible NULL
save_multi_panel_figure <- function(
    figure_id,
    plots,
    data,
    methods_text,
    technical_legend = NULL,
    output_dir,
    width = 6,
    height = 5,
    code_path = NULL
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Save each plot
 for (panel_name in names(plots)) {
    panel_id <- paste0(figure_id, panel_name)
    png_path <- file.path(output_dir, paste0(panel_id, ".png"))
    pdf_path <- file.path(output_dir, paste0(panel_id, ".pdf"))

    ggsave(png_path, plots[[panel_name]], width = width, height = height, dpi = 300, bg = "white")
    ggsave(pdf_path, plots[[panel_name]], width = width, height = height)
    cat(sprintf("  Saved: %s\n", basename(png_path)))
  }

  # Save combined XLSX
  xlsx_path <- file.path(output_dir, paste0(figure_id, "_data.xlsx"))
  save_xlsx(data, xlsx_path)
  cat(sprintf("  Saved: %s\n", basename(xlsx_path)))

  # Save methods
  methods_path <- file.path(output_dir, paste0(figure_id, "_methods.txt"))
  save_methods_txt(
    figure_id = figure_id,
    methods_text = methods_text,
    technical_legend = technical_legend,
    code_path = code_path,
    output_path = methods_path
  )
  cat(sprintf("  Saved: %s\n", basename(methods_path)))

  invisible(NULL)
}

#' Save data to XLSX with multiple sheets
#'
#' @param data Named list of data frames (names become sheet names)
#' @param output_path Full path to output xlsx file
#'
#' @return Invisible NULL
save_xlsx <- function(data, output_path) {
  wb <- createWorkbook()

  for (sheet_name in names(data)) {
    # Clean sheet name (max 31 chars, no special chars)
    clean_name <- substr(gsub("[^A-Za-z0-9_]", "_", sheet_name), 1, 31)
    addWorksheet(wb, clean_name)
    writeData(wb, clean_name, data[[sheet_name]])
  }

  saveWorkbook(wb, output_path, overwrite = TRUE)
  invisible(NULL)
}

#' Save methods text file
#'
#' @param figure_id Figure identifier
#' @param methods_text Methods description
#' @param technical_legend Technical legend text
#' @param code_path Path to generating code
#' @param output_path Full path to output txt file
#'
#' @return Invisible NULL
save_methods_txt <- function(
    figure_id,
    methods_text,
    technical_legend = NULL,
    code_path = NULL,
    output_path
) {
  lines <- c(
    paste0("# ", figure_id, " - Methods and Technical Details"),
    paste0("# Generated: ", Sys.time()),
    "",
    "## Brief Methods Description",
    methods_text,
    ""
  )

  if (!is.null(technical_legend)) {
    lines <- c(lines,
      "## Technical Figure Legend",
      technical_legend,
      ""
    )
  }

  if (!is.null(code_path)) {
    lines <- c(lines,
      "## Generating Code",
      paste0("File: ", code_path),
      ""
    )
  }

  writeLines(lines, output_path)
  invisible(NULL)
}

#' Create run output directory with timestamp
#'
#' @param base_dir Base results directory
#' @param method_name Name of selection method
#'
#' @return Character string of full path
create_run_directory <- function(base_dir, method_name) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_dir <- file.path(base_dir, paste0(method_name, "_", timestamp))
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  return(run_dir)
}

#' Generate README for a pipeline run
#'
#' @param run_dir Run output directory
#' @param selection_info Selection criteria info list
#' @param n_figures Number of figures generated
#'
#' @return Invisible NULL
generate_run_readme <- function(run_dir, selection_info, n_figures) {
  readme_lines <- c(
    "# SMC Analysis Pipeline Run",
    "",
    paste0("**Generated:** ", Sys.time()),
    "",
    "## Selection Criteria",
    paste0("- **Method:** ", selection_info$method_name),
    paste0("- **Description:** ", selection_info$description),
    paste0("- **SMC cells:** ", length(selection_info$smc)),
    paste0("- **PDGFRa cells:** ", length(selection_info$pdgfra)),
    paste0("- **ICC cells:** ", length(selection_info$icc)),
    "",
    "## Figures Generated",
    paste0("Total figures: ", n_figures),
    "",
    "## Directory Contents",
    "Each figure directory contains:",
    "- `FigXX.png` - Publication-quality figure (300 DPI)",
    "- `FigXX.pdf` - Vector format",
    "- `FigXX_data.xlsx` - Source data",
    "- `FigXX_methods.txt` - Methods description and technical legend",
    ""
  )

  readme_path <- file.path(run_dir, "README.md")
  writeLines(readme_lines, readme_path)
  invisible(NULL)
}
