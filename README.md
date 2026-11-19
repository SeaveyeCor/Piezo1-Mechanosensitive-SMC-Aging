# Maladaptive Piezo1 mechanotransduction drives smooth muscle aging in the gut: analysis code

Code for the single-cell RNA-sequencing and RNA velocity analyses in Joshi, Seavey, *et al.*, "Maladaptive Piezo1
Mechanotransduction Drives Smooth Muscle Aging in the Gut" (bioRxiv, https://doi.org/10.64898/2026.09.02.748183).

## Data
- Processed Seurat object for the six mice analysed (`seurat_paper_6samples.rds`, 15,514 cells, with the stored PCA/UMAP
  used for Fig 3A and Fig S3A) and a cell-level table linking every cell to the figures (`cell_annotations_paper.csv`):
  Zenodo, https://doi.org/10.5281/zenodo.24381866. Place the object at `data/seurat_paper_6samples.rds` before running.
- Raw FASTQ files and Cell Ranger filtered feature-barcode matrices: GEO, accession number upcoming.
- Bundled in `data/`: `annotation_template_REVIEWED.csv` (cluster to cell-type annotation), `doublet_cell_ids.txt`
  (scDblFinder doublets removed from the SMC set), `table_S*.csv` (descriptive text for Tables S0-S2 and the S4 gene
  list; the numbers in Table S1 are computed from the object when the tables are built).

Single-cell libraries (10x Chromium Single Cell 3' v3, mouse colonic muscularis, Cell Ranger with refdata-gex-mm10-2020-A):

| sample | group | library | Cell Ranger |
|---|---|---|---|
| WT_Young_1 | Young CTRL | GEX_WT | 7.2.0 |
| WT_Young_2 | Young CTRL | GEX_WT-2_Jan25 | 7.2.0 |
| WT_Aged_1 | Aged CTRL | GEX_Aged_WT_May21v | 9.0.0 |
| WT_Aged_2 | Aged CTRL | GEX_Aged-3_Apr09 | 7.2.0 |
| KO_Young_1 | Young P1KO | GEX_Piezo1-Young_May09 | 7.2.0 |
| KO_Aged | Aged P1KO | GEX_Piezo-Aged_May15v | 9.0.0 |

## Layout
All scripts run from the repository root, with inputs under `data/`. Outputs: `pipeline_results/` (one folder per figure module, internal IDs), `output/source_data/` (source data by paper panel), `output/tables/`, `output/figures/` and `output/cell_level/`.

| folder | contents |
|---|---|
| `00_preprocessing/` | Sample identity checks run on every load (`check_sample_identity.R`: no cell shared between samples, each sample matched to its own Cell Ranger file, stored UMAP), cell-type annotation from clusters (`11_`), SMC phenotype scores (`13_`) and scDblFinder doublet detection on SMCs (`60d_doublet_detection.R`). QC, normalization and clustering are described in the Methods. |
| `pipeline/` | Analysis modules (core steps, cell selection, figure modules, utilities) called by the entry scripts. `pipeline/core/gene_sets.R` holds every gene set. |
| `01_pipeline/` | Main analysis, entry point `run_complete_analysis.R`. Steps: remove non-muscularis cells; depth check (no library falls outside 0.5-2x the median of per-sample median UMIs, so none is filtered); select SMCs (clusters 1/5/15/16 minus doublets); classify phenotypes; figure modules; Tables S0-S12. Module IDs (`Fig3E`, `FigS3D`, ...) are internal names; the table below maps them to manuscript panels. |
| `02_source_data_postprocessing/` | Source-data workbooks named by paper panel (`output/source_data/Fig3A_data.xlsx`, ...), with SD, regression and cell-level sheets added, and the Fig S4B plot; run with `run_postprocessing.R`. |
| `03_velocity/` | RNA velocity: STARsolo alignment (`run_starsolo.sh`), SMC metadata export (`export_smc_metadata.R`), scVelo dynamical and stochastic models (`run_velocity.py`, `run_velocity_stochastic.py`), panels (`python -m velocity_pipeline.run_velocity_pipeline`) and source data (`build_velocity_source_data.py`). Some libraries were aligned against a STAR index built with `--genomeSAsparseD 2` (lower memory). |
| `04_cell_level/` | Per-cell analyses and panels: per-cell cache and scores (`01_`, `02_scores_singlecell.R`), cluster labels (`37_`), score panels and superplot inputs (`71_`, `72_`, `75_`), plots (`02_superplots.R`), mechanosensitive channel analyses (`46_`-`53c_`, `64_`, `68_`), Piezo1+ SMC fractions (`62_`), NFAT score robustness (`70_`), pseudotime module dynamics (`73_`-`74_`). |

## Figure panels and tables

Every figure script starts with a header giving the paper panel, what it plots, the cells used, the groups compared and
the statistical test. Internal module IDs (for example `Fig3E`) differ from the paper panel names; the header and the
table below give the mapping.

| panel | script |
|---|---|
| 3A UMAP of contractile and synthetic SMCs | `pipeline/figure_modules/Fig3E_umap_phenotype.R` |
| 3B phenotype proportions | `.../Fig3K_phenotype_proportions.R` |
| 3C Piezo1 vs contractile score | `.../Fig3F_piezo1_contractile_correlation.R` |
| 3D Piezo1 vs synthetic score | `04_cell_level/71_score_panels.R` |
| 3E / 3F contractile / synthetic score | `04_cell_level/75_score_superplot_inputs.R` + `02_superplots.R` |
| 3G / 3H DE of contractile / synthetic genes | `.../Fig3I_contractile_rescue.R`, `.../Fig3J_synthetic_rescue.R` |
| 3M ECM score | `04_cell_level/75_` + `02_superplots.R` (score: `.../supplementary/FigS3K_ecm_organization.R`) |
| 3N cytoskeletal remodeling score | `04_cell_level/75_` + `02_superplots.R` (score: `.../supplementary/FigS4H_cyto_remodeling.R`) |
| 4A-4C Piezo1 / NFAT correlations | `04_cell_level/71_score_panels.R` |
| 4D NFAT-associated remodeling score | `04_cell_level/71_score_panels.R` + `75_` + `02_superplots.R`; robustness `70_nfat_score_audit.R` |
| 4F gene-module dynamics on pseudotime | `04_cell_level/73_fig4F_cascade.R` + `74_fig4F_cascade_workbook.py` |
| 4G PHATE | `.../Fig4H_phate_trajectory.R` (PHATE implemented in R: kNN, Gaussian kernel, diffusion t = 5, MDS) |
| 4H-4J RNA velocity | `03_velocity/velocity_pipeline/figure_modules/Fig_velocity_summary.py`; source data `03_velocity/build_velocity_source_data.py` |
| S3A UMAP of annotated cell populations (Cl6, Cl13 labelled separately) | `.../supplementary/FigS3D_umap_celltypes.R` |
| S3B mechanosensitive channels with age | `04_cell_level/46_`-`48_`, `64_` |
| S3C / S3D Piezo1 change and % Piezo1+ SMCs | `.../supplementary/FigS3F_piezo1_log2fc.R`; `04_cell_level/62_piezo1_pct_S3D.R` |
| S3E % contractile SMCs detecting each channel | `04_cell_level/37_`, `52_`, `52b_`, `53_`, `53c_`, `68_` |
| S3F GO enrichment (GSEA) | `.../supplementary/FigS3I_go_enrichment.R`, `02_.../source_data_FigS3F_gsea.R` |
| S3J ECM gene heatmap | `.../supplementary/FigS3L_ecm_heatmap.R` (module `FigS4A`) |
| S4A Piezo1 vs representative genes | `.../supplementary/FigS3H_gene_correlations.R` |
| S4B GO programs, Piezo1+ vs Piezo1- CTRL synthetic SMCs | `.../Fig4J_p1pos_vs_p1neg_synthetic.R` + `02_.../plot_FigS4B_go_programs.R` |
| S4C % SMCs with positive Myocd velocity | `03_velocity/velocity_pipeline/figure_modules/Fig_velocity_summary.py` |
| Tables S0-S7 | `.../supplementary/SupplementalTables.R` with `PaperTables_S0_S7.R` (`Supplemental_Tables/Tables_S0_S7.xlsx`) |
| Tables S8-S12 (genome-wide DE, GSEA) | `01_pipeline/run_complete_analysis.R` (phase 2) |

The human analyses (Fig 5) are not part of this repository.

## Environment
`environment/` holds the R 4.4.2 `sessionInfo()` and the Python 3.12 package set for velocity
(`requirements_velocity.txt`, install with `pip install --no-deps`). STAR 2.7.11b (velocity alignment, GRCm39 / Ensembl 111).

## Running
From the repository root:

    bash run_all.sh

`run_all.sh` runs, in order: `01_pipeline/run_complete_analysis.R`, `02_source_data_postprocessing/run_postprocessing.R`,
the `04_cell_level/` scripts, and RNA velocity (`03_velocity/`, when STARsolo output is in `data/velocity/scvelo_ready/`).
Set `RSCRIPT` / `PYTHON` to point at specific interpreters if needed.

## Licence
MIT (see `LICENSE`). Copyright (c) 2026 the authors and the Mayo Clinic Enteric Neurosciences Program.
