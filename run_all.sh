#!/bin/bash
# Full analysis from the deposited Seurat object. Run from the repository root:
#   bash run_all.sh
# Outputs are written to pipeline_results/ and output/.
set -u
RS=${RSCRIPT:-Rscript}; PY=${PYTHON:-python3}
mkdir -p logs output/cell_level/panels
run() { echo "== $* $(date)"; "$@" > "logs/$(basename "$2").log" 2>&1 || echo "FAILED: $*"; }
run "$RS" 01_pipeline/run_complete_analysis.R
run "$RS" 02_source_data_postprocessing/run_postprocessing.R
for s in 01_build_cache.R 02_scores_singlecell.R 37_cluster_labels.R 46_conditional_positive_channels.R 47_conditional_figure.R \
         48_conditional_violins.R 52_channels_contractile_compartment.R 52b_channels_contractile_compartment_young.R \
         53_compartment_channels_figure.R 62_piezo1_pct_S3D.R 70_nfat_score_audit.R 71_score_panels.R 73_fig4F_cascade.R; do
  run "$RS" 04_cell_level/$s; done
run "$RS" 04_cell_level/75_score_superplot_inputs.R output/cell_level/superplots
run "$RS" 04_cell_level/02_superplots.R output/cell_level/superplots "3E,3F,3M,3N,4D_rec"
"$PY" 04_cell_level/53c_channel_detection_workbook.py output/cell_level/S3E_channel_detection.xlsx > logs/53c.log 2>&1 || echo "FAILED 53c"
"$PY" 04_cell_level/68_S3E_one_sheet.py output/cell_level/S3E_one_sheet.xlsx > logs/68.log 2>&1 || echo "FAILED 68"
"$PY" 04_cell_level/72_score_panels_workbook.py output/cell_level/score_panels_data.xlsx > logs/72.log 2>&1 || echo "FAILED 72"
"$PY" 04_cell_level/74_fig4F_cascade_workbook.py output/cell_level/Fig4F_cascade_data.xlsx > logs/74.log 2>&1 || echo "FAILED 74"
"$PY" 04_cell_level/64_S3B_expressing_cells_workbook.py output/cell_level/S3B_channels_with_age.xlsx > logs/64.log 2>&1 || echo "FAILED 64"
# RNA velocity (needs STARsolo output in data/velocity/scvelo_ready/, see 03_velocity/run_starsolo.sh)
if [ -d data/velocity/scvelo_ready ]; then
  run "$RS" 03_velocity/export_smc_metadata.R
  "$PY" 03_velocity/run_velocity.py > logs/vel_dyn.log 2>&1 || echo "FAILED velocity dynamical"
  "$PY" 03_velocity/run_velocity_stochastic.py > logs/vel_sto.log 2>&1 || echo "FAILED velocity stochastic"
  PYTHONPATH=03_velocity "$PY" -m velocity_pipeline.run_velocity_pipeline > logs/vel_panels.log 2>&1 || echo "FAILED velocity panels"
  "$PY" 03_velocity/build_velocity_source_data.py > logs/vel_sd.log 2>&1 || echo "FAILED velocity source data"
fi
echo "== DONE $(date)"
