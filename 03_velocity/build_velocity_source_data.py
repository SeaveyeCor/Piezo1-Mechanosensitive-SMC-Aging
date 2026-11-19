"""
Build a master Excel workbook with raw/summary data for all velocity figure panels.

Outputs:
  7_source_data/Velocity_Panels_source_data.xlsx

Worksheets (intuitive names):
  README                          — panel descriptions, gene sets, cell counts
  A_CCI_summary                   — Panel A summary stats (Contractile Velocity Index)
  A_CCI_stats                     — Panel A pairwise Mann-Whitney U tests
  A_CCI_cell_level                — Panel A per-cell CCI values
  B_Myocd_summary                 — Panel B summary stats (Myocd upregulation)
  B_Myocd_stats                   — Panel B pairwise tests
  B_Myocd_cell_level              — Panel B per-cell Myocd velocity indicator
  C_Piezo1_summary                — Panel C summary stats (Piezo1 upregulation)
  C_Piezo1_stats                  — Panel C pairwise tests
  C_Piezo1_cell_level             — Panel C per-cell Piezo1 velocity indicator
  D_Flow_summary                  — Panel D summary stats (contractile->synthetic)
  D_Flow_stats                    — Panel D pairwise tests
  D_Flow_cell_level               — Panel D per-cell transition flow (contractile only)
"""
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# -- Paths --
ROOT = Path('data/velocity')
H5AD_DYN = ROOT / 'velocity_output' / 'velocity_smc.h5ad'
H5AD_STO = ROOT / 'velocity_output' / 'velocity_smc_stochastic.h5ad'
PIPE_OUT = ROOT / 'velocity_pipeline' / 'output'
DEST = Path('output/source_data/velocity'); DEST.mkdir(parents=True, exist_ok=True)
OUT_FILE = DEST / 'Velocity_Panels_source_data.xlsx'

GROUP_ORDER = ['Young_Ctrl', 'Aged_Ctrl', 'Young_KO', 'Aged_KO']
CONTRACTILE_DYNAMICAL = ['Myh11', 'Cnn1', 'Lmod1', 'Actg2']


def get_vel(adata, gene):
    j = adata.var_names.get_loc(gene)
    v = adata.layers['velocity'][:, j]
    if hasattr(v, 'toarray'):
        v = v.toarray().flatten()
    return np.asarray(v).flatten()


def main():
    print('Loading h5ads...')
    ad_d = sc.read_h5ad(H5AD_DYN)
    ad_s = sc.read_h5ad(H5AD_STO)
    print(f'  Dynamical: {ad_d.n_obs} cells | Stochastic: {ad_s.n_obs} cells')

    groups_d = ad_d.obs['group'].astype(str).values
    groups_s = ad_s.obs['group'].astype(str).values
    phenotype_d = ad_d.obs['phenotype'].astype(str).values
    barcodes_d = ad_d.obs_names.values
    barcodes_s = ad_s.obs_names.values
    sample_d = ad_d.obs['sample'].astype(str).values if 'sample' in ad_d.obs.columns else ['NA'] * ad_d.n_obs
    sample_s = ad_s.obs['sample'].astype(str).values if 'sample' in ad_s.obs.columns else ['NA'] * ad_s.n_obs

    # -- Read existing summary/stats CSVs --
    def read_csv(panel_dir, suffix):
        return pd.read_csv(PIPE_OUT / panel_dir / f'{panel_dir}_{suffix}.csv')

    summ_a = read_csv('Fig4H_contractile_velocity_index', 'data')
    stat_a = read_csv('Fig4H_contractile_velocity_index', 'stats')
    summ_b = read_csv('FigS4C_myocd', 'data')
    stat_b = read_csv('FigS4C_myocd', 'stats')
    summ_c = read_csv('Fig4I_piezo1', 'data')
    stat_c = read_csv('Fig4I_piezo1', 'stats')
    summ_d = read_csv('Fig4J_contractile_to_synthetic', 'data')
    stat_d = read_csv('Fig4J_contractile_to_synthetic', 'stats')

    # -- Panel A: cell-level CCI --
    print('Building cell-level data for Panel A (CCI)...')
    vels_a = {}
    for g in CONTRACTILE_DYNAMICAL:
        vels_a[g] = get_vel(ad_d, g)
    frac_pos = np.zeros(ad_d.n_obs)
    for i in range(ad_d.n_obs):
        vals = [vels_a[g][i] for g in CONTRACTILE_DYNAMICAL]
        n_valid = sum(1 for v in vals if not np.isnan(v))
        n_pos = sum(1 for v in vals if not np.isnan(v) and v > 0)
        frac_pos[i] = n_pos / max(n_valid, 1) if n_valid > 0 else np.nan

    cell_a = pd.DataFrame({
        'barcode': barcodes_d,
        'sample': sample_d,
        'group': groups_d,
        'phenotype': phenotype_d,
    })
    for g in CONTRACTILE_DYNAMICAL:
        cell_a[f'velocity_{g}'] = vels_a[g]
    cell_a['CCI_fraction_positive'] = frac_pos
    cell_a['CCI_pct'] = frac_pos * 100

    # -- Panel B: cell-level Myocd (stochastic) --
    print('Building cell-level data for Panel B (Myocd)...')
    v_myocd = get_vel(ad_s, 'Myocd')
    cell_b = pd.DataFrame({
        'barcode': barcodes_s,
        'sample': sample_s,
        'group': groups_s,
        'Myocd_velocity': v_myocd,
        'Myocd_upregulating': np.where(np.isnan(v_myocd), np.nan, (v_myocd > 0).astype(float)),
    })

    # -- Panel C: cell-level Piezo1 (dynamical) --
    print('Building cell-level data for Panel C (Piezo1)...')
    v_piezo1 = get_vel(ad_d, 'Piezo1')
    cell_c = pd.DataFrame({
        'barcode': barcodes_d,
        'sample': sample_d,
        'group': groups_d,
        'phenotype': phenotype_d,
        'Piezo1_velocity': v_piezo1,
        'Piezo1_upregulating': np.where(np.isnan(v_piezo1), np.nan, (v_piezo1 > 0).astype(float)),
    })

    # -- Panel D: cell-level contractile->synthetic flow (dynamical) --
    print('Building cell-level data for Panel D (Flow)...')
    contr_mask = (phenotype_d == 'Contractile')
    synth_mask = (phenotype_d == 'Synthetic')
    vgraph = ad_d.uns['velocity_graph']
    trans = vgraph[contr_mask].toarray()
    f_t = trans[:, synth_mask].sum(axis=1)
    f_o = trans[:, ~synth_mask].sum(axis=1)
    tot = f_t + f_o
    valid = tot > 0
    flow = np.full(contr_mask.sum(), np.nan)
    flow[valid] = f_t[valid] / tot[valid]

    cell_d = pd.DataFrame({
        'barcode': barcodes_d[contr_mask],
        'sample': sample_d[contr_mask],
        'group': groups_d[contr_mask],
        'transition_flow_to_synthetic': flow,
        'transition_flow_pct': flow * 100,
    })

    # -- README sheet --
    readme_rows = [
        ['RNA velocity source data: Fig 4H, 4I, 4J and Fig S4C'],
        [''],
        ['This workbook contains raw data for all 4 velocity figure panels.'],
        ['Each panel has 3 sheets: _summary (group means), _stats (pairwise tests), _cell_level (per-cell values).'],
        [''],
        ['Panel A = Fig 4H: Dynamical Contractile Velocity Index (CCI)'],
        [f'  Genes: {", ".join(CONTRACTILE_DYNAMICAL)}'],
        ['  Model: scVelo dynamical mode'],
        [f'  Metric: fraction of genes with positive velocity per cell (0-1 scale; pct = x100)'],
        [f'  N cells: {ad_d.n_obs} ({", ".join(f"{g}={int((groups_d==g).sum())}" for g in GROUP_ORDER)})'],
        [''],
        ['Panel B = Fig S4C: Myocd Active Transcription'],
        ['  Gene: Myocd (sparse TF requires stochastic mode for velocity recovery)'],
        ['  Model: scVelo stochastic mode'],
        ['  Metric: binary indicator (1 = velocity > 0 = actively upregulating)'],
        [f'  N cells: {ad_s.n_obs}'],
        [''],
        ['Panel C = Fig 4I: Piezo1 Active Upregulation'],
        ['  Gene: Piezo1'],
        ['  Model: scVelo dynamical mode'],
        ['  Metric: binary indicator (1 = velocity > 0 = actively upregulating)'],
        [f'  N cells: {ad_d.n_obs}'],
        [''],
        ['Panel D = Fig 4J: Contractile-to-Synthetic Transition Flow'],
        ['  Source: contractile-phenotype cells only'],
        ['  Model: scVelo dynamical mode (velocity graph)'],
        ['  Metric: fraction of velocity-graph transition probability landing on synthetic cells'],
        [f'  N contractile cells: {int(contr_mask.sum())} ({", ".join(f"{g}={int(((groups_d==g) & contr_mask).sum())}" for g in GROUP_ORDER)})'],
        [''],
        ['Statistics: Mann-Whitney U, two-sided, Bonferroni-adjusted (m=4)'],
        ['  4 comparisons per panel (excludes diagonal age+genotype confounds)'],
        ['  sig column: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant'],
    ]
    readme_df = pd.DataFrame(readme_rows, columns=['Description'])

    # -- Write Excel --
    print(f'Writing {OUT_FILE}...')
    with pd.ExcelWriter(str(OUT_FILE), engine='openpyxl') as writer:
        readme_df.to_excel(writer, sheet_name='README', index=False)
        summ_a.to_excel(writer, sheet_name='A_CCI_summary', index=False)
        stat_a.to_excel(writer, sheet_name='A_CCI_stats', index=False)
        cell_a.to_excel(writer, sheet_name='A_CCI_cell_level', index=False)
        summ_b.to_excel(writer, sheet_name='B_Myocd_summary', index=False)
        stat_b.to_excel(writer, sheet_name='B_Myocd_stats', index=False)
        cell_b.to_excel(writer, sheet_name='B_Myocd_cell_level', index=False)
        summ_c.to_excel(writer, sheet_name='C_Piezo1_summary', index=False)
        stat_c.to_excel(writer, sheet_name='C_Piezo1_stats', index=False)
        cell_c.to_excel(writer, sheet_name='C_Piezo1_cell_level', index=False)
        summ_d.to_excel(writer, sheet_name='D_Flow_summary', index=False)
        stat_d.to_excel(writer, sheet_name='D_Flow_stats', index=False)
        cell_d.to_excel(writer, sheet_name='D_Flow_cell_level', index=False)

    print(f'Done! 13 sheets written to:\n  {OUT_FILE}')


if __name__ == '__main__':
    main()
