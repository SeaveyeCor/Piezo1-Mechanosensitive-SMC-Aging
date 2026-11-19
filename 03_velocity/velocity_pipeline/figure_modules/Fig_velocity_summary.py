"""
Velocity summary figures — 4 sub-panels of manuscript Figure 4F.

Manuscript mapping (Figure 4F, RNA velocity):
  Fig 4H   — Dynamical Contractile Velocity Index (4-gene)
  Fig S4C  — Myocd active transcription (stochastic)
  Fig 4I   — Piezo1 active upregulation (dynamical)
  Fig 4J   — Contractile→Synthetic transition flux

Pairwise statistics: 4 valid Mann-Whitney U comparisons per panel
(excluding YC×AK and YK×AC). Bonferroni-adjusted (m=4). ns brackets suppressed.

Outputs (per panel):
  output/<panel>.png  (300 DPI)
  output/<panel>.pdf  (vector)
  output/<panel>_data.csv  (group means + sd + sems + n)
  output/<panel>_stats.csv (4-comparison pairwise table with p_adjusted)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from velocity_pipeline.core.config import (
    H5AD_DYNAMICAL, H5AD_STOCHASTIC, OUTPUT_DIR,
    GROUP_ORDER, GROUP_LABELS, GROUP_COLORS,
    PAIRWISE_COMPARISONS, CONTRACTILE_DYNAMICAL,
)
from velocity_pipeline.core.metrics import (
    per_cell_fraction_positive, per_cell_gene_positive,
    transition_flow_to_target, group_summary, pairwise_stats,
)


# ── Publication style ────────────────────────────────────────────────────────
# Mirrors smc-analysis pipeline/utils/plotting_theme.R::theme_publication
# (base_size = 18; axis title 18; axis text 16; bold throughout).
def _apply_publication_style(ax, ylabel):
    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=18, fontweight='bold')
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_linewidth(0.8)
        ax.spines[spine].set_color('black')
    ax.tick_params(width=0.5, length=4, color='black')


def _bracket(ax, x1, x2, y, label, h_frac=0.025):
    yl = ax.get_ylim()
    h = (yl[1] - yl[0]) * h_frac
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.9, c='black')
    ax.text((x1 + x2) / 2, y + h * 1.1, label,
            ha='center', va='bottom', fontsize=14, fontweight='bold')


def _draw_panel(summary_df, stats_df, ylabel, comparisons, save_stem):
    """Render a single 4-group bar panel and write PNG + PDF + CSVs."""
    g2i = {g: i for i, g in enumerate(GROUP_ORDER)}
    ys = [summary_df.set_index('group').loc[g, 'pct_mean'] for g in GROUP_ORDER]
    es = [summary_df.set_index('group').loc[g, 'pct_sem']  for g in GROUP_ORDER]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(range(4), ys, yerr=es,
           color=[GROUP_COLORS[g] for g in GROUP_ORDER],
           edgecolor='black', linewidth=0.9, width=0.7,
           error_kw=dict(ecolor='black', elinewidth=1.4, capsize=6,
                          capthick=1.4, zorder=10))
    ax.set_xticks(range(4))
    ax.set_xticklabels([GROUP_LABELS[g] for g in GROUP_ORDER])
    _apply_publication_style(ax, ylabel)

    # Headroom for stacked brackets (4 valid comparisons after ns suppression)
    ymax = max(y + e for y, e in zip(ys, es) if not np.isnan(y))
    ax.set_ylim(0, ymax * 1.55)

    # Stack only significant brackets, sorted by span so wider ones nest above
    pair_rows = stats_df.copy()
    pair_rows = pair_rows[pair_rows['sig'] != 'ns']
    pair_rows['span'] = pair_rows.apply(lambda r: abs(g2i[r['b']] - g2i[r['a']]), axis=1)
    pair_rows['min_x'] = pair_rows.apply(lambda r: min(g2i[r['a']], g2i[r['b']]), axis=1)
    pair_rows = pair_rows.sort_values(['span', 'min_x'])
    bracket_y = ymax * 1.08
    step = ymax * 0.10
    for _, row in pair_rows.iterrows():
        x1 = g2i[row['a']]
        x2 = g2i[row['b']]
        _bracket(ax, min(x1, x2), max(x1, x2), bracket_y, row['sig'])
        bracket_y += step

    plt.tight_layout()
    # Each panel gets its own subfolder under OUTPUT_DIR (matches smc-analysis layout)
    panel_dir = OUTPUT_DIR / save_stem
    panel_dir.mkdir(parents=True, exist_ok=True)
    out_png = panel_dir / f'{save_stem}.png'
    out_pdf = panel_dir / f'{save_stem}.pdf'
    plt.savefig(str(out_png), dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(str(out_pdf), bbox_inches='tight', facecolor='white')
    plt.close(fig)

    summary_df.to_csv(panel_dir / f'{save_stem}_data.csv', index=False)
    stats_df.to_csv(panel_dir / f'{save_stem}_stats.csv', index=False)
    print(f'  Saved: {save_stem}/{{,.png,.pdf,_data.csv,_stats.csv}}')


def generate():
    print('Loading h5ads...')
    ad_d = sc.read_h5ad(H5AD_DYNAMICAL)
    ad_s = sc.read_h5ad(H5AD_STOCHASTIC)
    assert ad_d.n_obs == ad_s.n_obs, 'Cell count mismatch between h5ads'
    print(f'  dynamical:  {ad_d.n_obs} cells, {ad_d.n_vars} genes')
    print(f'  stochastic: {ad_s.n_obs} cells, {ad_s.n_vars} genes')

    groups_d = ad_d.obs['group'].astype(str).values
    groups_s = ad_s.obs['group'].astype(str).values
    phen_d = ad_d.obs['phenotype'].astype(str).values

    # ── Fig 4H: Dynamical Contractile Velocity Index ─────────────────────────
    print('\n[4H] Dynamical Contractile Velocity Index')
    cci, used_a = per_cell_fraction_positive(ad_d, CONTRACTILE_DYNAMICAL)
    summ_a  = group_summary(cci, groups_d)
    stats_a = pairwise_stats(cci, groups_d)
    _draw_panel(summ_a, stats_a,
                ylabel='Contractile Velocity Index (%)',
                comparisons=PAIRWISE_COMPARISONS,
                save_stem='Fig4H_contractile_velocity_index')

    # ── Fig S4C: Myocd alone (stochastic) ──────────────────────────────────
    print('\n[S4C] Myocd active transcription (stochastic)')
    myocd_pos = per_cell_gene_positive(ad_s, 'Myocd')
    summ_b  = group_summary(myocd_pos, groups_s)
    stats_b = pairwise_stats(myocd_pos, groups_s)
    _draw_panel(summ_b, stats_b,
                ylabel='% SMCs upregulating Myocd',
                comparisons=PAIRWISE_COMPARISONS,
                save_stem='FigS4C_myocd')

    # ── Fig 4I: Piezo1 upregulation (dynamical) ──────────────────────────
    print('\n[4I] Piezo1 active upregulation')
    p1_pos = per_cell_gene_positive(ad_d, 'Piezo1')
    summ_c  = group_summary(p1_pos, groups_d)
    stats_c = pairwise_stats(p1_pos, groups_d)
    _draw_panel(summ_c, stats_c,
                ylabel='% SMCs upregulating Piezo1',
                comparisons=PAIRWISE_COMPARISONS,
                save_stem='Fig4I_piezo1')

    # ── Fig 4J: contractile→synthetic flow (dynamical) ────────────────────
    print('\n[4J] Contractile cells: velocity flow → synthetic')
    contr = (phen_d == 'Contractile')
    synth = (phen_d == 'Synthetic')
    flow_full = np.full(ad_d.n_obs, np.nan)
    flows = transition_flow_to_target(ad_d, contr, synth)
    flow_full[contr] = flows
    summ_d  = group_summary(flow_full, groups_d)
    stats_d = pairwise_stats(flow_full, groups_d)
    _draw_panel(summ_d, stats_d,
                ylabel='% velocity flow → Synthetic',
                comparisons=PAIRWISE_COMPARISONS,
                save_stem='Fig4J_contractile_to_synthetic')

    # ── Velocity confidence per group ──────────────────────────────────────────
    # scVelo's velocity_confidence = per-cell correlation between velocity vector
    # and actual displacement in embedding. High = the velocity direction is
    # reliable for that cell.
    print('\n[Confidence] Per-group velocity confidence')
    conf_rows = []
    for label, ad, grp in [('dynamical', ad_d, groups_d),
                            ('stochastic', ad_s, groups_s)]:
        if 'velocity_confidence' not in ad.obs.columns:
            print(f'  WARNING: velocity_confidence not in {label} h5ad — skipping')
            continue
        vc = ad.obs['velocity_confidence'].values
        for g in GROUP_ORDER:
            mask = grp == g
            vals = vc[mask]
            valid = vals[~np.isnan(vals)]
            conf_rows.append({
                'model': label,
                'group': g,
                'n': len(valid),
                'mean_confidence': np.mean(valid) if len(valid) else np.nan,
                'sd_confidence': np.std(valid) if len(valid) else np.nan,
                'median_confidence': np.median(valid) if len(valid) else np.nan,
                'min_confidence': np.min(valid) if len(valid) else np.nan,
                'max_confidence': np.max(valid) if len(valid) else np.nan,
            })
            print(f'  {label} {g}: mean={np.mean(valid):.3f}, '
                  f'median={np.median(valid):.3f}, n={len(valid)}')
    conf_df = pd.DataFrame(conf_rows)
    conf_df.to_csv(OUTPUT_DIR / 'Fig_velocity_confidence_by_group.csv', index=False)
    print(f'  Saved: Fig_velocity_confidence_by_group.csv')

    # Combined source data + stats files for the README / supplementary table
    summ_a['panel'] = 'Fig4H_Contractile_Velocity_Index'
    summ_b['panel'] = 'FigS4C_Myocd'
    summ_c['panel'] = 'Fig4I_Piezo1_upregulation'
    summ_d['panel'] = 'Fig4J_Contractile_to_Synthetic'
    stats_a['panel'] = 'Fig4H_Contractile_Velocity_Index'
    stats_b['panel'] = 'FigS4C_Myocd'
    stats_c['panel'] = 'Fig4I_Piezo1_upregulation'
    stats_d['panel'] = 'Fig4J_Contractile_to_Synthetic'
    pd.concat([summ_a, summ_b, summ_c, summ_d], ignore_index=True).to_csv(
        OUTPUT_DIR / 'Fig_velocity_summary_groups.csv', index=False)
    pd.concat([stats_a, stats_b, stats_c, stats_d], ignore_index=True).to_csv(
        OUTPUT_DIR / 'Fig_velocity_summary_stats.csv', index=False)
    print(f'\nCombined CSVs:  Fig_velocity_summary_{{groups,stats}}.csv')

    return {'panels': ['A','B','C','D']}


if __name__ == '__main__':
    generate()
