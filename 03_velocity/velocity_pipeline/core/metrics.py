"""
Velocity pipeline metrics — per-cell scores and group statistics.
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .config import GROUP_ORDER, PAIRWISE_COMPARISONS


def get_velocity_vector(adata, gene):
    """Per-cell velocity values for a single gene (NaN for non-velocity genes)."""
    if gene not in adata.var_names:
        raise KeyError(f'{gene} not in adata.var_names')
    j = adata.var_names.get_loc(gene)
    v = adata.layers['velocity'][:, j]
    if hasattr(v, 'toarray'):
        v = v.toarray().flatten()
    return np.asarray(v).flatten()


def per_cell_fraction_positive(adata, genes):
    """For each cell, fraction of `genes` that have positive velocity (NaN-aware).
    Genes absent from var_names are silently skipped (warned)."""
    avail = [g for g in genes if g in adata.var_names]
    if len(avail) == 0:
        raise ValueError('No genes available for fraction-positive metric')
    if len(avail) < len(genes):
        missing = set(genes) - set(avail)
        print(f'  WARNING: skipping {len(missing)} absent genes: {sorted(missing)}')
    V = np.column_stack([get_velocity_vector(adata, g) for g in avail])
    n_avail = np.sum(~np.isnan(V), axis=1)
    n_pos = np.sum(V > 0, axis=1)
    return np.where(n_avail > 0, n_pos / np.maximum(n_avail, 1), np.nan), avail


def per_cell_gene_positive(adata, gene):
    """Binary per-cell indicator (1 if velocity > 0, else 0; NaN preserved)."""
    v = get_velocity_vector(adata, gene)
    out = np.where(np.isnan(v), np.nan, (v > 0).astype(float))
    return out


def transition_flow_to_target(adata, src_mask, tgt_mask):
    """For each source cell, fraction of velocity-graph transition probability
    mass that lands on target-mask cells. Returns one value per source cell."""
    vgraph = adata.uns['velocity_graph']
    trans = vgraph[src_mask].toarray()
    f_t = trans[:, tgt_mask].sum(axis=1)
    f_o = trans[:, ~tgt_mask].sum(axis=1)
    tot = f_t + f_o
    valid = tot > 0
    out = np.full(src_mask.sum(), np.nan)
    out[valid] = f_t[valid] / tot[valid]
    return out


def group_summary(per_cell_values, groups, order=None):
    """Mean ± SD ± SEM and n by group, for a per-cell scalar metric."""
    order = order or GROUP_ORDER
    df = pd.DataFrame({'value': per_cell_values, 'group': groups})
    rows = []
    for g in order:
        s = df[df['group'] == g]['value'].dropna()
        rows.append({
            'group': g,
            'n': len(s),
            'mean': s.mean() if len(s) else np.nan,
            'sd':   s.std()  if len(s) else np.nan,
            'sem':  s.sem()  if len(s) else np.nan,
            'pct_mean': 100 * s.mean() if len(s) else np.nan,
            'pct_sd':   100 * s.std()  if len(s) else np.nan,
            'pct_sem':  100 * s.sem()  if len(s) else np.nan,
        })
    return pd.DataFrame(rows)


def pairwise_stats(per_cell_values, groups, comparisons=None, adjust='bonferroni'):
    """Two-sided Mann-Whitney U for each (a, b, label) in `comparisons`.

    Defaults to PAIRWISE_COMPARISONS (excludes the two diagonal cross-group
    pairs that conflate age + genotype). P-values are Bonferroni-adjusted by
    default (m = number of comparisons), matching the smc-analysis pipeline
    helper at pipeline/utils/statistics.R::pairwise_group_tests.

    Returns a data frame with columns:
        comparison, a, b, n_a, n_b, mean_a, mean_b, U, p, p_adjusted, sig
    The `sig` column reflects `p_adjusted`, so display brackets are conservative.
    """
    comparisons = comparisons or PAIRWISE_COMPARISONS
    df = pd.DataFrame({'value': per_cell_values, 'group': groups})
    rows = []
    for a, b, label in comparisons:
        va = df[df['group'] == a]['value'].dropna()
        vb = df[df['group'] == b]['value'].dropna()
        if len(va) < 2 or len(vb) < 2:
            rows.append({
                'comparison': label, 'a': a, 'b': b,
                'n_a': len(va), 'n_b': len(vb),
                'U': np.nan, 'p': np.nan,
                'mean_a': va.mean() if len(va) else np.nan,
                'mean_b': vb.mean() if len(vb) else np.nan,
            })
            continue
        u, p = mannwhitneyu(vb, va, alternative='two-sided')
        rows.append({
            'comparison': label, 'a': a, 'b': b,
            'n_a': len(va), 'n_b': len(vb),
            'U': u, 'p': p,
            'mean_a': va.mean(), 'mean_b': vb.mean(),
        })
    out = pd.DataFrame(rows)

    # Bonferroni adjust (m = number of finite p-values)
    if adjust == 'bonferroni':
        m = int(out['p'].notna().sum())
        out['p_adjusted'] = (out['p'] * m).clip(upper=1.0)
    elif adjust in (None, 'none'):
        out['p_adjusted'] = out['p']
    else:
        raise ValueError(f'Unknown adjust method: {adjust}')

    out['sig'] = out['p_adjusted'].apply(sig_label)
    return out


def sig_label(p):
    if p is None or np.isnan(p):
        return ''
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

