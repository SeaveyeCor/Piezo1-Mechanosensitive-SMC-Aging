#!/usr/bin/env python3
"""
RNA Velocity Analysis — STOCHASTIC mode (hybrid companion to dynamical run)
============================================================================
Loads the existing depth-matched, phenotype-classified, moments-computed
adata from `velocity_smc.h5ad` (dynamical run) and re-computes velocity
using scVelo's stochastic model.

Why this exists:
  The dynamical model (run_velocity.py) filters out low-abundance
  transcription factors (Myocd, Srf, Mef2c, etc.) because their per-gene
  kinetic parameter fits collapse — too few cells with spliced reads to
  estimate alpha/beta/gamma reliably. Master SMC TFs are critical for any
  contractile-commitment story, so we run a parallel stochastic model on
  the same cell set, which estimates velocity from the spliced/unspliced
  moment correlation at steady state and is much more forgiving for sparse
  genes.

Hybrid usage:
  - dynamical h5ad: trajectory work (PHATE, latent_time, pseudotime cascade)
  - stochastic h5ad: gene-set scoring panels that need TF coverage
                    (Contractile Velocity Index, NFAT pathway flux, etc.)

Outputs:
  velocity_output/velocity_smc_stochastic.h5ad
"""
import os, sys, warnings
warnings.filterwarnings('ignore')

import numpy as np
import scanpy as sc
import scvelo as scv

DYNAMICAL_H5AD = 'data/velocity/velocity_output/velocity_smc.h5ad'
STOCHASTIC_H5AD = 'data/velocity/velocity_output/velocity_smc_stochastic.h5ad'

print('Loading dynamical h5ad...')
adata = sc.read_h5ad(DYNAMICAL_H5AD)
print(f'  cells: {adata.n_obs}  genes: {adata.n_vars}')

# Save dynamical-mode velocity to a separate layer key so we can compare per-gene
# coverage if needed, then overwrite the canonical 'velocity' layer with the
# stochastic estimate.
if 'velocity' in adata.layers:
    adata.layers['velocity_dynamical'] = adata.layers['velocity'].copy()
    print('  saved dynamical velocity to layer "velocity_dynamical"')

# Drop the dynamical-mode fit columns from var so they don't shadow the
# stochastic outputs (scVelo will write fresh columns).
drop_cols = [c for c in adata.var.columns if c.startswith('fit_') or
             c in ('velocity_genes', 'velocity_score', 'velocity_gamma',
                   'velocity_qreg_ratio', 'velocity_r2')]
if drop_cols:
    print(f'  dropping {len(drop_cols)} stale dynamical-fit var columns')
    adata.var.drop(columns=drop_cols, inplace=True)

# Stochastic-mode velocity. Requires Ms/Mu moments (already computed in the
# dynamical run) — does NOT require recover_dynamics.
print('Running scv.tl.velocity(mode="stochastic")...')
scv.tl.velocity(adata, mode='stochastic')

# Recompute velocity graph for downstream embedding/transition uses
print('Computing velocity graph...')
scv.tl.velocity_graph(adata, n_jobs=4)

# Diagnostics — focus on the genes we care about
print('\n=== Per-gene stochastic velocity coverage (target genes) ===')
target = ['Myh11','Acta2','Tagln','Cnn1','Mylk','Myocd','Smtn',
          'Lmod1','Actg2','Cald1','Tpm1','Tpm2','Des',
          'Srf','Mef2c','Klf4','Vim','Col1a1','Fn1','Spp1',
          'Mmp2','Mmp9','Timp1','Wnt5a',
          'Bgn','Col6a1','Dcn','Sparc','Timp2',
          'Nfatc1','Nfatc2','Nfatc3','Nfatc4','Nfat5','Calm1','Ppp3ca',
          'Piezo1']
V = adata.layers['velocity']
for g in target:
    if g in adata.var_names:
        j = adata.var_names.get_loc(g)
        v = np.asarray(V[:, j]).flatten()
        n_nan = int(np.isnan(v).sum())
        n_pos = int((v > 0).sum())
        n_neg = int((v < 0).sum())
        vg = bool(adata.var.loc[g, 'velocity_genes']) if 'velocity_genes' in adata.var.columns else None
        print(f'  {g:8s}  velocity_genes={vg}  non-NaN={adata.n_obs - n_nan}/{adata.n_obs}  pos={n_pos} neg={n_neg}')
    else:
        print(f'  {g:8s}  ABSENT from var_names')

# Save
print(f'\nSaving stochastic h5ad → {STOCHASTIC_H5AD}')
adata.write(STOCHASTIC_H5AD)
print('Done.')
