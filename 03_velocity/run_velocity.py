#!/usr/bin/env python3
"""
RNA Velocity Analysis — 7 smc-analysis samples
================================================
Loads STARsolo outputs, maps barcodes to Seurat identities,
selects the bona fide SMCs exported by export_smc_metadata.R,
and runs scVelo RNA velocity.

Folder → Seurat sample mapping:
  GEX_WT         → WT_Young_1
  GEX_WT-2_Jan25 → WT_Young_2
  Aged_WT        → WT_Aged_1
  Aged_3         → WT_Aged_2
  Piezo1_Young   → KO_Young_1
  Piezo_Aged         → KO_Aged
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── Configuration ─────────────────────────────────────────────────────────────
SCVELO_DIR = "data/velocity/scvelo_ready"
OUTPUT_DIR = "data/velocity/velocity_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Folder → Seurat sample mapping
FOLDER_MAP = {
    "GEX_WT": "WT_Young_1",
    "GEX_WT-2_Jan25": "WT_Young_2",
    "Aged_WT": "WT_Aged_1",
    "Aged_3": "WT_Aged_2",
    "Piezo1_Young": "KO_Young_1",
    "Piezo_Aged": "KO_Aged",
}

# Gene panels (from smc-analysis/pipeline/core/gene_sets.R)
CONTRACTILE_GENES = ["Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn"]
SYNTHETIC_GENES = ["Col1a1", "Fn1", "Spp1", "Klf4", "Vim"]

# Colors (lab palette)
GROUP_COLORS = {
    "Young_Ctrl": "#4DAF4A",
    "Aged_Ctrl": "#E41A1C",
    "Young_KO": "#7F7F7F",
    "Aged_KO": "#ff7fb2",
}
PHENOTYPE_COLORS = {
    "Contractile": "#E41A1C",
    "Synthetic": "#377EB8",
}

sc.settings.verbosity = 2
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo", dpi=150, dpi_save=300)
plt.rcParams.update({'figure.max_open_warning': 0})


# ── Load metadata from Seurat ────────────────────────────────────────────────
print("Loading Seurat metadata...")
meta = pd.read_csv(os.path.join(SCVELO_DIR, "velocity_metadata_corrected.csv"))
vel_meta = meta[meta["is_bonafide"] == True].copy()
print(f"  {len(meta)} cells, {len(vel_meta)} velocity SMCs")
print(f"  Samples: {vel_meta['sample'].unique().tolist()}")



# ── Load STARsolo RAW data and filter to Seurat barcodes ─────────────────────
def load_starsolo_raw(folder_name, keep_barcodes):
    """Load Gene + Velocyto RAW matrices, subset to keep_barcodes."""
    base = os.path.join(SCVELO_DIR, folder_name, "Solo.out")
    gene_dir = os.path.join(base, "Gene", "raw")
    vel_dir = os.path.join(base, "Velocyto", "raw")

    with open(os.path.join(gene_dir, "features.tsv")) as f:
        features = [line.strip().split("\t") for line in f]
    gene_names = [g[1] if len(g) > 1 else g[0] for g in features]

    with open(os.path.join(gene_dir, "barcodes.tsv")) as f:
        all_barcodes = [line.strip() for line in f]

    # Find indices of barcodes to keep
    bc_to_idx = {bc: i for i, bc in enumerate(all_barcodes)}
    keep_indices = [bc_to_idx[bc] for bc in keep_barcodes if bc in bc_to_idx]
    kept_barcodes = [all_barcodes[i] for i in keep_indices]

    if not keep_indices:
        return None

    # Load matrices and subset rows (cells)
    mat = scipy.io.mmread(os.path.join(gene_dir, "matrix.mtx")).T.tocsr()[keep_indices]
    spliced = scipy.io.mmread(os.path.join(vel_dir, "spliced.mtx")).T.tocsr()[keep_indices]
    unspliced = scipy.io.mmread(os.path.join(vel_dir, "unspliced.mtx")).T.tocsr()[keep_indices]
    ambiguous = scipy.io.mmread(os.path.join(vel_dir, "ambiguous.mtx")).T.tocsr()[keep_indices]

    adata = sc.AnnData(
        X=mat,
        obs=pd.DataFrame(index=kept_barcodes),
        var=pd.DataFrame(index=gene_names),
    )
    adata.layers["spliced"] = spliced
    adata.layers["unspliced"] = unspliced
    adata.layers["ambiguous"] = ambiguous
    adata.obs["folder"] = folder_name
    adata.var_names_make_unique()

    return adata


# ── Build per-folder barcode lists from velocity metadata ────────────────────
# Each folder is one library: barcode -> sample mapping is 1:1

FOLDER_LOADS = []  # (folder_name, {barcode: sample}, unique_prefix)

for folder_name, sample in FOLDER_MAP.items():
    if sample is not None:
        bcs = vel_meta[vel_meta["sample"] == sample]["barcode"].tolist()
        FOLDER_LOADS.append((folder_name, {bc: sample for bc in bcs}, sample))


# ── Load all folders ─────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("LOADING STARSOLO DATA (RAW + SEURAT BARCODE FILTER)")
print("=" * 70)

adata_parts = []

for folder_name, bc_map, prefix in FOLDER_LOADS:
    if not bc_map:
        print(f"  SKIP {folder_name}/{prefix}: no matching barcodes")
        continue

    keep_bcs = list(bc_map.keys())
    print(f"  Loading {folder_name} → {prefix} ({len(keep_bcs)} target barcodes)...")
    adata = load_starsolo_raw(folder_name, keep_bcs)

    if adata is None or adata.n_obs == 0:
        print(f"    WARNING: no cells loaded")
        continue

    # Assign sample identity
    sample_name = list(set(bc_map.values()))[0]
    adata.obs["sample"] = sample_name
    # Make cell IDs unique by prefixing with sample
    adata.obs["barcode"] = adata.obs.index.values
    adata.obs.index = [f"{prefix}_{bc}" for bc in adata.obs.index]
    print(f"    Loaded {adata.n_obs} cells")
    adata_parts.append(adata)

# Merge
print("\nMerging...")
adata = sc.concat(adata_parts, join="outer", fill_value=0)
del adata_parts
print(f"Merged: {adata.n_obs} cells, {adata.n_vars} genes")

# ── Annotate from metadata ───────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ANNOTATING VELOCITY SMCs")
print("=" * 70)

# Build sample+barcode lookup for metadata
vel_meta_keyed = vel_meta.set_index(["sample", "barcode"])

groups = []
clusters = []
for _, row in adata.obs.iterrows():
    bc = row["barcode"] if "barcode" in row.index else row.name
    key = (row["sample"], bc)
    if key in vel_meta_keyed.index:
        r = vel_meta_keyed.loc[key]
        groups.append(r["group"])
        clusters.append(str(r["cluster"]) if pd.notna(r["cluster"]) else "NA")
    else:
        groups.append("unknown")
        clusters.append("NA")

adata.obs["group"] = groups
adata.obs["cluster"] = clusters

print(f"Total velocity SMCs: {adata.n_obs}")

print("\nPer-sample counts:")
for s in sorted(adata.obs["sample"].unique()):
    n = (adata.obs["sample"] == s).sum()
    grp = adata.obs.loc[adata.obs["sample"] == s, "group"].iloc[0]
    print(f"  {s} ({grp}): {n}")

print(f"\nPer-group counts:")
for g in ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]:
    n = (adata.obs["group"] == g).sum()
    print(f"  {g}: {n}")


# ── Phenotype classification ─────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PHENOTYPE CLASSIFICATION")
print("=" * 70)

def score_genes(adata_obj, genes, layer="spliced"):
    """Mean expression of available genes."""
    avail = [g for g in genes if g in adata_obj.var_names]
    if not avail:
        return np.zeros(adata_obj.n_obs)
    mat = adata_obj[:, avail].layers[layer] if layer else adata_obj[:, avail].X
    if sp.issparse(mat):
        mat = np.array(mat.todense())
    return np.mean(mat, axis=1).A1 if hasattr(mat, 'A1') else np.mean(mat, axis=1)

# Normalize for scoring
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)

c_score = score_genes(adata_norm, CONTRACTILE_GENES, layer=None)
s_score = score_genes(adata_norm, SYNTHETIC_GENES, layer=None)

adata.obs["contractile_score"] = c_score
adata.obs["synthetic_score"] = s_score
adata.obs["phenotype"] = np.where(
    c_score > 0,
    np.where(c_score / (s_score + 0.01) > 1.5, "Contractile", "Synthetic"),
    "Synthetic"
)

del adata_norm

print("Phenotype distribution:")
print(adata.obs["phenotype"].value_counts())
print("\nPer-group:")
for g in ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]:
    mask = adata.obs["group"] == g
    n_contr = ((adata.obs["phenotype"] == "Contractile") & mask).sum()
    n_total = mask.sum()
    pct = 100 * n_contr / max(n_total, 1)
    print(f"  {g}: {pct:.1f}% contractile ({n_contr}/{n_total})")


# ── Preprocessing for velocity ───────────────────────────────────────────────
print("\n" + "=" * 70)
print("VELOCITY PREPROCESSING")
print("=" * 70)

# Use spliced as main matrix for scVelo
adata.X = adata.layers["spliced"].copy()

# Basic QC
adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()
adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()

print(f"Before QC: {adata.n_obs} cells")
print(f"  Median counts: {np.median(adata.obs['n_counts']):.0f}")
print(f"  Median genes: {np.median(adata.obs['n_genes']):.0f}")

# Filter low-quality cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")

# scVelo preprocessing
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
# Force-include Table S2 genes + Piezo1 in HVG selection
FORCE_GENES = ["Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn",
               "Col1a1", "Fn1", "Spp1", "Klf4", "Vim", "Mmp2", "Mmp9", "Timp1",
               "Piezo1"]
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="spliced")
# Force-include key genes even if not in top 2000 HVGs
for gene in FORCE_GENES:
    if gene in adata.var_names:
        adata.var.loc[gene, "highly_variable"] = True
n_forced = sum(1 for g in FORCE_GENES if g in adata.var_names and g not in
               adata.var_names[adata.var["highly_variable"]])
print(f"Force-included {sum(1 for g in FORCE_GENES if g in adata.var_names)} Table S2 genes")
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"After scVelo filter: {adata.n_obs} cells, {adata.n_vars} genes")

sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
sc.tl.umap(adata)

# Moments for velocity
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# ── Velocity computation ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("COMPUTING RNA VELOCITY")
print("=" * 70)

scv.tl.recover_dynamics(adata, n_jobs=1)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

# Velocity confidence
scv.tl.velocity_confidence(adata)
print(f"Mean velocity confidence: {adata.obs['velocity_confidence'].mean():.3f}")
print(f"Mean velocity length: {adata.obs['velocity_length'].mean():.3f}")


# ── Visualization ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

# --- Panel A: UMAP colored by group ---
fig, ax = plt.subplots(figsize=(8, 6))
for g in ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]:
    mask = adata.obs["group"] == g
    ax.scatter(
        adata.obsm["X_umap"][mask, 0],
        adata.obsm["X_umap"][mask, 1],
        c=GROUP_COLORS[g], label=g, s=10, alpha=0.7
    )
ax.legend(frameon=False, fontsize=10)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("SMC RNA Velocity — Groups")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "A_umap_groups.png"), dpi=300)
plt.close()
print("  A: UMAP by group")

# --- Panel B: UMAP colored by phenotype ---
fig, ax = plt.subplots(figsize=(8, 6))
for p in ["Contractile", "Synthetic"]:
    mask = adata.obs["phenotype"] == p
    ax.scatter(
        adata.obsm["X_umap"][mask, 0],
        adata.obsm["X_umap"][mask, 1],
        c=PHENOTYPE_COLORS[p], label=p, s=10, alpha=0.7
    )
ax.legend(frameon=False, fontsize=10)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.set_title("SMC RNA Velocity — Phenotype")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "B_umap_phenotype.png"), dpi=300)
plt.close()
print("  B: UMAP by phenotype")

# Set up categorical colors for scVelo compatibility
adata.obs["group"] = pd.Categorical(
    adata.obs["group"],
    categories=["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]
)
adata.uns["group_colors"] = [GROUP_COLORS[g] for g in adata.obs["group"].cat.categories]

adata.obs["phenotype"] = pd.Categorical(
    adata.obs["phenotype"],
    categories=["Contractile", "Synthetic"]
)
adata.uns["phenotype_colors"] = [PHENOTYPE_COLORS[p] for p in adata.obs["phenotype"].cat.categories]

# --- Panel C: Velocity streamplot colored by group ---
fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="group",
    ax=ax, legend_loc="right margin",
    title="RNA Velocity Streamplot — Groups", show=False
)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "C_velocity_stream_groups.png"), dpi=300)
plt.close()
print("  C: Velocity streamplot by group")

# --- Panel D: Velocity streamplot colored by phenotype ---
fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="phenotype",
    ax=ax, legend_loc="right margin",
    title="RNA Velocity Streamplot — Phenotype", show=False
)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "D_velocity_stream_phenotype.png"), dpi=300)
plt.close()
print("  D: Velocity streamplot by phenotype")

# --- Panel E: Velocity arrows (grid) ---
fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_grid(
    adata, basis="umap", color="group",
    ax=ax, legend_loc="right margin",
    title="RNA Velocity Grid — Groups", show=False,
    arrow_length=3, arrow_size=2
)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "E_velocity_grid_groups.png"), dpi=300)
plt.close()
print("  E: Velocity grid by group")

# --- Panel F: Velocity per sample ---
samples_order = [s for s in ["WT_Young_1", "WT_Young_2", "WT_Aged_1", "WT_Aged_2",
                              "KO_Young_1", "KO_Aged"] if s in set(adata.obs["sample"])]
ncol = 3; nrow = int(np.ceil(len(samples_order) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow), squeeze=False)
for i, s in enumerate(samples_order):
    ax = axes[i // ncol, i % ncol]
    mask = adata.obs["sample"] == s
    grp = adata.obs.loc[mask, "group"].iloc[0] if mask.sum() > 0 else "unknown"

    # Plot all cells in grey
    ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
               c="#EEEEEE", s=5, alpha=0.3)
    # Highlight this sample
    ax.scatter(adata.obsm["X_umap"][mask, 0], adata.obsm["X_umap"][mask, 1],
               c=GROUP_COLORS.get(grp, "#333333"), s=15, alpha=0.8)
    ax.set_title(f"{s}\n({grp}, n={mask.sum()})", fontsize=13)
    ax.set_xticks([])
    ax.set_yticks([])

# Hide unused subplots
for j in range(len(samples_order), nrow * ncol):
    axes[j // ncol, j % ncol].set_visible(False)
plt.suptitle("Per-Sample Distribution on Velocity UMAP", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "F_per_sample_umap.png"), dpi=300)
plt.close()
print("  F: Per-sample UMAP")

# --- Panel G: Key gene velocity ---
velocity_genes = ["Myh11", "Acta2", "Tagln", "Cnn1", "Piezo1",
                  "Col1a1", "Fn1", "Klf4", "Myocd", "Smtn"]
avail_genes = [g for g in velocity_genes if g in adata.var_names]
if avail_genes:
    try:
        # one row per gene: phase portrait, velocity and expression on the UMAP
        scv.pl.velocity(adata, var_names=avail_genes, ncols=1, figsize=(4.5, 4), show=False)
        plt.gcf().savefig(os.path.join(FIG_DIR, "G_gene_velocity.png"), dpi=200, bbox_inches="tight")
        plt.close()
        print(f"  G: Gene velocity phase plots ({len(avail_genes)} genes)")
    except Exception as e:
        print(f"  G: Gene velocity failed: {e}")

# --- Panel H: Velocity confidence ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="velocity_confidence", ax=axes[0], show=False,
           title="Velocity Confidence", color_map="coolwarm")
sc.pl.umap(adata, color="velocity_length", ax=axes[1], show=False,
           title="Velocity Length", color_map="coolwarm")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "H_velocity_confidence.png"), dpi=300)
plt.close()
print("  H: Velocity confidence & length")

# --- Panel I: Latent time ---
try:
    scv.tl.latent_time(adata)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="latent_time", ax=ax, show=False,
               title="Latent Time", color_map="gnuplot")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "I_latent_time.png"), dpi=300)
    plt.close()
    print("  I: Latent time")

    # Latent time by group
    fig, ax = plt.subplots(figsize=(8, 5))
    for g in ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]:
        mask = adata.obs["group"] == g
        vals = adata.obs.loc[mask, "latent_time"]
        ax.hist(vals, bins=30, alpha=0.5, label=g, color=GROUP_COLORS[g], density=True)
    ax.set_xlabel("Latent Time")
    ax.set_ylabel("Density")
    ax.set_title("Latent Time Distribution by Group")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "J_latent_time_by_group.png"), dpi=300)
    plt.close()
    print("  J: Latent time by group")
except Exception as e:
    print(f"  Latent time failed: {e}")

# --- Panel K: Top velocity genes ---
try:
    scv.tl.rank_velocity_genes(adata, groupby="group")
    top_genes = pd.DataFrame(adata.uns['rank_velocity_genes']['names']).head(10)
    top_genes.to_csv(os.path.join(OUTPUT_DIR, "top_velocity_genes_by_group.csv"), index=False)
    print("  K: Top velocity genes saved")
except Exception as e:
    print(f"  Top velocity genes failed: {e}")


# ── Save results ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SAVING")
print("=" * 70)

adata.write(os.path.join(OUTPUT_DIR, "velocity_smc.h5ad"))
print(f"Saved: {OUTPUT_DIR}/velocity_smc.h5ad")

# Summary stats
summary = {
    "total_cells": adata.n_obs,
    "total_genes": adata.n_vars,
    "mean_confidence": adata.obs["velocity_confidence"].mean(),
}
for g in ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]:
    mask = adata.obs["group"] == g
    summary[f"{g}_cells"] = int(mask.sum())
    summary[f"{g}_confidence"] = float(adata.obs.loc[mask, "velocity_confidence"].mean())

pd.Series(summary).to_csv(os.path.join(OUTPUT_DIR, "velocity_summary.csv"))
print("Saved: velocity_summary.csv")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
print(f"Figures: {FIG_DIR}/")
print(f"Data: {OUTPUT_DIR}/velocity_smc.h5ad")
