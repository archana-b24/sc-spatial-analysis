"""Visium analysis (DEVELOPMENT_LOG §3.8).

v1: Scanorama joint embedding + cosine "label transfer" assigned one cell type per 55 µm spot (S9);
three overlapping exploratory notebooks, none ran end-to-end; no spatial statistics (S10).

v2:
* Load Space Ranger output (official 10x 2.0.1 release, checksummed); also export a SpatialData zarr.
* **cell2location** (Kleshchevnikov 2022) — Bayesian deconvolution into per-spot cell-type *abundances*
  using reference signatures learned from the QC'd, CNV-annotated scRNA-seq (incl. a Malignant class).
  Runs in the isolated `c2l` pixi env (see `deconvolve_c2l` / scripts/run_c2l.py).
* **squidpy**: spatially variable genes (Moran's I), neighbourhood enrichment / co-occurrence of
  dominant cell types.
* **CellCharter** (Varrone 2024): spatial niches from neighbourhood-aggregated abundance profiles,
  with automatic selection of the number of niches by stability.
* n = 1 section → every spatial result is descriptive; stated in the report.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def load_visium(raw_dir: Path, sample_id: str) -> ad.AnnData:
    import squidpy as sq

    adata = sq.read.visium(raw_dir, counts_file="filtered_feature_bc_matrix.h5", library_id=sample_id)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    return adata


def write_spatialdata(raw_dir: Path, out_zarr: Path, sample_id: str) -> None:
    """Interoperable SpatialData/zarr export (images + shapes + table)."""
    try:
        import spatialdata_io

        sdata = spatialdata_io.visium(
            raw_dir, dataset_id=sample_id, counts_file="filtered_feature_bc_matrix.h5"
        )
        sdata.write(out_zarr, overwrite=True)
    except Exception as exc:  # noqa: BLE001 — export is a convenience artefact, never blocks analysis
        log.warning("SpatialData export skipped: %s", exc)


def qc_and_cluster(adata: ad.AnnData, scfg: dict, seed: int = 0) -> pd.DataFrame:
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=scfg["min_counts"])
    sc.pp.filter_cells(adata, min_genes=scfg["min_genes"])
    adata._inplace_subset_obs(adata.obs["pct_counts_mt"] < scfg["max_pct_mt"])
    sc.pp.filter_genes(adata, min_cells=scfg["min_cells_per_gene"])
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.pca(adata, random_state=seed)
    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata, key_added="expr_cluster", flavor="igraph", n_iterations=2, directed=False, random_state=seed
    )
    return pd.DataFrame(
        {
            "metric": ["spots_in", "spots_kept", "genes_kept", "median_counts", "median_genes"],
            "value": [
                n0,
                adata.n_obs,
                adata.n_vars,
                float(np.median(adata.obs.total_counts)),
                float(np.median(adata.obs.n_genes_by_counts)),
            ],
        }
    )


def spatial_stats(
    adata: ad.AnnData, scfg: dict, cell_types: list[str], n_jobs: int = 16
) -> dict[str, pd.DataFrame]:
    import squidpy as sq

    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=scfg["nhood_n_rings"])
    hv = adata.var_names[adata.var["highly_variable"]]
    sq.gr.spatial_autocorr(adata, mode="moran", genes=hv, n_perms=100, n_jobs=n_jobs, seed=0)
    out = {"svg": adata.uns["moranI"].rename_axis("gene").reset_index()}

    # dominant cell type per spot from cell2location proportions
    props = adata.obs[cell_types]
    adata.obs["dominant_type"] = pd.Categorical(props.idxmax(axis=1))
    sq.gr.nhood_enrichment(adata, cluster_key="dominant_type", seed=0, n_jobs=n_jobs)
    z = adata.uns["dominant_type_nhood_enrichment"]["zscore"]
    cats = adata.obs["dominant_type"].cat.categories
    out["nhood_enrichment"] = pd.DataFrame(z, index=cats, columns=cats).rename_axis("source").reset_index()
    # abundance-level co-localisation: Pearson correlation of proportions across spots
    out["colocalisation"] = props.corr().rename_axis("cell_type").reset_index()
    return out


def niches_cellcharter(adata: ad.AnnData, scfg: dict, cell_types: list[str], seed: int = 0) -> pd.DataFrame:
    import cellcharter as cc
    import squidpy as sq

    adata.obsm["X_c2l_prop"] = adata.obs[cell_types].to_numpy()
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
    cc.gr.aggregate_neighbors(
        adata, n_layers=scfg["cellcharter"]["n_layers"], use_rep="X_c2l_prop", out_key="X_cellcharter"
    )
    lo, hi = scfg["cellcharter"]["cluster_range"]
    model = cc.tl.ClusterAutoK(
        n_clusters=(lo, hi), max_runs=5, convergence_tol=1e-3, model_params={"random_state": seed}
    )
    model.fit(adata, use_rep="X_cellcharter")
    adata.obs["niche"] = pd.Categorical(model.predict(adata, use_rep="X_cellcharter").astype(str))
    comp = adata.obs.groupby("niche", observed=True)[cell_types].mean()
    return comp.rename_axis("niche").reset_index()
