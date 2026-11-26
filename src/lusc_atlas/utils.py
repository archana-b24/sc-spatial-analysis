"""Small shared helpers."""

from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

log = logging.getLogger(__name__)


def robust_hvg(
    adata: ad.AnnData,
    n_top: int,
    layer: str | None,
    batch_key: str | None,
    flavor: str = "seurat_v3",
    subset: bool = False,
) -> None:
    """Batch-aware HVG selection that survives small batches.

    seurat_v3 fits a LOESS mean–variance curve *per batch*; tiny batches make it singular
    ("near singularities"). We select HVGs batch-aware at the study level (large batches) and fall back to a
    non-batch fit only if that still fails — the fallback is logged, never silent.
    """
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top, flavor=flavor, layer=layer, batch_key=batch_key, subset=subset
        )
    except ValueError as exc:
        log.warning("batch-aware HVG failed (%s); falling back to batch_key=None", str(exc)[:120])
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=flavor, layer=layer, subset=subset)
