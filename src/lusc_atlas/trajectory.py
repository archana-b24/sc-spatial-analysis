"""Basal → malignant transcriptional progression (DEVELOPMENT_LOG §3.7).

v1: PAGA connectivity graphs only; `04_trajectory` never ran (S8).
v2: on the epithelial basal-lineage compartment (normal basal/suprabasal/hillock + CNV-called malignant):

* **Palantir** pseudotime (Setty 2019) rooted in normal-donor *Basal resting* cells.
* **CellRank 2** (Weiler 2024) combining a PseudotimeKernel (Palantir) with a CytoTRACEKernel
  (developmental potential, no RNA velocity needed — HLCA has no unspliced counts) → macrostates and
  fate probabilities toward the malignant state.

Interpretation caveat (logged in the report): malignant cells are patient-clonal; across donors this is a
*transcriptional-distance-from-normal-basal* axis, not a lineage tracing result.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import scanpy as sc

log = logging.getLogger(__name__)


def subset_lineage(adata: ad.AnnData, tcfg: dict, label_key: str) -> ad.AnnData:
    lab = adata.obs[label_key].astype(str)
    keep = lab.isin(tcfg["include_labels"]) | adata.obs.get("malignant", False)
    sub = adata[keep.values].copy()
    sub.obs["lineage_state"] = np.where(sub.obs["malignant"], "Malignant", lab[keep.values])
    sub.obs["lineage_state"] = sub.obs["lineage_state"].astype("category")
    return sub


def run_trajectory(sub: ad.AnnData, tcfg: dict, seed: int = 0) -> dict:
    import cellrank as cr
    import palantir

    sc.pp.neighbors(sub, use_rep="X_integrated", random_state=seed)
    sc.tl.umap(sub, random_state=seed)

    # Root: the normal-donor basal cell closest to the basal centroid in scVI space
    root_mask = (sub.obs["lineage_state"] == tcfg["root_label"]) & (sub.obs["condition"] == "healthy")
    if root_mask.sum() == 0:
        raise ValueError(f"no healthy-donor '{tcfg['root_label']}' cells to root the trajectory")
    Z = sub.obsm["X_integrated"]
    centroid = Z[root_mask.values].mean(0)
    cand = np.where(root_mask.values)[0]
    root = sub.obs_names[cand[np.argmin(((Z[cand] - centroid) ** 2).sum(1))]]

    sub.obsm["X_pca_for_palantir"] = sub.obsm["X_integrated"]
    dm = palantir.utils.run_diffusion_maps(
        sub, pca_key="X_pca_for_palantir", n_components=tcfg["n_components"], seed=seed
    )
    palantir.utils.determine_multiscale_space(sub)
    palantir.core.run_palantir(sub, root, num_waypoints=tcfg["num_waypoints"], seed=seed)
    del dm

    # CellRank 2: pseudotime kernel + CytoTRACE kernel
    # CytoTRACE aggregates *imputed* expression: kNN-smoothed log-expression over the integrated-space graph
    # (the scVelo-moments equivalent CellRank's tutorials use; HLCA/LuCA have no unspliced counts anyway)
    import scipy.sparse as sp

    W = sub.obsp["connectivities"] + sp.identity(sub.n_obs, format="csr")
    W = sp.diags(1.0 / np.asarray(W.sum(axis=1)).ravel()) @ W
    sub.layers["imputed"] = sp.csr_matrix(W @ sub.layers["lognorm"])
    ctk = (
        cr.kernels.CytoTRACEKernel(sub).compute_cytotrace().compute_transition_matrix(threshold_scheme="soft")
    )
    ptk = cr.kernels.PseudotimeKernel(sub, time_key="palantir_pseudotime").compute_transition_matrix()
    kernel = 0.5 * ctk + 0.5 * ptk
    g = cr.estimators.GPCCA(kernel)
    g.fit(n_states=tcfg["n_macrostates"], cluster_key="lineage_state")
    g.predict_terminal_states()
    try:
        g.compute_fate_probabilities()
    except Exception as exc:  # noqa: BLE001 — fate probs need ≥1 terminal state; report and continue
        log.warning("fate probabilities not computed: %s", exc)
    del sub.layers["imputed"]
    return {
        "root": root,
        "terminal_states": list(map(str, g.terminal_states.cat.categories))
        if g.terminal_states is not None
        else [],
    }
