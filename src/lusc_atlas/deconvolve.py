"""cell2location deconvolution — runs ONLY in the isolated `c2l` pixi environment.

cell2location 0.1.5 needs scvi-tools 1.3.x and pyro < 1.9.2, which conflict with the main stack
(scvi-tools 1.5). Isolating it keeps both sides on their tested versions (DEVELOPMENT_LOG §2.3).

Two stages (Kleshchevnikov 2022):
1. RegressionModel on scRNA-seq counts → per-cell-type reference expression signatures,
   with sample and study as batch/covariates so technical effects are not baked into signatures.
2. Cell2location on Visium counts → posterior cell abundance per spot: posterior mean (point estimate,
   used for proportions and niches) and q05 (conservative lower bound, used for presence calls).
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def prepare_reference_labels(
    labels: pd.Series,
    min_cells: int = 0,
    merge: dict[str, str] | None = None,
    exclude: list[str] | tuple[str, ...] = (),
) -> tuple[pd.Series, pd.DataFrame]:
    """Decide which labels enter the cell2location reference (§3.10 #3).

    Run 1 built a signature for every label, including `Rare epithelial` (66 healthy HLCA cells: ionocytes,
    neuroendocrine, brush cells). A signature estimated from so few cells is noisy enough to soak up
    unexplained epithelial signal in tumour spots: it took 66 % of the tumour niche. Rules, in order:

    1. ``exclude`` labels are removed (e.g. ``Unknown``: a mixed bucket has no coherent signature);
    2. a label with fewer than ``min_cells`` cells is merged into ``merge[label]`` if given, else dropped.

    Counts are taken on the full reference, before the per-label subsample. Returns the new labels
    (NaN = cell not used) and a decision table (label, n_cells, action, target).
    """
    merge = merge or {}
    lab = labels.astype(str)
    counts = lab.value_counts()
    rows, mapping = [], {}
    for label, n in counts.items():
        if label in exclude:
            action, target = "excluded", None
        elif n < min_cells and label in merge:
            action, target = "merged", merge[label]
            if target not in counts.index:
                raise ValueError(f"merge target {target!r} for {label!r} is not a reference label")
        elif n < min_cells:
            action, target = "dropped (< min_ref_cells)", None
        else:
            action, target = "kept", label
        mapping[label] = target
        rows.append({"label": label, "n_cells": int(n), "action": action, "target": target})
    decisions = pd.DataFrame(rows).sort_values("n_cells", ascending=False, ignore_index=True)
    for _, r in decisions[decisions.action != "kept"].iterrows():
        log.info("c2l reference: %s (%d cells) → %s", r.label, r.n_cells, r.action)
    new = lab.map(mapping)
    return new.where(new.notna(), np.nan), decisions


def reference_signatures(ref: ad.AnnData, c2l_cfg: dict, seed: int = 0) -> pd.DataFrame:
    import cell2location
    import scvi
    from cell2location.utils.filtering import filter_genes

    scvi.settings.seed = seed
    label = c2l_cfg["reference_label"]
    # balanced subsample per label (bounded CPU time; rare types kept whole)
    rng = np.random.default_rng(seed)
    keep = []
    for _, idx in ref.obs.groupby(label, observed=True).indices.items():
        k = c2l_cfg["max_ref_cells_per_label"]
        keep.append(idx if len(idx) <= k else rng.choice(idx, k, replace=False))
    ref = ref[np.sort(np.concatenate(keep))].copy()
    ref.X = ref.layers["counts"]
    sel = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    ref = ref[:, sel].copy()
    covs = ["study"] if "study" in ref.obs and ref.obs["study"].nunique() > 1 else None
    cell2location.models.RegressionModel.setup_anndata(
        ref, batch_key="sample", labels_key=label, categorical_covariate_keys=covs
    )
    mod = cell2location.models.RegressionModel(ref)
    mod.train(max_epochs=c2l_cfg["ref_max_epochs"], accelerator="cpu")
    ref = mod.export_posterior(ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})
    cols = [f"means_per_cluster_mu_fg_{c}" for c in ref.uns["mod"]["factor_names"]]
    sig = ref.varm["means_per_cluster_mu_fg"][cols].copy()
    sig.columns = ref.uns["mod"]["factor_names"]
    return sig


def map_spots(vis: ad.AnnData, sig: pd.DataFrame, c2l_cfg: dict, seed: int = 0) -> pd.DataFrame:
    import cell2location
    import scvi

    scvi.settings.seed = seed
    shared = vis.var_names.intersection(sig.index)
    log.info("cell2location: %d shared genes", len(shared))
    vis = vis[:, shared].copy()
    vis.X = vis.layers["counts"] if "counts" in vis.layers else vis.X
    cell2location.models.Cell2location.setup_anndata(vis)
    mod = cell2location.models.Cell2location(
        vis,
        cell_state_df=sig.loc[shared],
        N_cells_per_location=c2l_cfg["n_cells_per_location"],
        detection_alpha=c2l_cfg["detection_alpha"],
    )
    mod.train(max_epochs=c2l_cfg["max_epochs"], batch_size=None, train_size=1, accelerator="cpu")
    vis = mod.export_posterior(vis, sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs})
    props = posterior_proportions(mod, vis, c2l_cfg.get("proportion_samples", 500))
    # §3.12: posterior MEAN abundance is the point estimate used for proportions/niches; q05 (5 % quantile, a
    # conservative lower bound) is kept for "is this type confidently present" statements only.
    out = {}
    for key, prefix in (
        ("means_cell_abundance_w_sf", "meanscell_abundance_w_sf_"),
        ("q05_cell_abundance_w_sf", "q05cell_abundance_w_sf_"),
    ):
        ab = vis.obsm[key].copy()
        ab.columns = [c.replace(prefix, "") for c in ab.columns]
        out[key.split("_")[0]] = ab
    if props is not None:
        props.index = out["means"].index
        props.columns = list(out["means"].columns)
    return out["means"], out["q05"], props


def posterior_proportions(mod, vis: ad.AnnData, num_samples: int = 500) -> pd.DataFrame | None:
    """Per-spot cell-type proportions as the posterior mean of the RATIO (senior review item 13):
    proportions are
    computed in every posterior sample of the abundance ``w_sf`` and then averaged. A ratio of posterior
    means is
    not the posterior mean of the ratio. Returns None (caller falls back to ratio of means, logged) on
    failure."""
    try:
        smp = mod.sample_posterior(
            num_samples=num_samples, batch_size=vis.n_obs, return_samples=True, return_sites=["w_sf"]
        )
        w = np.asarray(smp["posterior_samples"]["w_sf"])  # samples × spots × types
        prop = (w / np.clip(w.sum(-1, keepdims=True), 1e-12, None)).mean(0)
        return pd.DataFrame(prop, index=vis.obs_names)
    except Exception as exc:  # noqa: BLE001
        log.warning("posterior proportions failed (%s); falling back to ratio of posterior means", exc)
        return None
