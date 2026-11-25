"""Cross-atlas integration: benchmark → guard against over-correction → choose (DEVELOPMENT_LOG §3.3 F–H).

v1 normalised/clustered tumour and control separately with no batch model (S5).
v2 integrates HLCA + LuCA cells in one latent space and *chooses* the method on evidence:

Candidates (all on the same HVGs, raw counts where the model needs counts):
* ``X_pca``      unintegrated PCA                                — the floor
* ``X_harmony``  Harmony on PCA, batch = sample (Korsunsky 2019) — classical baseline
* ``X_scVI``     scVI, batch = sample, covariates study + assay  (Lopez 2018)
* ``X_scANVI``   scANVI from scVI with harmonised coarse labels  (Xu 2021)
* ``X_sysVI``    sysVI, system = atlas, covariate = sample       (Hrovatin et al.; substantial batch effects)
* ``X_scArches`` transfer learning: scArches surgery onto the published HLCA scANVI model (Lotfollahi 2022)
* ``X_geneformer`` foundation model, zero-shot Geneformer-V2 (cancer CL) — benchmark subsample only

Selection:
1. scib-metrics (Luecken 2022) with batch = study, label = harmonised coarse label.
2. Over-correction guard: kNN label purity (fraction of a cell's neighbours sharing its harmonised label)
   must stay within `purity_tolerance` of the unintegrated PCA purity — integration may mix *batches*,
   never *cell types*.
3. Highest scib Total among methods that pass the guard. All scores are written out.

The chosen embedding is stored as ``X_integrated``; downstream neighbours/UMAP/Leiden use it.
Integration is for representation only — DE always uses raw counts (§3.3 H).
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def normalise(adata: ad.AnnData) -> None:
    """Log-normalised expression in .X and layers['lognorm']; counts stay in layers['counts']."""
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()


def select_hvgs(adata: ad.AnnData, icfg: dict) -> ad.AnnData:
    from .utils import robust_hvg

    robust_hvg(
        adata, icfg["n_hvg"], layer="counts", batch_key=icfg["hvg_batch_key"], flavor=icfg["hvg_flavor"]
    )
    return adata[:, adata.var["highly_variable"]].copy()


def _covariates(hv: ad.AnnData, keys: list[str]) -> list[str] | None:
    covs = [k for k in keys if k in hv.obs.columns and hv.obs[k].nunique() > 1]
    return covs or None


def embed_pca_harmony(adata: ad.AnnData, hv: ad.AnnData, icfg: dict, seed: int) -> None:
    tmp = hv.copy()
    tmp.X = tmp.layers["lognorm"]
    sc.pp.scale(tmp, max_value=10)
    sc.pp.pca(tmp, n_comps=icfg["n_latent"], random_state=seed)
    adata.obsm["X_pca"] = tmp.obsm["X_pca"]
    try:
        import harmonypy

        ho = harmonypy.run_harmony(
            tmp.obsm["X_pca"], tmp.obs, [icfg["batch_key"]], random_state=seed, max_iter_harmony=20
        )
        Z = np.asarray(ho.Z_corr)
        adata.obsm["X_harmony"] = Z.T if Z.shape[0] != adata.n_obs else Z
    except Exception as exc:  # noqa: BLE001 — baseline only; benchmark proceeds without it
        log.warning("Harmony failed: %s", exc)


def embed_scvi_scanvi(
    adata: ad.AnnData, hv: ad.AnnData, icfg: dict, label_key: str, heldout: np.ndarray | None = None
):
    import scvi

    covs = _covariates(hv, icfg["categorical_covariates"])
    scvi.model.SCVI.setup_anndata(
        hv, layer="counts", batch_key=icfg["batch_key"], categorical_covariate_keys=covs
    )
    vae = scvi.model.SCVI(hv, n_latent=icfg["n_latent"], n_layers=icfg["n_layers"], gene_likelihood="nb")
    vae.train(max_epochs=icfg["max_epochs"], early_stopping=icfg["early_stopping"])
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    hv.obs["_scanvi_label"] = hv.obs[label_key].astype(str).replace({"nan": "Unknown"})
    if heldout is not None:
        # §3.12: the benchmark cells' labels are hidden from scANVI, so its bio-conservation scores (and the
        # label-transfer check on the dashboard) are measured on cells whose labels it never saw. Runs 1–2
        # scored
        # scANVI on its own training labels — circular, and part of why it "won".
        lab = hv.obs["_scanvi_label"].to_numpy().copy()
        lab[heldout] = "Unknown"
        hv.obs["_scanvi_label"] = lab
    lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="Unknown", labels_key="_scanvi_label")
    lvae.train(max_epochs=icfg["scanvi_max_epochs"], n_samples_per_label=icfg["scanvi_samples_per_label"])
    adata.obsm["X_scANVI"] = lvae.get_latent_representation()
    adata.obs["scanvi_prediction"] = lvae.predict()
    return vae, lvae


def embed_sysvi(adata: ad.AnnData, hv: ad.AnnData, icfg: dict) -> None:
    from scvi.external import SysVI

    tmp = hv.copy()
    tmp.X = tmp.layers["lognorm"]  # sysVI models log-normalised expression (Gaussian likelihood)
    SysVI.setup_anndata(
        tmp, batch_key=icfg["sysvi_system_key"], categorical_covariate_keys=[icfg["batch_key"]]
    )
    m = SysVI(tmp, prior="vamp", n_latent=icfg["n_latent"])
    m.train(
        max_epochs=icfg["max_epochs"],
        plan_kwargs={"kl_weight": 1, "z_distance_cycle_weight": icfg["sysvi_cycle_weight"]},
    )
    adata.obsm["X_sysVI"] = m.get_latent_representation(adata=tmp)


def knn_label_purity(Z: np.ndarray, labels: pd.Series, k: int = 15) -> float:
    from sklearn.neighbors import NearestNeighbors

    lab = labels.astype(str).to_numpy()
    known = lab != "Unknown"
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Z)
    idx = nn.kneighbors(Z[known], return_distance=False)[:, 1:]
    return float((lab[idx] == lab[known][:, None]).mean())


BENCHMARK_ONLY = {"X_geneformer"}  # computed on the benchmark subsample only → never chosen downstream


def benchmark_cells(obs_names: pd.Index, n: int, seed: int) -> np.ndarray:
    """Deterministic benchmark subsample (positions). Shared by `fm_embed` and `integrate` so every
    candidate — including the CPU-expensive foundation model — is scored on the *same* cells."""
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(len(obs_names), size=min(n, len(obs_names)), replace=False))


def benchmark(
    adata: ad.AnnData, keys: list[str], label_key: str, batch_key: str, idx: np.ndarray, n_jobs: int
) -> pd.DataFrame:
    """scib-metrics + kNN label purity on the shared subsample (cells with known harmonised label)."""
    from scib_metrics.benchmark import Benchmarker

    sub = adata[idx].copy()
    sub = sub[sub.obs[label_key].astype(str) != "Unknown"].copy()
    bm = Benchmarker(
        sub,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=keys,
        pre_integrated_embedding_obsm_key="X_pca",
        n_jobs=n_jobs,
    )
    bm.benchmark()
    res = bm.get_results(min_max_scale=False)
    res = res.drop(index="Metric Type", errors="ignore").apply(pd.to_numeric, errors="coerce")
    res["knn_label_purity"] = [knn_label_purity(sub.obsm[k], sub.obs[label_key]) for k in res.index]
    res["benchmark_only"] = [k in BENCHMARK_ONLY for k in res.index]
    return res.rename_axis("embedding").reset_index()


def choose(scores: pd.DataFrame, tolerance: float) -> str:
    """Best scib Total among methods that pass the over-correction guard and can embed every cell."""
    base = float(scores.loc[scores.embedding == "X_pca", "knn_label_purity"].iloc[0])
    eligible = ~scores.embedding.isin(["X_pca", *BENCHMARK_ONLY])
    ok = scores[eligible & (scores.knn_label_purity >= base - tolerance)]
    if ok.empty:
        log.warning("every method failed the over-correction guard; falling back to scANVI")
        return "X_scANVI"
    return str(ok.sort_values("Total", ascending=False).embedding.iloc[0])


def run_integration(
    adata: ad.AnnData,
    icfg: dict,
    seed: int = 0,
    n_jobs: int = 16,
    repo=None,
    fm_embedding: tuple[np.ndarray, np.ndarray] | None = None,
) -> tuple[pd.DataFrame, str, object, dict]:
    import scvi
    import torch

    torch.set_num_threads(n_jobs)
    scvi.settings.seed = seed
    scvi.settings.num_threads = n_jobs
    label = icfg["harmonised_label_key"]
    extra: dict = {}

    hv = select_hvgs(adata, icfg)
    embed_pca_harmony(adata, hv, icfg, seed)
    idx = benchmark_cells(adata.obs_names, icfg["benchmark_n_cells"], seed)
    adata.obs["scanvi_heldout"] = False
    adata.obs.iloc[idx, adata.obs.columns.get_loc("scanvi_heldout")] = True
    vae, lvae = embed_scvi_scanvi(adata, hv, icfg, label, heldout=idx)
    if icfg["run_sysvi"]:
        try:
            embed_sysvi(adata, hv, icfg)
        except Exception as exc:  # noqa: BLE001
            log.warning("sysVI failed: %s", exc)
    if icfg["transfer"]["enabled"]:  # transfer learning: scArches surgery onto the HLCA reference
        from .transfer import map_to_reference

        try:
            extra["scarches"] = map_to_reference(adata, icfg["transfer"], repo, seed)
        except Exception as exc:  # noqa: BLE001
            log.warning("scArches mapping failed: %s", exc)
            extra["scarches_error"] = str(exc)

    if fm_embedding is not None:  # foundation model: zero-shot Geneformer, benchmark subsample only
        names, emb = fm_embedding
        if not np.array_equal(names, adata.obs_names[idx].to_numpy()):
            raise ValueError("foundation-model embedding was computed on a different subsample")
        full = np.full((adata.n_obs, emb.shape[1]), np.nan, dtype=np.float32)
        full[idx] = emb
        adata.obsm["X_geneformer"] = full

    keys = [
        k
        for k in ["X_pca", "X_harmony", "X_scVI", "X_scANVI", "X_sysVI", "X_scArches", "X_geneformer"]
        if k in adata.obsm
    ]
    scores = benchmark(adata, keys, label, icfg["benchmark_batch_key"], idx, n_jobs)
    chosen = choose(scores, icfg["purity_tolerance"])
    scores["chosen"] = scores.embedding == chosen
    log.info("integration benchmark → %s\n%s", chosen, scores.to_string())

    adata.obsm["X_integrated"] = adata.obsm[chosen]
    if "X_geneformer" in adata.obsm:  # NaN rows would break downstream neighbour graphs
        del adata.obsm["X_geneformer"]
    adata.uns["integration"] = {"chosen": chosen, "candidates": keys}
    sc.pp.neighbors(adata, use_rep="X_integrated", random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(
        adata,
        resolution=icfg["leiden_resolution"],
        random_state=seed,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return scores, chosen, lvae, extra


def atlas_mixing_per_label(
    adata: ad.AnnData,
    label_key: str,
    atlas_key: str = "atlas",
    k: int = 30,
    conditions: list[str] | None = None,
    scope: str = "all cells",
) -> pd.DataFrame:
    """Per harmonised cell type: fraction of kNN neighbours from the *other* atlas (high = well mixed).

    ``conditions`` restricts to cells of those conditions (§3.12). Over all cells the metric compares healthy
    HLCA cells with tumour-patient LuCA cells, so "below random" is expected biology, not under-integration.
    Restricted to tumour + adjacent tissue (HLCA's Thienpont patients vs LuCA patients) it compares like with
    like: the residual is the technical atlas effect.
    """
    from sklearn.neighbors import NearestNeighbors

    if conditions is not None:
        adata = adata[adata.obs["condition"].astype(str).isin(conditions).to_numpy()]
    Z = np.asarray(adata.obsm["X_integrated"])
    atlas = adata.obs[atlas_key].astype(str).to_numpy()
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Z)
    rows = []
    for lab, idx in adata.obs.groupby(label_key, observed=True).indices.items():
        if len(idx) < 50 or len(set(atlas[idx])) < 2:
            continue
        sample = idx if len(idx) <= 5000 else np.random.default_rng(0).choice(idx, 5000, replace=False)
        nbr = nn.kneighbors(Z[sample], return_distance=False)[:, 1:]
        other = (atlas[nbr] != atlas[sample][:, None]).mean()
        frac = pd.Series(atlas[idx]).value_counts(normalize=True)
        rows.append(
            {
                "label": lab,
                "n_cells": len(idx),
                "frac_other_atlas_neighbours": float(other),
                "expected_if_random": float(1 - (frac**2).sum()),
                "n_cells_hlca": int((atlas[idx] == "HLCA").sum()),
                "n_cells_luca": int((atlas[idx] == "LuCA").sum()),
                "scope": scope,
            }
        )
    out = pd.DataFrame(rows)
    return out.sort_values("n_cells", ascending=False) if len(out) else out
