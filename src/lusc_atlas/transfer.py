"""Transfer learning: map the cohort onto the published HLCA scANVI reference with scArches surgery.

Integration-benchmark candidate ``X_scArches`` (DEVELOPMENT_LOG §3.3 G).

Question it answers: for a *tumour* cohort, is it better to (a) integrate de novo (scVI/scANVI/sysVI) or
(b) re-use a reference model trained on healthy lung and fine-tune only the query-batch parameters
(Lotfollahi 2022, scArches)? (b) is the standard HLCA workflow and is fast, but tumour states are
out-of-distribution for a healthy reference — the benchmark + over-correction guard decide.

Mechanics
* The HLCA model (Zenodo 7599104; Sikkema 2023) was saved by scvi-tools < 0.15. scvi-tools' own converter
  unpickles without restrictions, so we re-implement it with :mod:`lusc_atlas.safe_load` (sha256-pinned,
  whitelisted types only) and ``registry_from_setup_dict``.
* Query = every cohort cell (raw counts, Ensembl IDs); query batch = our ``sample`` (new batch categories
  are added by surgery; reference weights stay frozen).
* Reference genes missing from the query are zero-filled by ``prepare_query_anndata``; the fraction is
  reported because heavy zero-filling degrades mapping.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .safe_load import load_pickle, load_torch_weights, verify

log = logging.getLogger(__name__)


def convert_reference(tcfg: dict, repo: Path) -> Path:
    """Legacy HLCA scANVI save → scvi-tools ≥ 1.0 format, using only restricted loaders."""
    import torch
    from scvi.data._compat import registry_from_setup_dict
    from scvi.model import SCANVI

    src = repo / tcfg["reference_dir"]
    out = repo / tcfg["converted_dir"]
    if (out / "model.pt").exists():
        return out
    sums = tcfg["sha256"]
    state_dict = load_torch_weights(src / "model_params.pt", sums["model_params.pt"])
    attr = load_pickle(src / "attr.pkl", sums["attr.pkl"])
    verify(src / "var_names.csv", sums["var_names.csv"])
    var_names = np.genfromtxt(src / "var_names.csv", delimiter=",", dtype=str)
    setup = attr.pop("scvi_setup_dict_")
    attr["registry_"] = registry_from_setup_dict(
        SCANVI, setup, unlabeled_category=attr.get("unlabeled_category_")
    )
    out.mkdir(parents=True, exist_ok=True)
    torch.save({"model_state_dict": state_dict, "var_names": var_names, "attr_dict": attr}, out / "model.pt")
    log.info("converted HLCA reference → %s (%d genes)", out, len(var_names))
    return out


def map_to_reference(adata: ad.AnnData, tcfg: dict, repo: Path, seed: int = 0) -> dict:
    """scArches surgery; writes ``X_scArches`` and ``hlca_ref_prediction`` into ``adata``."""
    import scvi
    from scvi.model import SCANVI

    scvi.settings.seed = seed
    ref_dir = convert_reference(tcfg, repo)
    q = ad.AnnData(
        X=adata.layers["counts"].copy(),
        obs=pd.DataFrame(
            {"dataset": adata.obs["sample"].astype(str).values, "scanvi_label": "unlabeled"},
            index=adata.obs_names,
        ),
        var=pd.DataFrame(index=adata.var["ensembl_id"].astype(str).values),
    )
    q.var_names_make_unique()
    SCANVI.prepare_query_anndata(q, str(ref_dir))
    n_ref = q.n_vars
    missing = int((np.asarray(q.X.sum(axis=0)).ravel() == 0).sum())
    model = SCANVI.load_query_data(q, str(ref_dir), freeze_dropout=True)
    model.train(
        max_epochs=tcfg["surgery_epochs"],
        plan_kwargs={"weight_decay": 0.0},
        early_stopping=True,
        check_val_every_n_epoch=1,
    )
    adata.obsm["X_scArches"] = model.get_latent_representation()
    adata.obs["hlca_ref_prediction"] = pd.Categorical(model.predict())
    return {
        "reference_genes": n_ref,
        "reference_genes_unobserved_in_query": missing,
        "frac_unobserved": round(missing / n_ref, 4),
        "surgery_epochs": tcfg["surgery_epochs"],
    }
