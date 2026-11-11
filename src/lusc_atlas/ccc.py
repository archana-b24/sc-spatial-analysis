"""Cell–cell communication (DEVELOPMENT_LOG §3.9) — LIANA+ (Dimitrov 2024).

v1: `rank_aggregate` on tumour cells only, pooled across donors (S11) → no way to say what *changes*.
v2:
1. `rank_aggregate.by_sample` in **every sample of both conditions** (consensus of CellPhoneDB,
   CellChat, NATMI, Connectome, log2FC, SingleCellSignalR).
2. **Differential CCC**: ligand/receptor-level tumour-vs-normal statistics from the pseudobulk DESeq2
   results (`li.mt.df_to_lr`) → interactions whose *both* partners change.
3. **Tensor-cell2cell** (Armingol 2022): sample × interaction × sender × receiver tensor decomposition →
   communication programmes and their association with condition.
4. **Spatial**: bivariate local co-expression of ligand–receptor pairs on Visium (`li.mt.bivariate`).
"""

from __future__ import annotations

import logging

import anndata as ad
import pandas as pd

log = logging.getLogger(__name__)


def by_sample(adata: ad.AnnData, ccfg: dict, sample_key: str, n_jobs: int = 16) -> pd.DataFrame:
    import liana as li

    li.mt.rank_aggregate.by_sample(
        adata,
        groupby=ccfg["groupby"],
        sample_key=sample_key,
        resource_name=ccfg["resource"],
        expr_prop=ccfg["expr_prop"],
        min_cells=ccfg["min_cells"],
        use_raw=False,
        layer="lognorm",
        n_perms=100,
        verbose=False,
        key_added="liana_res",
        n_jobs=n_jobs,
    )
    return adata.uns["liana_res"]


def differential(adata: ad.AnnData, de: pd.DataFrame, ccfg: dict) -> pd.DataFrame:
    """Combine per-cell-type DESeq2 stats into ligand–receptor differential scores."""
    import liana as li

    # One contrast only: each gene×group must appear once, and mixing designs (paired vs unpaired)
    # inside one ligand–receptor score would be statistically incoherent.
    if de.duplicated(["gene", "group"]).any():
        raise ValueError("DE table has duplicate gene×group rows — filter to a single contrast first")
    if de.empty:
        raise ValueError("no DE rows for the configured ccc.de_contrast")
    # LIANA 2.x: dea_df is indexed by gene (must match var_names) and carries the `groupby` column.
    dea_df = (
        de[["gene", "group", "stat", "pvalue", "padj", "log2FoldChange"]]
        .rename(columns={"group": ccfg["groupby"]})
        .set_index("gene")
    )
    dea_df = dea_df[dea_df[ccfg["groupby"]].isin(adata.obs[ccfg["groupby"]].unique())]
    lr = li.mt.df_to_lr(
        adata,
        dea_df=dea_df,
        resource_name=ccfg["resource"],
        expr_prop=ccfg["expr_prop"],
        min_cells=ccfg["min_cells"],
        groupby=ccfg["groupby"],
        stat_keys=["stat", "pvalue", "padj"],
        use_raw=False,
        layer="lognorm",
        complex_col="stat",
        verbose=False,
        return_all_lrs=False,
    )
    return lr.sort_values("interaction_stat", ascending=False)


def tensor_c2c(adata: ad.AnnData, liana_res: pd.DataFrame, sample_key: str, ccfg: dict, seed: int = 0):
    import cell2cell as c2c
    import liana as li

    # Scope: tumour-vs-normal samples only (healthy donors would dominate the factors), and drop
    # non-informative labels. how="outer" keeps cells AND LR pairs present in ≥ outer_fraction of
    # samples (missing → 0); "outer_cells" intersects LR pairs across all samples, which collapsed
    # the v2.0 run (215 samples → 1 LR pair).
    conds = ccfg["tensor_conditions"]
    keep = adata.obs.loc[adata.obs["condition"].isin(conds), sample_key].unique()
    excl = set(ccfg.get("tensor_exclude_labels", []))
    liana_res = liana_res[
        liana_res[sample_key].isin(keep) & ~liana_res["source"].isin(excl) & ~liana_res["target"].isin(excl)
    ]
    log.info("tensor-c2c: %d samples (%s)", liana_res[sample_key].nunique(), ", ".join(conds))
    tensor = li.ms.to_tensor_c2c(
        liana_res=liana_res,
        sample_key=sample_key,
        source_key="source",
        target_key="target",
        ligand_key="ligand_complex",
        receptor_key="receptor_complex",
        score_key="magnitude_rank",
        inverse_fn=lambda x: 1 - x,  # magnitude_rank ∈ [0,1], lower = stronger (v1 convention)
        how="outer",
        outer_fraction=float(ccfg["tensor_outer_fraction"]),
    )
    dims = dict(zip(["samples", "lr_pairs", "senders", "receivers"], tensor.tensor.shape, strict=True))
    log.info("tensor-c2c dims: %s", dims)
    if dims["lr_pairs"] < 50 or dims["senders"] < 3 or dims["samples"] < 10:
        raise ValueError(f"degenerate communication tensor {dims} — check sample scope / outer_fraction")
    meta = adata.obs[[sample_key, "condition"]].drop_duplicates().set_index(sample_key)
    ctx = c2c.tensor.generate_tensor_metadata(
        interaction_tensor=tensor,
        metadata_dicts=[meta["condition"].to_dict(), None, None, None],
        fill_with_order_elements=True,
    )
    rank = ccfg["tensor_rank"]
    tensor = c2c.analysis.run_tensor_cell2cell_pipeline(
        tensor,
        ctx,
        rank=None if rank == "auto" else int(rank),
        tf_optimization="regular",
        random_state=seed,
        device="cpu",
        elbow_metric="error",
        smooth_elbow=False,
        upper_rank=int(ccfg.get("tensor_upper_rank", 8)),
        tf_init="random",
        tf_svd="numpy_svd",
        cmaps=None,
        sample_col="Element",
        group_col="Category",
        output_fig=False,
    )
    return {k: v for k, v in tensor.factors.items()}


def spatial_bivariate(vis: ad.AnnData, ccfg: dict) -> pd.DataFrame:
    import liana as li

    li.pp.spatial_neighbors(  # LIANA 2.x: moved from li.ut → li.pp
        vis, bandwidth=ccfg["spatial_bandwidth"], cutoff=0.1, kernel="gaussian", set_diag=True
    )
    lrdata = li.mt.bivariate(
        vis,
        resource_name=ccfg["resource"],
        local_name="cosine",
        global_name="morans",
        n_perms=100,
        mask_negatives=True,
        add_categories=True,
        nz_prop=0.2,
        use_raw=False,
        layer="lognorm" if "lognorm" in vis.layers else None,
        verbose=False,
    )
    return lrdata.var.sort_values("morans", ascending=False).rename_axis("interaction").reset_index()
