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
    ckey = "tissue_type" if "tissue_type" in adata.obs else "condition"
    conds = ccfg["tensor_conditions"]
    keep = adata.obs.loc[adata.obs[ckey].astype(str).isin(conds), sample_key].unique()
    excl = set(ccfg.get("tensor_exclude_labels", []))
    liana_res = liana_res[
        liana_res[sample_key].isin(keep) & ~liana_res["source"].isin(excl) & ~liana_res["target"].isin(excl)
    ]
    liana_res = balanced_tensor_scope(liana_res, adata.obs, sample_key, ckey, conds, ccfg)
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
    meta = adata.obs[[sample_key, ckey]].drop_duplicates().set_index(sample_key)
    meta = meta[~meta.index.duplicated()]
    ctx = c2c.tensor.generate_tensor_metadata(
        interaction_tensor=tensor,
        metadata_dicts=[meta[ckey].astype(str).to_dict(), None, None, None],
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


def balanced_tensor_scope(
    liana_res: pd.DataFrame, obs: pd.DataFrame, sample_key: str, ckey: str, conds: list[str], ccfg: dict
) -> pd.DataFrame:
    """Keep only cell types present in BOTH conditions, then drop slices (not types) that lack one (§3.12).

    The tensor fills missing entries with 0. A type that exists only in tumour slices (tumour epithelium)
    would then separate the conditions *by construction*, and a type missing from a study's slices because
    of its
    sorting protocol (CD31-enriched Goveia lacks most immune types) would make the factors track protocol. A
    type is kept if it sends or receives in ≥ ``tensor_type_presence`` of the slices of every condition; a
    slice missing any kept type is dropped together with its paired slice (same donor)."""
    frac = float(ccfg.get("tensor_type_presence", 0.8))
    slice_cond = obs[[sample_key, ckey]].drop_duplicates().set_index(sample_key)[ckey].astype(str)
    slice_cond = slice_cond[~slice_cond.index.duplicated()]
    pres = (
        pd.concat(
            [
                liana_res[[sample_key, "source"]].set_axis([sample_key, "ct"], axis=1),
                liana_res[[sample_key, "target"]].set_axis([sample_key, "ct"], axis=1),
            ]
        )
        .drop_duplicates()
        .assign(present=True)
        .pivot_table(index=sample_key, columns="ct", values="present", aggfunc="any", fill_value=False)
    )
    by_cond = pres.groupby(slice_cond.reindex(pres.index).values).mean()
    types = [t for t in pres.columns if all(by_cond.loc[c, t] >= frac for c in conds if c in by_cond.index)]
    donor = pd.Series({s: str(s).rsplit("|", 1)[0] for s in pres.index})

    def complete_donors(ts: list[str]) -> set[str]:
        complete = pres[ts].all(axis=1)
        n_ok = donor[complete[complete].index].value_counts()
        return set(n_ok[n_ok == len(conds)].index)

    # Greedy trade-off between cell types and patients: while fewer than `tensor_min_donor_fraction` of the
    # paired patients have every kept type in both tissues, drop the type whose removal frees most patients.
    need = float(ccfg.get("tensor_min_donor_fraction", 0.75)) * donor.nunique()
    dropped = []
    while len(complete_donors(types)) < need and len(types) > 4:
        best = max(types, key=lambda t: len(complete_donors([x for x in types if x != t])))
        types = [x for x in types if x != best]
        dropped.append(best)
    if dropped:
        log.info("tensor scope: dropped %s to keep ≥ %.0f patients complete", dropped, need)
    ok_donors = complete_donors(types)
    ok_slices = set(donor[donor.isin(ok_donors)].index)
    lost = sorted(set(donor) - ok_donors)
    if lost:
        log.info("tensor scope: patients dropped (a kept cell type missing in one tissue): %s", lost)
    if len(ok_donors) < int(ccfg.get("tensor_min_donors", 5)):
        raise ValueError(f"tensor scope left {len(ok_donors)} complete patients — too few for a paired test")
    keep_slices = [s for s in ok_slices if donor[s] in ok_donors]
    log.info(
        "tensor scope: %d cell types present in ≥ %.0f %% of slices of every condition (%s); %d / %d slices "
        "kept "
        "(%d donors)",
        len(types),
        100 * frac,
        ", ".join(types),
        len(keep_slices),
        len(pres),
        len(ok_donors),
    )
    return liana_res[
        liana_res[sample_key].isin(keep_slices)
        & liana_res["source"].isin(types)
        & liana_res["target"].isin(types)
    ]


def paired_units(obs: pd.DataFrame, case: str, control: str, unit_key: str = "donor_condition") -> pd.Series:
    """Cells of donors with BOTH conditions → one unit per donor × condition (§3.12).

    Multi-region tumour samples of one patient are pooled into that patient's tumour unit, so each patient
    counts once per tissue (runs 1–2 counted Leader_Merad_706's 16 samples as 16 replicates). Returns the
    unit label per
    cell (NaN for cells outside the paired design)."""
    col = (
        "tissue_type" if "tissue_type" in obs else "condition"
    )  # §3.12: region-unknown tumour counts as tumour tissue
    cond = obs[col].astype(str)
    inc = cond.isin([case, control])
    both = obs[inc].groupby("donor_id", observed=True)[col].nunique()
    paired = inc & obs["donor_id"].isin(both[both == 2].index)
    unit = obs["donor_id"].astype(str) + "|" + cond
    return unit.where(paired)


def paired_loading_test(factors: dict, case: str, control: str) -> pd.DataFrame:
    """Per factor: Wilcoxon signed-rank test of tumour vs adjacent sample loadings WITHIN donor (§3.12),
    BH across
    factors. Run 1–2 used an unpaired Mann–Whitney on samples, which ignores that both slices come from
    one patient."""
    from scipy.stats import false_discovery_control, wilcoxon

    ctx = factors["Contexts"] if "Contexts" in factors else next(iter(factors.values()))
    idx = ctx.index.astype(str)
    donor = idx.str.rsplit("|", n=1).str[0]
    cond = idx.str.rsplit("|", n=1).str[1]
    rows = []
    for f in ctx.columns:
        v = pd.DataFrame({"donor": donor, "cond": cond, "x": ctx[f].to_numpy()}).pivot_table(
            index="donor", columns="cond", values="x"
        )
        v = v.dropna(subset=[case, control]) if {case, control} <= set(v.columns) else v.iloc[0:0]
        if len(v) < 5:
            continue
        d = v[case] - v[control]
        stat, p = wilcoxon(v[case], v[control])
        rows.append(
            {
                "factor": f,
                "n_donors": len(v),
                "median_tumour": float(v[case].median()),
                "median_adjacent": float(v[control].median()),
                "median_paired_difference": float(d.median()),
                "n_donors_higher_in_tumour": int((d > 0).sum()),
                "wilcoxon_stat": float(stat),
                "pvalue": float(p),
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out["padj"] = false_discovery_control(out["pvalue"].to_numpy())
    return out


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
