"""Differential expression, pathway/TF activity and compositional change (DEVELOPMENT_LOG §3.6).

v1 → v2:
* v1 tested tumour vs normal with cell-level Wilcoxon/t-tests and ran GSEA with cells as samples
  (S3, S4). Cells from one donor are not independent, so this is pseudoreplication: p-values are
  dramatically anti-conservative (Squair 2021, Nat Commun; Zimmerman 2021).
* v2 aggregates counts per **sample × cell type** (pseudobulk), fits **PyDESeq2** (negative binomial,
  donors as replicates), and infers pathway / TF activity with **decoupler 2.x** on the Wald statistic:
  PROGENy (pathways), CollecTRI (TF regulons), MSigDB Hallmark (gene-set enrichment, preranked).
* Study is added to the design when it is not collinear with condition (v1 ignored it — S5).
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def make_de_groups(adata: ad.AnnData, coarse_key: str) -> None:
    """Coarse cell-type groups for pseudobulk; malignant cells form their own group."""
    g = adata.obs[coarse_key].astype(str)
    if "malignant" in adata.obs:
        g = g.where(~adata.obs["malignant"].astype(bool), "Malignant")
    adata.obs["de_group"] = pd.Categorical(g)


def pseudobulk(adata: ad.AnnData, dcfg: dict) -> ad.AnnData:
    import decoupler as dc

    pdata = dc.pp.pseudobulk(
        adata,
        sample_col=dcfg["pseudobulk_sample_key"],
        groups_col=dcfg["group_key"],
        layer="counts",
        mode="sum",
    )
    dc.pp.filter_samples(pdata, min_cells=dcfg["min_cells_per_pseudobulk"], min_counts=1000)
    return pdata


def design_for(meta: pd.DataFrame, cond: str, paired: bool) -> str:
    """DESeq2 design for one contrast.

    paired   → ``~donor_id + condition``: every donor is its own control, so donor, study, assay and
               atlas effects cancel exactly (they are constant within donor).
    unpaired → ``~[study +] [assay +] condition``, adding a covariate only when it is *not* collinear
               with condition (a covariate constant within each condition makes the model
               unidentifiable — exactly v1's study/condition confounding, S5).
    """
    if paired:
        return f"~donor_id + {cond}"
    terms = []
    for cov in ("study", "assay"):
        if cov in meta and meta[cov].nunique() > 1:
            if (meta.groupby(cov, observed=True)[cond].nunique() > 1).any():
                terms.append(cov)
    return "~" + " + ".join(terms + [cond])


def _design(meta: pd.DataFrame, dcfg: dict) -> str:  # thin wrapper used by tests
    return design_for(meta, dcfg["condition_key"], paired=False)


def run_contrast(pdata: ad.AnnData, dcfg: dict, name: str, spec: dict, n_jobs: int = 8) -> pd.DataFrame:
    """One contrast (case vs control) in every cell-type group with enough replicates."""
    import decoupler as dc
    from pydeseq2.dds import DefaultInference, DeseqDataSet
    from pydeseq2.ds import DeseqStats

    cond, gkey = dcfg["condition_key"], dcfg["group_key"]
    case, ctrl, paired = spec["case"], spec["control"], spec.get("paired", False)
    inference = DefaultInference(n_cpus=n_jobs)
    out = []
    for grp in pdata.obs[gkey].cat.categories:
        if grp == "Malignant":
            continue  # malignant cells have no counterpart in control tissue → own contrast
        sub = pdata[(pdata.obs[gkey] == grp) & pdata.obs[cond].isin([case, ctrl])].copy()
        if paired:  # keep donors contributing both conditions for this cell type
            both = sub.obs.groupby("donor_id", observed=True)[cond].nunique()
            sub = sub[sub.obs["donor_id"].isin(both[both == 2].index)].copy()
        n = sub.obs[cond].value_counts()
        mn = dcfg["min_samples_per_condition"]
        if n.get(case, 0) < mn or n.get(ctrl, 0) < mn:
            log.info("[%s] skip %s (replicates: %s)", name, grp, dict(n))
            continue
        sub.obs[cond] = pd.Categorical(sub.obs[cond].astype(str), categories=[ctrl, case])
        for c in ("donor_id", "study", "assay"):
            if c in sub.obs:
                sub.obs[c] = sub.obs[c].astype(str).astype("category")
        dc.pp.filter_by_expr(sub, group=cond, min_count=10, min_total_count=15)
        if sub.n_vars < 100:
            continue
        design = design_for(sub.obs, cond, paired)
        try:
            dds = DeseqDataSet(adata=sub, design=design, refit_cooks=True, inference=inference, quiet=True)
            dds.deseq2()
            st = DeseqStats(dds, contrast=[cond, case, ctrl], inference=inference, quiet=True)
            st.summary()
        except Exception as exc:  # noqa: BLE001 — e.g. rank-deficient design in one group
            log.warning("[%s] %s failed: %s", name, grp, exc)
            continue
        res = st.results_df.rename_axis("gene").reset_index()
        out.append(res.assign(contrast=name, group=grp, design=design, n_case=int(n.get(case, 0)),
                              n_control=int(n.get(ctrl, 0)), n_donors=int(sub.obs["donor_id"].nunique())))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def malignant_vs_healthy_basal(pdata: ad.AnnData, dcfg: dict, basal_label: str = "Basal", n_jobs: int = 8) -> pd.DataFrame:
    """The 'progression' contrast: malignant pseudobulks vs healthy-donor basal pseudobulks.

    Unpaired by necessity (healthy donors have no tumour); study/assay added when not collinear.
    """
    import decoupler as dc
    from pydeseq2.dds import DefaultInference, DeseqDataSet
    from pydeseq2.ds import DeseqStats

    g, cond = dcfg["group_key"], dcfg["condition_key"]
    mask = (pdata.obs[g] == "Malignant") | ((pdata.obs[g] == basal_label) & (pdata.obs[cond] == "healthy"))
    sub = pdata[mask].copy()
    sub.obs["epi_state"] = pd.Categorical(
        np.where(sub.obs[g] == "Malignant", "malignant", "healthy_basal"), categories=["healthy_basal", "malignant"]
    )
    n = sub.obs["epi_state"].value_counts()
    if (n == 0).any() or n.min() < dcfg["min_samples_per_condition"]:
        log.warning("malignant vs healthy basal: insufficient replicates %s", dict(n))
        return pd.DataFrame()
    for c in ("study", "assay"):
        if c in sub.obs:
            sub.obs[c] = sub.obs[c].astype(str).astype("category")
    dc.pp.filter_by_expr(sub, group="epi_state", min_count=10, min_total_count=15)
    design = design_for(sub.obs, "epi_state", paired=False)
    inference = DefaultInference(n_cpus=n_jobs)
    dds = DeseqDataSet(adata=sub, design=design, inference=inference, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["epi_state", "malignant", "healthy_basal"], inference=inference, quiet=True)
    st.summary()
    return st.results_df.rename_axis("gene").reset_index().assign(
        contrast="malignant_vs_healthy_basal", group="Malignant", design=design,
        n_case=int(n["malignant"]), n_control=int(n["healthy_basal"]))


def pathway_activities(de: pd.DataFrame, pcfg: dict) -> dict[str, pd.DataFrame]:
    """PROGENy / CollecTRI (ULM) and Hallmark (preranked GSEA) on the DESeq2 Wald statistic."""
    import decoupler as dc

    de = de.assign(key=de["contrast"] + "|" + de["group"])
    mat = de.pivot_table(index="key", columns="gene", values="stat").fillna(0)
    out = {}
    progeny = dc.op.progeny(organism="human", top=pcfg["progeny_top"])
    es, pv = dc.mt.ulm(data=mat, net=progeny)
    out["progeny"] = _long(es, pv)
    if pcfg["collectri"]:
        collectri = dc.op.collectri(organism="human")
        es, pv = dc.mt.ulm(data=mat, net=collectri)
        out["collectri"] = _long(es, pv)
    if pcfg["hallmark"]:
        hallmark = dc.op.hallmark(organism="human")
        es, pv = dc.mt.gsea(data=mat, net=hallmark)
        out["hallmark"] = _long(es, pv)
    return out


def _long(es: pd.DataFrame, pv: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import false_discovery_control

    long = es.stack().rename("score").to_frame().join(pv.stack().rename("pval"))
    long = long.rename_axis(["key", "source"]).reset_index()
    long[["contrast", "group"]] = long["key"].str.split("|", n=1, expand=True)
    long["padj"] = long.groupby("key")["pval"].transform(lambda p: false_discovery_control(p.fillna(1)))
    return long


def composition_sccoda(
    adata: ad.AnnData, cell_type_key: str, sample_key: str, condition_key: str, seed: int = 0
):
    """Compositional differential abundance (scCODA via pertpy; Büttner 2021)."""
    import pertpy as pt

    sccoda = pt.tl.Sccoda()
    mdata = sccoda.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=cell_type_key,
        sample_identifier=sample_key,
        covariate_obs=[condition_key],
    )
    mdata = sccoda.prepare(mdata, modality_key="coda", formula=condition_key, reference_cell_type="automatic")
    sccoda.run_nuts(mdata, modality_key="coda", rng_key=seed)
    eff = sccoda.get_effect_df(mdata, modality_key="coda")
    return eff.reset_index()


def milo_da(
    adata: ad.AnnData, acfg: dict, sample_key: str, condition_key: str, seed: int = 0
) -> pd.DataFrame:
    """Neighbourhood-level differential abundance (Milo via pertpy; Dann 2022)."""
    import pertpy as pt
    import scanpy as sc

    milo = pt.tl.Milo()
    mdata = milo.load(adata)
    sc.pp.neighbors(mdata["rna"], use_rep="X_integrated", n_neighbors=acfg["milo_k"], random_state=seed)
    milo.make_nhoods(mdata["rna"], prop=acfg["milo_prop"])
    mdata = milo.count_nhoods(mdata, sample_col=sample_key)
    milo.da_nhoods(mdata, design=f"~{condition_key}")
    milo.annotate_nhoods(mdata, anno_col="de_group")
    return mdata["milo"].var.reset_index()
