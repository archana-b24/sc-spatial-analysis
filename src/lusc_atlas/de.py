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


def make_de_groups(
    adata: ad.AnnData,
    coarse_key: str,
    cnv_low_label: str | None = None,
    author_label: str = "malignant cell",
    tumour_conditions: tuple[str, ...] | list[str] = ("tumour", "tumour_region_unknown"),
    author_key: str = "cell_type",
    unknown_label: str = "Unknown",
) -> None:
    """Coarse cell-type groups for pseudobulk, CCC, scCODA and the cell2location reference.

    * CNV-called malignant cells → ``Malignant`` (definition unchanged since run 1).
    * If ``cnv_low_label`` is given (§3.10 #2): cells the LuCA authors labelled ``author_label`` that sit in a
      tumour sample but did **not** pass the CNV threshold → ``cnv_low_label``. In run 1 these ~25k cells fell
      into ``Unknown`` (the crosswalk maps the author label to ``Unknown`` on purpose, so it never steers
      integration), which then mixed un-called tumour cells with truly unmapped cells.
    * Everything else keeps its harmonised label; ``Unknown`` now means "no usable label" only.
    """
    g = adata.obs[coarse_key].astype(str)
    if cnv_low_label:
        author = adata.obs[author_key].astype(str) == author_label if author_key in adata.obs else False
        cnv_low = (
            (g == unknown_label)
            & author
            & adata.obs["condition"].astype(str).isin(list(tumour_conditions))
            & ~adata.obs.get("malignant", pd.Series(False, index=adata.obs_names)).astype(bool)
        )
        g = g.where(~cnv_low, cnv_low_label)
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
    tumour_only = set(dcfg.get("tumour_only_groups", ["Malignant"]))
    for grp in pdata.obs[gkey].cat.categories:
        if grp in tumour_only:
            continue  # no counterpart in control tissue → compared with healthy basal instead
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
        out.append(
            res.assign(
                contrast=name,
                group=grp,
                design=design,
                n_case=int(n.get(case, 0)),
                n_control=int(n.get(ctrl, 0)),
                n_donors=int(sub.obs["donor_id"].nunique()),
            )
        )
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def tumour_epithelium_vs_healthy_basal(
    pdata: ad.AnnData,
    dcfg: dict,
    group: str = "Malignant",
    basal_label: str = "Basal",
    contrast_name: str = "malignant_vs_healthy_basal",
    n_jobs: int = 8,
) -> pd.DataFrame:
    """Tumour-epithelium group (``Malignant`` or the CNV-low group) vs healthy-donor basal pseudobulks.

    Unpaired by necessity (healthy donors have no tumour); study/assay added when not collinear.
    """
    import decoupler as dc
    from pydeseq2.dds import DefaultInference, DeseqDataSet
    from pydeseq2.ds import DeseqStats

    g, cond = dcfg["group_key"], dcfg["condition_key"]
    mask = (pdata.obs[g] == group) | ((pdata.obs[g] == basal_label) & (pdata.obs[cond] == "healthy"))
    sub = pdata[mask].copy()
    sub.obs["epi_state"] = pd.Categorical(
        np.where(sub.obs[g] == group, "tumour_epi", "healthy_basal"),
        categories=["healthy_basal", "tumour_epi"],
    )
    n = sub.obs["epi_state"].value_counts()
    if (n == 0).any() or n.min() < dcfg["min_samples_per_condition"]:
        log.warning("%s: insufficient replicates %s", contrast_name, dict(n))
        return pd.DataFrame()
    for c in ("study", "assay"):
        if c in sub.obs:
            sub.obs[c] = sub.obs[c].astype(str).astype("category")
    dc.pp.filter_by_expr(sub, group="epi_state", min_count=10, min_total_count=15)
    design = design_for(sub.obs, "epi_state", paired=False)
    inference = DefaultInference(n_cpus=n_jobs)
    dds = DeseqDataSet(adata=sub, design=design, inference=inference, quiet=True)
    dds.deseq2()
    st = DeseqStats(
        dds, contrast=["epi_state", "tumour_epi", "healthy_basal"], inference=inference, quiet=True
    )
    st.summary()
    return (
        st.results_df.rename_axis("gene")
        .reset_index()
        .assign(
            contrast=contrast_name,
            group=group,
            design=design,
            n_case=int(n["tumour_epi"]),
            n_control=int(n["healthy_basal"]),
        )
    )


def malignant_vs_healthy_basal(
    pdata: ad.AnnData, dcfg: dict, basal_label: str = "Basal", n_jobs: int = 8
) -> pd.DataFrame:
    """The 'progression' contrast (CNV-called malignant vs healthy basal); kept for callers/tests."""
    return tumour_epithelium_vs_healthy_basal(
        pdata, dcfg, "Malignant", basal_label, "malignant_vs_healthy_basal", n_jobs
    )


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
    adata: ad.AnnData,
    cell_type_key: str,
    sample_key: str,
    condition_key: str,
    seed: int = 0,
    reference: str = "automatic",
    num_samples: int = 10000,
    num_warmup: int = 1000,
):
    """Compositional differential abundance (scCODA via pertpy; Büttner 2021).

    ``reference`` is set explicitly in config (§3.10 #2): run 1's "automatic" choice was `Unknown`, so every
    effect was relative to a bucket of un-called tumour cells.
    """
    import pertpy as pt

    if reference != "automatic" and reference not in set(adata.obs[cell_type_key].astype(str)):
        raise ValueError(f"scCODA reference {reference!r} is not a level of {cell_type_key}")
    sccoda = pt.tl.Sccoda()
    mdata = sccoda.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=cell_type_key,
        sample_identifier=sample_key,
        covariate_obs=[condition_key],
    )
    mdata = sccoda.prepare(mdata, modality_key="coda", formula=condition_key, reference_cell_type=reference)
    sccoda.run_nuts(mdata, modality_key="coda", rng_key=seed, num_samples=num_samples, num_warmup=num_warmup)
    eff = sccoda.get_effect_df(mdata, modality_key="coda")
    return eff.reset_index().assign(reference=reference)


def milo_design(obs: pd.DataFrame, condition_key: str, spec: dict) -> tuple[pd.Series, str, str]:
    """Which cells Milo tests, the design formula and the explicit contrast (§3.10 #1).

    Run 1 called ``da_nhoods(design="~condition")`` without ``model_contrasts``; pertpy then tests the *last*
    level (``tumour_region_unknown``, 2 HLCA donors) against the first. Here only case + control samples are
    kept; with ``paired`` only donors that have both, and donor enters the design, so donor/study effects
    cancel exactly as in the primary pseudobulk DE.

    Returns (boolean cell mask, design, model_contrasts).
    """
    case, ctrl = spec["case"], spec["control"]
    cond = obs[condition_key].astype(str)
    keep = cond.isin([case, ctrl])
    if spec.get("paired", False):
        both = obs[keep].groupby("donor_id", observed=True)[condition_key].nunique()
        keep &= obs["donor_id"].isin(both[both == 2].index)
        design = f"~donor_id + {condition_key}"
    else:
        design = f"~{condition_key}"
    if not keep.any():
        raise ValueError(f"Milo: no cells for {case} vs {ctrl}")
    return keep, design, f"{condition_key}{case} - {condition_key}{ctrl}"


def milo_da(
    adata: ad.AnnData, acfg: dict, sample_key: str, condition_key: str, seed: int = 0
) -> pd.DataFrame:
    """Neighbourhood-level differential abundance (Milo via pertpy; Dann 2022), explicit contrast."""
    import pertpy as pt
    import scanpy as sc

    spec = acfg["milo_contrast"]
    keep, design, contrast = milo_design(adata.obs, condition_key, spec)
    sub = adata[keep.values].copy()
    sub.obs[condition_key] = pd.Categorical(
        sub.obs[condition_key].astype(str), categories=[spec["control"], spec["case"]]
    )
    sub.obs["donor_id"] = sub.obs["donor_id"].astype(str).astype("category")
    log.info(
        "Milo: %d cells, %d samples, %d donors; design %s; contrast %s",
        sub.n_obs,
        sub.obs[sample_key].nunique(),
        sub.obs["donor_id"].nunique(),
        design,
        contrast,
    )
    milo = pt.tl.Milo()
    mdata = milo.load(sub)
    sc.pp.neighbors(mdata["rna"], use_rep="X_integrated", n_neighbors=acfg["milo_k"], random_state=seed)
    milo.make_nhoods(mdata["rna"], prop=acfg["milo_prop"])
    mdata = milo.count_nhoods(mdata, sample_col=sample_key)
    milo.da_nhoods(mdata, design=design, model_contrasts=contrast)
    milo.annotate_nhoods(mdata, anno_col="de_group")
    return (
        mdata["milo"]
        .var.reset_index()
        .assign(
            design=design,
            contrast=contrast,
            n_samples=int(sub.obs[sample_key].nunique()),
            n_donors=int(sub.obs["donor_id"].nunique()),
        )
    )
