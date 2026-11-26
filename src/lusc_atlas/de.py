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


TUMOUR_EPITHELIUM = "Tumour epithelium"
CNV_UNRESOLVED = "Tumour epithelial (CNV unresolved)"


def make_de_groups(
    adata: ad.AnnData,
    coarse_key: str,
    cnv_low_label: str | None = None,
    author_label: str = "malignant cell",
    tumour_conditions: tuple[str, ...] | list[str] = ("tumour", "tumour_region_unknown"),
    author_key: str = "cell_type",
    unknown_label: str = "Unknown",
    refine: dict | None = None,
    fine_key: str = "fine_label",
) -> None:
    """Coarse cell-type groups for pseudobulk, CCC, composition and the cell2location reference.

    * CNV-called malignant cells → ``Malignant`` (definition unchanged since run 1).
    * If ``cnv_low_label`` is given (§3.10 #2): cells the LuCA authors labelled ``author_label`` that sit in a
      tumour sample but did **not** pass the CNV threshold → ``cnv_low_label``. In run 1 these ~25k cells fell
      into ``Unknown`` (the crosswalk maps the author label to ``Unknown`` on purpose, so it never steers
      integration), which then mixed un-called tumour cells with truly unmapped cells.
    * Everything else keeps its harmonised label; ``Unknown`` now means "no usable label" only.
    """
    g = adata.obs[coarse_key].astype(str)
    if refine and fine_key in adata.obs:
        # §3.12: split coarse labels whose subtype MIX differs between tumour and adjacent tissue (alveolar vs
        # monocyte-derived macrophages; regulatory vs conventional CD4 T). Without the split, paired DE in
        # "Macrophage" mostly measures the subtype shift, not a change of expression state.
        fine = adata.obs[fine_key].astype(str)
        for parent, subs in refine.items():
            for child, fines in subs.items():
                g = g.where(~((g == parent) & fine.isin(fines)), child)
    if cnv_low_label:
        author = adata.obs[author_key].astype(str) == author_label if author_key in adata.obs else False
        cnv_low = (
            (g == unknown_label)
            & author
            & adata.obs["condition"].astype(str).isin(list(tumour_conditions))
            & ~adata.obs.get("malignant", pd.Series(False, index=adata.obs_names)).astype(bool)
        )
        if "cnv_status" in adata.obs:
            # §3.12 (senior review B2): only cells that were actually scored and not called are "CNV-low";
            # author-malignant cells in samples too small to score get their own label
            status = adata.obs["cnv_status"].astype(str)
            unresolved = cnv_low & status.str.startswith("unresolved")
            cnv_low = cnv_low & (status == "not called")
            g = g.where(~unresolved, CNV_UNRESOLVED)
        g = g.where(~cnv_low, cnv_low_label)
    if "malignant" in adata.obs:
        g = g.where(~adata.obs["malignant"].astype(bool), "Malignant")
    adata.obs["de_group"] = pd.Categorical(g)
    # §3.12: for spatial deconvolution, CCC and composition the two tumour-epithelium groups are merged — they
    # differ by a CNV threshold, not by an expression programme a signature/communication model could
    # separate.
    merged = g.where(~g.isin(["Malignant", cnv_low_label or "Malignant", CNV_UNRESOLVED]), TUMOUR_EPITHELIUM)
    adata.obs["analysis_group"] = pd.Categorical(merged)


def pseudobulk(adata: ad.AnnData, dcfg: dict) -> ad.AnnData:
    import decoupler as dc

    pdata = dc.pp.pseudobulk(
        adata,
        sample_col=dcfg["pseudobulk_sample_key"],
        groups_col=dcfg["group_key"],
        layer="counts",
        mode="sum",
    )
    # decoupler drops obs columns that vary within a unit; donor × tissue units of two Leader_Merad patients
    # pool 10x 3' v2 and 5' v1 libraries, so `assay` vanished and unpaired designs / the cross-study null
    # silently lost it (run 3 bug). Restore design covariates as the sorted "+"-joined set of values.
    key = dcfg["pseudobulk_sample_key"]
    for col in ("assay", "study", "sex", "atlas"):
        if col in adata.obs and col not in pdata.obs:
            joined = (
                adata.obs[[key, col]].astype(str).groupby(key)[col].agg(lambda v: " + ".join(sorted(set(v))))
            )
            pdata.obs[col] = pdata.obs[key].astype(str).map(joined).to_numpy()
    dc.pp.filter_samples(pdata, min_cells=dcfg["min_cells_per_pseudobulk"], min_counts=1000)
    return pdata


def design_for(meta: pd.DataFrame, cond: str, paired: bool) -> str:
    """DESeq2 design for one contrast.

    paired   → ``~donor_id + condition``: every donor is its own control, so what is constant within a
               donor cancels (genotype, study, atlas, chemistry). Differences in how the two tissues were
               handled (dissociation, blood content, ambient RNA) do not (§3.13).
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


def add_pseudobulk_unit(adata: ad.AnnData, key: str = "donor_condition") -> None:
    """One pseudobulk unit per donor × condition (§3.12).

    Runs 1–2 pseudobulked per *sample*: a patient with 8 tumour regions contributed 8 "replicates" of the
    tumour condition (Leader_Merad_706 = 16 of 48 paired samples), inflating degrees of freedom and
    over-weighting that
    patient. Summing counts over a patient's samples of one tissue makes the patient the unit of
    replication."""
    col = (
        "tissue_type" if "tissue_type" in adata.obs else "condition"
    )  # tissue merges region-unknown tumour into tumour
    adata.obs[key] = pd.Categorical(adata.obs["donor_id"].astype(str) + "|" + adata.obs[col].astype(str))


def _fit(
    sub: ad.AnnData,
    cond: str,
    case: str,
    ctrl: str,
    design: str,
    n_jobs: int,
    lfc_null: float | None = None,
    shrink: bool = False,
) -> pd.DataFrame:
    """PyDESeq2 Wald test of case vs control. Adds 95 % CI (MLE ± 1.96 SE); optionally a thresholded test
    (H0: |log2FC| ≤ ``lfc_null``; columns ``pvalue_thr``/``padj_thr``) and apeGLM-style shrunken log2FC."""
    from pydeseq2.dds import DefaultInference, DeseqDataSet
    from pydeseq2.ds import DeseqStats

    inference = DefaultInference(n_cpus=n_jobs)
    dds = DeseqDataSet(adata=sub, design=design, refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=[cond, case, ctrl], inference=inference, quiet=True)
    st.summary()
    res = st.results_df.copy()
    res["ci_low"] = res["log2FoldChange"] - 1.96 * res["lfcSE"]
    res["ci_high"] = res["log2FoldChange"] + 1.96 * res["lfcSE"]
    if lfc_null:
        st2 = DeseqStats(
            dds,
            contrast=[cond, case, ctrl],
            inference=inference,
            quiet=True,
            lfc_null=lfc_null,
            alt_hypothesis="greaterAbs",
        )
        st2.summary()
        res["pvalue_thr"] = st2.results_df["pvalue"]
        res["padj_thr"] = st2.results_df["padj"]
    if shrink:
        coeff = [c for c in dds.varm["LFC"].columns if c.startswith(f"{cond}[")]
        if coeff:
            try:
                st.lfc_shrink(coeff=coeff[0])
                res["log2FoldChange_shrunk"] = st.results_df["log2FoldChange"]
            except Exception as exc:  # noqa: BLE001 — shrinkage is a display aid; never fails the contrast
                log.warning("lfc_shrink failed: %s", exc)
    disp = dds.var["dispersions"] if "dispersions" in dds.var else dds.varm["dispersions"]
    res["dispersion"] = np.asarray(disp).ravel()
    return res.rename_axis("gene").reset_index()


def _subset(pdata: ad.AnnData, dcfg: dict, grp: str, case: str, ctrl: str, paired: bool) -> ad.AnnData:
    cond, gkey = dcfg["condition_key"], dcfg["group_key"]
    sub = pdata[(pdata.obs[gkey] == grp) & pdata.obs[cond].isin([case, ctrl])].copy()
    if paired:  # keep donors contributing both conditions for this cell type
        both = sub.obs.groupby("donor_id", observed=True)[cond].nunique()
        sub = sub[sub.obs["donor_id"].isin(both[both == 2].index)].copy()
    return sub


def paired_eligibility(sub: ad.AnnData, dcfg: dict) -> tuple[str | None, str]:
    """§3.12 (senior review): paired groups need ≥ ``min_paired_donors`` pairs AND ≥
    ``min_paired_studies`` studies
    with ≥ 2 pairs each → "primary"; ≥ ``min_descriptive_donors`` pairs → "descriptive" (fitted and shown,
    never headlined; e.g. a type found in one study only); fewer → not fitted. Returns (tier or None,
    reason)."""
    nd = sub.obs["donor_id"].nunique()
    per_study = (
        sub.obs.groupby("study", observed=True)["donor_id"].nunique() if "study" in sub.obs else pd.Series()
    )
    n_studies_2 = int((per_study >= 2).sum())
    per = ", ".join(f"{k} {int(v)}" for k, v in per_study.items())
    why = f"{nd} pairs; studies with ≥ 2 pairs: {n_studies_2} ({per})"
    if nd >= dcfg.get("min_paired_donors", 6) and n_studies_2 >= dcfg.get("min_paired_studies", 2):
        return "primary", why
    if nd >= dcfg.get("min_descriptive_donors", 3):
        return "descriptive", why
    return None, why


def _prepare(sub: ad.AnnData, cond: str, case: str, ctrl: str) -> ad.AnnData:
    import decoupler as dc

    sub.obs[cond] = pd.Categorical(sub.obs[cond].astype(str), categories=[ctrl, case])
    for c in ("donor_id", "study", "assay"):
        if c in sub.obs:
            sub.obs[c] = sub.obs[c].astype(str).astype("category")
    dc.pp.filter_by_expr(sub, group=cond, min_count=10, min_total_count=15)
    return sub


def run_contrast(pdata: ad.AnnData, dcfg: dict, name: str, spec: dict, n_jobs: int = 8) -> pd.DataFrame:
    """One contrast (case vs control) in every cell-type group with enough replicates."""
    cond, gkey = dcfg["condition_key"], dcfg["group_key"]
    case, ctrl, paired = spec["case"], spec["control"], spec.get("paired", False)
    out = []
    skip = set(dcfg.get("tumour_only_groups", ["Malignant"])) | set(dcfg.get("exclude_groups", []))
    for grp in pdata.obs[gkey].cat.categories:
        if grp in skip:
            continue  # tumour-only (no counterpart in control tissue) or a mixed bucket (`Unknown`)
        sub = _subset(pdata, dcfg, grp, case, ctrl, paired)
        if paired:
            tier, why = paired_eligibility(sub, dcfg)
            inference = spec.get("inference", "primary") if tier == "primary" else f"descriptive ({why})"
        else:
            n = sub.obs[cond].value_counts()
            mn = dcfg["min_samples_per_condition"]
            tier = "ok" if (n.get(case, 0) >= mn and n.get(ctrl, 0) >= mn) else None
            why, inference = str(dict(n)), spec.get("inference", "exploratory")
        if tier is None:
            log.info("[%s] skip %s (%s)", name, grp, why)
            continue
        n = sub.obs[cond].value_counts()
        sub = _prepare(sub, cond, case, ctrl)
        if sub.n_vars < 100:
            continue
        design = design_for(sub.obs, cond, paired)
        try:
            res = _fit(
                sub,
                cond,
                case,
                ctrl,
                design,
                n_jobs,
                lfc_null=dcfg.get("lfc_threshold") if paired else None,
                shrink=paired,
            )
        except Exception as exc:  # noqa: BLE001 — e.g. rank-deficient design in one group
            log.warning("[%s] %s failed: %s", name, grp, exc)
            continue
        out.append(
            res.assign(
                contrast=name,
                group=grp,
                design=design,
                inference=inference,
                tier=tier,
                n_case=int(n.get(case, 0)),
                n_control=int(n.get(ctrl, 0)),
                n_donors=int(sub.obs["donor_id"].nunique()),
                n_studies=int(sub.obs["study"].nunique()) if "study" in sub.obs else np.nan,
                studies=";".join(sorted(sub.obs["study"].astype(str).unique())) if "study" in sub.obs else "",
            )
        )
    if not out:
        return pd.DataFrame()
    res = pd.concat(out, ignore_index=True)
    # §3.12: BH is per cell type (DESeq2 default, `padj`); also a GLOBAL BH over every gene × cell-type test
    # of the
    # contrast (`padj_global`) — the stricter scope, used for headline counts of the primary contrast.
    from scipy.stats import false_discovery_control

    # global scope over the primary tier only; also for the thresholded test that headlines use
    for pcol, out_col in (("pvalue", "padj_global"), ("pvalue_thr", "padj_thr_global")):
        res[out_col] = np.nan
        if pcol not in res:
            continue
        ok = res[pcol].notna() & (res["tier"].isin(["primary", "ok"]))
        if ok.any():
            res.loc[ok, out_col] = false_discovery_control(res.loc[ok, pcol].to_numpy())
    return res


def _study_design(obs: pd.DataFrame, cond: str, case: str) -> pd.DataFrame:
    """Design matrix ``~ donor + study:condition`` built explicitly (one tumour-effect column per study; donor
    intercepts absorb study main effects). Shared dispersions across studies — the joint model the senior
    review asked for instead of separate per-study fits on 1–3 donors."""
    donors = pd.get_dummies(obs["donor_id"].astype(str), drop_first=True, dtype=float)
    donors.columns = [f"donor_{c}" for c in donors.columns]
    is_case = (obs[cond].astype(str) == case).astype(float)
    st = pd.get_dummies(obs["study"].astype(str), dtype=float)
    eff = st.mul(is_case, axis=0)
    eff.columns = [f"tumour_in_{c}" for c in eff.columns]
    X = pd.concat([pd.Series(1.0, index=obs.index, name="Intercept"), donors, eff], axis=1)
    return X


def per_study_effects(sub: ad.AnnData, cond: str, case: str, n_jobs: int) -> pd.DataFrame:
    """Per-study log2FC ± SE from the joint model ``~ donor + study:condition`` (PyDESeq2, explicit
    design)."""
    from pydeseq2.dds import DefaultInference, DeseqDataSet
    from pydeseq2.ds import DeseqStats

    X = _study_design(sub.obs, cond, case)
    # §3.12 (senior review B3): rank-deficient, or < 2 residual degrees of freedom → not estimable
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or X.shape[0] - X.shape[1] < 2:
        log.info(
            "per-study model skipped: %d units, %d parameters (rank-deficient or no residual df)", *X.shape
        )
        return pd.DataFrame()
    inference = DefaultInference(n_cpus=n_jobs)
    dds = DeseqDataSet(adata=sub, design=X, refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    cols = list(dds.varm["LFC"].columns)
    nd = sub.obs.groupby("study", observed=True)["donor_id"].nunique()
    rows = []
    for j, c in enumerate(cols):
        if not c.startswith("tumour_in_"):
            continue
        vec = np.zeros(len(cols))
        vec[j] = 1.0
        st = DeseqStats(dds, contrast=vec, inference=inference, quiet=True, independent_filter=False)
        st.summary()
        r = st.results_df[["log2FoldChange", "lfcSE", "pvalue"]].rename_axis("gene").reset_index()
        study = c.removeprefix("tumour_in_")
        rows.append(r.assign(study=study, n_donors=int(nd.get(study, 0))))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def meta_summary(per_study: pd.DataFrame, pooled: pd.DataFrame, min_donors: int = 2) -> pd.DataFrame:
    """Fixed-effect (inverse-variance) meta-analysis over studies with ≥ ``min_donors`` pairs, Cochran's Q
    and I²
    (descriptive; unreliable with K = 2–4), concordance with the pooled sign (graded, k of K), and whether
    any such
    study contradicts the pooled direction with a 95 % CI excluding 0. 1-pair studies are point estimates
    only."""
    from scipy.stats import chi2

    ps = per_study[per_study["n_donors"] >= min_donors].dropna(subset=["log2FoldChange", "lfcSE"])
    ps = ps[ps["lfcSE"] > 0].merge(
        pooled[["gene", "group", "log2FoldChange"]].rename(columns={"log2FoldChange": "lfc_pooled"}),
        on=["gene", "group"],
    )
    ps["w"] = 1 / ps["lfcSE"] ** 2
    ps["wb"] = ps["w"] * ps["log2FoldChange"]
    ps["agree"] = np.sign(ps["log2FoldChange"]) == np.sign(ps["lfc_pooled"])
    ps["contra"] = (~ps["agree"]) & ((ps["log2FoldChange"].abs() - 1.96 * ps["lfcSE"]) > 0)
    g = ps.groupby(["gene", "group"])
    out = g.agg(
        k_studies=("study", "nunique"),
        k_concordant=("agree", "sum"),
        any_study_contradicts=("contra", "any"),
        sw=("w", "sum"),
        swb=("wb", "sum"),
    )
    out["meta_log2fc"] = out["swb"] / out["sw"]
    out["meta_se"] = 1 / np.sqrt(out["sw"])
    ps = ps.join(out["meta_log2fc"], on=["gene", "group"])
    ps["q"] = ps["w"] * (ps["log2FoldChange"] - ps["meta_log2fc"]) ** 2
    out["cochran_q"] = ps.groupby(["gene", "group"])["q"].sum()
    df = (out["k_studies"] - 1).clip(lower=1)
    out["q_pvalue"] = chi2.sf(out["cochran_q"], df)
    out["i2"] = ((out["cochran_q"] - (out["k_studies"] - 1)) / out["cochran_q"]).clip(lower=0).fillna(0)
    return out.drop(columns=["sw", "swb"]).reset_index()


def paired_robustness(
    pdata: ad.AnnData, dcfg: dict, name: str, spec: dict, main: pd.DataFrame, n_jobs: int = 8
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Study-grounded checks for the paired contrast (§3.12; DEVELOPMENT_LOG §3.3 K; senior review items 2–3).

    1. **Per-study effects** from ONE joint model ``~ donor + study:condition`` (shared dispersions) →
       fixed-effect meta-analysis, Cochran's Q / I², concordance k of K (studies with ≥ 2 pairs).
    2. **Leave-one-study-out (LOSO)** refits of the pooled model, for every study.
    3. ``robust`` = significant in the thresholded test (|log2FC| > log2 1.5, padj_thr < alpha) AND no study
    with ≥ 2 pairs contradicting the direction with a 95 % CI excluding 0 AND, for every *major* study (≥
       ``loso_major_fraction`` of the pairs), the refit without it keeps the sign with padj < ``loso_alpha``.
       Inference is conditional on these studies (2–4 per cell type: too few for a random study effect).
    """
    cond = dcfg["condition_key"]
    case, ctrl = spec["case"], spec["control"]
    loso, per_study = [], []
    groups = main.loc[
        main["tier"] == "primary", "group"
    ].unique()  # robustness only where inference is claimed
    for grp in groups:
        full = _subset(pdata, dcfg, grp, case, ctrl, True)
        full = _prepare(full, cond, case, ctrl)
        try:
            ps = per_study_effects(full, cond, case, n_jobs)
            if len(ps):
                per_study.append(ps.assign(group=grp))
        except Exception as exc:  # noqa: BLE001
            log.warning("[%s per-study] %s failed: %s", name, grp, exc)
        n_pairs = full.obs.groupby("study", observed=True)["donor_id"].nunique()
        for st, k in n_pairs.items():
            sub = full[full.obs["study"].astype(str) != st].copy()
            if sub.obs["donor_id"].nunique() < 3:
                continue
            try:
                r = _fit(sub, cond, case, ctrl, design_for(sub.obs, cond, True), n_jobs)
            except Exception as exc:  # noqa: BLE001
                log.warning("[%s LOSO] %s without %s failed: %s", name, grp, st, exc)
                continue
            loso.append(
                r[["gene", "log2FoldChange", "padj"]].assign(
                    group=grp,
                    left_out=st,
                    n_pairs_left_out=int(k),
                    major=bool(k / n_pairs.sum() >= dcfg.get("loso_major_fraction", 0.25)),
                    n_donors=sub.obs["donor_id"].nunique(),
                )
            )
    loso = (
        pd.concat(loso, ignore_index=True)
        if loso
        else pd.DataFrame(
            columns=[
                "gene",
                "log2FoldChange",
                "padj",
                "group",
                "left_out",
                "n_pairs_left_out",
                "major",
                "n_donors",
            ]
        )
    )
    per_study = (
        pd.concat(per_study, ignore_index=True)
        if per_study
        else pd.DataFrame(columns=["gene", "log2FoldChange", "lfcSE", "pvalue", "study", "n_donors", "group"])
    )
    meta = meta_summary(per_study, main) if len(per_study) else pd.DataFrame(columns=["gene", "group"])
    lo = loso.merge(
        main[["gene", "group", "log2FoldChange"]].rename(columns={"log2FoldChange": "lfc"}),
        on=["gene", "group"],
    )
    lo["ok"] = (np.sign(lo["log2FoldChange"]) == np.sign(lo["lfc"])) & (
        lo["padj"] < dcfg.get("loso_alpha", 0.1)
    )
    lo_major = (
        lo[lo["major"]]
        .groupby(["gene", "group"])
        .agg(n_major_loso=("left_out", "nunique"), major_loso_ok=("ok", "all"))
    )
    lo_all = lo.groupby(["gene", "group"]).agg(
        n_loso_fits=("left_out", "nunique"), loso_max_padj=("padj", "max")
    )
    summ = (
        main[["gene", "group", "padj_thr", "tier"]]
        .set_index(["gene", "group"])
        .join(meta.set_index(["gene", "group"]))
        .join(lo_all)
        .join(lo_major)
    )
    for c in (
        "k_studies",
        "k_concordant",
        "any_study_contradicts",
        "meta_log2fc",
        "meta_se",
        "cochran_q",
        "i2",
        "n_loso_fits",
        "loso_max_padj",
        "n_major_loso",
        "major_loso_ok",
    ):
        if c not in summ:
            summ[c] = np.nan
    summ["per_study_available"] = summ["k_studies"].fillna(0) >= 1
    summ["robust"] = (
        (summ["tier"] == "primary")
        & summ["per_study_available"]
        & (summ["padj_thr"] < dcfg["alpha"])
        & ~summ["any_study_contradicts"].fillna(False).astype(bool)
        & summ["major_loso_ok"].fillna(False).astype(bool)
        & (summ["n_major_loso"].fillna(0) > 0)
    )
    summ = summ.drop(columns=["padj_thr", "tier"]).reset_index()
    return loso.assign(contrast=name), per_study.assign(contrast=name), summ


def supported_by_paired(exploratory: pd.DataFrame, paired: pd.DataFrame, alpha: float) -> pd.Series:
    """For an atlas-confounded contrast: is the gene × group also significant (thresholded test, > 1.5-fold),
    same sign, in the paired within-donor contrast? (§3.12: the paired contrast is the only one protected from
    study/atlas effects.)"""
    col = "padj_thr" if "padj_thr" in paired else "padj"  # same evidence bar as the paired headlines
    p = paired[["gene", "group", "log2FoldChange", col]].rename(
        columns={"log2FoldChange": "lfc_p", col: "padj_p"}
    )
    m = exploratory[["gene", "group", "log2FoldChange"]].merge(p, on=["gene", "group"], how="left")
    ok = (m["padj_p"] < alpha) & (np.sign(m["lfc_p"]) == np.sign(m["log2FoldChange"]))
    return ok.fillna(False).to_numpy()


def ambient_flags(pdata: ad.AnnData, res: pd.DataFrame, dcfg: dict, ratio: float = 10.0) -> pd.DataFrame:
    """Ambient-RNA / contamination flag for up-regulated hits (§3.12, senior review item 6; no raw
    droplets, so
    CellBender/SoupX are impossible). For gene g up in type X: ``ambient_ratio`` = max over other types S of
    mean CPM_S / mean CPM_X in tumour units, ``ambient_source`` = that S. ``ambient_suspect`` when the ratio
    is ≥
    ``ratio`` (g is a marker of another type, e.g. Ig genes in non-B cells). ``up_in_frac_groups`` =
    fraction of the
    other tested types where g is also up (p < 0.05): a pan-type rise is the other signature of ambient
    RNA."""
    cond, gkey = dcfg["condition_key"], dcfg["group_key"]
    tum = pdata[pdata.obs[cond].astype(str) == "tumour"]
    X = tum.X.toarray() if hasattr(tum.X, "toarray") else np.asarray(tum.X)
    cpm = X / np.clip(X.sum(1, keepdims=True), 1, None) * 1e6
    mean = pd.DataFrame(cpm, columns=tum.var_names).groupby(tum.obs[gkey].astype(str).values).mean()
    up = res[(res["log2FoldChange"] > 0)][["gene", "group"]]
    rows = []
    sig_up = res[(res["log2FoldChange"] > 0) & (res["pvalue"] < 0.05)].groupby("gene")["group"].nunique()
    n_groups = res["group"].nunique()
    for grp, d in up.groupby("group"):
        if grp not in mean.index:
            continue
        genes = [g for g in d["gene"] if g in mean.columns]
        own = mean.loc[grp, genes].clip(lower=1e-3)
        others = mean.drop(index=grp)[genes]
        r = others.div(own, axis=1)
        rows.append(
            pd.DataFrame(
                {
                    "gene": genes,
                    "group": grp,
                    "ambient_ratio": r.max(0).values,
                    "ambient_source": r.idxmax(0).values,
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=["gene", "group", "ambient_ratio", "ambient_source", "ambient_suspect"])
    out = pd.concat(rows, ignore_index=True)
    out["up_in_frac_groups"] = (out["gene"].map(sig_up).fillna(0) - 1).clip(lower=0) / max(n_groups - 1, 1)
    out["ambient_suspect"] = out["ambient_ratio"] >= ratio
    return out


def power_table(res: pd.DataFrame, n_pairs: dict, alpha: float = 0.05, power: float = 0.8) -> pd.DataFrame:
    """Approximate minimum detectable |log2FC| at the given power for each tested type (Wald approximation
    for a
    paired NB design: SE ≈ sqrt(2 (1/μ + φ) / n) / ln 2, with μ the median normalised count and φ the median
    fitted dispersion of the type's genes). A planning aid, not a replacement for simulation."""
    from scipy.stats import norm

    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    rows = []
    for grp, d in res.groupby("group"):
        n = n_pairs.get(grp)
        if not n:
            continue
        for label, q in (("median-expressed gene", 0.5), ("lower-quartile gene", 0.25)):
            mu = float(d["baseMean"].quantile(q))
            phi = float(d["dispersion"].median())
            se = np.sqrt(2 * (1 / max(mu, 1e-6) + phi) / n) / np.log(2)
            rows.append(
                {
                    "group": grp,
                    "n_pairs": n,
                    "gene_class": label,
                    "base_mean": mu,
                    "median_dispersion": phi,
                    "mde_log2fc": z * se,
                }
            )
    return pd.DataFrame(rows)


def sex_check(pdata: ad.AnnData, cond: str) -> pd.DataFrame:
    """Sample-swap check (§3.12): per donor × tissue unit (all cell types summed), XIST vs RPS4Y1 CPM →
    inferred sex; compared between a patient's tumour and adjacent units and with the recorded sex."""
    genes = [g for g in ("XIST", "RPS4Y1") if g in pdata.var_names]
    if len(genes) < 2:
        return pd.DataFrame()
    unit = pdata.obs["donor_id"].astype(str) + "|" + pdata.obs[cond].astype(str)
    X = pdata.X.toarray() if hasattr(pdata.X, "toarray") else np.asarray(pdata.X)
    tot = pd.Series(X.sum(1), index=pdata.obs_names).groupby(unit.values).sum()
    sub = (
        pd.DataFrame(X[:, [pdata.var_names.get_loc(g) for g in genes]], columns=genes)
        .groupby(unit.values)
        .sum()
    )
    cpm = sub.div(tot, axis=0) * 1e6
    out = cpm.reset_index(names="unit")
    out[["donor_id", "tissue_type"]] = out["unit"].str.rsplit("|", n=1, expand=True)
    out["inferred_sex"] = np.where(
        (out["RPS4Y1"] == 0) & (out["XIST"] == 0),
        "undetermined",  # both markers absent (run 3 review: the tie rule called these "female")
        np.where(out["RPS4Y1"] > out["XIST"], "male", "female"),
    )
    if "sex" in pdata.obs:
        rec = pdata.obs.groupby("donor_id", observed=True)["sex"].agg(lambda x: str(x.iloc[0]))
        out["recorded_sex"] = out["donor_id"].map(rec)
    consistent = out.groupby("donor_id")["inferred_sex"].nunique() == 1
    out["consistent_within_donor"] = out["donor_id"].map(consistent)
    return out


def cross_study_null(pdata: ad.AnnData, dcfg: dict, n_jobs: int = 8) -> pd.DataFrame:
    """A real technical null (§3.12, senior review item 9): healthy HLCA donors from studies of ONE chemistry,
    studies split into two sets; set B vs set A per cell type. Healthy = lung tissue from cancer-free donors
    only (§3.13: lavage and tumour-adjacent surgery tissue removed), so no biology should differ; the share of
    genes with p < 0.01 is the cross-study DE floor (one atlas and chemistry: a lower bound)."""
    cond = dcfg["condition_key"]
    h = pdata[pdata.obs[cond].astype(str) == "healthy"].copy()
    if "assay" not in h.obs:
        return pd.DataFrame()
    assay = h.obs.groupby("assay", observed=True)["study"].nunique().idxmax()
    h = h[h.obs["assay"].astype(str) == str(assay)].copy()
    studies = h.obs.groupby("study", observed=True)["donor_id"].nunique().sort_values(ascending=False)
    set_b = set(studies.index[1::2])
    h.obs["null_set"] = np.where(h.obs["study"].isin(set_b), "B", "A")
    res = run_contrast(
        h,
        {**dcfg, "condition_key": "null_set", "tumour_only_groups": [], "lfc_threshold": None},
        "healthy_study_null",
        {"case": "B", "control": "A", "paired": False, "inference": "technical null"},
        n_jobs,
    )
    if res.empty:
        return res
    summ = (
        res.groupby("group")
        .agg(
            n_genes=("gene", "size"),
            frac_p_lt_001=("pvalue", lambda p: float((p < 0.01).mean())),
            n_padj_lt_005=("padj", lambda p: int((p < 0.05).sum())),
            n_case=("n_case", "first"),
            n_control=("n_control", "first"),
        )
        .reset_index()
    )
    summ["assay"] = str(assay)
    summ["set_a"] = ";".join(sorted(set(studies.index) - set_b))
    summ["set_b"] = ";".join(sorted(set_b))
    return summ


def within_study_null(pdata: ad.AnnData, dcfg: dict, n_jobs: int = 8, n_splits: int = 3) -> pd.DataFrame:
    """Calibration check for the unpaired pipeline (§3.12, final review).

    Healthy donors of the ONE study with most healthy donors, split at random into two halves (``n_splits``
    seeds). Nothing differs by design, not even study, so the share of genes with p < 0.01 should be ≈ 1 %.
    Read next to ``cross_study_null``: if this is near 1 % and that is far above it, the cross-study excess
    is a study effect, not a miscalibrated test."""
    cond = dcfg["condition_key"]
    h = pdata[pdata.obs[cond].astype(str) == "healthy"].copy()
    study = h.obs.groupby("study", observed=True)["donor_id"].nunique().idxmax()
    h = h[h.obs["study"].astype(str) == str(study)].copy()
    donors = sorted(h.obs["donor_id"].astype(str).unique())
    rows = []
    for seed in range(n_splits):
        rng = np.random.default_rng(seed)
        half = set(rng.choice(donors, len(donors) // 2, replace=False))
        h.obs["null_set"] = np.where(h.obs["donor_id"].astype(str).isin(half), "B", "A")
        res = run_contrast(
            h,
            {**dcfg, "condition_key": "null_set", "tumour_only_groups": [], "lfc_threshold": None},
            "within_study_null",
            {"case": "B", "control": "A", "paired": False, "inference": "calibration"},
            n_jobs,
        )
        if res.empty:
            continue
        g = (
            res.groupby("group")
            .agg(
                n_genes=("gene", "size"),
                frac_p_lt_001=("pvalue", lambda p: float((p < 0.01).mean())),
                n_case=("n_case", "first"),
                n_control=("n_control", "first"),
            )
            .reset_index()
        )
        g["seed"], g["study"] = seed, str(study)
        rows.append(g)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def paired_sign_flip_null(
    pdata: ad.AnnData, dcfg: dict, spec: dict, n_jobs: int = 8, seeds: int = 3
) -> pd.DataFrame:
    """Calibration check for the PAIRED test (§3.13, decision 5a).

    For each primary cell type: swap the tumour and adjacent labels in exactly half of the patients (with an
    odd number, one random patient is left out so the halves balance), refit the same paired model and record
    the share of genes with p < 0.01. If the test is calibrated this is ≈ 1 % or less; real tumour effects do
    not vanish in the swap but become noise (up in half the patients, down in the other half), which inflates
    the dispersion and pushes the share *down*, so a pass is weaker evidence than a fail."""
    cond = dcfg["condition_key"]
    case, ctrl = spec["case"], spec["control"]
    skip = set(dcfg.get("tumour_only_groups", [])) | set(dcfg.get("exclude_groups", []))
    rows = []
    for grp in pdata.obs[dcfg["group_key"]].cat.categories:
        if grp in skip:
            continue
        base = _subset(pdata, dcfg, grp, case, ctrl, paired=True)
        tier, _ = paired_eligibility(base, dcfg)
        if tier != "primary":
            continue
        for seed in range(seeds):
            rng = np.random.default_rng(seed)
            donors = sorted(base.obs["donor_id"].astype(str).unique())
            if len(donors) % 2:
                donors.remove(str(rng.choice(donors)))
            flip = set(rng.choice(donors, len(donors) // 2, replace=False))
            sub = base[base.obs["donor_id"].astype(str).isin(donors)].copy()
            lab = sub.obs[cond].astype(str)
            swap = sub.obs["donor_id"].astype(str).isin(flip)
            sub.obs[cond] = np.where(swap, lab.map({case: ctrl, ctrl: case}), lab)
            sub = _prepare(sub, cond, case, ctrl)
            try:
                res = _fit(sub, cond, case, ctrl, design_for(sub.obs, cond, True), n_jobs)
            except Exception as exc:  # noqa: BLE001
                log.warning("[sign-flip null] %s seed %d failed: %s", grp, seed, exc)
                continue
            rows.append(
                {
                    "group": grp,
                    "seed": seed,
                    "n_pairs": len(donors),
                    "n_flipped": len(flip),
                    "n_genes": int(res["pvalue"].notna().sum()),
                    "frac_p_lt_001": float((res["pvalue"].dropna() < 0.01).mean()),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        flag = float(dcfg.get("sign_flip_flag_frac", 0.03))
        mean = out.groupby("group")["frac_p_lt_001"].transform("mean")
        out["mean_frac_p_lt_001"] = mean
        out["flagged"] = mean > flag
    return out


AIRWAY_EPI = ("Basal", "Multiciliated", "Secretory")
ALVEOLAR_EPI = ("AT1", "AT2")
BLOOD_MYELOID = ("Monocyte classical", "Monocyte non-classical", "Neutrophil")


def adjacent_makeup(obs: pd.DataFrame, cond: str, case: str, ctrl: str, group_key: str) -> pd.DataFrame:
    """What each paired patient's adjacent ("control") sample actually contains (§3.13, decision 5b).

    Descriptive only; changes no result. Per paired patient: compartment shares, airway vs alveolar share of
    labelled epithelium, blood-type myeloid share of immune cells (monocytes + neutrophils, which also reach
    tissue from blood), and the CNV call rate in adjacent epithelium (cells that look tumour-like)."""
    both = obs.groupby("donor_id", observed=True)[cond].agg(lambda s: {case, ctrl} <= set(s.astype(str)))
    adj = obs[obs["donor_id"].isin(both[both].index) & (obs[cond].astype(str) == ctrl)]
    rows = []
    for (study, donor), o in adj.groupby(["study", "donor_id"], observed=True):
        comp = o["compartment"].astype(str).value_counts(normalize=True)
        epi = o.loc[o["compartment"].astype(str) == "Epithelial", "harmonised_label"].astype(str)
        n_air, n_alv = int(epi.isin(AIRWAY_EPI).sum()), int(epi.isin(ALVEOLAR_EPI).sum())
        imm = o.loc[o["compartment"].astype(str) == "Immune", group_key].astype(str)
        cs = o.loc[o["compartment"].astype(str) == "Epithelial", "cnv_status"].astype(str)
        called, tested = (
            int(cs.str.startswith("CNV-high").sum()),
            int(cs.str.startswith(("CNV-high", "not called")).sum()),
        )
        rows.append(
            {
                "study": str(study),
                "donor_id": str(donor),
                "n_cells": len(o),
                **{
                    f"share_{k.lower()}": float(comp.get(k, 0.0))
                    for k in ("Immune", "Epithelial", "Endothelial", "Stroma", "Unknown")
                },
                "n_epithelial_labelled": n_air + n_alv,
                "share_airway_of_epithelium": n_air / (n_air + n_alv) if n_air + n_alv else np.nan,
                "share_blood_myeloid_of_immune": float(imm.isin(BLOOD_MYELOID).mean())
                if len(imm)
                else np.nan,
                "n_epithelial_cnv_tested": tested,
                "cnv_call_rate": called / tested if tested else np.nan,
            }
        )
    return pd.DataFrame(rows)


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
            inference="hypothesis-generating (unpaired; atlas-confounded)",
            n_case=int(n["tumour_epi"]),
            n_control=int(n["healthy_basal"]),
            n_donors=int(sub.obs["donor_id"].nunique()),
            n_studies=int(sub.obs["study"].nunique()) if "study" in sub.obs else np.nan,
            studies=";".join(sorted(sub.obs["study"].astype(str).unique())) if "study" in sub.obs else "",
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


def composition_counts(
    obs: pd.DataFrame, acfg: dict, group_key: str = "analysis_group", unit_key: str = "donor_condition"
) -> pd.DataFrame:
    """Long table: unit (donor × tissue) × compartment × cell type → cell count and share within compartment.
    All units of all tissues (used for the descriptive views); paired testing selects from it."""
    o = obs.copy()
    o["_unit"] = o[unit_key].astype(str)
    excl = set(acfg.get("exclude_groups", []))
    o = o[~o[group_key].astype(str).isin(excl)]
    cnt = (
        o.groupby(["_unit", "donor_id", "study", "tissue_type", "compartment", group_key], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={"_unit": "unit", group_key: "cell_type"})
    )
    tot = cnt.groupby(["unit", "compartment"])["n_cells"].transform("sum")
    cnt["compartment_total"] = tot
    cnt["share_in_compartment"] = cnt["n_cells"] / tot
    cnt["unit_total"] = cnt.groupby("unit")["n_cells"].transform("sum")
    cnt["share_of_unit"] = cnt["n_cells"] / cnt["unit_total"]
    return cnt


def composition_paired(cnt: pd.DataFrame, acfg: dict) -> pd.DataFrame:
    """Paired compositional test within one compartment (§3.12; senior review item 7).

    Units = donor × tissue of patients with both tumour and adjacent tissue; proportions are WITHIN the
    compartment (e.g. each immune type as a share of immune cells), which neutralises study differences in
    CD45/CD31 sorting.
    Units with < ``min_compartment_cells`` cells in the compartment are dropped with their partner.

    Three views per cell type, BH across types within each:
    1. **NB GLM** (primary): count ~ donor + tissue, offset log(compartment total) — the edgeR-style model
       behind Milo/propeller; donor fixed effects make it paired. statsmodels NegativeBinomial (MLE
       dispersion); falls back to quasi-Poisson if it does not converge.
    2. **propeller-style**: arcsin √proportion, paired t-test (Phipson 2022).
    3. **CLR**: centred log-ratio with pseudocount 0.5, paired t-test (compositional).
    ``consistent`` = all three padj < alpha with the same sign. ``k_studies_same_direction`` counts studies
    whose mean within-donor change has the GLM's sign (per-study replication).
    """
    import statsmodels.api as sm
    from scipy.stats import false_discovery_control, ttest_rel

    comp, case, ctrl = acfg["compartment"], acfg["case"], acfg["control"]
    d = cnt[(cnt["compartment"] == comp) & cnt["tissue_type"].isin([case, ctrl])].copy()
    # plain strings: categorical columns make get_dummies/groupby emit every category, including donors
    # outside this test → all-zero design columns → singular NB fit (run 3 bug, silently quasi-Poisson)
    for c in ("unit", "donor_id", "study", "tissue_type", "cell_type"):
        d[c] = d[c].astype(str)
    tot = (
        d.groupby(["unit", "donor_id", "study", "tissue_type"])["n_cells"].sum().rename("total").reset_index()
    )
    tot = tot[tot["total"] >= acfg["min_compartment_cells"]]
    both = tot.groupby("donor_id")["tissue_type"].nunique()
    tot = tot[tot["donor_id"].isin(both[both == 2].index)]
    wide = (
        d[d["unit"].isin(tot["unit"])]
        .pivot_table(index="unit", columns="cell_type", values="n_cells", aggfunc="sum", fill_value=0)
        .reindex(tot["unit"])
    )
    tot = tot.set_index("unit").loc[wide.index]
    types = [t for t in wide.columns if (wide[t] > 0).mean() >= acfg.get("min_type_presence", 0.5)]
    wide = wide[types]
    is_case = (tot["tissue_type"] == case).astype(float).to_numpy()
    donors = pd.get_dummies(tot["donor_id"], drop_first=True, dtype=float).to_numpy()
    X = np.column_stack([np.ones(len(tot)), donors, is_case])
    off = np.log(tot["total"].to_numpy(dtype=float))
    prop = wide.div(tot["total"], axis=0)
    clr = np.log(wide + 0.5)
    clr = clr.sub(clr.mean(axis=1), axis=0)
    pair = tot.reset_index().set_index(["donor_id", "tissue_type"])["unit"]
    donors_ok = sorted(set(tot["donor_id"]))
    u_case = [pair[(dn, case)] for dn in donors_ok]
    u_ctrl = [pair[(dn, ctrl)] for dn in donors_ok]
    study_of = tot.reset_index().drop_duplicates("donor_id").set_index("donor_id")["study"]
    rows = []
    for t in types:
        y = wide[t].to_numpy(dtype=float)
        method = "NB GLM (MLE dispersion)"
        try:
            # Newton from statsmodels' default start returns NaN for several types (run 3); BFGS from Poisson
            # start values converges and matches Newton wherever Newton works
            p0 = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit().params
            fit = sm.NegativeBinomial(y, X, offset=off).fit(
                start_params=np.append(p0, 0.3), method="bfgs", disp=0, maxiter=5000
            )
            if not fit.mle_retvals.get("converged", True) or not np.isfinite(fit.bse[X.shape[1] - 1]):
                raise RuntimeError("not converged")
            b, se, p = fit.params[X.shape[1] - 1], fit.bse[X.shape[1] - 1], fit.pvalues[X.shape[1] - 1]
        except Exception:  # noqa: BLE001
            fit = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(scale="X2")
            b, se, p = fit.params[-1], fit.bse[-1], fit.pvalues[-1]
            method = "quasi-Poisson GLM (NB did not converge)"
        a_case, a_ctrl = np.arcsin(np.sqrt(prop.loc[u_case, t])), np.arcsin(np.sqrt(prop.loc[u_ctrl, t]))
        t1 = ttest_rel(a_case.to_numpy(), a_ctrl.to_numpy())
        c_case, c_ctrl = clr.loc[u_case, t].to_numpy(), clr.loc[u_ctrl, t].to_numpy()
        t2 = ttest_rel(c_case, c_ctrl)
        dprop = pd.Series(prop.loc[u_case, t].to_numpy() - prop.loc[u_ctrl, t].to_numpy(), index=donors_ok)
        # plain strings: a categorical study column would make groupby return every category (run 3 bug)
        per_study = dprop.groupby(study_of.reindex(donors_ok).astype(str).to_numpy()).mean().dropna()
        rows.append(
            {
                "compartment": comp,
                "cell_type": t,
                "n_donors": len(donors_ok),
                "n_studies": int(study_of.reindex(donors_ok).nunique()),
                "mean_share_adjacent": float(prop.loc[u_ctrl, t].mean()),
                "mean_share_tumour": float(prop.loc[u_case, t].mean()),
                "n_donors_higher_in_tumour": int((dprop > 0).sum()),
                "glm_log2fc": float(b / np.log(2)),
                "glm_ci_low": float((b - 1.96 * se) / np.log(2)),
                "glm_ci_high": float((b + 1.96 * se) / np.log(2)),
                "glm_pvalue": float(p),
                "glm_method": method,
                "propeller_t": float(t1.statistic),
                "propeller_pvalue": float(t1.pvalue),
                "clr_mean_difference": float(np.mean(c_case - c_ctrl)),
                "clr_pvalue": float(t2.pvalue),
                "k_studies": int(per_study.size),
                "k_studies_same_direction": int((np.sign(per_study) == np.sign(b)).sum()),
            }
        )
    out = pd.DataFrame(rows)
    alpha = acfg.get("alpha", 0.05)
    for c in ("glm", "propeller", "clr"):
        p = out[f"{c}_pvalue"].fillna(1).to_numpy()
        out[f"{c}_padj"] = false_discovery_control(p)
    sig = (out["glm_padj"] < alpha) & (out["propeller_padj"] < alpha) & (out["clr_padj"] < alpha)
    same = (np.sign(out["glm_log2fc"]) == np.sign(out["propeller_t"])) & (
        np.sign(out["glm_log2fc"]) == np.sign(out["clr_mean_difference"])
    )
    out["consistent"] = np.where(
        sig & same,
        np.where(out["glm_log2fc"] > 0, "higher in tumour", "lower in tumour"),
        "not consistent across models",
    )
    return out.sort_values("glm_log2fc")


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
    sub.obs[sample_key] = sub.obs[sample_key].astype(str).astype("category")
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
    milo.annotate_nhoods(mdata, anno_col=acfg.get("group_key", "de_group"))
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
