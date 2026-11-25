"""Tests for the run-3 overhaul (DEVELOPMENT_LOG §3.12)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# pipeline tests: run in the default env (`pixi run pytest`); the dashboard env has no analysis stack
pytest.importorskip("scanpy")
pytest.importorskip("statsmodels")


# ------------------------------------------------------------------ §3.10 #9 mitochondrial fraction
def test_mito_fraction_from_full_matrix():
    from lusc_atlas.select import mito_fraction

    var = pd.DataFrame({"feature_name": ["MT-CO1", "MT-ND1", "ACTB", "GAPDH"]}, index=list("abcd"))
    X = sp.csr_matrix(np.array([[10, 10, 60, 20], [0, 0, 50, 50], [0, 0, 0, 0]], dtype=np.float32))
    pct, src = mito_fraction(X, var, pd.DataFrame(index=range(3)))
    assert src.startswith("raw matrix (2")
    assert pct[0] == pytest.approx(20.0) and pct[1] == 0.0 and np.isnan(pct[2])


def test_mito_fraction_falls_back_to_author_column():
    from lusc_atlas.select import mito_fraction

    var = pd.DataFrame({"feature_name": ["ACTB", "GAPDH"]}, index=list("ab"))
    obs = pd.DataFrame({"pct_counts_mito": [1.5, 12.0]})
    pct, src = mito_fraction(sp.csr_matrix((2, 2)), var, obs, "pct_counts_mito")
    assert "pct_counts_mito" in src and pct.tolist() == [1.5, 12.0]
    pct, src = mito_fraction(sp.csr_matrix((2, 2)), var, obs, None)
    assert src == "none" and np.isnan(pct).all()


def test_mito_rule_is_upper_tail_and_skips_mad_zero_samples():
    from lusc_atlas.qc import mito_flags

    rng = np.random.default_rng(0)
    v = pd.Series(np.r_[rng.normal(5, 1, 500).clip(0.5), [0.0, 30.0]])
    f, rule = mito_flags(v, 3, 20)
    assert rule == "upper MAD + ceiling"
    assert f.iloc[-1] and not f.iloc[-2]  # very low mito is NOT a QC failure
    zeros = pd.Series(np.r_[np.zeros(300), [0.5, 25.0]])
    f, rule = mito_flags(zeros, 3, 20)
    assert rule == "ceiling only (MAD 0)"
    assert f.tolist()[-2:] == [False, True] and f.sum() == 1
    f, rule = mito_flags(pd.Series([np.nan] * 5), 3, 20)
    assert rule.startswith("none") and not f.any()


def test_compute_qc_metrics_refuses_missing_mito():
    import anndata as ad

    from lusc_atlas.qc import compute_qc_metrics

    a = ad.AnnData(
        X=sp.csr_matrix(np.ones((5, 25), dtype=np.float32)),
        var=pd.DataFrame(index=["ACTB", "RPS3"] + [f"G{i}" for i in range(23)]),
    )
    with pytest.raises(ValueError, match="mitochondrial"):
        compute_qc_metrics(a)
    a.obs["pct_counts_mt"] = [1.0, 2, 3, 4, 50]
    compute_qc_metrics(a)
    assert a.obs["pct_counts_mt"].tolist() == [1.0, 2, 3, 4, 50]  # not overwritten with 0 %


# ------------------------------------------------------------------ §3.12 age parsing (row-wise fallback)
def test_age_years_falls_back_row_by_row():
    from lusc_atlas.select import age_years

    obs = pd.DataFrame(
        {
            "age_or_mean_of_age_range": [60.0, np.nan, np.nan],
            "development_stage": ["60-year-old stage", "24-year-old stage", "unknown"],
        }
    )
    out = age_years(obs)
    assert out.iloc[0] == 60 and out.iloc[1] == 24 and np.isnan(out.iloc[2])


# ------------------------------------------------------------------ §3.12 Scrublet threshold guard
def test_scrublet_guard_replaces_failed_thresholds_only():
    from lusc_atlas.qc import guard_scrublet

    rng = np.random.default_rng(0)
    n = 200
    obs = pd.DataFrame(
        {
            "sample": ["ok"] * n + ["bad"] * n,
            "predicted_doublet": [False] * (n - 4) + [True] * 4 + [True] * (n - 10) + [False] * 10,
            "doublet_score": np.r_[rng.uniform(0, 0.2, n), rng.uniform(0, 1, n)],
        }
    )
    calls, rules = guard_scrublet(obs, "sample", 0.25, 0.008, 0.1)
    assert calls[:n].sum() == 4 and rules["ok"].startswith("scrublet")
    # 200 cells → expected rate 0.16 % → top 1 cell only
    assert calls[n:].sum() == 1 and rules["bad"].startswith("fallback")


# ------------------------------------------------------------------ §3.12 DE: joint per-study design, meta,
# eligibility
def test_study_design_has_one_tumour_effect_per_study_and_full_rank():
    from lusc_atlas.de import _study_design

    obs = pd.DataFrame(
        {
            "donor_id": ["a", "a", "b", "b", "c", "c"],
            "study": ["S1", "S1", "S1", "S1", "S2", "S2"],
            "tissue_type": ["tumour", "adjacent_normal"] * 3,
        }
    )
    X = _study_design(obs, "tissue_type", "tumour")
    assert {"tumour_in_S1", "tumour_in_S2"} <= set(X.columns)
    assert np.linalg.matrix_rank(X.to_numpy()) == X.shape[1]


def test_meta_summary_concordance_and_contradiction():
    from lusc_atlas.de import meta_summary

    per = pd.DataFrame(
        {
            "gene": ["g"] * 3,
            "group": ["T"] * 3,
            "study": ["A", "B", "C"],
            "n_donors": [5, 4, 1],
            "log2FoldChange": [1.0, -2.0, 3.0],
            "lfcSE": [0.2, 0.3, 1.0],
        }
    )
    pooled = pd.DataFrame({"gene": ["g"], "group": ["T"], "log2FoldChange": [0.5]})
    m = meta_summary(per, pooled).iloc[0]
    assert m.k_studies == 2  # the 1-pair study is a point estimate only
    assert m.k_concordant == 1 and bool(m.any_study_contradicts)
    assert m.i2 > 0.5


def test_paired_eligibility_tiers():
    import anndata as ad

    from lusc_atlas.de import paired_eligibility

    def make(studies):
        rows = [{"donor_id": f"{s}{i}", "study": s} for s, k in studies.items() for i in range(k)]
        return ad.AnnData(obs=pd.DataFrame(rows, index=[str(i) for i in range(len(rows))]))

    cfg = {"min_paired_donors": 6, "min_paired_studies": 2, "min_descriptive_donors": 3}
    assert paired_eligibility(make({"L": 5, "U": 5}), cfg)[0] == "primary"
    assert paired_eligibility(make({"G": 3, "T": 1, "L": 1, "U": 1}), cfg)[0] == "descriptive"
    assert paired_eligibility(make({"U": 2}), cfg)[0] is None


# ------------------------------------------------------------------ §3.12 composition: paired, within
# compartment
def test_composition_paired_detects_a_real_shift_and_is_paired():
    from lusc_atlas.de import composition_paired

    rng = np.random.default_rng(1)
    rows = []
    for d in range(10):
        depth = rng.integers(500, 3000)
        for tissue, shift in (("adjacent_normal", 1.0), ("tumour", 3.0)):
            p = np.array([0.1 * shift, 0.3, 0.3, 0.3])
            p = p / p.sum()
            n = rng.multinomial(depth, p)
            for t, k in zip(["Treg", "NK", "B", "Mono"], n, strict=True):
                rows.append(
                    {
                        "unit": f"d{d}|{tissue}",
                        "donor_id": f"d{d}",
                        "study": "S1" if d < 5 else "S2",
                        "tissue_type": tissue,
                        "compartment": "Immune",
                        "cell_type": t,
                        "n_cells": int(k),
                    }
                )
    cnt = pd.DataFrame(rows)
    out = composition_paired(
        cnt,
        {
            "compartment": "Immune",
            "case": "tumour",
            "control": "adjacent_normal",
            "min_compartment_cells": 200,
        },
    )
    treg = out.set_index("cell_type").loc["Treg"]
    assert treg.glm_log2fc > 1 and treg.glm_padj < 0.01 and treg.consistent == "higher in tumour"
    assert treg.k_studies_same_direction == 2


# ------------------------------------------------------------------ §3.12 CNV helpers
def test_cross_fit_consensus_never_uses_own_half():
    from lusc_atlas.cnv import _cross_fit_corr

    rng = np.random.default_rng(0)
    prof = rng.normal(size=50)
    tumour = prof + rng.normal(scale=0.3, size=(60, 50))
    normal = rng.normal(scale=0.3, size=(60, 50))
    X = np.vstack([tumour, normal]).astype(np.float32)
    score = (X**2).mean(1)
    corr = _cross_fit_corr(X, score, np.random.default_rng(1), 0.1)
    assert np.nanmean(corr[:60]) > 0.8 and abs(np.nanmean(corr[60:])) < 0.2


def test_tissue_counts_region_unknown_as_tumour():
    from lusc_atlas.cnv import tissue_of

    t = tissue_of(pd.Series(["tumour", "tumour_region_unknown", "adjacent_normal", "healthy"]))
    assert t.tolist() == ["tumour", "tumour", "adjacent_normal", "healthy"]


# ------------------------------------------------------------------ §3.12 tensor scope
def test_tensor_scope_drops_tumour_only_types_and_incomplete_pairs():
    from lusc_atlas.ccc import balanced_tensor_scope

    rows = []
    for d in ["d1", "d2", "d3"]:
        for tis in ["tumour", "adjacent_normal"]:
            types = ["T", "NK", "Mac"] + (["TumourEpi"] if tis == "tumour" else [])
            if d == "d3" and tis == "adjacent_normal":
                types = ["T", "Mac"]  # sorting protocol lost NK in this slice
            for a in types:
                for b in types:
                    rows.append({"unit": f"{d}|{tis}", "source": a, "target": b})
    res = pd.DataFrame(rows)
    obs = pd.DataFrame({"unit": res["unit"].unique()})
    obs["tissue_type"] = obs["unit"].str.split("|").str[1]
    out = balanced_tensor_scope(
        res,
        obs,
        "unit",
        "tissue_type",
        ["tumour", "adjacent_normal"],
        {"tensor_type_presence": 0.6, "tensor_min_donors": 2},
    )
    assert "TumourEpi" not in set(out["source"])
    assert set(out["unit"].str.split("|").str[0]) == {"d1", "d2"}  # d3 dropped as a pair


def test_studies_with_all_zero_mito_become_unavailable():
    from lusc_atlas.select import mask_studies_without_mito

    pct = np.array([0, 0, 0, 1.5, 0, 12.0], dtype=np.float32)
    study = pd.Series(["Wu", "Wu", "Wu", "Leader", "Leader", "Leader"])
    out, src = mask_studies_without_mito(pct, "author obs column", study)
    assert np.isnan(out[:3]).all() and out[3] == 1.5 and out[4] == 0.0
    assert src[0].startswith("unavailable") and src[3] == "author obs column"


def test_unscored_author_malignant_cells_are_not_called_cnv_low():
    import anndata as ad

    from lusc_atlas.de import CNV_UNRESOLVED, TUMOUR_EPITHELIUM, make_de_groups

    obs = pd.DataFrame(
        {
            "harmonised_label": ["Unknown"] * 3 + ["NK"],
            "cell_type": ["malignant cell"] * 3 + ["natural killer cell"],
            "condition": ["tumour"] * 4,
            "malignant": [True, False, False, False],
            "cnv_status": ["CNV-high", "not called", "unresolved (< 50 epithelial cells)", "not run"],
        },
        index=list("abcd"),
    )
    a = ad.AnnData(obs=obs)
    make_de_groups(a, "harmonised_label", cnv_low_label="Tumour epithelial (CNV-low)")
    assert a.obs["de_group"].astype(str).tolist() == [
        "Malignant",
        "Tumour epithelial (CNV-low)",
        CNV_UNRESOLVED,
        "NK",
    ]
    assert (a.obs["analysis_group"].astype(str)[:3] == TUMOUR_EPITHELIUM).all()


def test_composition_paired_categorical_inputs_fit_nb():
    """Categorical donor/study columns (as in the pipeline) must not add empty design columns."""
    from lusc_atlas.de import composition_paired

    rng = np.random.default_rng(0)
    rows = []
    for i in range(8):
        for tis in ("tumour", "adjacent_normal"):
            for ct, base in (("A", 100), ("B", 60), ("C", 40)):
                lam = base * (2.0 if (ct == "A" and tis == "tumour") else 1.0)
                rows.append(
                    {
                        "unit": f"d{i}|{tis}",
                        "donor_id": f"d{i}",
                        "study": f"s{i % 2}",
                        "tissue_type": tis,
                        "compartment": "Immune",
                        "cell_type": ct,
                        "n_cells": int(rng.poisson(lam)),
                    }
                )
    cnt = pd.DataFrame(rows)
    extra = ["d_unused", "s_unused"]
    cnt["donor_id"] = pd.Categorical(cnt["donor_id"], categories=sorted(cnt["donor_id"].unique()) + extra[:1])
    cnt["study"] = pd.Categorical(cnt["study"], categories=["s0", "s1", extra[1]])
    acfg = {
        "compartment": "Immune",
        "case": "tumour",
        "control": "adjacent_normal",
        "min_compartment_cells": 50,
    }
    out = composition_paired(cnt, acfg).set_index("cell_type")
    assert (out["k_studies"] == 2).all()
    assert out.loc["A", "glm_log2fc"] > 0


def test_pseudobulk_restores_mixed_assay_covariate():
    """A unit pooling two chemistries keeps an `assay` value (decoupler drops varying columns)."""
    import anndata as ad

    from lusc_atlas.de import pseudobulk

    rng = np.random.default_rng(1)
    n = 120
    obs = pd.DataFrame(
        {
            "unit": ["u1"] * 60 + ["u2"] * 60,
            "grp": ["T"] * n,
            "donor_id": ["d1"] * 60 + ["d2"] * 60,
            "assay": ["v2"] * 30 + ["v3"] * 30 + ["v2"] * 60,
            "study": ["s"] * n,
        },
        index=[f"c{i}" for i in range(n)],
    )
    X = sp.csr_matrix(rng.poisson(5, size=(n, 50)).astype(np.float32))
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(50)]))
    a.layers["counts"] = X.copy()
    p = pseudobulk(a, {"pseudobulk_sample_key": "unit", "group_key": "grp", "min_cells_per_pseudobulk": 10})
    got = dict(zip(p.obs["unit"].astype(str), p.obs["assay"].astype(str), strict=True))
    assert got == {"u1": "v2 + v3", "u2": "v2"}


def test_sex_check_both_markers_absent_is_undetermined():
    import anndata as ad

    from lusc_atlas.de import sex_check

    X = np.array([[0, 0, 50], [10, 0, 40], [0, 8, 30]], dtype=float)
    obs = pd.DataFrame(
        {"donor_id": ["a", "b", "c"], "tissue_type": ["tumour"] * 3, "sex": ["male", "female", "male"]},
        index=["u1", "u2", "u3"],
    )
    p = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=["XIST", "RPS4Y1", "ACTB"]))
    out = sex_check(p, "tissue_type").set_index("donor_id")
    assert out.loc["a", "inferred_sex"] == "undetermined"
    assert out.loc["b", "inferred_sex"] == "female" and out.loc["c", "inferred_sex"] == "male"
