"""Unit tests for the §3.10 fixes of run 2 (Milo contrast, CNV-low tumour label, c2l reference labels).

Synthetic data only; no model fitting.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from lusc_atlas.config import load_config

CNV_LOW = "Tumour epithelial (CNV-low)"


def _obs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "harmonised_label": [
                "Unknown",
                "Unknown",
                "Unknown",
                "Unknown",
                "Basal",
                "Unknown",
                "Fibroblast",
            ],
            "cell_type": [
                "malignant cell",  # tumour, CNV-high          → Malignant
                "malignant cell",  # tumour, CNV-low           → CNV-low group
                "malignant cell",  # adjacent tissue, CNV-low → stays Unknown (no tumour label there)
                "unknown",  # unannotated HLCA cell       → stays Unknown
                "respiratory basal cell",
                "malignant cell",  # tumour_region_unknown, CNV-low → CNV-low group
                "fibroblast",
            ],
            "condition": [
                "tumour",
                "tumour",
                "adjacent_normal",
                "healthy",
                "healthy",
                "tumour_region_unknown",
                "tumour",
            ],
            "malignant": [True, False, False, False, False, False, False],
        },
        index=[f"c{i}" for i in range(7)],
    )


# ------------------------------------------------------------------ #2 de_group split
def test_cnv_low_tumour_cells_leave_unknown():
    from lusc_atlas.de import make_de_groups

    a = ad.AnnData(obs=_obs())
    make_de_groups(a, "harmonised_label", cnv_low_label=CNV_LOW)
    assert list(a.obs["de_group"].astype(str)) == [
        "Malignant",
        CNV_LOW,
        "Unknown",
        "Unknown",
        "Basal",
        CNV_LOW,
        "Fibroblast",
    ]


def test_de_groups_unchanged_without_cnv_low_label():
    """Backwards compatible: run-1 behaviour when the option is off (Malignant definition unchanged)."""
    from lusc_atlas.de import make_de_groups

    a = ad.AnnData(obs=_obs())
    make_de_groups(a, "harmonised_label")
    assert list(a.obs["de_group"].astype(str))[:4] == ["Malignant", "Unknown", "Unknown", "Unknown"]


def test_config_declares_cnv_low_group_tumour_only():
    cfg = load_config()
    label = cfg["cnv"]["cnv_low_label"]
    assert label in cfg["de"]["tumour_only_groups"] and "Malignant" in cfg["de"]["tumour_only_groups"]
    assert label in cfg["de"]["vs_healthy_basal"].values()
    # run 3 (§3.12): scCODA replaced by paired within-compartment tests; the mixed bucket is never tested
    assert "Unknown" in cfg["abundance"]["exclude_groups"] and "Unknown" in cfg["de"]["exclude_groups"]


# ------------------------------------------------------------------ #1 Milo contrast
def _milo_obs() -> pd.DataFrame:
    rows = []
    for donor, conds in {
        "d1": ["tumour", "adjacent_normal"],
        "d2": ["tumour", "adjacent_normal"],
        "d3": ["tumour"],  # unpaired → excluded when paired
        "h1": ["healthy"],
        "x1": ["tumour_region_unknown"],
    }.items():
        for cond in conds:
            rows += [{"donor_id": donor, "condition": cond, "sample": f"{donor}_{cond}"}] * 3
    return pd.DataFrame(rows)


def test_milo_design_paired_explicit_contrast():
    from lusc_atlas.de import milo_design

    obs = _milo_obs()
    keep, design, contrast = milo_design(
        obs, "condition", {"case": "tumour", "control": "adjacent_normal", "paired": True}
    )
    assert design == "~donor_id + condition"
    assert contrast == "conditiontumour - conditionadjacent_normal"
    assert set(obs.loc[keep, "donor_id"]) == {"d1", "d2"}
    assert set(obs.loc[keep, "condition"]) == {"tumour", "adjacent_normal"}


def test_milo_design_unpaired_keeps_only_case_and_control():
    from lusc_atlas.de import milo_design

    obs = _milo_obs()
    keep, design, _ = milo_design(obs, "condition", {"case": "tumour", "control": "adjacent_normal"})
    assert design == "~condition"
    assert set(obs.loc[keep, "condition"]) == {"tumour", "adjacent_normal"}
    assert "d3" in set(obs.loc[keep, "donor_id"])


def test_milo_contrast_matches_pertpy_coefficients():
    """The contrast string must resolve to the tumour coefficient of the donor-blocked design in the
    installed pertpy (guards against the run-1 failure: an implicit 'last level' test)."""
    pytest.importorskip("pertpy")
    from formulaic import model_matrix
    from pertpy.tools._milo import _contrast_vector

    from lusc_atlas.de import milo_design

    obs = _milo_obs()
    keep, design, contrast = milo_design(
        obs, "condition", {"case": "tumour", "control": "adjacent_normal", "paired": True}
    )
    df = obs[keep].drop_duplicates("sample").set_index("sample")[["donor_id", "condition"]]
    df["condition"] = pd.Categorical(df["condition"], categories=["adjacent_normal", "tumour"])
    cols = list(model_matrix(design + " + 0", df).columns)
    w = _contrast_vector(cols, contrast, {"donor_idd1", "conditionadjacent_normal"})
    assert w.sum() == 1 and w[[i for i, c in enumerate(cols) if "tumour" in c][0]] == 1


# ------------------------------------------------------------------ #3 c2l reference labels
def test_reference_labels_drop_small_and_excluded():
    from lusc_atlas.deconvolve import prepare_reference_labels

    lab = pd.Series(["A"] * 300 + ["Rare"] * 66 + ["Unknown"] * 500 + ["B"] * 250)
    new, dec = prepare_reference_labels(lab, min_cells=200, exclude=["Unknown"])
    assert set(new.dropna()) == {"A", "B"}
    assert new.isna().sum() == 66 + 500
    actions = dec.set_index("label")["action"].to_dict()
    assert actions == {"Unknown": "excluded", "A": "kept", "B": "kept", "Rare": "dropped (< min_ref_cells)"}


def test_reference_labels_merge_into_parent():
    from lusc_atlas.deconvolve import prepare_reference_labels

    lab = pd.Series(["A"] * 300 + ["Rare"] * 66)
    new, dec = prepare_reference_labels(lab, min_cells=200, merge={"Rare": "A"})
    assert (new == "A").sum() == 366
    assert dec.set_index("label").loc["Rare", "target"] == "A"
    with pytest.raises(ValueError):
        prepare_reference_labels(lab, min_cells=200, merge={"Rare": "Missing"})


def test_reference_labels_noop_by_default():
    from lusc_atlas.deconvolve import prepare_reference_labels

    lab = pd.Series(np.repeat(["A", "B"], [5, 3]))
    new, dec = prepare_reference_labels(lab)
    assert new.tolist() == lab.tolist() and set(dec.action) == {"kept"}
