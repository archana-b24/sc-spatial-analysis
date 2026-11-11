"""Cross-atlas merge logic (DEVELOPMENT_LOG §3.3 A–E) on synthetic metadata."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from lusc_atlas.de import design_for
from lusc_atlas.select import assign_condition, is_count_matrix, keys_overlap, stratified_cap, study_key

CCFG = {
    "tumour_disease_ontology": "MONDO:0005097",
    "normal_disease": "normal",
    "origin_map": {
        "tumour": ["tumor_primary"],
        "adjacent_normal": ["normal_adjacent"],
        "exclude": ["tumor_metastasis"],
    },
}


def test_same_deposit_named_differently_is_detected():
    assert keys_overlap(study_key("Thienpont_2018"), study_key("Lambrechts_Thienpont_2018_6653"))
    assert not keys_overlap(
        study_key("Thienpont_2018"), study_key("Thienpont_2021")
    )  # same author, other year
    assert not keys_overlap(study_key("Travaglini_2020"), study_key("Kim_Lee_2020"))


def test_condition_uses_origin_and_imports_missing_origin_by_sample():
    obs = pd.DataFrame(
        {
            "disease_ontology_term_id": ["MONDO:0005097"] * 4 + ["PATO:0000461"],
            "disease": ["squamous cell lung carcinoma"] * 4 + ["normal"],
            "study": ["Thienpont_2018"] * 5,
            "sample": ["BT1001", "BT1002", "BT1003", "BT1004", "h1"],
        }
    )
    k = study_key("Lambrechts_Thienpont_2018_6653")  # same deposit, other atlas' naming
    lookup = {
        (k, "1001"): "tumor_primary",
        (k, "1002"): "normal_adjacent",
        (k, "1003"): "tumor_metastasis",
        (study_key("Other_2018"), "1004"): "tumor_primary",
    }  # same id, different study → must NOT match
    acfg = {"origin_col": None, "sample_col": "sample", "study_col": "study"}
    cond = assign_condition(obs, acfg, CCFG, lookup)
    assert cond.tolist() == ["tumour", "adjacent_normal", "excluded", "tumour_region_unknown", "healthy"]


def test_count_check_rejects_normalised_data():
    counts = sp.csr_matrix(np.random.default_rng(0).poisson(1, (50, 20)).astype(np.float32))
    lognorm = counts.copy()
    lognorm.data = np.log1p(lognorm.data)
    assert is_count_matrix(counts)
    assert not is_count_matrix(lognorm)


def test_stratified_cap_respects_cap_and_keeps_rare_types():
    obs = pd.DataFrame({"donor_id": ["d1"] * 1000, "ct": ["common"] * 995 + ["rare"] * 5})
    idx = stratified_cap(obs, 100, "ct", seed=0)
    assert len(idx) <= 101
    assert "rare" in set(obs.iloc[idx]["ct"])


def test_paired_design_blocks_on_donor():
    meta = pd.DataFrame({"study": ["a"] * 4, "condition": ["tumour", "adjacent_normal"] * 2})
    assert design_for(meta, "condition", paired=True) == "~donor_id + condition"


def test_sample_core_id_matches_same_sample_only():
    from lusc_atlas.select import sample_core_id as f

    assert f("BT1247") == f("Lambrechts_Thienpont_2018_6149v1_1247") == "1247"
    assert f("scrBT1425") == f("Lambrechts_Thienpont_2018_6653_scrBT1425_hg19") == "1425"
    assert f("Lambrechts_Thienpont_2018_6149v1_BT1A") == ""  # batch number in prefix must not leak
