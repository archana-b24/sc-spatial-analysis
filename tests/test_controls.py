"""§3.13 control-group decisions: healthy allow-list and the adjacent-sample description."""

import pandas as pd

from lusc_atlas.de import adjacent_makeup
from lusc_atlas.select import healthy_source_ok

REQ = {"lung_condition": ["Healthy"], "tissue_sampling_method": ["donor_lung", "lung_explant", "autopsy"]}


def test_healthy_allow_list_keeps_only_cancer_free_lung_tissue():
    obs = pd.DataFrame(
        {
            "lung_condition": ["Healthy", "Healthy", "Healthy (tumor adjacent)", "Healthy"],
            "tissue_sampling_method": ["donor_lung", "balf", "surgical_resection", "autopsy"],
        }
    )
    assert healthy_source_ok(obs, REQ).tolist() == [True, False, False, True]


def test_healthy_allow_list_missing_column_never_qualifies():
    obs = pd.DataFrame({"lung_condition": ["Healthy", "Healthy"]})
    assert not healthy_source_ok(obs, REQ).any()


def test_adjacent_makeup_only_paired_patients_and_shares():
    rows = []
    # patient p1: tumour + adjacent (paired); patient p2: tumour only (not paired)
    for d, tissue, comp, lab, grp, cnv in [
        ("p1", "adjacent_normal", "Immune", "Monocyte classical", "Monocyte classical", "not run"),
        ("p1", "adjacent_normal", "Immune", "CD8 T", "CD8 T", "not run"),
        ("p1", "adjacent_normal", "Epithelial", "AT2", "AT2", "not called"),
        (
            "p1",
            "adjacent_normal",
            "Epithelial",
            "Basal",
            "Basal",
            "CNV-high (adjacent: false-positive check)",
        ),
        ("p1", "tumour", "Immune", "CD8 T", "CD8 T", "not run"),
        ("p2", "tumour", "Immune", "CD8 T", "CD8 T", "not run"),
    ]:
        rows.append(
            {
                "study": "S",
                "donor_id": d,
                "tissue_type": tissue,
                "compartment": comp,
                "harmonised_label": lab,
                "de_group": grp,
                "cnv_status": cnv,
            }
        )
    out = adjacent_makeup(pd.DataFrame(rows), "tissue_type", "tumour", "adjacent_normal", "de_group")
    assert out["donor_id"].tolist() == ["p1"]
    r = out.iloc[0]
    assert r.n_cells == 4 and r.share_immune == 0.5 and r.share_epithelial == 0.5
    assert r.share_airway_of_epithelium == 0.5
    assert r.share_blood_myeloid_of_immune == 0.5
    assert r.cnv_call_rate == 0.5
