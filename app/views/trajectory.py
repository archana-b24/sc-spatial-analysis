"""Trajectory: Palantir pseudotime from healthy basal cells, per-patient summaries, label-free scVI
sensitivity; descriptive, atlas-confounded (DEVELOPMENT_LOG §3.7, §3.12)."""

import json

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import has_table, meta, q
from lib.explain import guide
from lib.plots import umap_categorical, umap_continuous
from lib.ui import (
    ATLAS_COLOR,
    BASAL_COLOR,
    MALIGNANT_COLOR,
    TISSUE_COLOR,
    caveats,
    cnv_gate_passed,
    decisions,
    figure_legend,
    footer,
    inference_badge,
    methods_box,
    page_header,
    show_table,
    style,
)

page_header(
    "Trajectory",
    "How far, in expression space, are malignant cells from the healthy basal cells that LUSC is thought to "
    "arise from?",
    part="context",
)
guide("trajectory")
caveats("trajectory")
inference_badge("hypothesis" if cnv_gate_passed() else "provisional")
if not cnv_gate_passed():
    st.error(
        "The endpoint of this trajectory is the set of CNV-called malignant cells, which are provisional: the "
        "pre-registered false-positive check failed (issue 16, Malignant cells page)."
    )
st.info(
    "**Distance, not lineage.** Tumour cells are clonal within each patient and there is no time course, so "
    "pseudotime measures how far a cell's expression is from healthy basal cells. The root (healthy basal) "
    "cells all come from HLCA and the malignant cells from LuCA, so part of the distance is technical."
)
methods_box(
    "Palantir + CellRank 2",
    "Palantir orders cells by diffusion distance from a root cell (here the healthy basal cell closest to the "
    "basal centroid); CellRank 2 combines that pseudotime with CytoTRACE (differentiation potential from the "
    "number of expressed genes) to estimate end states.",
    "Works without RNA velocity (the atlases have no unspliced counts) and gives a continuous measure.",
    "PAGA (graph only, used in v1); Slingshot/Monocle 3 (similar assumptions).",
    "Setty et al. 2019 Nat Biotechnol; Weiler et al. 2024 Nat Methods; Gulati et al. 2020 Science.",
)
decisions("trajectory")

lc = q("SELECT * FROM lineage_cells")
lc["harmonised_label"] = lc.lineage_state
state_color = {"Basal": BASAL_COLOR, "Malignant": MALIGNANT_COLOR}
c1, c2 = st.columns(2)
with c1:
    by = st.radio("Colour cells by", ["lineage state", "atlas", "tissue"], horizontal=True)
    if by == "atlas":
        fig = umap_categorical(lc, "atlas", colors=ATLAS_COLOR, height=440)
    elif by == "tissue":
        lc["tissue"] = lc.condition.replace({"tumour_region_unknown": "tumour"})
        fig = umap_categorical(lc, "tissue", colors=TISSUE_COLOR, height=440)
    else:
        fig = umap_categorical(lc, "lineage_state", colors=state_color, height=440)
    st.plotly_chart(fig)
with c2:
    st.radio("Shown", ["Palantir pseudotime"], horizontal=True, disabled=True)
    st.plotly_chart(umap_continuous(lc, lc.palantir_pseudotime, "pseudotime", height=440))
figure_legend(
    1,
    "Basal and malignant cells on their own map.",
    "UMAP of basal and CNV-called malignant cells (left, coloured by state, atlas or tissue) and their pseudotime "
    "(right; 0 = healthy basal root).",
    "colour by atlas to see how much of the separation coincides with HLCA vs LuCA.",
    f"{len(lc):,} cells.",
    "neighbour graph on the integrated embedding; Palantir pseudotime.",
    "malignant cells sit at higher pseudotime, but atlas and state coincide.",
)

st.subheader("Patient by patient")
if has_table("trajectory_per_donor"):
    pdn = q("SELECT * FROM trajectory_per_donor WHERE n_cells >= 20")
    fig = go.Figure()
    for s, color in state_color.items():
        d = pdn[pdn.lineage_state == s]
        fig.add_trace(
            go.Box(
                y=d.median_pseudotime,
                name=f"{s} ({d.donor_id.nunique()} donors)",
                marker_color=color,
                boxpoints="all",
                jitter=0.4,
                pointpos=0,
                text=d.donor_id + " · " + d.study,
                hovertemplate="%{text}<br>median %{y:.2f}<extra></extra>",
            )
        )
    fig.update_layout(showlegend=False)
    st.plotly_chart(style(fig, 360, x="cell state (donors with ≥ 20 cells)", y="median pseudotime per donor"))
    figure_legend(
        2,
        "Pseudotime with the donor as the unit.",
        "each dot is one donor's median pseudotime for basal or malignant cells (≥ 20 cells).",
        "if the malignant shift holds patient by patient, the malignant dots sit above the basal ones with little "
        "overlap; few malignant dots = few patients behind the result.",
        f"{pdn[pdn.lineage_state == 'Malignant'].donor_id.nunique()} malignant donors, "
        f"{pdn[pdn.lineage_state == 'Basal'].donor_id.nunique()} basal donors.",
        "Palantir pseudotime summarised per donor.",
        "descriptive: basal and malignant donors are different people from different atlases.",
    )
info = meta().get("info_trajectory_info")
if info:
    j = json.loads(info)
    rho = j.get("sensitivity_spearman")
    if rho is None:
        st.caption("Sensitivity analysis unavailable in this snapshot.")
    else:
        verdict = (
            " — the ordering is largely robust to the choice of map."
            if rho > 0.9
            else " — the ordering is largely, but not entirely, robust to the choice of map."
            if rho > 0.75
            else " — the ordering depends substantially on the map. Read the trajectory as a rough picture of "
            "which cells are closer to or further from healthy basal cells, not as a measurement."
        )
        st.markdown(
            "**Sensitivity to the map.** The main map (scANVI) uses cell-type labels during training. Recomputing "
            "pseudotime from the same root cell on the label-free scVI map gives a Spearman rank correlation of "
            f"ρ = **{rho:.2f}** with the main result (1 = identical order){verdict}"
        )
        st.caption(
            "Pseudotime is scaled so that the most distant cell is 1, so its values are relative: compare cells "
            "and patients within this snapshot, not numbers between pipeline runs."
        )
    ts = j.get("terminal_states") or []
    if ts:
        st.caption("CellRank end states: " + ", ".join(ts) + " (macrostates named after the dominant state).")
fig = go.Figure()
for s, d in lc.groupby("lineage_state"):
    fig.add_trace(
        go.Violin(
            x=d.palantir_pseudotime,
            name=f"{s} ({len(d):,} cells)",
            line_color=state_color.get(s, "#8a8984"),
            orientation="h",
            side="positive",
            meanline_visible=True,
            points=False,
        )
    )
st.plotly_chart(style(fig, 300, x="pseudotime (0 = healthy basal root)", y=None))
figure_legend(
    3,
    "Cell-level pseudotime distributions.",
    "density of pseudotime per state (cells, not donors).",
    "overlap means pseudotime separates the states only partly.",
    f"{len(lc):,} cells.",
    "Palantir.",
    "use Figure 2 for the donor-level view.",
)
show_table(q("SELECT * FROM trajectory_summary"), "Pseudotime summary per state (cells).")
_ = np
footer()
