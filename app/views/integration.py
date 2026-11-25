"""Integration: 7-method scIB benchmark (scANVI scored on held-out labels), the shared map, atlas mixing
over all cells and over tumour-patient tissue only, held-out label transfer (DEVELOPMENT_LOG §3.3 F–H, §3.12)."""

import plotly.graph_objects as go
import streamlit as st
from lib.db import config, q
from lib.explain import guide
from lib.plots import heatmap, umap_categorical
from lib.ui import (
    ATLAS_COLOR,
    COMPARTMENT_COLOR,
    CONDITION_COLOR,
    GREY,
    INK_2,
    caveats,
    decisions,
    figure_legend,
    footer,
    methods_box,
    page_header,
    show_table,
    style,
)

page_header(
    "Integration",
    "Cells from 20 studies (19 after QC) and two atlases differ for technical reasons. Which method removes those batch "
    "effects without also removing the biology we want to study?",
    part="built",
)
guide("integration")
caveats("integration")
icfg = config().get("integration", {})

st.subheader("What integration is for — and what it is not used for")
st.markdown(
    "Integration builds a shared low-dimensional **map** in which the same cell type from different studies "
    "sits together. The map is used for neighbours, UMAP pictures, clustering, label transfer, Milo "
    "neighbourhoods and the trajectory. It is **never** used to compare gene expression: corrected values are "
    "not counts, so every DE test uses the raw counts with patient in the model. Integration removes variation "
    "between **samples** (the batch key, the finest technical unit) with study and assay as covariates; the "
    "benchmark then judges how well **studies** mix."
)
c1, c2 = st.columns(2)
with c1:
    methods_box(
        "the candidates",
        "PCA (no correction); Harmony (adjusts principal components per batch); scVI (variational autoencoder on "
        "counts); scANVI (scVI + cell-type labels); sysVI (scVI for large system differences); scArches (maps "
        "cells onto the published HLCA model); Geneformer (foundation model, zero-shot, benchmark subsample only).",
        "Covers the classical baseline, the deep generative standard, label-guided, cross-system, transfer-learning "
        "and foundation-model approaches.",
        "Scanorama, BBKNN, LIGER (lower in the scIB benchmark for atlas-scale data).",
        "Korsunsky 2019; Lopez 2018; Xu 2021; Hrovatin 2024; Lotfollahi 2022; Theodoris 2023.",
    )
with c2:
    methods_box(
        "the benchmark (scIB) and its guard",
        "On a shared 20,000-cell subsample: bio-conservation (do labelled cell types stay apart?) and batch "
        "correction (do studies mix within a type?), combined as 0.6 × bio + 0.4 × batch; plus a kNN label-purity "
        "guard — a method may not blur cell types more than no correction does.",
        "scANVI uses cell-type labels during training; since run 3 the labels of the 20,000 benchmark cells are hidden "
        "from it, so it is scored on cells it never saw labelled (runs 1–2 were circular).",
        "visual inspection only; batch metrics alone (reward over-correction).",
        "Luecken et al. 2022 Nat Methods.",
    )
decisions("integration")

bm = q("SELECT * FROM integration_benchmark ORDER BY Total DESC")
bm["method"] = bm.embedding.str.replace("X_", "", regex=False)
chosen = bm.loc[bm.chosen, "method"].iloc[0]
st.success(f"**Chosen: {chosen}** — best scIB total among methods that pass the label-purity guard.")

c1, c2 = st.columns(2)
with c1:
    colors = [
        "#1c5cab" if m == chosen else (GREY if b else "#86b6ef")
        for m, b in zip(bm.method, bm.benchmark_only, strict=True)
    ]
    fig = go.Figure(
        go.Bar(
            x=bm.Total,
            y=bm.method,
            orientation="h",
            marker_color=colors,
            text=bm.Total.round(3),
            textposition="outside",
        )
    )
    fig.update_layout(yaxis=dict(autorange="reversed"), xaxis=dict(range=[0, float(bm.Total.max()) * 1.2]))
    st.plotly_chart(style(fig, 340, x="scIB total score (0.6 bio + 0.4 batch)"))
with c2:
    fig = go.Figure()
    for m, bx, by_ in zip(bm.method, bm["Batch correction"], bm["Bio conservation"], strict=True):
        fig.add_trace(
            go.Scatter(
                x=[bx],
                y=[by_],
                mode="markers+text",
                text=[m],
                textposition="top center",
                marker=dict(
                    size=14 if m == chosen else 10,
                    color="#1c5cab" if m == chosen else "#86b6ef",
                    line=dict(width=2, color="white"),
                ),
                showlegend=False,
            )
        )
    st.plotly_chart(style(fig, 340, x="batch correction score", y="bio-conservation score"))
figure_legend(
    1,
    f"{chosen} balances keeping cell types apart with mixing studies best.",
    "left: overall scIB score per method (grey = benchmark-only Geneformer); right: the two components.",
    "top-right in the scatter is best; a method far right but low has over-corrected (mixed away biology).",
    f"{icfg.get('benchmark_n_cells', 20000):,} cells, same subsample for every method; batch = study.",
    "scib-metrics; kNN label purity guard (tolerance 0.05 vs unintegrated PCA).",
    "the chosen map is used for neighbourhood-based analyses only, never for DE.",
)
show_table(bm.drop(columns=["method"]).round(3), "All scIB metrics.")

st.subheader("The shared map")
cells = q("""SELECT umap_1, umap_2, atlas, study, condition, harmonised_label, compartment
             FROM cells WHERE in_subsample""")
by = st.radio("Colour by", ["atlas", "condition", "compartment", "one cell type"], horizontal=True)
if by == "one cell type":
    labs = sorted(cells.harmonised_label.unique())
    ct = st.selectbox("Cell type", labs, index=labs.index("Basal") if "Basal" in labs else 0)
    fig = umap_categorical(cells, "harmonised_label", highlight=ct)
else:
    cmap = {"atlas": ATLAS_COLOR, "condition": CONDITION_COLOR, "compartment": COMPARTMENT_COLOR}[by]
    fig = umap_categorical(cells, by, colors=cmap)
st.plotly_chart(fig)
n_types = cells.harmonised_label.nunique()
figure_legend(
    2,
    "Cells from both atlases on one map.",
    f"UMAP of the chosen embedding; a stratified subsample of {len(cells):,} cells (every cell type × condition "
    "represented).",
    "nearby cells have similar expression; distances between far-apart clusters are not meaningful. Colouring "
    "by atlas mostly also colours by condition, because healthy donors are all HLCA.",
    f"{n_types} harmonised cell types (too many to colour at once — pick one to highlight).",
    f"UMAP on {chosen} latent space (30 dimensions).",
    "same cell types from both atlases should overlap.",
)

st.subheader("Do the two atlases mix within each cell type?")
mx = q("SELECT * FROM atlas_mixing_per_label")
scopes = mx["scope"].unique().tolist() if "scope" in mx else ["all cells"]
sc = st.radio("Cells compared", scopes, horizontal=True, index=len(scopes) - 1)
m = mx[mx["scope"] == sc] if "scope" in mx else mx
m = m.sort_values("frac_other_atlas_neighbours")
fig = go.Figure()
fig.add_trace(
    go.Bar(
        y=m.label, x=m.frac_other_atlas_neighbours, orientation="h", name="observed", marker_color="#86b6ef"
    )
)
fig.add_trace(
    go.Scatter(
        y=m.label,
        x=m.expected_if_random,
        mode="markers",
        name="expected if perfectly mixed",
        marker=dict(symbol="line-ns", size=14, line=dict(width=3, color=INK_2)),
    )
)
st.plotly_chart(
    style(fig, 120 + 24 * len(m), x="fraction of a cell's 30 nearest neighbours from the other atlas")
)
figure_legend(
    3,
    "The two atlases stay partly separated within every cell type.",
    "for each harmonised cell type, the share of neighbours from the other atlas (bar) vs the share expected if "
    "the atlases were perfectly mixed (tick).",
    "'all cells' compares healthy HLCA cells with tumour-patient LuCA cells, so bars below the ticks partly "
    "reflect real biology. 'Tumour-patient tissue only' compares like with like (HLCA's Thienpont patients vs "
    "LuCA patients); the gap there is closer to the technical atlas effect, but HLCA contributes few cells.",
    (
        f"{int(m.n_cells_hlca.sum()):,} HLCA and {int(m.n_cells_luca.sum()):,} LuCA cells in this view."
        if "n_cells_hlca" in m
        else ""
    ),
    "kNN (k = 30) on the chosen embedding; types with < 50 cells or one atlas only omitted.",
    "the atlases are not perfectly mixed; cross-atlas comparisons are treated as confounded throughout.",
)

st.subheader("Label transfer on cells scANVI never saw labelled")
cm = q("SELECT * FROM scanvi_label_confusion").set_index("harmonised_label")
acc = float(sum(cm.loc[i, i] for i in cm.index if i in cm.columns) / cm.to_numpy().sum())
st.plotly_chart(
    heatmap(cm.div(cm.sum(axis=1), axis=0), "row fraction", diverging=False, height=640, fmt=".2f")
)
figure_legend(
    4,
    f"scANVI recovers held-out labels for {acc:.0%} of cells.",
    "rows = harmonised (author-derived) label, columns = scANVI's prediction, for the benchmark cells whose labels "
    "were hidden during training; each row sums to 1.",
    "a strong diagonal = labels are recoverable from expression; off-diagonal blocks show which types are "
    "confused (often related subtypes).",
    f"{int(cm.to_numpy().sum()):,} held-out cells with a known label.",
    "scANVI trained with these cells' labels masked as unknown.",
    "an honest label-transfer accuracy (runs 1–2 checked scANVI against its own training labels).",
)
footer()
