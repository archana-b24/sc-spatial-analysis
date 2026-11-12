import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.plots import COMPARTMENT_COLOR, heatmap, umap_categorical
from lib.ui import (
    ATLAS_COLOR,
    CONDITION_COLOR,
    GREY,
    SERIES,
    caveats,
    decisions,
    footer,
    page_header,
    show_table,
    style,
)

page_header(
    "3 · Integration",
    "Cells from 20 studies and two atlases differ for technical reasons. Which method removes those batch "
    "effects without also removing the tumour biology we want to study?",
)
guide("integration")
decisions("integration")
caveats("integration")

bm = q("SELECT * FROM integration_benchmark ORDER BY Total DESC")
bm["method"] = bm.embedding.str.replace("X_", "", regex=False)
chosen = bm.loc[bm.chosen, "method"].iloc[0]
st.success(f"**Chosen: {chosen}.** Best scIB total score among methods that pass the kNN label-purity guard.")

c1, c2 = st.columns(2)
with c1:
    st.subheader("scIB total score")
    colors = [
        SERIES[0] if m == chosen else (GREY if b else "#86b6ef")
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
            hovertemplate="%{y}: %{x:.3f}<extra></extra>",
        )
    )
    fig.update_layout(
        yaxis=dict(autorange="reversed"), xaxis=dict(range=[0, 0.8], title="Total (0.6 bio + 0.4 batch)")
    )
    st.plotly_chart(style(fig, 360))
    st.caption("Grey = benchmark-only candidate (zero-shot Geneformer on a 20k-cell subsample).")
with c2:
    st.subheader("Biology kept vs batch removed")
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
                    color=SERIES[0] if m == chosen else "#86b6ef",
                    line=dict(width=2, color="white"),
                ),
                showlegend=False,
                hovertemplate=f"{m}<br>batch %{{x:.3f}}<br>bio %{{y:.3f}}<extra></extra>",
            )
        )
    fig.update_xaxes(title="Batch correction")
    fig.update_yaxes(title="Bio conservation")
    st.plotly_chart(style(fig, 360))
show_table(bm.drop(columns=["method"]), "scib-metrics on a stratified subsample; batch = study")

st.subheader("The integrated atlas (UMAP)")
cells = q("""SELECT umap_1, umap_2, atlas, study, condition, harmonised_label, compartment
             FROM cells WHERE in_subsample""")
by = st.radio("Colour by", ["atlas", "condition", "compartment", "one cell type"], horizontal=True)
if by == "one cell type":
    ct = st.selectbox(
        "Cell type",
        sorted(cells.harmonised_label.unique()),
        index=sorted(cells.harmonised_label.unique()).index("Basal"),
    )
    fig = umap_categorical(cells, "harmonised_label", highlight=ct)
else:
    cmap = {"atlas": ATLAS_COLOR, "condition": CONDITION_COLOR, "compartment": COMPARTMENT_COLOR}[by]
    fig = umap_categorical(cells, by, colors=cmap)
st.plotly_chart(fig)
st.caption(
    f"{len(cells):,}-cell stratified subsample (every cell type × condition represented). 27 cell types are too "
    "many to tell apart by colour, so pick one to highlight; the Cell & gene explorer has more options."
)

st.subheader("Do the two atlases mix within each cell type?")
mx = q("SELECT * FROM atlas_mixing_per_label ORDER BY frac_other_atlas_neighbours")
fig = go.Figure()
fig.add_trace(
    go.Bar(
        y=mx.label, x=mx.frac_other_atlas_neighbours, orientation="h", name="observed", marker_color=SERIES[0]
    )
)
fig.add_trace(
    go.Scatter(
        y=mx.label,
        x=mx.expected_if_random,
        mode="markers",
        name="expected if random",
        marker=dict(symbol="line-ns", size=14, line=dict(width=3, color=SERIES[1])),
    )
)
fig.update_xaxes(title="fraction of neighbours from the other atlas")
st.plotly_chart(style(fig, 620))

st.subheader("Label transfer check (scANVI prediction vs harmonised label)")
cm = q("SELECT * FROM scanvi_label_confusion").set_index("harmonised_label")
st.plotly_chart(
    heatmap(cm.div(cm.sum(axis=1), axis=0), "row fraction", diverging=False, height=640, fmt=".2f")
)
footer()
