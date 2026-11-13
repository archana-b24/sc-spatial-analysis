import base64
import io

import plotly.graph_objects as go
import streamlit as st
from lib.db import has_table, q
from lib.explain import guide
from lib.plots import heatmap
from lib.ui import (
    GREY,
    SERIES,
    caveats,
    decisions,
    footer,
    metric_row,
    page_header,
    series_color,
    show_table,
    style,
)
from PIL import Image

page_header(
    "8 · Spatial",
    "Where do the cell types found in the single-cell atlas sit in a real LUSC tissue section, and which "
    "neighbourhoods (niches) do they form?",
)
guide("spatial")
decisions("spatial")
caveats("spatial")
st.info(
    "**One tissue section (n = 1).** Everything on this page is descriptive; no comparison between patients is possible."
)

img = q("SELECT png, scalef, spot_diameter_fullres FROM visium_images WHERE resolution = 'hires'").iloc[0]
spots = q("SELECT * FROM visium_spots")
types = [c for c in q("SELECT * FROM c2l_abundance LIMIT 0").columns if c != "spot"]
metric_row(
    [
        ("Spots", f"{len(spots):,}", None),
        ("Median UMIs / spot", f"{spots.total_counts.median():,.0f}", None),
        ("Cell types deconvolved", f"{len(types)}", "cell2location, signatures from the single-cell atlas"),
        ("Niches", f"{spots.niche.nunique()}", "CellCharter"),
    ]
)


@st.cache_data
def _uri(png: bytes) -> str:
    return "data:image/png;base64," + base64.b64encode(png).decode()


w, h = Image.open(io.BytesIO(img.png)).size  # PNG bytes from the DB  noqa: source-of-truth
x, y = spots.x_fullres * img.scalef, spots.y_fullres * img.scalef

c1, c2 = st.columns([1, 3])
with c1:
    what = st.radio("Colour spots by", ["niche", "cell-type abundance", "gene"])
    if what == "cell-type abundance":
        ct = st.selectbox("Cell type", types, index=types.index("Fibroblast"))
        val, title = spots[ct], f"{ct} (abundance)"
    elif what == "gene":
        genes = q("SELECT DISTINCT gene FROM spot_expr ORDER BY gene").gene.tolist()
        g = st.selectbox("Gene", genes, index=genes.index("KRT5") if "KRT5" in genes else 0)
        e = q("SELECT spot, value FROM spot_expr WHERE gene = ?", (g,)).set_index("spot").value
        val, title = spots.spot.map(e).fillna(0), f"{g} (log-norm)"
    opacity = st.slider("Spot opacity", 0.0, 1.0, 0.8, 0.1)
    st.caption(
        "Gene list: top 200 spatially variable genes (Moran's I), spatial ligand–receptor genes and markers."
    )
fig = go.Figure()
fig.add_layout_image(
    source=_uri(img.png), x=0, y=0, sizex=w, sizey=h, xref="x", yref="y", sizing="stretch", layer="below"
)
size = max(3, img.spot_diameter_fullres * img.scalef * 0.9)
if what == "niche":
    for i, (n, d) in enumerate(spots.assign(x=x, y=y).groupby("niche")):
        fig.add_trace(
            go.Scattergl(
                x=d.x,
                y=d.y,
                mode="markers",
                name=f"niche {n} ({len(d):,})",
                marker=dict(size=size, color=series_color(i), opacity=opacity),
            )
        )
else:
    fig.add_trace(
        go.Scattergl(
            x=x,
            y=y,
            mode="markers",
            showlegend=False,
            text=spots.spot,
            marker=dict(
                size=size,
                color=val,
                opacity=opacity,
                colorscale=[[0, "#f0efec"], [1, "#0d366b"]],
                colorbar=dict(title=title, thickness=12),
            ),
            hovertemplate="%{text}<br>%{marker.color:.3f}<extra></extra>",
        )
    )
fig.update_xaxes(visible=False, range=[0, w])
fig.update_yaxes(visible=False, range=[h, 0], scaleanchor="x")
with c2:
    st.plotly_chart(style(fig, 640))

st.subheader("What is in each niche")
ni = q("SELECT * FROM spatial_niches").set_index("niche")
ni.index = [f"niche {i}" for i in ni.index]
fig = go.Figure()
top = ni.max().sort_values(ascending=False).index[:7]
for i, ct in enumerate(top):
    fig.add_trace(go.Bar(y=ni.index, x=ni[ct], name=ct, orientation="h", marker_color=series_color(i)))
fig.add_trace(
    go.Bar(y=ni.index, x=ni.drop(columns=top).sum(axis=1), name="other", orientation="h", marker_color=GREY)
)
fig.update_layout(barmode="stack", xaxis=dict(title="mean fraction of spot", tickformat=".0%"))
st.plotly_chart(style(fig, 300))
show_table(ni.reset_index(names="niche"))
if has_table("c2l_reference_labels"):
    rl = q("SELECT label, n_cells, action, target FROM c2l_reference_labels ORDER BY n_cells DESC")
    off = rl[rl.action != "kept"]
    st.caption(
        f"cell2location reference: {int((rl.action == 'kept').sum())} cell types with a signature. Not used: "
        + ", ".join(f"*{r.label}* ({int(r.n_cells):,} cells, {r.action})" for r in off.itertuples())
        + ". Tiny reference types give noisy signatures that soak up unexplained signal (in run 1 a 66-cell type "
        "took 66 % of the tumour niche), and a mixed bucket has no coherent signature (DEVELOPMENT_LOG §3.11)."
    )
    show_table(rl, "Reference labels: kept, merged, dropped or excluded")

c1, c2 = st.columns(2)
with c1:
    st.subheader("Co-localisation of cell types")
    co = q("SELECT * FROM spatial_colocalisation").set_index("cell_type")
    st.plotly_chart(heatmap(co, "correlation", height=640))
with c2:
    st.subheader("Top spatially variable genes")
    svg = q('SELECT gene, I AS "Moran I", pval_norm_fdr_bh AS fdr FROM spatial_svg ORDER BY I DESC LIMIT 25')
    fig = go.Figure(go.Bar(y=svg.gene, x=svg["Moran I"], orientation="h", marker_color=SERIES[0]))
    fig.update_layout(yaxis=dict(autorange="reversed"), xaxis=dict(title="Moran's I"))
    st.plotly_chart(style(fig, 640))
footer()
