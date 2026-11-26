import base64
import io

import plotly.graph_objects as go
import streamlit as st
from lib.db import has_table, q
from lib.explain import guide
from lib.plots import heatmap
from lib.ui import (
    GREY,
    INK_2,
    MALIGNANT_COLOR,
    caveats,
    decisions,
    figure_legend,
    footer,
    inference_badge,
    methods_box,
    metric_row,
    page_header,
    series_color,
    show_table,
    style,
)
from PIL import Image

page_header(
    "Spatial",
    "Where do the cell types found in the single-cell atlas sit in a real LUSC tissue section, and which "
    "neighbourhoods (niches) do they form?",
    part="context",
)
guide("spatial")
caveats("spatial")
inference_badge("descriptive")
# the one-section caveat is known issue 15, shown above by caveats(); not repeated here

# ------------------------------------------------------------------------------------------------ identity
st.subheader("Is this section really a squamous cell carcinoma?")
st.markdown(
    "Runs 1 and 2 analysed 10x Genomics' 'CytAssist 11 mm FFPE Human Lung Cancer' section, which turned out to be "
    "a **neuroendocrine** carcinoma (chromogranin A was about 10 % of every spot's RNA; squamous markers were "
    "absent). Since run 3 the site uses 10x's lung **squamous cell carcinoma** section, and the pipeline now refuses any section "
    "in which the squamous lineage markers KRT5, TP63 and SOX2 are not detected in at least 30 % of spots."
)
if has_table("visium_identity"):
    vi = q("SELECT * FROM visium_identity").sort_values("spot_fraction")
    fig = go.Figure(
        go.Bar(
            y=vi.gene,
            x=100 * vi.spot_fraction,
            orientation="h",
            marker_color=["#1c5cab" if r else "#b9b8b2" for r in vi.required],
            text=[f"{v:.1f} per 10k" for v in vi.mean_per_10k],
            textposition="outside",
        )
    )
    fig.add_vline(x=30, line_dash="dash", line_color=INK_2, annotation_text="gate: 30 %")
    st.plotly_chart(style(fig, 420, x="% of spots where the gene is detected (≥ 1 UMI)", y="marker"))
    figure_legend(
        1,
        "Squamous identity confirmed.",
        "share of spots detecting each marker; dark = required squamous markers, grey = context markers "
        "(adenocarcinoma: NKX2-1, NAPSA; alveolar: SFTPC; neuroendocrine: CHGA, SYP, INSM1; immune: PTPRC; "
        "fibroblast: COL1A1); labels = mean counts per 10,000.",
        "squamous markers widespread and neuroendocrine markers near zero = the right tumour type.",
        "all spots of the section.",
        "computed in `spatial_prep` before any analysis (the step fails if the gate is not met).",
        "the spatial results below describe a LUSC tumour.",
    )

methods_box(
    "cell2location",
    "Learns an expression signature per cell type from the single-cell atlas (with sample and study as "
    "covariates), then models each spot's counts as a weighted mix of signatures; reports each spot's posterior "
    "abundance per type. Proportions are computed in every posterior sample and averaged.",
    "A 55 µm Visium spot holds several cells, so one label per spot is wrong; cell2location performed best in "
    "deconvolution benchmarks and models technical differences between the atlas and the slide.",
    "RCTD, Tangram, Stereoscope, CARD; label transfer (v1 approach, one label per spot).",
    "Kleshchevnikov et al. 2022 Nat Biotechnol; Li et al. 2023 (benchmark).",
)
decisions("spatial")

img = q("SELECT png, scalef, spot_diameter_fullres FROM visium_images WHERE resolution = 'hires'").iloc[0]
spots = q("SELECT * FROM visium_spots")
types = [c for c in q("SELECT * FROM c2l_abundance LIMIT 0").columns if c != "spot"]
metric_row(
    [
        ("Spots", f"{len(spots):,}", None),
        ("Median UMIs / spot", f"{spots.total_counts.median():,.0f}", None),
        (
            "Cell types deconvolved",
            f"{len(types)}",
            "cell2location; CNV-high and CNV-low tumour cells share one 'Tumour epithelium' signature",
        ),
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
    what = st.radio("Colour spots by", ["niche", "cell-type proportion", "gene"])
    if what == "cell-type proportion":
        dflt = "Tumour epithelium" if "Tumour epithelium" in types else types[0]
        ct = st.selectbox("Cell type", types, index=types.index(dflt))
        val, title = spots[ct] if ct in spots else 0, f"{ct} (share of spot)"
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
figure_legend(
    2,
    "The tissue section with deconvolved spots.",
    "the H&E image with each Visium spot coloured by niche, by a cell type's estimated share of the spot, or by "
    "a gene's log-normalised expression.",
    "compare tumour-epithelium-rich regions with the H&E morphology; niches are groups of spots with similar "
    "neighbourhoods.",
    f"{len(spots):,} spots after QC.",
    "cell2location posterior proportions; CellCharter niches; gene list = top spatially variable genes, "
    "ligand–receptor genes and markers.",
    "descriptive only (one section).",
)

st.subheader("What is in each niche")
ni = q("SELECT * FROM spatial_niches").set_index("niche")
names = {n: " + ".join(row.sort_values(ascending=False).head(2).index) for n, row in ni.iterrows()}
ni.index = [f"niche {i}: {names[i]}" for i in ni.index]
fig = go.Figure()
top = ni.max().sort_values(ascending=False).index[:7]
others = [c for c in top if c != "Tumour epithelium"]
for (
    ct
) in top:  # tumour epithelium keeps the site-wide malignant red; the rest take generic hues in fixed order
    color = MALIGNANT_COLOR if ct == "Tumour epithelium" else series_color(others.index(ct))
    fig.add_trace(go.Bar(y=ni.index, x=ni[ct], name=ct, orientation="h", marker_color=color))
fig.add_trace(
    go.Bar(y=ni.index, x=ni.drop(columns=top).sum(axis=1), name="other", orientation="h", marker_color=GREY)
)
fig.update_layout(barmode="stack", xaxis=dict(tickformat=".0%"))
st.plotly_chart(
    style(fig, 120 + 34 * len(ni), x="mean share of the spot", y="niche (named by its two largest types)")
)
figure_legend(
    3,
    "Most niches are dominated by tumour epithelium; the others are plasma/fibroblast, alveolar, airway and "
    "smooth-muscle regions.",
    "mean cell-type shares of the spots in each niche; the seven largest types coloured, the rest grey.",
    "niches are named by their two largest cell types; similar bars mean CellCharter's automatic choice of the "
    "number of niches split one tissue region into several.",
    f"{len(ni)} niches.",
    "CellCharter (neighbour-aggregated proportions, 3 layers; number of niches chosen by stability).",
    "where tumour epithelium meets immune and stromal cells.",
)
if {"Tumour epithelium", "Basal"} <= set(spots.columns):
    st.caption(
        "Squamous tumour spots can split their signal between the tumour-epithelium signature and the healthy "
        f"basal signature (both squamous-lineage): their summed share averages "
        f"{(spots['Tumour epithelium'] + spots['Basal']).mean():.0%} per spot."
    )
show_table(ni.reset_index(names="niche"))
if has_table("c2l_reference_labels"):
    rl = q("SELECT label, n_cells, action, target FROM c2l_reference_labels ORDER BY n_cells DESC")
    off = rl[rl.action != "kept"]
    st.caption(
        f"cell2location reference: {int((rl.action == 'kept').sum())} cell types with a signature. Not used: "
        + ", ".join(f"*{r.label}* ({int(r.n_cells):,} cells, {r.action})" for r in off.itertuples())
        + ". Tiny reference types give noisy signatures that soak up unexplained signal, and a mixed bucket has no "
        "coherent signature (DEVELOPMENT_LOG §3.11)."
    )
    show_table(rl, "Reference labels: kept, merged, dropped or excluded")

c1, c2 = st.columns(2)
with c1:
    st.subheader("Co-localisation of cell types")
    co = q("SELECT * FROM spatial_colocalisation").set_index("cell_type")
    st.plotly_chart(heatmap(co, "Pearson r", height=640))
    figure_legend(
        4,
        "Which cell types share spots.",
        "Pearson correlation of two cell types' shares across all spots.",
        "red = found together, blue = mutually exclusive.",
        f"{len(spots):,} spots.",
        "correlation of cell2location proportions (no test; spots are not independent).",
        "descriptive tissue architecture.",
    )
with c2:
    st.subheader("Top spatially variable genes")
    svg = q('SELECT gene, I AS "Moran I", pval_norm_fdr_bh AS fdr FROM spatial_svg ORDER BY I DESC LIMIT 25')
    fig = go.Figure(go.Bar(y=svg.gene, x=svg["Moran I"], orientation="h", marker_color="#8a8984"))
    fig.update_layout(yaxis=dict(autorange="reversed"))
    st.plotly_chart(style(fig, 640, x="Moran's I (spatial autocorrelation)"))
    figure_legend(
        5,
        "Most spatially patterned genes.",
        "the 25 genes with the highest Moran's I.",
        "I near 1 = expression forms spatial patches.",
        f"{len(spots):,} spots.",
        "squidpy Moran's I with 100 permutations (smallest attainable p ≈ 0.01); all shown have BH FDR < 0.05.",
        "genes ranked by I (nearly all genes are 'significant' on Visium).",
    )
footer()
