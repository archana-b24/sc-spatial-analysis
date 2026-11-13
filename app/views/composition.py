import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.ui import (
    CONDITION_COLOR,
    CONDITION_LABEL,
    GREY,
    caveats,
    decisions,
    footer,
    page_header,
    show_table,
    style,
)

page_header(
    "6 · Composition",
    "Do tumours contain different proportions of cell types than normal lung? Proportions must sum to 1, "
    "so a rise in one type forces others down; the model has to account for that.",
)
guide("composition")
decisions("composition")
caveats("composition")
caveats("de")  # issue 2 (Unknown / scCODA reference) lives under DE but changes this page too

sc = q("SELECT * FROM composition_sccoda")
sc["cond"] = sc.Covariate.str.extract(r"T\.(.+)\]")[0]
sc["credible"] = sc["Inclusion probability"] >= 0.8
cond = st.radio(
    "Compared with adjacent-normal tissue:",
    [c for c in ["tumour", "healthy", "tumour_region_unknown"] if c in set(sc.cond)],
    format_func=CONDITION_LABEL.get,
    horizontal=True,
)
d = sc[sc.cond == cond].sort_values("log2-fold change")
st.subheader("scCODA: credible changes in proportion")
fig = go.Figure(
    go.Bar(
        y=d["Cell Type"],
        x=d["log2-fold change"],
        orientation="h",
        marker_color=[CONDITION_COLOR[cond] if c else GREY for c in d.credible],
        customdata=d["Inclusion probability"],
        hovertemplate="%{y}<br>log2FC %{x:.2f}<br>inclusion prob %{customdata:.2f}<extra></extra>",
    )
)
fig.add_vline(x=0, line_color=GREY)
fig.update_xaxes(title=f"log2 fold change, {CONDITION_LABEL[cond]} vs adjacent normal")
st.plotly_chart(style(fig, 640))
ref = str(sc["reference"].iloc[0]) if "reference" in sc.columns else "automatic"
st.caption(
    "Coloured = credible (inclusion probability ≥ 0.8); grey = not credible. scCODA estimates every change "
    f"relative to one reference cell type assumed not to change: here **{ref}** (set in config; chosen because it "
    "is present in most samples of every condition with a near-constant mean proportion, DEVELOPMENT_LOG §3.11). "
    "A type with log2FC 0 is the reference itself."
)
show_table(d.drop(columns=["cond"]))

st.subheader("Raw proportions per sample (descriptive)")
types = sorted(sc["Cell Type"].unique())
ct = st.selectbox("Cell type", types, index=types.index("CD8 T") if "CD8 T" in types else 0)
prop = q(
    """SELECT sample, condition, avg(CASE WHEN de_group = ? THEN 1 ELSE 0 END) frac, count(*) n
            FROM cells GROUP BY 1, 2 HAVING count(*) >= 50""",
    (ct,),
)
fig = go.Figure()
for c, color in CONDITION_COLOR.items():
    s = prop[prop.condition == c]
    if len(s):
        fig.add_trace(
            go.Box(
                y=s.frac,
                name=CONDITION_LABEL[c],
                marker_color=color,
                boxpoints="all",
                jitter=0.4,
                pointpos=0,
                text=s["sample"],
                hovertemplate="%{text}<br>%{y:.1%}<extra></extra>",
            )
        )
fig.update_yaxes(title=f"fraction {ct}", tickformat=".0%")
st.plotly_chart(style(fig, 380))
st.caption("Samples with ≥ 50 cells. Raw fractions are for orientation only; inference is scCODA above.")

st.subheader("Milo: neighbourhood-level changes, tumour vs adjacent normal")
mi = q("SELECT logFC, SpatialFDR, nhood_annotation, design, contrast, n_samples, n_donors FROM milo_nhoods")
mi = mi.dropna(subset=["logFC"])
mi["sig"] = mi.SpatialFDR < 0.1
r0 = mi.iloc[0]
st.markdown(
    f"Milo groups cells into overlapping neighbourhoods on the integrated graph and tests each for a change in "
    f"abundance. Design `{r0.design}`, contrast `{r0.contrast}`: only tumour and adjacent-normal samples from the "
    f"**{int(r0.n_donors)}** patients who gave both (**{int(r0.n_samples)}** samples), so every patient is "
    f"compared with themself. **{int(mi.sig.sum()):,}** of {len(mi):,} neighbourhoods change at SpatialFDR < 0.1 "
    f"({int((mi.sig & (mi.logFC > 0)).sum()):,} more abundant in tumour, "
    f"{int((mi.sig & (mi.logFC < 0)).sum()):,} less)."
)
order = mi.groupby("nhood_annotation")["logFC"].median().sort_values().index.tolist()
fig = go.Figure()
for sig, color, name in [
    (False, GREY, "not significant"),
    (True, CONDITION_COLOR["tumour"], "SpatialFDR < 0.1"),
]:
    d = mi[mi.sig == sig]
    fig.add_trace(
        go.Box(
            y=d.nhood_annotation,
            x=d.logFC,
            orientation="h",
            name=name,
            marker_color=color,
            boxpoints="all",
            jitter=0.6,
            pointpos=0,
            line_width=0,
            fillcolor="rgba(0,0,0,0)",
            marker=dict(size=3, opacity=0.5),
            hovertemplate="%{y}<br>log2FC %{x:.2f}<extra></extra>",
        )
    )
fig.add_vline(x=0, line_color=GREY)
fig.update_yaxes(categoryorder="array", categoryarray=order)
fig.update_xaxes(title="log2 fold change in neighbourhood abundance, tumour vs adjacent normal")
fig.update_layout(boxmode="overlay")
st.plotly_chart(style(fig, 720))
st.caption(
    "One point per neighbourhood, labelled by its most common cell type. Right of zero = more abundant in tumour. "
    "Neighbourhoods are tested on the integrated embedding; cell types with no adjacent-normal counterpart (e.g. "
    "malignant cells) can only appear on the tumour side."
)
show_table(
    mi.assign(direction=mi.logFC.gt(0).map({True: "up in tumour", False: "down in tumour"}))
    .loc[mi.sig]
    .groupby(["nhood_annotation", "direction"])
    .size()
    .unstack(fill_value=0)
    .reset_index(),
    "Significant neighbourhoods per cell type",
)
footer()
