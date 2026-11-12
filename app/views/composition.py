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
caveats("de")  # issue 2 (Unknown as scCODA reference) lives under DE but changes this page too

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
st.caption(
    "Coloured = credible (inclusion probability ≥ 0.8); grey = not credible. The reference type was chosen "
    "automatically (`Unknown`), so effects are relative to that bucket."
)
show_table(d.drop(columns=["cond"]))

st.subheader("Raw proportions per sample (descriptive)")
ct = st.selectbox(
    "Cell type", sorted(sc["Cell Type"].unique()), index=sorted(sc["Cell Type"].unique()).index("CD8 T")
)
prop = q(
    """SELECT sample, condition, avg(CASE WHEN harmonised_label = ? THEN 1 ELSE 0 END) frac, count(*) n
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

st.subheader("Milo neighbourhoods")
mi = q("SELECT logFC, SpatialFDR, nhood_annotation FROM milo_nhoods")
st.markdown(
    f"{len(mi):,} neighbourhoods, **{(mi.SpatialFDR < 0.1).sum()}** with SpatialFDR < 0.1. "
    "This result is not interpretable until known issue 1 (wrong contrast) is fixed, so no plot is shown."
)
footer()
