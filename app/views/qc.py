import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.ui import ATLAS_COLOR, SERIES, caveats, decisions, footer, metric_row, page_header, show_table, style

page_header(
    "2 · Quality control",
    "Which cells are damaged, empty or doublets? Thresholds are set per sample, because samples sequenced "
    "at different depths cannot share one cutoff.",
)
guide("qc")
decisions("qc")
caveats("qc")
qc = q("""SELECT qc.*, split_part(sample, ':', 1) AS atlas FROM qc_per_sample qc""")
tot_in, tot_kept = qc.n_cells.sum(), qc.n_kept.sum()
metric_row(
    [
        ("Cells in", f"{tot_in:,}", None),
        ("Cells kept", f"{tot_kept:,}", f"{tot_kept / tot_in:.1%}"),
        ("MAD outliers", f"{qc.n_outlier.sum():,}", "median-absolute-deviation outliers, per sample"),
        ("Consensus doublets", f"{qc.n_doublet.sum():,}", "flagged by both Scrublet and SOLO"),
    ]
)

c1, c2 = st.columns(2)
with c1:
    st.subheader("Sequencing depth per sample")
    fig = go.Figure()
    for atlas, d in qc.groupby("atlas"):
        fig.add_trace(
            go.Scatter(
                x=d.median_counts,
                y=d.median_genes,
                mode="markers",
                name=atlas,
                marker=dict(color=ATLAS_COLOR.get(atlas), size=8, line=dict(width=1, color="white")),
                text=d["sample"],
                hovertemplate="%{text}<br>median UMIs %{x:,.0f}<br>median genes %{y:,.0f}<extra></extra>",
            )
        )
    fig.update_xaxes(type="log", title="median UMIs per cell")
    fig.update_yaxes(type="log", title="median genes per cell")
    st.plotly_chart(style(fig, 400))
    st.caption(
        "A ten-fold spread in depth is why a single global threshold (v1) was replaced by per-sample MADs."
    )
with c2:
    st.subheader("Fraction of cells kept per sample")
    fig = go.Figure(go.Histogram(x=qc.n_kept / qc.n_cells, nbinsx=40, marker_color=SERIES[0]))
    fig.update_xaxes(title="fraction kept", tickformat=".0%")
    fig.update_yaxes(title="samples")
    st.plotly_chart(style(fig, 400))

st.subheader("Doublet callers disagree")
dd = {
    "Scrublet": qc.n_doublet_scrublet.sum(),
    "scVI-SOLO": qc.n_doublet_solo.sum(),
    "Both (removed)": qc.n_doublet.sum(),
}
fig = go.Figure(
    go.Bar(
        x=list(dd.values()),
        y=list(dd),
        orientation="h",
        marker_color=[SERIES[0], SERIES[1], SERIES[2]],
        text=[f"{v:,}" for v in dd.values()],
        textposition="outside",
    )
)
st.plotly_chart(style(fig, 240))
st.caption(
    "Cells are removed only when both callers agree, so removal is effectively Scrublet-only (known issue 5)."
)
show_table(qc.drop(columns="atlas"))
footer()
