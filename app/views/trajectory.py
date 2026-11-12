import plotly.graph_objects as go
import streamlit as st
from lib.db import meta, q
from lib.explain import guide
from lib.plots import umap_categorical, umap_continuous
from lib.ui import SERIES, caveats, decisions, footer, page_header, show_table, style

page_header(
    "7 · Trajectory",
    "How far, in expression space, are malignant cells from the healthy basal cells that LUSC is thought to "
    "arise from?",
)
guide("trajectory")
decisions("trajectory")
caveats("trajectory")
st.info(
    "**This is not lineage tracing.** Tumour cells are clonal within each patient, so pseudotime measures "
    "transcriptional distance from healthy basal cells, not a developmental path."
)

lc = q("SELECT * FROM lineage_cells")
lc["harmonised_label"], lc["study"] = lc.lineage_state, lc.atlas
c1, c2 = st.columns(2)
with c1:
    st.subheader("Cell state")
    st.plotly_chart(
        umap_categorical(lc, "lineage_state", colors={"Basal": SERIES[0], "Malignant": SERIES[1]}, height=440)
    )
with c2:
    st.subheader("Palantir pseudotime")
    st.plotly_chart(umap_continuous(lc, lc.palantir_pseudotime, "pseudotime", height=440))

st.subheader("Pseudotime by state")
fig = go.Figure()
for i, (s, d) in enumerate(lc.groupby("lineage_state")):
    fig.add_trace(
        go.Violin(
            x=d.palantir_pseudotime,
            name=f"{s} ({len(d):,})",
            line_color=SERIES[i],
            orientation="h",
            side="positive",
            meanline_visible=True,
            points=False,
        )
    )
fig.update_xaxes(title="pseudotime (0 = healthy basal root)", range=[0, 1])
st.plotly_chart(style(fig, 320))
show_table(q("SELECT * FROM trajectory_summary"), "pseudotime summary per state")
info = meta().get("info_trajectory_info")
if info:
    st.caption(f"Root cell and CellRank terminal states: `{info}`")
footer()
