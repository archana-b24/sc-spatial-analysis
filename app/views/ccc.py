import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.plots import heatmap
from lib.ui import (
    CONDITION_COLOR,
    CONDITION_LABEL,
    SERIES,
    caveats,
    decisions,
    footer,
    page_header,
    show_table,
    style,
)

page_header(
    "9 · Cell–cell communication",
    "Which ligand–receptor signals between cell types change in tumours, and which communication programmes "
    "vary across samples?",
)
guide("ccc")
decisions("ccc")
caveats("ccc")
tab1, tab2, tab3, tab4 = st.tabs(
    ["Differential (tumour vs adjacent)", "Tensor programmes", "Spatial co-expression", "Per-sample browser"]
)

with tab1:
    st.markdown(
        "Interactions scored from the **paired** pseudobulk DE: a high score means ligand (in the sender) and "
        "receptor (in the receiver) both go up in tumour."
    )
    n = st.slider("Show top", 10, 50, 20, 5)
    d = q(
        """SELECT source || ' → ' || target AS pair, ligand_complex || ' → ' || receptor_complex AS lr,
                    interaction_stat, interaction_padj, ligand_stat, receptor_stat, source, target
             FROM ccc_differential ORDER BY interaction_stat DESC LIMIT ?""",
        (n,),
    )
    d["label"] = d.lr + "  (" + d.pair + ")"
    fig = go.Figure(
        go.Bar(
            y=d.label,
            x=d.interaction_stat,
            orientation="h",
            marker_color="#c23b3a",
            customdata=d.interaction_padj,
            hovertemplate="%{y}<br>stat %{x:.2f}<br>padj %{customdata:.1e}<extra></extra>",
        )
    )
    fig.update_layout(
        yaxis=dict(autorange="reversed"), xaxis=dict(title="interaction statistic (up in tumour)")
    )
    st.plotly_chart(style(fig, 28 * n + 80))
    show_table(q("SELECT * FROM ccc_differential ORDER BY interaction_stat DESC LIMIT 500"), "top 500 rows")

with tab2:
    st.markdown(
        "Tensor-cell2cell decomposes samples × LR pairs × senders × receivers into **5 factors** (elbow rank). "
        "Each factor is a communication programme with loadings on every dimension."
    )
    k = st.selectbox("Factor", [f"Factor {i}" for i in range(1, 6)])
    c1, c2 = st.columns(2)
    for col, tbl, title in [
        (c1, "ccc_tensor_sender_cells", "Senders"),
        (c2, "ccc_tensor_receiver_cells", "Receivers"),
    ]:
        t = q(f'SELECT element, "{k}" AS v FROM {tbl} ORDER BY v DESC')
        fig = go.Figure(go.Bar(y=t.element, x=t.v, orientation="h", marker_color=SERIES[0]))
        fig.update_layout(yaxis=dict(autorange="reversed"), title=title)
        col.plotly_chart(style(fig, 360))
    lr = q(f'SELECT element, "{k}" AS v FROM ccc_tensor_ligand_receptor_pairs ORDER BY v DESC LIMIT 15')
    fig = go.Figure(go.Bar(y=lr.element, x=lr.v, orientation="h", marker_color=SERIES[0]))
    fig.update_layout(yaxis=dict(autorange="reversed"), title="Top ligand–receptor pairs")
    st.plotly_chart(style(fig, 440))
    ctx = q(f"""SELECT t.element AS sample, t."{k}" AS v, any_value(c.condition) AS condition
                FROM ccc_tensor_contexts t JOIN cells c ON c.sample = t.element GROUP BY 1, 2""")
    fig = go.Figure()
    for cnd, color in CONDITION_COLOR.items():
        s = ctx[ctx.condition == cnd]
        if len(s):
            fig.add_trace(
                go.Box(
                    y=s.v,
                    name=CONDITION_LABEL[cnd],
                    marker_color=color,
                    boxpoints="all",
                    pointpos=0,
                    jitter=0.4,
                )
            )
    fig.update_layout(title="Sample loadings by condition")
    st.plotly_chart(style(fig, 360))
    st.caption(
        "No factor separates tumour from adjacent (all p ≥ 0.31, known issue 4): read them as sample-level programmes."
    )
    st.subheader("All factors at a glance (senders)")
    st.plotly_chart(
        heatmap(
            q("SELECT * FROM ccc_tensor_sender_cells").set_index("element"),
            "loading",
            diverging=False,
            height=380,
        )
    )

with tab3:
    st.markdown(
        "Ligand–receptor co-expression across neighbouring Visium spots (LIANA+ bivariate Moran's R)."
    )
    b = q("SELECT interaction, morans, mean, std FROM ccc_spatial_bivariate ORDER BY morans DESC LIMIT 25")
    fig = go.Figure(go.Bar(y=b.interaction, x=b.morans, orientation="h", marker_color=SERIES[0]))
    fig.update_layout(yaxis=dict(autorange="reversed"), xaxis=dict(title="global bivariate Moran's R"))
    st.plotly_chart(style(fig, 640))

with tab4:
    st.markdown("Per-sample LIANA+ results (shipped rows: consensus magnitude rank ≤ 0.01).")
    labels = q("SELECT DISTINCT source FROM ccc_by_sample ORDER BY 1").source.tolist()
    c1, c2, c3 = st.columns(3)
    src = c1.selectbox("Sender", labels, index=labels.index("Macrophage") if "Macrophage" in labels else 0)
    tgt = c2.selectbox("Receiver", labels, index=labels.index("CD8 T") if "CD8 T" in labels else 0)
    cond = c3.multiselect(
        "Conditions",
        list(CONDITION_COLOR),
        default=["tumour", "adjacent_normal"],
        format_func=CONDITION_LABEL.get,
    )
    r = q(
        f"""WITH sc AS (SELECT sample, any_value(condition) condition FROM cells GROUP BY 1)
              SELECT b.ligand_complex || ' → ' || b.receptor_complex AS lr, sc.condition, count(*) n_samples,
                     median(b.magnitude_rank) median_rank
              FROM ccc_by_sample b JOIN sc USING (sample)
              WHERE b.source = ? AND b.target = ? AND sc.condition IN ({",".join("?" * len(cond)) or "''"})
              GROUP BY 1, 2 ORDER BY n_samples DESC LIMIT 60""",
        (src, tgt, *cond),
    )
    st.dataframe(r, hide_index=True)
    st.caption("`n_samples` = samples in which this interaction is in the top 1 % by magnitude.")
footer()
