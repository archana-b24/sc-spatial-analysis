import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.ui import (
    CONDITION_COLOR,
    SERIES,
    caveats,
    decisions,
    footer,
    metric_row,
    page_header,
    show_table,
    style,
)

page_header(
    "4 · Malignant cells",
    "Which epithelial cells are tumour cells? Tumour cells carry copy-number changes (gains and losses of "
    "chromosome segments) that can be inferred from expression; normal cells do not.",
)
guide("malignant")
decisions("malignant")
caveats("malignant")
cs = q("SELECT * FROM cnv_summary")
thr = float(cs.threshold.iloc[0])
va = q("SELECT * FROM cnv_vs_author_malignant")
mal_cnv_author = int(va.loc[va.author_malignant.astype(str) == "True", "True"].iloc[0])
author_total = int(va.loc[va.author_malignant.astype(str) == "True", ["False", "True"]].sum(axis=1).iloc[0])
called_total = int(va["True"].sum())
metric_row(
    [
        ("Malignant cells", f"{int(cs.loc[cs.cnv_group == 'malignant', 'n_cells'].iloc[0]):,}", None),
        ("Threshold", f"{thr:.4f}", "95th percentile of CNV score in healthy-donor epithelium"),
        (
            "Precision vs author labels",
            f"{mal_cnv_author / called_total:.1%}",
            "CNV calls that authors also labelled malignant",
        ),
        (
            "Sensitivity vs author labels",
            f"{mal_cnv_author / author_total:.0%}",
            "author-malignant cells that pass the threshold",
        ),
    ]
)
st.info(
    "**Reading:** the calls are trustworthy (few false positives) but very conservative. Many tumour cells fall below "
    "the threshold; they show up as *CNV-low tumour-sample epithelium* below."
)

c1, c2 = st.columns(2)
order = [
    "reference (immune/endothelial/stroma)",
    "healthy epithelium",
    "adjacent-normal epithelium",
    "tumour-sample epithelium, CNV-low",
    "malignant",
]
cs = cs.set_index("cnv_group").loc[[o for o in order if o in set(cs.cnv_group)]].reset_index()
with c1:
    st.subheader("Positive control: chr3q gain")
    fig = go.Figure(
        go.Bar(
            y=cs.cnv_group,
            x=cs.mean_chr3q,
            orientation="h",
            marker_color=SERIES[0],
            text=cs.mean_chr3q.round(3),
            textposition="outside",
        )
    )
    fig.update_xaxes(title="mean chr3q CNV score")
    st.plotly_chart(style(fig, 340))
    st.caption(
        "chr3q gain (SOX2, TP63, PIK3CA) is the hallmark LUSC event: ~9× higher in called malignant cells than in healthy epithelium."
    )
with c2:
    st.subheader("Overall CNV score")
    fig = go.Figure(go.Bar(y=cs.cnv_group, x=cs.mean_cnv_score, orientation="h", marker_color=SERIES[0]))
    fig.add_vline(x=thr, line_dash="dash", line_color=SERIES[7], annotation_text="threshold")
    fig.update_xaxes(title="mean CNV score")
    st.plotly_chart(style(fig, 340))

st.subheader("CNV score distribution in epithelial cells")
ep = q("""SELECT condition, cnv_score FROM cells WHERE in_subsample AND compartment = 'Epithelial'""")
fig = go.Figure()
for cond, color in CONDITION_COLOR.items():
    d = ep[ep.condition == cond]
    if len(d):
        fig.add_trace(
            go.Violin(
                x=d.cnv_score,
                name=cond,
                line_color=color,
                orientation="h",
                side="positive",
                points=False,
                meanline_visible=True,
            )
        )
fig.add_vline(x=thr, line_dash="dash", line_color=SERIES[7], annotation_text="threshold")
fig.update_xaxes(title="CNV score")
st.plotly_chart(style(fig, 360))

st.subheader("Malignant cells per tumour donor")
pd_ = q("""SELECT donor_id, condition, "sum" AS malignant, size AS cells FROM cnv_malignant_per_donor
           WHERE condition != 'healthy' ORDER BY malignant DESC""")
tum = pd_[pd_.condition == "tumour"]
st.markdown(
    f"**{(tum.malignant >= 1).sum()} of {len(tum)}** tumour donors have ≥ 1 malignant call and **{(tum.malignant >= 10).sum()}** have ≥ 10, "
    "so malignant-cell analyses rest on few patients."
)
fig = go.Figure(
    go.Bar(
        x=tum.donor_id.str.split(":").str[-1].str[-12:],
        y=tum.malignant,
        marker_color=SERIES[0],
        hovertext=tum.donor_id,
        hovertemplate="%{hovertext}<br>%{y} malignant cells<extra></extra>",
    )
)
fig.update_yaxes(type="log", title="malignant cells (log)")
st.plotly_chart(style(fig, 320))
show_table(pd_)
st.subheader("Agreement with author labels")
st.dataframe(va.rename(columns={"False": "CNV: not malignant", "True": "CNV: malignant"}), hide_index=True)
footer()
