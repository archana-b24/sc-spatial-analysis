import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.plots import heatmap
from lib.ui import GREY, SERIES, caveats, decisions, footer, page_header, show_table, style

page_header(
    "5 · Differential expression & pathways",
    "Within each cell type, which genes and pathways change in LUSC tumour tissue?",
)
guide("de")
decisions("de")
caveats("de")

CONTRASTS = {
    "paired_tumour_vs_adjacent": "Tumour vs adjacent normal, same patients (primary, ~ donor + condition)",
    "tumour_vs_healthy": "Tumour vs healthy donors (sensitivity; confounded with study by design)",
    "malignant_vs_healthy_basal": "Malignant cells vs healthy basal cells",
}
summ = q("""SELECT contrast, "group", count(*) FILTER (WHERE padj < 0.05) n_sig, max(n_donors) n_donors,
                   max(n_case) n_case, max(n_control) n_control
            FROM de_pseudobulk GROUP BY 1, 2 ORDER BY 1, 3 DESC""")
contrast = st.selectbox("Contrast", list(CONTRASTS), format_func=CONTRASTS.get)
cs = summ[summ.contrast == contrast]

st.subheader("Significant genes per cell type (padj < 0.05)")
fig = go.Figure(
    go.Bar(
        y=cs.group,
        x=cs.n_sig,
        orientation="h",
        marker_color=SERIES[0],
        customdata=np.stack([cs.n_case, cs.n_control], axis=1),
        hovertemplate="%{y}: %{x:,} genes<br>%{customdata[0]} case vs %{customdata[1]} control samples<extra></extra>",
    )
)
fig.update_layout(yaxis=dict(autorange="reversed"))
st.plotly_chart(style(fig, max(220, 28 * len(cs))))
if contrast == "paired_tumour_vs_adjacent":
    st.caption(
        "No epithelial group passes the minimum of 10 cells in ≥ 3 samples per condition; adjacent tissue holds "
        "few epithelial cells. The epithelial question is answered by *malignant vs healthy basal*."
    )

st.subheader("Volcano plot")
c1, c2 = st.columns([1, 3])
with c1:
    group = st.selectbox("Cell type", cs.group.tolist())
    fdr = st.select_slider("padj cut", [0.001, 0.01, 0.05, 0.1], value=0.05)
    lfc = st.slider("|log2 fold change| cut", 0.0, 3.0, 1.0, 0.25)
d = q(
    """SELECT gene, log2FoldChange, padj, baseMean, stat FROM de_pseudobulk
         WHERE contrast = ? AND "group" = ? AND padj IS NOT NULL""",
    (contrast, group),
)
d["nlp"] = -np.log10(d.padj.clip(lower=1e-300))
up = (d.padj < fdr) & (d.log2FoldChange >= lfc)
dn = (d.padj < fdr) & (d.log2FoldChange <= -lfc)
fig = go.Figure()
for mask, name, color in [
    (~(up | dn), "not significant", GREY),
    (dn, f"down ({dn.sum()})", SERIES[0]),
    (up, f"up ({up.sum()})", "#c23b3a"),
]:
    s = d[mask]
    fig.add_trace(
        go.Scattergl(
            x=s.log2FoldChange,
            y=s.nlp,
            mode="markers",
            name=name,
            text=s.gene,
            marker=dict(color=color, size=5 if color == GREY else 7, opacity=0.6 if color == GREY else 0.9),
            hovertemplate="%{text}<br>log2FC %{x:.2f}<br>-log10 padj %{y:.1f}<extra></extra>",
        )
    )
top = d[up | dn].nlargest(12, "nlp")
fig.add_trace(
    go.Scatter(
        x=top.log2FoldChange,
        y=top.nlp,
        mode="text",
        text=top.gene,
        textposition="top center",
        showlegend=False,
        hoverinfo="skip",
        textfont=dict(size=11),
    )
)
fig.update_xaxes(title="log2 fold change (tumour / control)")
fig.update_yaxes(title="-log10 adjusted p")
with c2:
    st.plotly_chart(style(fig, 480))
show_table(d[up | dn].sort_values("padj").drop(columns="nlp"), f"{group}: significant genes")

st.subheader("Pathway activity (PROGENy)")
pg = q('SELECT "group", source, score, padj FROM pathways_progeny WHERE contrast = ?', (contrast,))
if not pg.empty:
    mat = pg.pivot(index="group", columns="source", values="score")
    st.plotly_chart(heatmap(mat, "score", height=max(260, 26 * len(mat) + 120)))
    st.caption(
        "Positive = pathway more active in the case group. Scores come from the DESeq2 Wald statistic (decoupler)."
    )

c1, c2 = st.columns(2)
for col, tbl, title in [
    (c1, "pathways_hallmark", "Hallmark gene sets"),
    (c2, "pathways_collectri", "Transcription factors (CollecTRI)"),
]:
    with col:
        st.subheader(title)
        t = q(
            f"""SELECT source, score, padj FROM {tbl} WHERE contrast = ? AND "group" = ? AND padj < 0.05
                  ORDER BY abs(score) DESC LIMIT 20""",
            (contrast, group),
        )
        if t.empty:
            st.caption(f"No significant entries for {group} in this contrast.")
            continue
        fig = go.Figure(
            go.Bar(
                y=t.source,
                x=t.score,
                orientation="h",
                marker_color=[("#c23b3a" if v > 0 else SERIES[0]) for v in t.score],
            )
        )
        fig.update_layout(yaxis=dict(autorange="reversed"))
        st.plotly_chart(style(fig, 520))

st.subheader("Look up a gene across cell types")
gene = st.text_input("Gene symbol", "SPP1").strip().upper()
g = q(
    """SELECT contrast, "group", log2FoldChange, padj FROM de_pseudobulk WHERE gene = ? ORDER BY contrast, "group" """,
    (gene,),
)
if g.empty:
    st.caption("Gene not tested (not expressed enough in any pseudobulk group).")
else:
    fig = go.Figure()
    for i, (c, s) in enumerate(g.groupby("contrast")):
        fig.add_trace(
            go.Scatter(
                x=s.log2FoldChange,
                y=s.group,
                mode="markers",
                name=c,
                marker=dict(
                    color=SERIES[i], size=np.where(s.padj < 0.05, 13, 7), line=dict(width=1, color="white")
                ),
                hovertemplate="%{y}<br>log2FC %{x:.2f}<extra>" + c + "</extra>",
            )
        )
    fig.add_vline(x=0, line_color=GREY)
    fig.update_xaxes(title=f"{gene} log2 fold change (large marker = padj < 0.05)")
    st.plotly_chart(style(fig, max(300, 22 * g.group.nunique() + 120)))
footer()
