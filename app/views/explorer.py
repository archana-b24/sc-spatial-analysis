import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.plots import COMPARTMENT_COLOR, umap_categorical, umap_continuous  # noqa: F401
from lib.ui import ATLAS_COLOR, CONDITION_COLOR, CONDITION_LABEL, footer, style

st.title("Cell & gene explorer")
st.caption(
    "A free-exploration view: no statistics are computed here. Colouring by atlas largely also colours by tissue, "
    "because healthy donors all come from HLCA — differences between those colours are not tumour effects. For "
    "tested comparisons use pages 5 (genes), 6 (composition) and 9 (communication)."
)
st.markdown(
    "Colour the integrated atlas by metadata or by a marker gene, and compare a gene list across cell "
    "types and conditions."
)

cells = q("""SELECT cell_idx, umap_1, umap_2, atlas, study, condition, harmonised_label, compartment,
                    cnv_score, total_counts, pct_counts_mt FROM cells WHERE in_subsample""")
panel = q("SELECT gene FROM genes WHERE in_panel ORDER BY gene").gene.tolist()
c1, c2 = st.columns([1, 3])
with c1:
    mode = st.radio(
        "Colour by", ["gene", "cell type (highlight)", "condition", "atlas", "compartment", "QC metric"]
    )
    cond_filter = st.multiselect(
        "Show conditions",
        list(CONDITION_COLOR),
        default=list(CONDITION_COLOR),
        format_func=CONDITION_LABEL.get,
    )
    view = cells[cells.condition.isin(cond_filter)]
    if mode == "gene":
        gene = st.selectbox(
            f"Gene ({len(panel):,} in panel)", panel, index=panel.index("KRT5") if "KRT5" in panel else 0
        )
        e = q("SELECT cell_idx, value FROM cell_expr WHERE gene = ?", (gene,)).set_index("cell_idx").value
        fig = umap_continuous(view, view.cell_idx.map(e).fillna(0.0), gene)
    elif mode == "cell type (highlight)":
        labs = sorted(cells.harmonised_label.unique())
        fig = umap_categorical(
            view, "harmonised_label", highlight=st.selectbox("Cell type", labs, index=labs.index("Basal"))
        )
    elif mode == "QC metric":
        m = st.selectbox("Metric", ["total_counts", "pct_counts_mt", "cnv_score"])
        v = np.log10(view[m] + 1) if m == "total_counts" else view[m]
        fig = umap_continuous(view, v, "log10 UMIs" if m == "total_counts" else m)
    else:
        fig = umap_categorical(
            view,
            mode,
            colors={"condition": CONDITION_COLOR, "atlas": ATLAS_COLOR, "compartment": COMPARTMENT_COLOR}[
                mode
            ],
        )
    st.caption(
        f"{len(view):,} cells shown (stratified subsample). Per-cell expression covers a panel of canonical "
        "markers, top DE genes and top variable genes; every gene is available as pseudobulk below."
    )
with c2:
    st.plotly_chart(fig)

st.subheader("Dot plot: any gene, every cell type")
all_genes = q("SELECT gene FROM genes ORDER BY gene").gene.tolist()
sel = st.multiselect(
    "Genes",
    all_genes,
    default=[g for g in ["KRT5", "TP63", "SOX2", "SPP1", "CD8A", "PDCD1", "TIGIT"] if g in all_genes],
)
cnd = st.selectbox("Condition", list(CONDITION_COLOR), index=2, format_func=CONDITION_LABEL.get)
if sel:
    d = q(
        f"""SELECT gene, cell_type, mean_lognorm, frac_expressing, n_cells FROM expr_pseudobulk
              WHERE condition = ? AND gene IN ({",".join("?" * len(sel))})""",
        (cnd, *sel),
    )
    fig = go.Figure(
        go.Scatter(
            x=d.gene,
            y=d.cell_type,
            mode="markers",
            marker=dict(
                size=4 + 22 * d.frac_expressing,
                color=d.mean_lognorm,
                colorscale=[[0, "#e8e7e3"], [1, "#0d366b"]],
                colorbar=dict(title="mean log-expr", thickness=12),
                line=dict(width=0.5, color="white"),
            ),
            customdata=np.stack([d.frac_expressing, d.n_cells], axis=1),
            hovertemplate="%{x} in %{y}<br>mean %{marker.color:.2f}<br>%{customdata[0]:.0%} of %{customdata[1]:,} cells<extra></extra>",
        )
    )
    fig.update_xaxes(categoryorder="array", categoryarray=sel)
    fig.update_yaxes(categoryorder="category descending")
    st.plotly_chart(style(fig, 760))
    st.caption(
        "Dot size = fraction of cells expressing; colour = mean log-normalised expression. All cells, all genes; pooled over donors (see page 5 for per-patient tests)."
    )
footer()
