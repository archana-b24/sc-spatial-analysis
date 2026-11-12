import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.plots import condition_bar
from lib.ui import SERIES, caveats, decisions, footer, metric_row, page_header, show_table, style

page_header(
    "1 · Cohort & merge",
    "Which patients and cells can fairly compare LUSC tumour with normal lung, and how do two independently "
    "curated atlases become one dataset without counting anyone twice?",
)
guide("cohort")
decisions("cohort")
caveats("cohort")

gs = q("SELECT * FROM merge_gene_space")
donors = q("SELECT * FROM cohort_donors")
paired = q("""SELECT donor_id, atlas FROM cohort_donors GROUP BY donor_id, atlas
              HAVING bool_or(condition='tumour') AND bool_or(condition='adjacent_normal')""")
lusc = donors[donors.condition != "healthy"].donor_id.nunique()
metric_row(
    [
        ("Donors", f"{donors.donor_id.nunique()}", None),
        ("LUSC patients", f"{lusc}", "donors with any tumour or adjacent sample"),
        ("Paired patients", f"{len(paired)} / {lusc}", "tumour + adjacent from the same person"),
        ("Shared genes", f"{int(gs.n_shared.iloc[0]):,}", "Ensembl-ID intersection of the two atlases"),
    ]
)

st.subheader("Cells per atlas and condition")
cells = q("""SELECT atlas AS grp, condition, sum(n_cells) n FROM cohort_donors GROUP BY 1, 2""")
st.plotly_chart(condition_bar(cells, "grp", "n", height=280))
st.caption("Cells selected per atlas (before QC). HLCA alone has only 3 LUSC donors from one study.")

st.subheader("Study crosswalk: one copy of every study")
cw = q("SELECT * FROM merge_study_crosswalk ORDER BY atlas, n_cells DESC")
fig = go.Figure()
for i, cond in enumerate(["n_healthy", "n_adjacent_normal", "n_tumour", "n_tumour_region_unknown"]):
    fig.add_trace(
        go.Bar(
            y=cw["atlas"] + " · " + cw["study"],
            x=cw[cond],
            name=cond[2:].replace("_", " "),
            orientation="h",
            marker_color=SERIES[i],
        )
    )
fig.update_layout(barmode="stack", yaxis=dict(autorange="reversed"))
st.plotly_chart(style(fig, 560))
st.caption(
    "Each study is taken from exactly one atlas. `decision` and `counts_ok` (integer-count check) are in the table."
)
show_table(cw)

st.subheader("What was excluded, and why")
ex = q(
    "SELECT reason, atlas, sum(n_cells) n_cells, count(*) n_studies FROM merge_exclusions GROUP BY 1, 2 ORDER BY 3 DESC"
)
fig = go.Figure(
    go.Bar(
        y=ex["atlas"] + " · " + ex["reason"],
        x=ex["n_cells"],
        orientation="h",
        marker_color=SERIES[0],
        hovertemplate="%{y}<br>%{x:,} cells<extra></extra>",
    )
)
fig.update_layout(yaxis=dict(autorange="reversed"))
st.plotly_chart(style(fig, 360))
show_table(q("SELECT * FROM merge_exclusions ORDER BY n_cells DESC"))

st.subheader("Paired patients")
st.markdown(
    "The primary differential-expression design compares tumour with adjacent-normal tissue **within the same patient**."
)
show_table(
    donors[donors.donor_id.isin(paired.donor_id)].sort_values(["donor_id", "condition"]), "Paired donors"
)
footer()
