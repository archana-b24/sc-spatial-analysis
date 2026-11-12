import streamlit as st
from lib.db import q
from lib.plots import condition_bar
from lib.ui import footer, metric_row

st.title("Lung squamous cell carcinoma: a single-cell and spatial atlas")
st.markdown(
    "This site shows the results of a reproducible analysis of **lung squamous cell carcinoma (LUSC)**. "
    "It asks what changes in the cells of LUSC tumours compared with normal lung, using public single-cell "
    "RNA-seq atlases and one spatial transcriptomics section. Each page covers one pipeline stage and gives "
    "**the question**, **the decision we made and why**, **the evidence**, and **the caveats**."
)

n = (
    q("""SELECT count(*) cells, count(DISTINCT donor_id) donors, count(DISTINCT study) studies,
               count(DISTINCT sample) samples, sum(CASE WHEN malignant::VARCHAR IN ('True','true','1') THEN 1 ELSE 0 END) malignant
        FROM cells""")
    .iloc[0]
    .astype(int)
)
paired = q("""SELECT count(*) n FROM (SELECT donor_id FROM cohort_donors GROUP BY donor_id
              HAVING bool_or(condition='tumour') AND bool_or(condition='adjacent_normal'))""").iloc[0, 0]
spots = q("SELECT count(*) FROM visium_spots").iloc[0, 0]
open_issues = q("SELECT count(*) FROM known_issues WHERE status = 'open'").iloc[0, 0]
metric_row(
    [
        ("Cells after QC", f"{n.cells:,}", "HLCA + LuCA, after per-sample QC and doublet removal"),
        ("Donors", f"{n.donors:,}", f"{n.samples} samples from {n.studies} studies"),
        (
            "Paired LUSC patients",
            f"{paired}",
            "tumour and adjacent-normal tissue from the same patient: the primary DE design",
        ),
        ("Malignant cells (CNV)", f"{n.malignant:,}", "called from inferred copy-number changes"),
        ("Visium spots", f"{spots:,}", "one FFPE LUSC section (10x CytAssist)"),
    ]
)
if open_issues:
    st.warning(
        f"**{open_issues} open issues** change how some results should be read. They are flagged on each page."
    )
    st.page_link("views/issues.py", label="See all known issues", icon="⚠️")

c1, c2 = st.columns([3, 2])
with c1:
    st.subheader("Who is in the cohort")
    df = q("""SELECT atlas, condition, count(DISTINCT donor_id) donors, count(*) cells
              FROM cells GROUP BY 1, 2 ORDER BY 1, 2""")
    df["group"] = df["atlas"]
    st.plotly_chart(condition_bar(df, "group", "donors", height=300))
    st.caption(
        "Donors per atlas and condition. HLCA supplies the healthy reference; LuCA supplies the tumours."
    )
with c2:
    st.subheader("The pipeline")
    st.graphviz_chart(
        """digraph {
        rankdir=TB; node [shape=box, style="rounded", fontname="Helvetica", fontsize=11]; edge [color="#888"];
        fetch [label="Fetch HLCA + LuCA + Visium\\n(pinned, checksummed)"];
        select [label="1 · Cohort & merge"]; qc [label="2 · QC & doublets"];
        integ [label="3 · Integration benchmark"]; cnv [label="4 · Malignant cells (CNV)"];
        de [label="5 · DE & pathways"]; ab [label="6 · Composition"]; traj [label="7 · Trajectory"];
        sp [label="8 · Spatial (cell2location, niches)"]; ccc [label="9 · Cell–cell communication"];
        exp [label="export → this dashboard", style="rounded,dashed"];
        fetch -> select -> qc -> integ -> cnv; cnv -> de; cnv -> ab; cnv -> traj; cnv -> sp; de -> ccc; sp -> ccc;
        ccc -> exp }""",
        width="stretch",
    )

st.subheader("Where to go next")
cols = st.columns(3)
links = [
    ("views/cohort.py", "1 · Cohort & merge", "why two atlases, and how duplicates were removed"),
    ("views/integration.py", "3 · Integration", "7 methods benchmarked, one chosen"),
    ("views/malignant.py", "4 · Malignant cells", "copy-number calling and its sensitivity"),
    ("views/de.py", "5 · DE & pathways", "paired tumour vs adjacent, per cell type"),
    ("views/spatial.py", "8 · Spatial", "deconvolved H&E section and niches"),
    ("views/explorer.py", "Cell & gene explorer", "colour the atlas by any marker gene"),
]
for i, (p, label, sub) in enumerate(links):
    with cols[i % 3]:
        st.page_link(p, label=label)
        st.caption(sub)
footer()
