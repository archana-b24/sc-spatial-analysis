"""Start here: what this site is, the analysis workflow as a map, and one computed headline per stage."""

import streamlit as st
from lib.db import meta, q
from lib.explain import findings
from lib.ui import footer

PAGES = {
    "cohort": "views/cohort.py",
    "qc": "views/qc.py",
    "integration": "views/integration.py",
    "malignant": "views/malignant.py",
    "de": "views/de.py",
    "composition": "views/composition.py",
    "trajectory": "views/trajectory.py",
    "spatial": "views/spatial.py",
    "ccc": "views/ccc.py",
}

st.title("Start here")
st.markdown(
    "This site walks through a reproducible analysis of **lung squamous cell carcinoma (LUSC)**: what changes "
    "in the cells of LUSC tumours compared with normal lung, using two public single-cell atlases and one "
    "spatial transcriptomics section. You do not need to code to use it."
)

c1, c2 = st.columns(2)
with c1, st.container(border=True):
    st.markdown("**1 · Follow the workflow**")
    st.caption(
        "Each stage below has its own page: the question, the decision and why, the evidence, and caveats."
    )
with c2, st.container(border=True):
    st.markdown("**2 · Read the guide on each page**")
    st.caption(
        "*What you're looking at · How to read it · What it means · Caveats.* Interpretations are "
        "computed from the data, not written by hand."
    )

st.subheader("The workflow")
st.graphviz_chart(
    """digraph {
    rankdir=LR; node [shape=box, style="rounded", fontname="Helvetica", fontsize=11, margin="0.15,0.06"];
    edge [color="#888"];
    data [label="Public data\\nHLCA + LuCA + Visium", style="rounded,dashed"];
    cohort [label="1 · Cohort\\n& merge"]; qc [label="2 · QC"]; integ [label="3 · Integration"];
    mal [label="4 · Malignant\\ncells (CNV)"]; de [label="5 · DE &\\npathways"]; comp [label="6 · Composition"];
    traj [label="7 · Trajectory"]; sp [label="8 · Spatial"]; ccc [label="9 · Cell–cell\\ncommunication"];
    db [label="DuckDB snapshot\\n→ this site", style="rounded,dashed"];
    data -> cohort -> qc -> integ -> mal; mal -> de; mal -> comp; mal -> traj; data -> sp; integ -> sp;
    de -> ccc; sp -> ccc; ccc -> db }""",
    width="stretch",
)
st.caption(
    "Stages 5–8 all start from the integrated, malignant-annotated atlas; spatial analysis also uses the "
    "Visium section directly. The dashed box is the read-only database everything on this site reads."
)

open_issues = q("SELECT count(*) n FROM known_issues WHERE status = 'open'").n.iloc[0]
if open_issues:
    st.warning(f"**{open_issues} open known issues** affect some stages; they are flagged where they apply.")

st.subheader("Stage by stage")
guide = q('SELECT * FROM page_guide ORDER BY "order"')
for r in guide.itertuples():
    with st.container(border=True):
        a, b = st.columns([2, 3], gap="large")
        with a:
            st.page_link(PAGES[r.stage], label=f"{r.order} · {r.title}")
            st.caption(r.one_line)
        with b:
            f = findings(r.stage)
            if f:
                st.markdown(f[0])

st.subheader("Try it yourself")
m = meta()
st.markdown(
    "Everything on this site comes from one file, `lusc_atlas.duckdb` "
    f"(schema v{m.get('schema_version', '?')}, built {m.get('built_at', '?')}). You can query it directly:"
)
st.code(
    "pip install duckdb\n"
    "duckdb lusc_atlas.duckdb -readonly\n"
    "  SELECT * FROM tables_info;                      -- what is in it\n"
    "  SELECT gene, log2FoldChange, padj FROM de_pseudobulk\n"
    "   WHERE contrast = 'paired_tumour_vs_adjacent' AND \"group\" = 'Macrophage'\n"
    "   ORDER BY padj LIMIT 10;",
    language="sql",
)
st.caption("Download links and local setup are in the repository's docs/design/C_dashboard_hosting.md.")
cols = st.columns(3)
cols[0].page_link("views/glossary.py", label="Glossary", icon="📚")
cols[1].page_link("views/issues.py", label="Caveats & limitations", icon="⚠️")
cols[2].page_link("views/provenance.py", label="Provenance", icon="🧾")
footer()
