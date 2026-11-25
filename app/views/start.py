"""Start here: what the project is and why it exists, the workflow as a map, what / why / how for every stage,
and how to reproduce it (data, website, full pipeline)."""

import streamlit as st
from lib.db import meta, q, setting
from lib.explain import findings
from lib.ui import PART, REPO_BRANCH, REPO_URL, STAGE_PART, confound_banner, footer, part_tag, two_parts

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
TREE = f"{REPO_URL}/tree/{REPO_BRANCH}"
BLOB = f"{REPO_URL}/blob/{REPO_BRANCH}"

# ---------------------------------------------------------------- what this is
st.title("Start here")
st.markdown(
    "**Lung squamous cell carcinoma (LUSC)** is one of the two main types of non-small-cell lung cancer and is "
    "closely linked to smoking. Unlike lung adenocarcinoma, it rarely carries the driver mutations that targeted "
    "drugs exist for, so treatment relies mostly on chemotherapy and immunotherapy. Those treatments act on the "
    "tumour **and on the immune and structural cells around it**, so it matters what all of those cells are doing."
)
st.markdown(
    "This project asks: **what changes in each type of cell when lung tissue becomes LUSC, and where are those "
    "cells in the tissue?** It answers that by re-analysing public single-cell and spatial data with one "
    "reproducible pipeline. Every number on this site is read from that pipeline's outputs."
)

c1, c2, c3 = st.columns(3)
with c1, st.container(border=True):
    st.markdown("**Why single-cell?**")
    st.caption(
        "Bulk sequencing averages every cell in a sample, so a change in a few immune cells is lost. Single-cell "
        "RNA sequencing measures each cell separately, which shows changes per cell type."
    )
with c2, st.container(border=True):
    st.markdown("**Why spatial?**")
    st.caption(
        "Single-cell data loses where each cell was. Spatial transcriptomics measures gene expression on a tissue "
        "slide, which shows which cells sit next to the tumour and which stay away."
    )
with c3, st.container(border=True):
    st.markdown("**What this project adds**")
    st.caption(
        "Two large public atlases merged without double-counting patients, statistics done per patient (not per "
        "cell), tumour versus normal compared within the same patients, and every decision logged with its reason."
    )

# ---------------------------------------------------------------- the data
st.subheader("The data, in plain terms")
st.caption("All public. Nothing new was sequenced; fixed releases are downloaded and checksum-verified.")
d1, d2, d3 = st.columns(3)
with d1, st.container(border=True):
    st.markdown("**HLCA · Human Lung Cell Atlas**")
    st.caption(
        "A curated combination of many single-cell studies of **healthy lung** (Sikkema et al. 2023, *Nature "
        "Medicine*). Here it is the reference for what normal lung cells look like."
    )
with d2, st.container(border=True):
    st.markdown("**LuCA · Lung Cancer Atlas**")
    st.caption(
        "A single-cell atlas of **lung tumours** (Salcher et al. 2022, *Cancer Cell*). It supplies the LUSC "
        "tumour samples and, for some patients, **adjacent normal** tissue from the same person, which allows "
        "within-patient comparisons."
    )
with d3, st.container(border=True):
    st.markdown("**10x Visium section**")
    st.caption(
        "One lung squamous cell carcinoma slide from 10x Genomics (one patient), with gene expression measured "
        "in a few thousand small spots (55 µm) on top of an H&E image. It shows **where** cell types sit in the "
        "tissue. Its identity is checked with squamous markers before use."
    )

two_parts()
confound_banner()

# ---------------------------------------------------------------- workflow
st.subheader("The workflow")
st.markdown(
    "The pages fall into three groups. **How the data were built** makes one clean, comparable atlas. "
    "**Evidence** compares tumour with adjacent normal tissue of the same patients; the claims of this site "
    "come from there. **Context** describes the cancer cells, their trajectory, one tissue section and "
    "predicted cell–cell signals; it orients the reader but is not evidence. Each box is a page on this site."
)
st.graphviz_chart(
    """digraph {
    rankdir=TB; newrank=true; nodesep=0.25; ranksep=0.35;
    node [shape=box, style="rounded", fontname="Helvetica", fontsize=10, margin="0.14,0.06"];
    edge [color="#888", fontname="Helvetica", fontsize=9, fontcolor="#888"];
    graph [fontname="Helvetica", fontsize=11, fontcolor="#888", color="#bbb", style="rounded,dashed",
           labeljust=l, margin=10];

    subgraph cluster_in { label="INPUTS · public data";
        hlca [label="HLCA\\nHuman Lung Cell Atlas\\nhealthy lung reference"];
        luca [label="LuCA\\nLung Cancer Atlas\\nLUSC tumours + adjacent normal"];
        vis  [label="10x Visium\\none LUSC tissue section\\nexpression + H&E image"]; }

    subgraph cluster_a { label="HOW THE DATA WERE BUILT";
        cohort [label="Cohort & merge\\ncombine atlases, keep shared studies once,\\nshared genes, one cell-type vocabulary"];
        qc     [label="Quality control\\nremove damaged cells, empty droplets\\nand doublets, per-sample cut-offs"];
        integ  [label="Integration\\nremove lab/study effects from the map;\\n7 methods benchmarked"]; }

    subgraph cluster_e { label="EVIDENCE · same-patient comparisons (tumour vs adjacent normal, immune cells)";
        color="#1c5cab"; fontcolor="#1c5cab";
        de   [label="DE & pathways\\ngenes, pathways and transcription\\nfactors per cell type"];
        comp [label="Composition\\nmore or fewer of each\\nimmune cell type?"]; }

    subgraph cluster_c { label="CONTEXT · descriptive (tumour-only patients, healthy donors, one section)";
        color="#a3501f"; fontcolor="#a3501f";
        mal  [label="Malignant cells\\ncopy-number changes inferred from expression"];
        traj [label="Trajectory\\nhow far cancer cells are\\nfrom healthy basal cells"];
        sp   [label="Spatial\\ncell-type mix per tissue spot,\\ntissue niches"];
        ccc  [label="Cell–cell communication\\npredicted ligand–receptor signals"]; }

    db [label="DuckDB snapshot\\none read-only file → this website", style="rounded,dashed"];

    hlca -> cohort; luca -> cohort; cohort -> qc -> integ;
    integ -> de; integ -> comp; integ -> mal -> traj;
    vis -> sp; integ -> sp [label=" annotated cells\\n as reference"];
    de -> ccc [label=" differential scores\\n from the paired DE"]; sp -> ccc [label=" spatial check"];
    ccc -> db; de -> db [style=dotted]; comp -> db [style=dotted]; traj -> db [style=dotted];
    sp -> db [style=dotted];
    }""",
    width="stretch",
)

# ---------------------------------------------------------------- stage by stage
st.subheader("Stage by stage: what, why and how")
st.caption(
    "Open a stage to see its evidence, decisions and caveats. The last line of each card is computed from the data."
)
guide = q('SELECT * FROM page_guide ORDER BY "order"').set_index("stage", drop=False)
for part in PART:
    stages = [st_ for st_ in guide.index if STAGE_PART.get(st_) == part]
    if not stages:
        continue
    part_tag(part)
    for stage in stages:
        r = guide.loc[stage]
        with st.container(border=True):
            st.page_link(PAGES[stage], label=f"**{r.title}**: {r.one_line}")
            w, y, h = st.columns(3, gap="medium")
            w.markdown(f"**What**  \n{r.stage_what}")
            y.markdown(f"**Why**  \n{r.stage_why}")
            h.markdown(f"**How**  \n{r.stage_how}")
            f = findings(stage)
            if f:
                st.caption(f"**Headline result:** {f[0]}")

# ---------------------------------------------------------------- reproduce
st.subheader("Reproduce it")
m = meta()
hf_repo = setting("LUSC_DB_HF_REPO")
hf_rev = setting("LUSC_DB_HF_REVISION") or "main"
st.markdown(
    "Everything here can be rebuilt from public data with the code in the "
    f"[GitHub repository]({TREE}). Three levels, from quickest to most complete:"
)
t1, t2, t3 = st.tabs(
    ["1 · Query the results (minutes)", "2 · Run this website locally", "3 · Rerun the whole analysis"]
)
with t1:
    if hf_repo:
        st.markdown(
            f"Download the snapshot this site is serving, `lusc_atlas.duckdb` "
            f"(schema v{m.get('schema_version', '?')}, built {m.get('built_at', '?')}), from the public "
            f"[Hugging Face Dataset](https://huggingface.co/datasets/{hf_repo}) and query it with DuckDB from "
            "Python, R or the command line. No account is needed."
        )
    else:
        st.markdown(
            f"Download the snapshot (`lusc_atlas.duckdb`, schema v{m.get('schema_version', '?')}) from the "
            f"project's public Hugging Face Dataset (linked in the [README]({BLOB}/README.md)) and query it with "
            "DuckDB from Python, R or the command line."
        )
    st.code(
        "pip install duckdb huggingface_hub\n\n"
        "from huggingface_hub import hf_hub_download\n"
        "import duckdb\n\n"
        f'path = hf_hub_download("{hf_repo or "<hf-user>/lusc-atlas-db"}", "lusc_atlas.duckdb",\n'
        f'                       repo_type="dataset", revision="{hf_rev}")\n'
        "con = duckdb.connect(path, read_only=True)\n\n"
        'con.sql("SELECT * FROM tables_info").show()     # every table, its size and the pipeline output it came from\n'
        'con.sql("SELECT * FROM known_issues").show()    # read these before quoting a number\n'
        'con.sql("""SELECT gene, log2FoldChange, padj FROM de_pseudobulk\n'
        "           WHERE contrast = 'paired_tumour_vs_adjacent'\n"
        "             AND \"group\" = 'Macrophage (monocyte-derived)'\n"
        '           ORDER BY padj LIMIT 10""").show()',
        language="python",
    )
    st.caption(
        "Each table is a direct copy or summary of a pipeline output, and `tables_info` says which. The database "
        "is also exposed as an MCP server for AI assistants (`pixi run -e dashboard mcp`)."
    )
with t2:
    st.markdown(
        "Clone the repository and start the app with [pixi](https://pixi.sh), which installs the exact "
        "pinned environment. Point it at the downloaded snapshot, or build the snapshot yourself from pipeline "
        "outputs."
    )
    st.code(
        f"git clone {REPO_URL} && cd sc-spatial-analysis\n"
        f"git checkout {REPO_BRANCH}\n"
        "pixi install -e dashboard\n\n"
        "export LUSC_DB=/path/to/lusc_atlas.duckdb    # the file from step 1 …\n"
        "# pixi run -e dashboard export               # … or rebuild it from your own pipeline run\n"
        "pixi run -e dashboard app                    # → http://localhost:8501",
        language="bash",
    )
    st.caption(
        f"Hosting your own public copy for free (Streamlit Community Cloud + Hugging Face Dataset): "
        f"[docs/design/C_dashboard_hosting.md]({BLOB}/docs/design/C_dashboard_hosting.md)."
    )
with t3:
    st.markdown(
        "The full pipeline downloads the raw data, runs all nine stages and writes every table, figure and report. "
        "It runs on CPU only (developed on 16 cores and 125 GB RAM, no GPU), needs about 40 GB for raw data and "
        "takes one to two days. Steps are resumable, so an interrupted run continues where it stopped."
    )
    st.code(
        f"git clone {REPO_URL} && cd sc-spatial-analysis\n"
        f"git checkout {REPO_BRANCH}\n"
        "pixi install --all                           # every environment, pinned in pixi.lock\n"
        "pixi run python scripts/run_all.py --dry-run # show the plan\n"
        "pixi run run-all                             # run every step (or: pixi run -e workflow smk)\n"
        "pixi run -e dashboard export                 # pipeline outputs → results/dashboard/lusc_atlas.duckdb\n"
        "pixi run -e dashboard app                    # browse your own results",
        language="bash",
    )
    st.caption(
        "Inputs, versions and settings are in "
        f"[config/config.yaml]({BLOB}/config/config.yaml); every decision and its reason is in "
        f"[DEVELOPMENT_LOG.md]({BLOB}/DEVELOPMENT_LOG.md); the step-by-step guide is the "
        f"[README Quick start]({BLOB}/README.md#quick-start)."
    )

cols = st.columns(2)
cols[0].page_link("views/glossary.py", label="Glossary of terms", icon="📚")
cols[1].page_link("views/provenance.py", label="Provenance of every number", icon="🧾")
footer()
