"""LUSC atlas results dashboard — entry point.

    pixi run -e dashboard app            # local, reads results/dashboard/lusc_atlas.duckdb

One page per pipeline stage; each answers question → decision → evidence → caveats
(docs/ROADMAP_v2.1.md §6), with a reader guide (what / how to read / what it means / caveats; v3,
docs/v3/README.md). The app is read-only: the DuckDB snapshot that `lusc-atlas export` builds from pipeline
outputs is its single serving source of truth.
"""

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).parent))
# make <repo>/src importable for a plain `pip install -r app/requirements.txt` checkout (e.g. Streamlit
# Community Cloud); harmless when pixi has installed the package
for extra in (Path(__file__).parents[1] / "src", Path(__file__).parents[1]):
    if (extra / "lusc_atlas").is_dir() and str(extra) not in sys.path:
        sys.path.append(str(extra))

st.set_page_config(page_title="LUSC atlas", page_icon="🫁", layout="wide")

V = "views"  # page paths are resolved relative to this entrypoint by Streamlit
pages = {
    "Start here": [
        st.Page(f"{V}/start.py", title="Start here", icon="🧭", default=True),
        st.Page(f"{V}/overview.py", title="Overview", icon="🏠"),
        st.Page(f"{V}/issues.py", title="Known issues", icon="⚠️"),
        st.Page(f"{V}/glossary.py", title="Glossary", icon="📚"),
    ],
    "Pipeline stages": [
        st.Page(f"{V}/cohort.py", title="1 · Cohort & merge", icon="🧬"),
        st.Page(f"{V}/qc.py", title="2 · Quality control", icon="🧹"),
        st.Page(f"{V}/integration.py", title="3 · Integration", icon="🧩"),
        st.Page(f"{V}/malignant.py", title="4 · Malignant cells", icon="🎯"),
        st.Page(f"{V}/de.py", title="5 · DE & pathways", icon="📊"),
        st.Page(f"{V}/composition.py", title="6 · Composition", icon="🥧"),
        st.Page(f"{V}/trajectory.py", title="7 · Trajectory", icon="🧭"),
        st.Page(f"{V}/spatial.py", title="8 · Spatial", icon="🗺️"),
        st.Page(f"{V}/ccc.py", title="9 · Cell–cell communication", icon="📡"),
    ],
    "Explore": [
        st.Page(f"{V}/explorer.py", title="Cell & gene explorer", icon="🔎"),
        st.Page(f"{V}/provenance.py", title="Provenance", icon="🧾"),
    ],
}
st.navigation(pages).run()
