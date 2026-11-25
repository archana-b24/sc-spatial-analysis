"""LUSC atlas results dashboard — entry point.

    pixi run -e dashboard app            # local, reads results/dashboard/lusc_atlas.duckdb

One page per pipeline stage, grouped into how the data were built, Evidence (same-patient comparisons)
and Context (descriptive) (DEVELOPMENT_LOG §3.13); each answers question → decision → evidence → caveats
(docs/ROADMAP.md §6), with a reader guide (what / how to read / what it means / caveats;
docs/design/README.md). The app is read-only: the DuckDB snapshot that `lusc-atlas export` builds from
pipeline outputs is its single serving source of truth.
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
        st.Page(f"{V}/glossary.py", title="Glossary", icon="📚"),
    ],
    "How the data were built": [
        st.Page(f"{V}/cohort.py", title="Cohort & merge", icon="🧬"),
        st.Page(f"{V}/qc.py", title="Quality control", icon="🧹"),
        st.Page(f"{V}/integration.py", title="Integration", icon="🧩"),
    ],
    # §3.13: the claims come from the paired, within-patient comparisons; everything else is descriptive
    "Evidence · same-patient comparisons": [
        st.Page(f"{V}/de.py", title="DE & pathways", icon="📊"),
        st.Page(f"{V}/composition.py", title="Composition", icon="🥧"),
    ],
    "Context · descriptive": [
        st.Page(f"{V}/malignant.py", title="Malignant cells", icon="🎯"),
        st.Page(f"{V}/trajectory.py", title="Trajectory", icon="🧭"),
        st.Page(f"{V}/spatial.py", title="Spatial", icon="🗺️"),
        st.Page(f"{V}/ccc.py", title="Cell–cell communication", icon="📡"),
    ],
    "Explore": [
        st.Page(f"{V}/explorer.py", title="Cell & gene explorer", icon="🔎"),
        st.Page(f"{V}/provenance.py", title="Provenance", icon="🧾"),
    ],
}
st.navigation(pages).run()
