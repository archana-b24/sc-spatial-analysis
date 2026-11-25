import json

import plotly.graph_objects as go
import streamlit as st
from lib.db import meta, q
from lib.ui import decisions, footer, style

st.title("Provenance")
st.markdown(
    "Where every number on this site comes from. The database is rebuilt from pipeline outputs by the "
    "`export` step, so it cannot drift from the code that produced the results."
)
decisions("provenance")
m = meta()
c = st.columns(4)
c[0].metric("Commit", m.get("git_sha", "?")[:7])
c[1].metric("Branch", m.get("git_branch", "?"))
c[2].metric("Uncommitted changes at build", "yes" if m.get("git_dirty") == "True" else "no")
c[3].metric("Snapshot schema", f"v{m.get('schema_version', '?')}")

st.subheader("Step runtimes (last full run)")
runs = q("SELECT * FROM step_runs ORDER BY finished")
fig = go.Figure(
    go.Bar(
        y=runs.step,
        x=runs.runtime_s / 3600,
        orientation="h",
        marker_color="#1c5cab",
        customdata=runs.finished,
        hovertemplate="%{y}: %{x:.2f} h<br>finished %{customdata}<extra></extra>",
    )
)
fig.update_layout(yaxis=dict(autorange="reversed"), xaxis=dict(title="hours (16 CPU cores, no GPU)"))
st.plotly_chart(style(fig, 460))
st.caption(f"Total ≈ {runs.runtime_s.sum() / 3600:.1f} h. cell2location deconvolution dominates (CPU-only).")

st.subheader("Package versions per step")
step = st.selectbox("Step", runs.step.tolist())
pk = json.loads(runs.set_index("step").loc[step, "packages"])
st.dataframe([{"package": k, "version": v} for k, v in pk.items()], hide_index=True)

st.subheader("One serving source of truth")
st.markdown(
    "The pipeline's AnnData/SpatialData objects and result tables are the **computational** source of truth. "
    "This DuckDB file is the **serving** source of truth: every page on this site "
    "reads only this file, and it is rebuilt from scratch by `export`, never edited."
)
st.dataframe(
    [
        {"field": "config.yaml sha256", "value": m.get("config_sha256", "?")},
        {"field": "dashboard_content.yaml sha256", "value": m.get("content_sha256", "?")},
        {"field": "built at", "value": m.get("built_at", "?")},
    ],
    hide_index=True,
)
stamps = q("SELECT * FROM step_stamps ORDER BY stamped_at")
if not stamps.empty:
    st.caption(f"{len(stamps)} pipeline steps had a clean-completion stamp when this snapshot was built:")
    st.dataframe(stamps, hide_index=True)

st.subheader("What this database contains")
st.dataframe(q("SELECT * FROM tables_info ORDER BY table_name"), hide_index=True)
with st.expander("Schema contract (tables and columns the site relies on)"):
    st.dataframe(
        q(
            "SELECT table_name, string_agg(column_name, ', ') AS columns FROM schema_contract GROUP BY 1 ORDER BY 1"
        ),
        hide_index=True,
    )
with st.expander("config/config.yaml used for this snapshot"):
    st.code(m.get("config_yaml", ""), language="yaml")
with st.expander("Auto-generated results summary (results/reports/summary.md)"):
    st.markdown(m.get("summary_md", ""))
footer()
