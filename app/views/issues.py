import streamlit as st
from lib.db import q
from lib.ui import SEVERITY_ICON, footer

st.title("Known issues")
st.markdown(
    "The pipeline ran end to end, but reading the outputs (not only the exit codes) showed the problems below. "
    "They are listed so that no result is quoted without its caveat. Open issues are planned fixes; documented "
    "ones are limitations that stay and are stated."
)
df = q(
    "SELECT * FROM known_issues ORDER BY CASE severity WHEN 'high' THEN 0 WHEN 'medium' THEN 1 ELSE 2 END, id"
)
for r in df.itertuples():
    with st.container(border=True):
        st.markdown(f"#### {SEVERITY_ICON.get(r.severity, '')} {r.title}")
        st.markdown(r.detail)
        st.caption(
            f"Severity: {r.severity} · Status: {r.status} · Stage: {r.stage} · DEVELOPMENT_LOG {r.log_ref}"
        )
footer()
