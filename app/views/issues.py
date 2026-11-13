import streamlit as st
from lib.db import q
from lib.ui import SEVERITY_ICON, footer

st.title("Caveats & limitations")
st.markdown(
    "What to keep in mind before quoting a result. Each item below is a known limitation of the data or the "
    "methods; it is not a bug waiting for a fix, so it is also flagged on the stage page it affects. "
    "Problems found after run 1 and fixed in run 2 are listed at the bottom for the record."
)

df = q(
    "SELECT * FROM known_issues ORDER BY CASE severity WHEN 'high' THEN 0 WHEN 'medium' THEN 1 ELSE 2 END, id"
)
resolved = df["status"].astype(str).str.startswith("resolved")
current, fixed = df[~resolved], df[resolved]


def card(r, icon: str) -> None:
    with st.container(border=True):
        st.markdown(f"#### {icon} {r.title}")
        st.markdown(r.detail)
        st.caption(
            f"Severity: {r.severity} · Status: {r.status} · Stage: {r.stage} · DEVELOPMENT_LOG {r.log_ref}"
        )


if current.empty:
    st.info("No current limitations are recorded.")
for r in current.itertuples():
    card(r, SEVERITY_ICON.get(r.severity, "⚠️"))

if not fixed.empty:
    st.subheader("Fixed in run 2")
    with st.expander(f"{len(fixed)} problem{'s' if len(fixed) > 1 else ''} found in run 1 and fixed"):
        for r in fixed.itertuples():
            card(r, "✅")
footer()
