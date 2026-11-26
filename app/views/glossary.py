import streamlit as st
from lib.db import q
from lib.ui import footer

st.title("Glossary")
st.markdown("Plain-language definitions of the terms used on this site.")
g = q("SELECT term, definition FROM glossary ORDER BY lower(term)")
needle = st.text_input("Filter", placeholder="e.g. padj, pseudobulk, niche").strip().lower()
if needle:
    g = g[
        g.term.str.lower().str.contains(needle, regex=False)
        | g.definition.str.lower().str.contains(needle, regex=False)
    ]
st.caption(f"{len(g)} terms")
for r in g.itertuples():
    st.markdown(f"**{r.term}** — {r.definition}")
footer()
