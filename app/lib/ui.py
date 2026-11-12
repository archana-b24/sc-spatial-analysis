"""Shared page furniture: palette, the question → decision → evidence → caveats panel, provenance footer.

Colour rules (dataviz method): categorical hues in fixed order and never cycled; condition always has
the same colour on every page; sequential = one hue; diverging = blue ↔ red around a grey midpoint.
27 cell types exceed any safe categorical palette, so cell-type views use compartment colour or
"highlight one, grey the rest" rather than 27 hues.
"""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from .db import meta, q

# validated reference palette, light-mode steps (dataviz references/palette.md)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
GREY = "#b9b8b2"
CONDITION_ORDER = ["healthy", "adjacent_normal", "tumour", "tumour_region_unknown"]
CONDITION_COLOR = dict(zip(CONDITION_ORDER, SERIES[:4], strict=True))
CONDITION_LABEL = {
    "healthy": "Healthy donor",
    "adjacent_normal": "Adjacent normal (LUSC patient)",
    "tumour": "LUSC tumour",
    "tumour_region_unknown": "Tumour, region unknown",
}
ATLAS_COLOR = {"HLCA": SERIES[0], "LuCA": SERIES[1]}
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
DIVERGING = [[0, "#1c5cab"], [0.5, "#f0efec"], [1, "#c23b3a"]]
SEVERITY_ICON = {"high": "🔴", "medium": "🟠", "low": "🟡"}

LOG_URL = "https://github.com/archana-b24/sc-spatial-analysis/blob/main/DEVELOPMENT_LOG.md"


def page_header(title: str, question: str) -> None:
    st.title(title)
    st.markdown(f"**The question.** {question}")


def decisions(stage: str) -> None:
    """The 'why' panel: each decision for this stage, as question → decision → why → evidence."""
    df = q("SELECT * FROM decisions WHERE stage = ?", (stage,))
    if df.empty:
        return
    with st.expander(
        f"Why it was done this way ({len(df)} decision{'s' if len(df) > 1 else ''})", expanded=False
    ):
        for r in df.itertuples():
            st.markdown(
                f"**{r.question}**  \n"
                f"*Decision:* {r.decision}  \n"
                f"*Why:* {r.why}  \n"
                f"<span style='opacity:.7'>Evidence: `{r.evidence}` · DEVELOPMENT_LOG {r.log_ref}</span>",
                unsafe_allow_html=True,
            )


def caveats(stage: str) -> None:
    """Known issues for this stage, shown above the results so nobody quotes a number without them."""
    df = q("SELECT * FROM known_issues WHERE stage = ? ORDER BY id", (stage,))
    for r in df.itertuples():
        icon = SEVERITY_ICON.get(r.severity, "⚠️")
        head = f"{icon} **Known issue ({r.severity}, {r.status}): {r.title}.**"
        body = f"{head} {r.detail} *(DEVELOPMENT_LOG {r.log_ref})*"
        (st.warning if r.severity == "high" else st.info)(body)


def footer() -> None:
    m = meta()
    sha = m.get("git_sha", "unknown")[:7]
    dirty = " + uncommitted changes" if m.get("git_dirty") == "True" else ""
    st.divider()
    st.caption(
        f"Snapshot built {m.get('built_at', '?')} from commit `{sha}`{dirty} "
        f"(branch `{m.get('git_branch', '?')}`), schema v{m.get('schema_version', '?')}. "
        "Every number is read from pipeline outputs. See the Provenance page."
    )


def style(fig: go.Figure, height: int = 420) -> go.Figure:
    """Recessive grid/axes, legend on top, no chart junk."""
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, title_text=""),
        hoverlabel=dict(namelength=-1),
        font=dict(size=13),
    )
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(gridcolor="rgba(128,128,128,0.15)", zeroline=False)
    return fig


def metric_row(items: list[tuple[str, str, str | None]]) -> None:
    cols = st.columns(len(items))
    for c, (label, value, help_) in zip(cols, items, strict=True):
        c.metric(label, value, help=help_)


def show_table(df, caption: str | None = None, **kw) -> None:
    """Every chart gets a table view (accessibility + the numbers people actually quote)."""
    with st.expander("Show as table"):
        if caption:
            st.caption(caption)
        st.dataframe(df, hide_index=True, **kw)
