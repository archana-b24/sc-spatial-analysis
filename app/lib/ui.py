"""Shared page furniture: palette, the question → decision → evidence → caveats panel, provenance footer.

Colour rules (dataviz method): categorical hues in fixed order and never cycled; condition always has
the same colour on every page; sequential = one hue; diverging = blue ↔ red around a grey midpoint.
More than 25 cell types exceed any safe categorical palette, so cell-type views use compartment colour or
"highlight one, grey the rest" rather than one hue per type.
"""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from .db import meta, q

# validated reference palette, light-mode steps (dataviz references/palette.md)
# run 3 (figure review): generic categories (niches, technologies, top-k cell types) must never borrow a
# hue that means something site-wide (tissue, atlas, up/down), so SERIES excludes those; validated with the
# dataviz checker (light mode: all checks pass; two light hues need labels/table views, which every chart has)
SERIES = ["#eda100", "#b8338a", "#5c8a00", "#e87ba4", "#b0651a", "#7a3db8", "#d6336c", "#8f7a00"]
# extra hues for categorical sets above 8; beyond 12 colours cycle
SERIES_EXTRA = ["#8c5a2b", "#5c6f1e", "#3d7a9e", "#9b2f6f"]
GREY = "#b9b8b2"


def series_color(i: int) -> str:
    """Categorical colour i: the 8 validated SERIES first, then SERIES_EXTRA, then cycle."""
    pal = SERIES + SERIES_EXTRA
    return pal[i % len(pal)]


CONDITION_ORDER = ["healthy", "adjacent_normal", "tumour", "tumour_region_unknown"]
# §3.12 semantic colour system: one meaning per colour across every page (validated with the dataviz
# palette checker,
# all pairs, light mode). Tissue = aqua / blue / orange (region-unknown tumour IS tumour tissue → same
# orange);
# atlas = violet / teal; compartments use hues that never stand for a tissue or an atlas.
TISSUE_COLOR = {"healthy": "#1baf7a", "adjacent_normal": "#2a78d6", "tumour": "#eb6834"}
# region-unknown tumour IS tumour tissue: same hue family, lighter tint, and always its own legend label
CONDITION_COLOR = {**TISSUE_COLOR, "tumour_region_unknown": "#f4a57f"}
CONDITION_LABEL = {
    "healthy": "Healthy donor (HLCA)",
    "adjacent_normal": "Adjacent normal (LUSC patient)",
    "tumour": "LUSC tumour",
    "tumour_region_unknown": "LUSC tumour (region not recorded)",
}
# atlas hues differ from the tissue hues so "HLCA" never reads as "healthy" on the same page
ATLAS_COLOR = {"HLCA": "#6250d6", "LuCA": "#0f9fb3"}
COMPARTMENT_COLOR = {
    "Epithelial": "#d55181",
    "Immune": "#008300",
    "Endothelial": "#c98500",
    "Stroma": "#9085e9",
    "Unknown": "#b9b8b2",
}
# direction of change (tumour vs adjacent): diverging poles, never the tissue hues
UP_COLOR, DOWN_COLOR, NS_COLOR = "#c23b3a", "#1c5cab", "#b9b8b2"
# identities that must not reuse the direction poles: malignant (CNV-called) keeps the red it has had on
# every page; basal (the trajectory root) gets a non-directional brown
MALIGNANT_COLOR, BASAL_COLOR = "#c23b3a", "#8c5a2b"
INK, INK_2 = "#0b0b0b", "#52514e"
SEQ = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
DIVERGING = [[0, "#1c5cab"], [0.5, "#f0efec"], [1, "#c23b3a"]]
SEVERITY_ICON = {"high": "🔴", "medium": "🟠", "low": "🟡"}

REPO_URL = "https://github.com/archana-b24/sc-spatial-analysis"
REPO_BRANCH = "v2-dev"  # `main` still holds the 2024 v1 notebooks; v2 lives on v2-dev until it is merged
LOG_URL = f"{REPO_URL}/blob/{REPO_BRANCH}/DEVELOPMENT_LOG.md"


def short_donor(donor_id: str) -> str:
    """Readable patient label: drop the atlas prefix and HLCA's long collection prefix.

    'HLCA:homosapiens_…_483747lambrechts_2' → 'Thienpont lambrechts_2'; 'LuCA:UKIM-V-2_P8' → 'UKIM-V-2_P8'."""
    d = str(donor_id).split(":", 1)[-1]
    if "lambrechts_" in d:
        return "Thienpont " + d[d.index("lambrechts_") :]
    if "sikkemalisa" in d:
        return "HLCA " + d.rsplit("_", 1)[-1][-12:]
    return d


# §3.13: the site has two parts. Evidence = tumour vs adjacent normal of the same patients (immune cells);
# every claim comes from there. Context = descriptive (tumour-only patients and/or healthy donors from other
# studies). "built" = the methods both parts rest on. One line under each page title names a page's part.
PART = {
    "built": (
        "How the data were built",
        "#52514e",
        "The methods behind both parts of this site: who is in the data, which cells were kept, and the "
        "shared map of cell types.",
    ),
    "evidence": (
        "Evidence · same-patient comparisons",
        "#1c5cab",
        "Tumour vs adjacent normal tissue of the same patients, immune cells. The claims of this site come "
        "from here.",
    ),
    "context": (
        "Context · descriptive",
        "#a3501f",
        "Descriptive. Uses tumour-only patients and/or healthy donors from other studies. For orientation, "
        "not evidence.",
    ),
}
STAGE_PART = {
    "cohort": "built",
    "qc": "built",
    "integration": "built",
    "de": "evidence",
    "composition": "evidence",
    "malignant": "context",
    "trajectory": "context",
    "spatial": "context",
    "ccc": "context",
}


def part_tag(part: str, text: str | None = None) -> None:
    """The coloured part tag (always with words, never colour alone) and its one-line meaning."""
    name, color, line = PART[part]
    st.markdown(
        f"<span style='border:1.5px solid {color};color:{color};border-radius:4px;padding:1px 8px;"
        f"font-size:0.85rem;font-weight:600'>{name}</span> "
        f"<span style='font-size:0.9rem;opacity:.85'>{text or line}</span>",
        unsafe_allow_html=True,
    )


def page_header(title: str, question: str, part: str | None = None, part_text: str | None = None) -> None:
    st.title(title)
    if part:
        part_tag(part, part_text)
    st.markdown(f"**The question.** {question}")


def two_parts(heading: str = "This site has two parts") -> None:
    """§3.13 cohort rule: every sample has one job. Counts are computed from the database (after QC)."""
    n = q("""WITH d AS (SELECT donor_id, any_value(study) study,
                              bool_or(condition LIKE 'tumour%') t, bool_or(condition = 'adjacent_normal') a,
                              bool_or(condition = 'healthy') h FROM cells GROUP BY 1)
             SELECT count(*) FILTER (WHERE t AND a) paired,
                    count(DISTINCT study) FILTER (WHERE t AND a) paired_studies,
                    count(*) FILTER (WHERE t AND NOT a) tumour_only,
                    count(DISTINCT study) FILTER (WHERE t AND NOT a) tumour_only_studies,
                    count(*) FILTER (WHERE h) healthy,
                    count(DISTINCT study) FILTER (WHERE h) healthy_studies FROM d""").iloc[0]
    st.subheader(heading)
    st.markdown(
        "Every sample has one job, and the job decides how much weight a result can bear. Counts are after "
        "quality control."
    )
    st.markdown(
        "| Group | Who | Its job | Part of the site |\n|---|---|---|---|\n"
        f"| LUSC patients with tumour **and** adjacent normal tissue | {int(n.paired)} patients, "
        f"{int(n.paired_studies)} studies | The patient is their own control | **Evidence**: DE & pathways, "
        "Composition |\n"
        f"| LUSC patients with tumour tissue only | {int(n.tumour_only)} patients, "
        f"{int(n.tumour_only_studies)} studies | Tumour cells for malignant-cell calling, the trajectory and "
        "the spatial reference | **Context** |\n"
        f"| Healthy lung (organ-donor, explant and autopsy lungs of people without cancer) | "
        f"{int(n.healthy)} donors, {int(n.healthy_studies)} studies | What normal lung looks like: labels, "
        "shared map, copy-number baseline, trajectory start | **Context** (reference) |"
    )
    st.caption(
        "Why only the first group counts as evidence: tumour-only patients have no control tissue of their "
        "own, and no study that gives LUSC tumours also gives healthy donors, so any tumour-vs-healthy "
        "difference is also a difference between studies. Comparing two "
        "tissues of the same patient cancels that. It measures tumour relative to the nearby lung, not "
        "relative to healthy lung (DEVELOPMENT_LOG §3.13)."
    )


def decisions(stage: str) -> None:
    """Key decisions of this stage, VISIBLE on the page (§3.12: decisions are not hidden in expanders):
    one line
    each — question → decision → why; the evidence/log reference sits in a small caption under each."""
    df = q("SELECT * FROM decisions WHERE stage = ?", (stage,))
    if df.empty:
        return
    st.subheader("Key decisions at this stage")
    for r in df.itertuples():
        with st.container(border=True):
            st.markdown(f"**{r.question}**  \n**Decision:** {r.decision}  \n**Why:** {r.why}")
            st.caption(f"Data: `{r.evidence}` · DEVELOPMENT_LOG {r.log_ref}")


def figure_legend(
    n: int | str, title: str, shown: str, read: str, sample: str, method: str, takeaway: str
) -> None:
    """Journal-style figure legend under every chart: **Figure n | finding.** What is shown · How to read it ·
    n (patients / studies / cells) · Method and statistics · Takeaway."""
    st.markdown(
        f"<div style='font-size:0.88rem;line-height:1.45;margin:-0.4rem 0 1.2rem 0;opacity:.92'>"
        f"<b>Figure {n} | {title}</b> "
        f"<b>Shown:</b> {shown} <b>How to read:</b> {read} <b>n:</b> {sample} "
        f"<b>Method:</b> {method} <b>Takeaway:</b> {takeaway}</div>",
        unsafe_allow_html=True,
    )


def methods_box(name: str, what: str, why: str, alternatives: str, cite: str) -> None:
    """A visible Methods box for a named method: what it is, why it was chosen over alternatives, citation."""
    with st.container(border=True):
        st.markdown(f"**Method · {name}**")
        a, b, c = st.columns(3, gap="medium")
        a.markdown(f"*What it is.* {what}")
        b.markdown(f"*Why this one.* {why}")
        c.markdown(f"*Alternatives considered.* {alternatives}")
        st.caption(f"Reference: {cite}")


INFERENCE_BADGE = {
    "primary": ("Primary · paired within patient", "#1c5cab"),
    "descriptive": ("Descriptive · too few patients/studies", "#52514e"),
    "exploratory": ("Exploratory · confounded by atlas and study", "#a3501f"),
    "hypothesis": ("Hypothesis-generating · confounded by atlas", "#a3501f"),
    "decomposition": ("Decomposition term · atlas + study + field effect", "#a3501f"),
    "provisional": ("Hypothesis-generating · atlas-confounded · provisional (CNV check failed)", "#a3501f"),
}


def cnv_gate_passed() -> bool:
    """The pre-registered CNV acceptance gate (cnv_summary.acceptance_passed); True if not recorded."""
    try:
        cs = q("SELECT acceptance_passed FROM cnv_summary LIMIT 1")
    except Exception:  # noqa: BLE001 — older DBs
        return True
    return cs.empty or bool(cs.acceptance_passed.iloc[0])


def inference_badge(label: str) -> None:
    """How much weight a result can bear, as a coloured tag (always paired with text, never colour alone)."""
    key = next((k for k in INFERENCE_BADGE if str(label).lower().startswith(k)), None)
    text, color = INFERENCE_BADGE.get(key, (str(label), "#52514e"))
    st.markdown(
        f"<span style='border:1.5px solid {color};color:{color};border-radius:4px;padding:1px 8px;"
        f"font-size:0.85rem;font-weight:600'>{text}</span>",
        unsafe_allow_html=True,
    )


def n_note(
    patients: int | None = None, studies: int | None = None, cells: int | None = None, extra: str = ""
) -> str:
    """'n = 15 patients · 4 studies · 12,345 cells' — put next to every result."""
    parts = []
    if patients is not None:
        parts.append(f"{int(patients)} patient{'s' if int(patients) != 1 else ''}")
    if studies is not None:
        parts.append(f"{int(studies)} stud{'ies' if int(studies) != 1 else 'y'}")
    if cells is not None:
        parts.append(f"{int(cells):,} cells")
    return "n = " + " · ".join(parts) + (f" · {extra}" if extra else "")


def confound_banner() -> None:
    """The one design fact every reader must know before any healthy-vs-tumour number (§3.12, §3.13)."""
    st.warning(
        "**Read before comparing healthy with tumour tissue.** Every healthy donor comes from one atlas "
        "(HLCA) and almost every tumour from the other (LuCA), and no study contains both. A "
        "healthy-vs-tumour difference is therefore also an atlas, study, chemistry and tissue-handling "
        "difference. Only the **tumour vs adjacent-normal comparison within the same patient** is protected "
        "from this. That comparison is the **Evidence** part of this site (DE & pathways, Composition); "
        "everything that involves healthy donors or tumour-only patients is **Context**.",
        icon="⚠️",
    )


# issues whose consequences reach beyond their own stage (run 3 review: the failed CNV gate must show wherever
# malignant calls are used)
ALSO_SHOWN_ON = {"de": [16], "trajectory": [16]}


def stage_issues(stage: str, cols: str = "*"):
    """Known issues of a stage plus those that also apply to it (ALSO_SHOWN_ON)."""
    extra = ALSO_SHOWN_ON.get(stage, [])
    ids = ",".join(str(int(i)) for i in extra) or "NULL"
    return q(f"SELECT {cols} FROM known_issues WHERE stage = ? OR id IN ({ids}) ORDER BY id", (stage,))


def fixed_in(status: str) -> str:
    """'resolved (run 3)' → 'Fixed in run 3' (never a hard-coded run number)."""
    import re

    m = re.search(r"run (\d+)", str(status))
    return f"Fixed in run {m.group(1)}" if m else "Fixed"


def caveats(stage: str) -> None:
    """Known issues for this stage, shown above the results so nobody quotes a number without them."""
    df = stage_issues(stage)
    for r in df.itertuples():
        icon = SEVERITY_ICON.get(r.severity, "⚠️")
        head = f"{icon} **Caveat ({r.severity}): {r.title}.**"
        body = f"{head} {r.detail} *(DEVELOPMENT_LOG {r.log_ref})*"
        if str(r.status).startswith("resolved"):
            st.caption(f"✅ {fixed_in(r.status)}: {r.title} (issue {r.id}, DEVELOPMENT_LOG {r.log_ref}).")
        else:
            (st.warning if r.severity == "high" else st.info)(body)


def footer() -> None:
    m = meta()
    sha = m.get("git_sha", "unknown")[:7]
    dirty = " + uncommitted changes" if m.get("git_dirty") == "True" else ""
    run = m.get("run_label") or "run not recorded"
    st.divider()
    st.caption(
        f"**Results: {run}.** Snapshot built {m.get('built_at', '?')} from commit `{sha}`{dirty} "
        f"(branch `{m.get('git_branch', '?')}`), database schema {m.get('schema_version', '?')}. "
        "Every number on this page is read from the pipeline outputs; see the Provenance page."
    )


def style(fig: go.Figure, height: int = 420, x: str | None = None, y: str | None = None) -> go.Figure:
    """Publication style: recessive grid, outside ticks, axis titles (pass x/y), legend on top, no chart
    junk."""
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=10, t=40, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, title_text=""),
        hoverlabel=dict(namelength=-1),
        font=dict(size=13, family="Helvetica, Arial, sans-serif"),
    )
    fig.update_xaxes(
        showgrid=False, zeroline=False, ticks="outside", showline=True, linecolor="rgba(128,128,128,.6)"
    )
    fig.update_yaxes(gridcolor="rgba(128,128,128,0.15)", zeroline=False, ticks="outside")
    if x is not None:
        fig.update_xaxes(title_text=x)
    if y is not None:
        fig.update_yaxes(title_text=y)
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
