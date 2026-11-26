"""Overview: the question and the study design, the key numbers with what each one means, who is in the cohort,
the atlas at a glance, and what was found (computed), with the method used for each question."""

import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import findings
from lib.ui import (
    ATLAS_COLOR,
    COMPARTMENT_COLOR,
    CONDITION_COLOR,
    CONDITION_LABEL,
    CONDITION_ORDER,
    GREY,
    UP_COLOR,
    confound_banner,
    figure_legend,
    footer,
    show_table,
    style,
    two_parts,
)

st.title("Overview")
st.markdown(
    "**The question.** When lung tissue becomes **lung squamous cell carcinoma (LUSC)**, what changes in each "
    "type of cell (the cancer cells and the immune, vascular and structural cells around them), and where do "
    "those cells sit in the tumour?"
)
st.markdown(
    "**The approach.** Public single-cell RNA-seq data from two atlases are merged into one atlas of healthy "
    "lung and LUSC. Cancer cells are identified from inferred copy-number changes. Tumour tissue is then "
    "compared with **adjacent normal tissue from the same patients**, cell type by cell type, using "
    "patients (not cells) as the statistical unit, and checked study by study. One spatial transcriptomics "
    "section places the cell types in real tissue."
)
two_parts()
confound_banner()

# ---------------------------------------------------------------- the design question
st.subheader("Is merging donors from many studies the right approach?")
xt = q("""SELECT study, any_value(atlas) atlas,
                 count(DISTINCT donor_id) FILTER (WHERE condition = 'healthy') healthy,
                 count(DISTINCT donor_id) FILTER (WHERE condition LIKE 'tumour%') tumour,
                 count(DISTINCT donor_id) FILTER (WHERE condition = 'adjacent_normal') adjacent
          FROM cells GROUP BY 1""")
pair = q("""SELECT study, count(*) n FROM (SELECT donor_id, any_value(study) study FROM cells GROUP BY 1
            HAVING bool_or(condition LIKE 'tumour%') AND bool_or(condition = 'adjacent_normal')) GROUP BY 1""")
xt["paired"] = xt.study.map(pair.set_index("study").n).fillna(0).astype(int)
both = int(((xt.healthy > 0) & (xt.tumour > 0)).sum())
st.markdown(
    f"The data come from **{len(xt)} studies**. Healthy donors come from **{int((xt.healthy > 0).sum())}** "
    f"studies (all HLCA), LUSC patients from **{int((xt.tumour > 0).sum())}**, and **{both}** studies contain "
    "both. So pooling everyone and comparing tumour with healthy would compare studies as much as tissues. The "
    "approach taken here is therefore: **pool for representation, infer within the patient.**"
)
c1, c2 = st.columns(2)
with c1, st.container(border=True, height="stretch"):
    st.markdown("**What the pooled atlas is used for**")
    st.markdown(
        "- one shared map of cell states (integration)\n"
        "- consistent cell-type labels across studies\n"
        "- the copy-number-normal reference cells for malignant calling\n"
        "- the cell-type signatures that deconvolve the tissue section"
    )
with c2, st.container(border=True, height="stretch"):
    st.markdown("**What every statistical claim rests on**")
    st.markdown(
        f"- tumour vs adjacent normal **within the same patient** ({int(xt.paired.sum())} patients from "
        f"{int((xt.paired > 0).sum())} studies)\n"
        "- replication: each study's own effect, and refits leaving each study out\n"
        "- two checks of the paired test itself: a label-swap (sign-flip) null, and what each adjacent "
        "sample contains\n"
        "- comparisons with healthy donors are **Context** (descriptive), shown next to a "
        "healthy-vs-healthy cross-study null"
    )
show_table(
    xt.sort_values(["paired", "tumour"], ascending=False).rename(columns={"paired": "paired patients"}),
    "Donors per study and tissue (after QC). Only the studies with paired patients enter the primary analysis.",
)
fig = go.Figure()
xs = xt[xt.paired > 0].sort_values("paired")
fig.add_trace(
    go.Bar(
        y=xs.study,
        x=xs.paired,
        orientation="h",
        marker_color="#1c5cab",
        text=xs.paired,
        textposition="outside",
    )
)
st.plotly_chart(style(fig, 80 + 40 * len(xs), x="patients with both tumour and adjacent tissue", y="study"))
figure_legend(
    1,
    "The primary design rests on few studies.",
    "patients with both tumour and adjacent-normal tissue, per study.",
    "every within-patient result is conditional on these studies; the DE page shows each study's own effect.",
    f"{int(xs.paired.sum())} patients, {len(xs)} studies.",
    "counted after QC; a tumour sample without a recorded region counts as tumour tissue.",
    "replication across these studies, not the size of the pooled atlas, is what supports a finding.",
)

# ---------------------------------------------------------------- key numbers
n = q("""SELECT count(*) cells, count(DISTINCT donor_id) donors, count(DISTINCT study) studies,
                count(DISTINCT sample) samples,
                count(*) FILTER (WHERE de_group = 'Malignant') malignant,
                count(*) FILTER (WHERE de_group = 'Tumour epithelial (CNV-low)') cnv_low,
                count(DISTINCT harmonised_label) FILTER (WHERE harmonised_label <> 'Unknown') labels
         FROM cells""").iloc[0]
paired = q("""SELECT count(*) n FROM (SELECT donor_id FROM cells GROUP BY donor_id
              HAVING bool_or(condition LIKE 'tumour%') AND bool_or(condition='adjacent_normal'))""").n.iloc[0]
spots = q("SELECT count(*) n FROM visium_spots").n.iloc[0]

st.subheader("Key numbers")
tiles = [
    (
        "Cells analysed",
        f"{int(n.cells):,}",
        "Individual cells whose RNA was measured, after removing damaged cells and doublets. Each cell is "
        f"assigned one of {int(n.labels)} cell types (e.g. macrophage, CD8 T cell, basal cell).",
    ),
    (
        "Donors · samples · studies",
        f"{int(n.donors)} · {int(n.samples)} · {int(n.studies)}",
        "Donors are people; a sample is one piece of tissue from one donor (some donors gave several); "
        "studies are the original publications the data come from.",
    ),
    (
        "Paired LUSC patients",
        f"{int(paired)}",
        "Patients who gave both a tumour sample and an adjacent-normal sample. Comparing within these "
        "patients removes person-to-person and study differences: the main tumour-versus-normal design.",
    ),
    (
        "Malignant cells",
        f"{int(n.malignant):,}",
        "Cells with copy-number evidence of cancer, called patient by patient against that patient's own "
        f"normal cells. Another {int(n.cnv_low):,} cells labelled as tumour by the original authors but not "
        "called are kept as a separate 'Tumour epithelial (CNV-low)' group.",
    ),
    (
        "Spatial spots",
        f"{int(spots):,}",
        "Measurement spots on one LUSC (squamous) tissue slide (10x Visium), one patient. Each 55 µm spot "
        "captures RNA from several neighbouring cells, so spots are deconvolved into cell-type mixtures.",
    ),
]
for row in (tiles[:3], tiles[3:]):
    for col, (label, value, text) in zip(st.columns(3), row, strict=False):
        with col, st.container(border=True):
            st.metric(label, value)
            st.caption(text)

# ---------------------------------------------------------------- study design
st.subheader("Study design: the four sample groups")
st.markdown(
    "Every sample belongs to one group, based on the donor's diagnosis and where the tissue was taken. The "
    "colours below are used for these groups on every page."
)
grp = q("""SELECT condition, count(DISTINCT donor_id) donors, count(DISTINCT sample) samples, count(*) cells,
                  string_agg(DISTINCT atlas, ' + ' ORDER BY atlas) atlases
           FROM cells GROUP BY 1""").set_index("condition")
MEANING = {
    "healthy": (
        "Lung tissue from people without cancer (organ-donor, explant and autopsy lungs), from the Human Lung "
        "Cell Atlas (HLCA).",
        "Reference for normal cell states (e.g. healthy basal cells start the trajectory). Every healthy "
        "comparison is Context (descriptive), because healthy donors come from different studies than the "
        "patients.",
    ),
    "adjacent_normal": (
        "Non-tumour lung taken from a LUSC patient, away from the tumour.",
        "The control in the main paired comparison. Not the same as healthy lung: it comes from a lung that "
        "also carries a tumour.",
    ),
    "tumour": (
        "Tissue from the LUSC tumour itself (mostly from the Lung Cancer Atlas, LuCA).",
        "The case in the main paired comparison; source of the malignant cells.",
    ),
    "tumour_region_unknown": (
        "Tissue from a LUSC patient where the original study did not record whether it came from the tumour "
        "or from adjacent lung.",
        "Small group (2 HLCA patients). Counted as tumour tissue in the paired tests, because the samples "
        "are from the tumour-bearing lung and LuCA's matching annotation places them in the tumour.",
    ),
}
rows = []
for c in CONDITION_ORDER:
    if c not in grp.index:
        continue
    g = grp.loc[c]
    rows.append(
        f"| <span style='color:{CONDITION_COLOR[c]}'>■</span> **{CONDITION_LABEL[c]}** | {MEANING[c][0]} | "
        f"{MEANING[c][1]} | {int(g.donors)} | {int(g.samples)} | {int(g.cells):,} | {g.atlases} |"
    )
st.markdown(
    "| Group | What it is | Role in the analysis | Donors | Samples | Cells | Source |\n"
    "|---|---|---|--:|--:|--:|---|\n" + "\n".join(rows),
    unsafe_allow_html=True,
)
st.caption(
    "Counts are after quality control. A paired patient appears in both the tumour and the adjacent-normal "
    "row, so donor counts do not add up to the total."
)

# ---------------------------------------------------------------- cohort chart
st.subheader("Who is in the cohort")
unit = st.radio(
    "Count", ["donors", "samples", "cells"], horizontal=True, key="ov_unit", label_visibility="collapsed"
)
df = q("""SELECT atlas, condition, count(DISTINCT donor_id) donors, count(DISTINCT sample) samples,
                  count(*) cells FROM cells GROUP BY 1, 2""")
order = [c for c in CONDITION_ORDER if c in set(df.condition)]
fig = go.Figure()
for atlas, color in ATLAS_COLOR.items():
    d = df[df.atlas == atlas].set_index("condition").reindex(order)
    fig.add_trace(
        go.Bar(
            y=[CONDITION_LABEL[c] for c in order],
            x=d[unit].fillna(0),
            orientation="h",
            name=f"from {atlas}",
            marker_color=color,
            hovertemplate="%{y}<br>" + atlas + ": %{x:,} " + unit + "<extra></extra>",
        )
    )
fig.update_layout(barmode="stack", legend_title_text="Source atlas")
fig.update_yaxes(autorange="reversed", title_text="Sample group")
fig.update_xaxes(title_text=f"Number of {unit}", showgrid=True, gridcolor="rgba(128,128,128,0.15)")
style(fig, 300)
fig.update_layout(legend=dict(title_text="Source atlas:  "))
st.plotly_chart(fig, width="stretch")
st.caption(
    f"**Figure 2.** Each bar is one sample group; its length is the number of {unit}, split by the atlas the data came from "
    "(violet = HLCA, teal = LuCA). **All** healthy lung comes from HLCA and almost all tumour and "
    "adjacent-normal tissue from LuCA: that is why the two atlases were merged, and why healthy-vs-tumour "
    "comparisons cannot be separated from atlas differences."
)

# ---------------------------------------------------------------- atlas map
st.subheader("The atlas at a glance")
BROAD = {
    "Malignant (CNV)": UP_COLOR,
    "Tumour epithelial (CNV-low)": "#f2a07b",
    "Epithelial (non-tumour)": COMPARTMENT_COLOR["Epithelial"],
    "Immune": COMPARTMENT_COLOR["Immune"],
    "Endothelial": COMPARTMENT_COLOR["Endothelial"],
    "Stroma": COMPARTMENT_COLOR["Stroma"],
    "Unknown": GREY,
}
view = st.radio("Colour cells by", ["Cell category", "Sample group"], horizontal=True, key="ov_umap")
cells = q("""SELECT umap_1, umap_2, condition, harmonised_label,
                    CASE WHEN de_group = 'Malignant' THEN 'Malignant (CNV)'
                         WHEN de_group = 'Tumour epithelial (CNV-low)' THEN de_group
                         WHEN compartment = 'Epithelial' THEN 'Epithelial (non-tumour)'
                         WHEN compartment = 'Unknown' THEN 'Unknown'
                         ELSE compartment END AS broad
             FROM cells WHERE in_subsample""")
col, colors, names = (
    ("broad", BROAD, {k: k for k in BROAD})
    if view == "Cell category"
    else ("condition", CONDITION_COLOR, CONDITION_LABEL)
)
fig = go.Figure()
for cat, color in colors.items():
    d = cells[cells[col] == cat]
    if d.empty:
        continue
    fig.add_trace(
        go.Scattergl(
            x=d.umap_1,
            y=d.umap_2,
            mode="markers",
            name=f"{names[cat]} ({len(d):,})",
            marker=dict(color=color, size=3, opacity=0.7),
            customdata=d[["harmonised_label"]].to_numpy(),
            hovertemplate="%{customdata[0]}<extra></extra>",
        )
    )
style(fig, 560)
fig.update_xaxes(title_text="UMAP 1", showticklabels=False)
fig.update_yaxes(title_text="UMAP 2", showticklabels=False, scaleanchor="x", showgrid=False)
fig.update_layout(legend=dict(itemsizing="constant"))
st.plotly_chart(fig, width="stretch")
st.caption(
    f"**Figure 3.** Each dot is one cell ({len(cells):,} shown, a stratified subsample of all {int(n.cells):,}). Cells with "
    "similar gene expression sit close together, so cells of the same type tend to cluster. The axes (UMAP 1 and 2) have no "
    "units: only closeness matters, not the distance between far-apart clusters. 'Unknown' cells had no "
    "confident cell-type label in the source atlases. Click a legend entry to hide or show it."
)

# ---------------------------------------------------------------- findings
st.subheader("What was found")
st.caption(
    "One line per question, computed from the data each time the page loads. Open the stage for the evidence "
    "and caveats."
)
QUESTIONS = [
    (
        "malignant",
        "views/malignant.py",
        "Which cells are cancer cells?",
        "inferCNV per patient, two-axis call against a same-patient null",
    ),
    (
        "de",
        "views/de.py",
        "Which genes and pathways change in tumours, per cell type?",
        "Pseudobulk PyDESeq2, paired; decoupler",
    ),
    (
        "composition",
        "views/composition.py",
        "Do tumours contain different proportions of cell types?",
        "paired NB GLM within the immune compartment, Milo",
    ),
    (
        "trajectory",
        "views/trajectory.py",
        "How far have cancer cells moved from healthy basal cells?",
        "Palantir, CellRank 2",
    ),
    (
        "spatial",
        "views/spatial.py",
        "Where are the cell types in tumour tissue?",
        "cell2location, CellCharter",
    ),
    ("ccc", "views/ccc.py", "Which cell-to-cell signals change in tumours?", "LIANA+, Tensor-cell2cell"),
]
for stage, page, question, method in QUESTIONS:
    f = findings(stage)
    with st.container(border=True):
        a, b = st.columns([2, 5], gap="medium")
        with a:
            st.page_link(page, label=f"**{question}**")
            st.caption(f"Method: {method}")
        with b:
            st.markdown("\n".join(f"- {x}" for x in f[:2]) if f else "_No result available._")
st.caption(
    "Earlier stages (cohort, quality control, integration) build the atlas these answers rest on. They are "
    "explained on the Start here page and on their own pages."
)
footer()
