"""Cohort & merge: how the donors were assembled from two atlases, who they are, and how the merge
was made safe. Every number is computed from the DB; thresholds come from the config stored in `meta`."""

import plotly.graph_objects as go
import streamlit as st
from lib.db import config, q
from lib.explain import guide
from lib.ui import (
    ATLAS_COLOR,
    CONDITION_COLOR,
    CONDITION_LABEL,
    CONDITION_ORDER,
    GREY,
    caveats,
    decisions,
    footer,
    page_header,
    show_table,
    style,
)

page_header(
    "Cohort & merge",
    "Which patients and cells can fairly compare LUSC tumour with normal lung, and how do two independently "
    "built public atlases become one dataset without counting anyone twice?",
    part="built",
)
guide("cohort")
caveats("cohort")

cfg = config()
coh = cfg.get("data", {}).get("cohort", {})
min_age = coh.get("healthy_min_age", "?")
cap_h = coh.get("healthy_max_cells_per_donor", "?")
cap_d = coh.get("disease_max_cells_per_donor", "?")
min_cells = cfg.get("qc", {}).get("min_cells_per_sample", "?")

donors = q("SELECT * FROM cohort_donors")
cw = q("SELECT * FROM merge_study_crosswalk")
ex = q("SELECT * FROM merge_exclusions")
gs = q("SELECT * FROM merge_gene_space").set_index("atlas")
post = q("SELECT count(DISTINCT donor_id) d, count(*) c, count(DISTINCT study) s FROM cells").iloc[0]
lost = q("""SELECT condition, count(DISTINCT donor_id) n FROM cohort_donors
            WHERE donor_id NOT IN (SELECT DISTINCT donor_id FROM cells) GROUP BY 1""")
n_lost_h = int(lost.loc[lost.condition == "healthy", "n"].sum())
n_lost_l = int(lost.loc[lost.condition != "healthy", "n"].sum())


def nd(mask) -> int:
    return donors.loc[mask, "donor_id"].nunique()


healthy = donors.condition == "healthy"
n_h, n_l = nd(healthy), nd(~healthy)
n_l_hlca, n_l_luca = nd(~healthy & (donors.atlas == "HLCA")), nd(~healthy & (donors.atlas == "LuCA"))
studies = {a: cw.loc[cw.atlas == a, "study"].nunique() for a in ("HLCA", "LuCA")}
lusc_studies_hlca = cw.loc[(cw.atlas == "HLCA") & (cw.n_healthy < cw.n_cells), "study"].tolist()

# LUSC patients: which samples did each give?
pt = (
    donors[~healthy]
    .pivot_table(index=["donor_id", "atlas", "study"], columns="condition", values="n_cells", aggfunc="sum")
    .reset_index()
)
for c in ("tumour", "adjacent_normal", "tumour_region_unknown"):
    if c not in pt:
        pt[c] = float("nan")
has_t, has_a, has_u = pt.tumour.notna(), pt.adjacent_normal.notna(), pt.tumour_region_unknown.notna()
# §3.12: a tumour sample whose region was not recorded is still tumour tissue → it can pair with adjacent tissue
n_paired = int(((has_t | has_u) & has_a).sum())
n_t_only = int((has_t & ~has_a).sum())
n_u_only = int((has_u & ~has_t & ~has_a).sum())

# ---------------------------------------------------------------- cohort in one picture
st.subheader("The cohort in one picture")
st.markdown(
    f"**{donors.donor_id.nunique()} donors were selected**: **{n_h} healthy donors** and **{n_l} LUSC "
    f"patients**. Quality control then removed {n_lost_h + n_lost_l} donors whose samples were too small "
    f"({n_lost_h} healthy, {n_lost_l} LUSC), leaving **{int(post.d)} donors after QC** ({n_h - n_lost_h} "
    "healthy). Counts on other pages are after QC unless they say 'selected'. The diagram shows where each "
    "group comes from."
)
st.graphviz_chart(
    f"""digraph {{
    rankdir=LR; nodesep=0.3; ranksep=0.45;
    node [shape=box, style="rounded", fontname="Helvetica", fontsize=10, margin="0.14,0.07"];
    edge [color="#888", fontname="Helvetica", fontsize=9, fontcolor="#888"];
    hlca [label="HLCA\\nHuman Lung Cell Atlas\\n{studies["HLCA"]} studies kept", color="{ATLAS_COLOR["HLCA"]}", penwidth=2];
    luca [label="LuCA\\nLung Cancer Atlas\\n{studies["LuCA"]} studies kept", color="{ATLAS_COLOR["LuCA"]}", penwidth=2];
    h [label="{n_h} healthy donors selected\\nlung tissue, no cancer, aged {min_age}+", color="{CONDITION_COLOR["healthy"]}", penwidth=2];
    l [label="{n_l} LUSC patients\\n{n_paired} paired (tumour + adjacent)\\n{n_t_only} tumour only · {n_u_only} region unknown",
       color="{CONDITION_COLOR["tumour"]}", penwidth=2];
    sel [label="Selected cohort\\n{donors.donor_id.nunique()} donors\\n{int(donors.n_cells.sum()):,} cells"];
    qc [label="After quality control\\n{int(post.d)} donors · {int(post.s)} studies\\n{int(post.c):,} cells"];
    hlca -> h [label=" {n_h}"]; hlca -> l [label=" {n_l_hlca}"]; luca -> l [label=" {n_l_luca}"];
    h -> sel; l -> sel;
    sel -> qc [label=" −{n_lost_h + n_lost_l} donors ({n_lost_h} healthy, {n_lost_l} LUSC):\\n every sample < {min_cells} cells"];
    }}""",
    width="stretch",
)

# ---------------------------------------------------------------- why two atlases
st.subheader("Why two atlases?")
a, b, c = st.columns(3, gap="medium")
with a, st.container(border=True, height="stretch"):
    st.markdown(
        f"<span style='color:{ATLAS_COLOR['HLCA']}'>■</span> **HLCA · Human Lung Cell Atlas**",
        unsafe_allow_html=True,
    )
    st.caption("Sikkema et al. 2023, *Nature Medicine*")
    st.markdown(
        "A curated integration of many single-cell studies, built as the reference for **healthy human "
        f"lung**. Here it supplies **{n_h} healthy donors** (selected; {n_h - n_lost_h} after QC) from "
        f"{studies['HLCA'] - len(lusc_studies_hlca)} studies. It contains only **{n_l_hlca} LUSC patients**, all from one study "
        f"({', '.join(lusc_studies_hlca)}), far too few for a tumour comparison and fully confounded with "
        "that study."
    )
with b, st.container(border=True, height="stretch"):
    st.markdown(
        f"<span style='color:{ATLAS_COLOR['LuCA']}'>■</span> **LuCA · Lung Cancer Atlas**",
        unsafe_allow_html=True,
    )
    st.caption("Salcher et al. 2022, *Cancer Cell*")
    st.markdown(
        "A single-cell atlas of **lung cancer** (several subtypes) with tumour and adjacent tissue. Only its "
        f"LUSC samples are used: **{n_l_luca} LUSC patients** from {studies['LuCA']} studies, "
        f"{nd(~healthy & (donors.atlas == 'LuCA') & (donors.condition == 'adjacent_normal'))} of them with "
        "adjacent-normal tissue from the same person. Its healthy-lung samples are copies of studies already "
        "in HLCA, so they were removed as duplicates."
    )
with c, st.container(border=True, height="stretch"):
    st.markdown("**Together**")
    st.caption("one cohort, used in two ways")
    st.markdown(
        f"**{n_l} LUSC patients** ({n_paired} with both tumour and adjacent tissue) and **{n_h} healthy donors** (selected) across "
        f"{studies['HLCA'] + studies['LuCA']} studies. Because healthy donors come from one atlas and tumours "
        "from the other, healthy-versus-tumour differences are partly atlas differences. So the **main "
        f"comparison is within the {n_paired} paired patients**, and healthy lung is used as a reference "
        "(normal cell states, copy-number baseline, trajectory start)."
    )

# ---------------------------------------------------------------- the merge
st.subheader("How the two atlases were merged: six problems and what was done")
st.markdown(
    "Merging two curated atlases is not a simple concatenation. Each is itself an integration of many studies "
    "with its own processing, gene annotation and cell-type names. These are the problems that would bias the "
    "results if ignored, and the fix for each."
)
dup = ex[ex.reason.str.startswith("not primary")]
dup_cells = dup.groupby("atlas").n_cells.sum()
lusc_dup = dup[dup.disease.str.contains("carcinoma")]
smart = ex[ex.reason.str.startswith("non-UMI")]
young = ex[ex.reason.str.startswith("healthy donor age")]
nonlung = ex[ex.reason.str.startswith("tissue not")]
notlung_h = ex[ex.reason.str.startswith("healthy sample not lung tissue")]
n_labels = int(
    q("SELECT count(DISTINCT harmonised_label) n FROM cells WHERE harmonised_label <> 'Unknown'").n.iloc[0]
)
n_unknown = int(q("SELECT count(*) n FROM cells WHERE harmonised_label = 'Unknown'").n.iloc[0])
avail = int(donors.n_cells_available.sum())

PROBLEMS = [
    (
        "1 · The same patients appear in both atlases",
        "LuCA re-distributes published studies that HLCA also contains. Patient IDs differ between atlases, so "
        "matching IDs finds nothing. Duplicates would count a patient twice (inflating the statistics) and "
        "make integration look better than it is, since identical cells mix perfectly.",
        "Keep only cells CELLxGENE marks as primary data (the original copy), then build a study-by-study "
        "crosswalk and check that no study contributes donors from both atlases.",
        f"Non-primary (duplicate) cells removed: {int(dup_cells.get('HLCA', 0)):,} HLCA and {int(dup_cells.get('LuCA', 0)):,} LuCA "
        "cells. Example: LuCA's copy of "
        + (
            f"{lusc_dup.study.iloc[0]} ({int(lusc_dup.n_cells.iloc[0]):,} LUSC cells)"
            if len(lusc_dup)
            else "a LUSC study"
        )
        + f" is the same data as HLCA's {', '.join(lusc_studies_hlca)}. After the fix, "
        f"{int(cw.overlaps_with.notna().sum())} studies appear in both atlases.",
    ),
    (
        "2 · Measurements must be the same kind of number",
        "Most studies count unique molecules (UMIs); full-length methods such as Smart-seq2 count reads, "
        "which scale with gene length. The count models used later (scVI, DESeq2, cell2location) assume "
        "integer UMI counts.",
        "Drop non-UMI assays, and check per study that the matrix is non-negative integer raw counts.",
        f"Removed: {', '.join(smart.study)} ({int(smart.n_cells.sum()):,} cells, Smart-seq2). "
        f"Raw-count check passed for {int(cw.counts_ok.sum())} of {len(cw)} studies.",
    ),
    (
        "3 · Genes are annotated differently",
        "Gene symbols get renamed between annotation releases, and an atlas may simply not report a gene, "
        "which must not be read as zero expression.",
        "Join on Ensembl gene IDs (version suffix removed), never symbols, and keep only genes measured in "
        "both atlases.",
        f"HLCA reports {int(gs.loc['HLCA', 'n_genes']):,} genes and LuCA {int(gs.loc['LuCA', 'n_genes']):,}; "
        f"**{int(gs.loc['HLCA', 'n_shared']):,} genes** are shared and used. LuCA's smaller gene list is the "
        "limit.",
    ),
    (
        "4 · The two atlases name cell types differently",
        "HLCA and LuCA built their cell-type labels independently, and LuCA has tumour-specific states that "
        "HLCA lacks. Comparing their fine labels directly is meaningless.",
        "A reviewed crosswalk (config/label_crosswalk.csv) maps both atlases' Cell Ontology terms to one "
        "shared set of cell types in four compartments. Cells with no confident mapping are labelled "
        "Unknown; fine labels are kept per atlas but never compared across atlases.",
        f"**{n_labels} shared cell types**. {n_unknown:,} cells are Unknown after QC, many of them LuCA "
        "tumour cells (stage 4 separates those out).",
    ),
    (
        "5 · 'Normal' tissue means different things",
        "Adjacent tissue from a cancer patient is not healthy lung, and some samples labelled 'healthy' are not "
        "healthy lung tissue either (airway lavage, or tumour-free tissue from cancer surgery). Many healthy "
        "HLCA donors are also much younger than LUSC patients. Lumping these together would mix up cancer, "
        "tissue, smoking and age effects.",
        f"Four sample groups instead of tumour/normal (see Overview). 'Healthy' means lung tissue from a "
        "person without cancer: organ-donor, explant or autopsy lungs only (an allow-list on HLCA's own "
        f"metadata). Healthy donors restricted to age {min_age}+ (donors without a recorded age are kept). "
        "Only lung tissue: nose, airway, lymph-node and metastasis samples are dropped.",
        f"Removed: {int(notlung_h.n_cells.sum()):,} 'healthy' cells that are lavage or cancer-surgery tissue "
        f"({', '.join(sorted(notlung_h.study.unique())) or 'none'}), {int(young.n_cells.sum()):,} cells from "
        f"healthy donors under {min_age} ({young.study.nunique()} studies) and {int(nonlung.n_cells.sum()):,} "
        "cells from non-lung tissue.",
    ),
    (
        "6 · A huge healthy reference would drown out the tumours",
        "HLCA holds hundreds of thousands of healthy cells. Taken in full, healthy lung would dominate the "
        "integration, the cell-type map and the computing cost.",
        f"Cap the cells taken per donor: at most {cap_h:,} per healthy donor and {cap_d:,} per LUSC patient, "
        "sampled at random within each donor and stratified by cell type, so each donor's cell-type mix is kept.",
        f"{avail:,} cells were available after the filters above; **{int(donors.n_cells.sum()):,}** were "
        f"selected, of which {int(donors.loc[~healthy, 'n_cells'].sum()) / int(donors.n_cells.sum()):.0%} come "
        "from LUSC patients.",
    ),
]
for i in range(0, len(PROBLEMS), 2):
    cols = st.columns(2, gap="medium")
    for col, (title, problem, fix, evidence) in zip(cols, PROBLEMS[i : i + 2], strict=False):
        with col, st.container(border=True, height="stretch"):
            st.markdown(f"**{title}**")
            st.markdown(f"*Problem.* {problem}")
            st.markdown(f"*What was done.* {fix}")
            st.caption(f"Data: {evidence}")
st.info(
    "**What the merge does not fix.** Study and technology differences (batch effects) remain in the data. "
    "They are handled in stage 3 (integration, with each sample as a batch) and, for statistics, by comparing "
    "tumour with adjacent tissue within the same patients."
)
decisions("cohort")

# ---------------------------------------------------------------- who the donors are
st.subheader("Who the donors are")
tab_s, tab_p, tab_a, tab_x, tab_t = st.tabs(
    ["Studies", "LUSC patients' samples", "Age", "Sex", "Sequencing technology"]
)

with tab_s:
    by = donors.groupby(["study", "atlas", "condition"]).donor_id.nunique().reset_index(name="donors")
    tot = by.groupby("study").donors.max().sort_values()
    ylab = {s: f"{s} ({by.loc[by.study == s, 'atlas'].iloc[0]})" for s in tot.index}
    fig = go.Figure()
    for cond in CONDITION_ORDER:
        d = by[by.condition == cond].set_index("study").reindex(tot.index)
        fig.add_trace(
            go.Bar(
                y=[ylab[s] for s in tot.index],
                x=d.donors.fillna(0),
                orientation="h",
                name=CONDITION_LABEL[cond],
                marker_color=CONDITION_COLOR[cond],
                hovertemplate="%{y}<br>" + CONDITION_LABEL[cond] + ": %{x} donors<extra></extra>",
            )
        )
    fig.update_layout(barmode="group", bargap=0.2)
    fig.update_xaxes(title_text="Number of donors", showgrid=True, gridcolor="rgba(128,128,128,0.15)")
    fig.update_yaxes(title_text="Study (atlas)")
    st.plotly_chart(style(fig, 620), width="stretch")
    st.caption(
        "**Figure 1.** One row per original study, with the atlas it was taken from. Bars count donors per sample group, so "
        "a paired patient is counted once under tumour and once under adjacent normal. Healthy donors come "
        "from many HLCA studies; LUSC patients from a few LuCA studies."
    )

with tab_p:
    m = pt.sort_values(["atlas", "study", "tumour"], ascending=[True, True, False]).reset_index(drop=True)
    m["label"] = m.study + " · patient " + (m.groupby("study").cumcount() + 1).astype(str)
    cols_ = ["tumour", "adjacent_normal", "tumour_region_unknown"]
    fig = go.Figure()
    for cond in cols_:
        d = m[m[cond].notna()]
        fig.add_trace(
            go.Scatter(
                x=[CONDITION_LABEL[cond]] * len(d),
                y=d.label,
                mode="markers",
                name=CONDITION_LABEL[cond],
                marker=dict(
                    color=CONDITION_COLOR[cond],
                    size=(d[cond] ** 0.5 / (m[cols_].max().max() ** 0.5) * 22 + 5),
                    line=dict(width=0),
                ),
                customdata=d[cond],
                hovertemplate="%{y}<br>%{x}: %{customdata:,} cells<extra></extra>",
            )
        )
    fig.update_yaxes(
        autorange="reversed", title_text="LUSC patient (study)", categoryorder="array", categoryarray=m.label
    )
    fig.update_xaxes(title_text="Sample type given", side="top")
    style(fig, max(420, 17 * len(m)))
    fig.update_layout(showlegend=False)
    st.plotly_chart(fig, width="stretch")
    st.caption(
        f"**Figure 2.** One row per LUSC patient ({len(m)}); a dot means the patient gave that sample type, and its size is "
        f"the number of cells selected. **{n_paired} patients have dots in both the tumour and adjacent-normal "
        f"columns: these are the paired patients.** {n_t_only} gave tumour only, and {n_u_only} gave tissue "
        "whose region was not recorded."
    )

with tab_a:
    ag = donors.dropna(subset=["age_years"]).drop_duplicates("donor_id")
    fig = go.Figure()
    for cond in ["healthy", "adjacent_normal", "tumour"]:
        d = ag[ag.condition == cond]
        if d.empty:
            continue
        fig.add_trace(
            go.Box(
                y=[CONDITION_LABEL[cond]] * len(d),
                x=d.age_years,
                name=CONDITION_LABEL[cond],
                orientation="h",
                boxpoints="all",
                jitter=0.4,
                pointpos=0,
                marker=dict(color=CONDITION_COLOR[cond], size=6, opacity=0.7),
                line=dict(color=CONDITION_COLOR[cond]),
                fillcolor="rgba(0,0,0,0)",
                hovertemplate="%{x} years<extra></extra>",
            )
        )
    fig.add_vline(
        x=min_age,
        line_dash="dash",
        line_color=GREY,
        annotation_text=f"healthy cut-off: {min_age}",
        annotation_position="top",
    )
    fig.update_xaxes(title_text="Age at sampling (years)", showgrid=True, gridcolor="rgba(128,128,128,0.15)")
    fig.update_yaxes(title_text="Sample group")
    fig.update_layout(showlegend=False)
    st.plotly_chart(style(fig, 320), width="stretch")
    known = ag.groupby("condition").donor_id.nunique()
    AGE_NOUN = {
        "healthy": "healthy donors",
        "adjacent_normal": "patients with adjacent tissue",
        "tumour": "patients with tumour tissue",
    }
    total = donors.groupby("condition").donor_id.nunique()
    st.caption(
        "**Figure 3.** Each dot is one donor; the box spans the middle half of ages and the line inside is the median. Age is "
        "recorded for "
        + ", ".join(
            f"{int(known.get(c, 0))} of {int(total[c])} {AGE_NOUN[c]}"
            for c in ["healthy", "adjacent_normal", "tumour"]
        )
        + " (selected donors, before QC). Healthy donors younger than the cut-off were excluded so the reference is closer in age to "
        "the patients; healthy donors with no recorded age are kept."
    )

with tab_x:
    sx = q("""SELECT condition, sex, count(DISTINCT donor_id) n FROM
              (SELECT donor_id, any_value(sex) sex, condition FROM cells GROUP BY donor_id, condition) GROUP BY 1, 2""")
    sex_color = {"female": "#7a3db8", "male": "#b0651a", "unknown": GREY}  # no tissue/direction hues
    order = [c for c in CONDITION_ORDER if c in set(sx.condition)]
    fig = go.Figure()
    for sex, color in sex_color.items():
        d = sx[sx.sex == sex].set_index("condition").reindex(order)
        fig.add_trace(
            go.Bar(
                y=[CONDITION_LABEL[c] for c in order],
                x=d.n.fillna(0),
                orientation="h",
                name=sex,
                marker_color=color,
                hovertemplate="%{y}<br>" + sex + ": %{x} donors<extra></extra>",
            )
        )
    fig.update_layout(barmode="stack", legend_title_text="Sex")
    fig.update_yaxes(autorange="reversed", title_text="Sample group")
    fig.update_xaxes(
        title_text="Number of donors (after QC)", showgrid=True, gridcolor="rgba(128,128,128,0.15)"
    )
    st.plotly_chart(style(fig, 300), width="stretch")
    st.caption(
        "**Figure 4.** Donors per sample group, split by recorded sex. LUSC patients are mostly male, as in LUSC "
        "epidemiology; many healthy HLCA donors have no recorded sex. Sex is not a covariate in the "
        "paired comparison because each patient is compared with themselves."
    )

with tab_t:
    tc = q(
        """SELECT assay, atlas, count(DISTINCT donor_id) donors, count(*) cells FROM cells GROUP BY 1, 2"""
    )
    tot = tc.groupby("assay").donors.sum().sort_values()
    fig = go.Figure()
    for atlas, color in ATLAS_COLOR.items():
        d = tc[tc.atlas == atlas].set_index("assay").reindex(tot.index)
        fig.add_trace(
            go.Bar(
                y=tot.index,
                x=d.donors.fillna(0),
                orientation="h",
                name=f"from {atlas}",
                marker_color=color,
                customdata=d.cells.fillna(0),
                hovertemplate="%{y}<br>" + atlas + ": %{x} donors, %{customdata:,} cells<extra></extra>",
            )
        )
    fig.update_layout(barmode="stack", legend_title_text="Source atlas")
    fig.update_xaxes(
        title_text="Number of donors (after QC)", showgrid=True, gridcolor="rgba(128,128,128,0.15)"
    )
    fig.update_yaxes(title_text="Sequencing technology")
    st.plotly_chart(style(fig, 380), width="stretch")
    st.caption(
        "**Figure 5.** Which single-cell technology each donor's cells were measured with (a donor can appear under two "
        "if its samples used different chemistries). Many technologies are mixed here, which is why "
        "integration (stage 3) treats every sample as its own batch."
    )

# ---------------------------------------------------------------- what was removed
st.subheader("What was removed, and why")
PLAIN = {
    "not primary data": "Duplicate of a study held elsewhere",
    "tissue not primary lung": "Not lung tissue (nose, airway, lymph node …)",
    "healthy donor age unknown": "Healthy donor with no recorded age (policy: exclude)",
    "healthy donor age": f"Healthy donor younger than {min_age}",
    "non-UMI assay": "Smart-seq2 (read counts, not UMIs)",
    "healthy sample not lung tissue": "'Healthy' sample that is lavage or cancer-surgery tissue",
}
ex["why"] = ex.reason.map(lambda r: next((v for k, v in PLAIN.items() if r.startswith(k)), r))
agg = (
    ex.groupby(["why", "atlas"]).agg(n_cells=("n_cells", "sum"), n_studies=("study", "nunique")).reset_index()
)
order = agg.groupby("why").n_cells.sum().sort_values().index
fig = go.Figure()
for atlas, color in ATLAS_COLOR.items():
    d = agg[agg.atlas == atlas].set_index("why").reindex(order)
    fig.add_trace(
        go.Bar(
            y=list(order),
            x=d.n_cells.fillna(0),
            orientation="h",
            name=f"from {atlas}",
            marker_color=color,
            customdata=d.n_studies.fillna(0),
            hovertemplate="%{y}<br>" + atlas + ": %{x:,} cells from %{customdata} studies<extra></extra>",
        )
    )
fig.update_layout(barmode="stack", legend_title_text="Source atlas")
fig.update_xaxes(title_text="Cells removed", showgrid=True, gridcolor="rgba(128,128,128,0.15)")
fig.update_yaxes(title_text="Reason")
st.plotly_chart(style(fig, 300), width="stretch")
st.caption(
    "**Figure 6.** Cells removed from the healthy and LUSC samples of each atlas before selection, by reason. Samples of "
    "other diseases (e.g. lung adenocarcinoma, fibrosis) were never selected and are not counted here."
)
show_table(ex.drop(columns="why").sort_values("n_cells", ascending=False), "Every exclusion, per study")

# ---------------------------------------------------------------- crosswalk
st.subheader("Study crosswalk: every study, taken from exactly one atlas")
cwo = cw.sort_values("n_cells")
fig = go.Figure()
for cond in CONDITION_ORDER:
    fig.add_trace(
        go.Bar(
            y=cwo.study + " (" + cwo.atlas + ")",
            x=cwo[f"n_{cond}"],
            orientation="h",
            name=CONDITION_LABEL[cond],
            marker_color=CONDITION_COLOR[cond],
            hovertemplate="%{y}<br>" + CONDITION_LABEL[cond] + ": %{x:,} cells<extra></extra>",
        )
    )
fig.update_layout(barmode="stack")
fig.update_xaxes(
    title_text="Cells available (after filters, before the per-donor cap)",
    showgrid=True,
    gridcolor="rgba(128,128,128,0.15)",
)
fig.update_yaxes(title_text="Study (atlas)")
st.plotly_chart(style(fig, 560), width="stretch")
st.caption(
    "**Figure 7.** One bar per study that was kept, coloured by sample group. The table also shows the overlap check "
    "(`overlaps_with`: empty means the study is not in the other atlas) and the raw-count check (`counts_ok`)."
)
show_table(cw.sort_values(["atlas", "n_cells"], ascending=[True, False]), "Study crosswalk")
show_table(
    donors[donors.donor_id.isin(pt.loc[(has_t | has_u) & has_a, "donor_id"])].sort_values(
        ["donor_id", "condition"]
    ),
    "Paired patients: one row per patient and sample type",
)
footer()
