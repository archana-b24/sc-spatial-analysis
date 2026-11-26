"""Composition: paired, within-compartment tests of immune cell-type shares (NB GLM + propeller-style +
CLR), per-patient slopes, per-study direction, Milo neighbourhoods; healthy shown descriptively only
(DEVELOPMENT_LOG §3.12)."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import q
from lib.explain import guide
from lib.ui import (
    DOWN_COLOR,
    GREY,
    INK,
    INK_2,
    NS_COLOR,
    UP_COLOR,
    caveats,
    decisions,
    figure_legend,
    footer,
    inference_badge,
    methods_box,
    n_note,
    page_header,
    short_donor,
    show_table,
    style,
)

page_header(
    "Composition",
    "Do LUSC tumours contain different proportions of immune cell types than adjacent normal lung of the same "
    "patients?",
    part="evidence",
)
guide("composition")
caveats("composition")

st.subheader("Why proportions need care here")
st.markdown(
    "Two things make cell-type proportions easy to misread in a merged atlas. **They are compositional:** "
    "shares add up to 100 %, so if tumour cells expand, every other share falls even when nothing else changed. "
    "**They depend on how each study prepared its cells:** some studies sorted for immune cells (CD45+), one "
    "for endothelial cells (CD31+), others sequenced everything, so the share of a cell type among *all* cells "
    "mostly reflects the protocol. Both problems shrink when we (1) compare tumour with adjacent tissue of the "
    "**same patient** (same protocol), and (2) express each immune type as a share of **immune cells only**."
)
st.page_link(
    "views/de.py",
    label="What each patient's adjacent-normal (control) sample actually contains: DE & pathways page, "
    "'Two checks of the paired comparison itself'",
    icon="🔍",
)
methods_box(
    "three paired models that must agree",
    "(1) Negative-binomial GLM on counts: count ~ patient + tissue, offset = log(immune cells in that "
    "patient-tissue); (2) propeller-style: arcsin-√ of each share, paired t-test; (3) centred log-ratio (CLR, "
    "pseudocount 0.5), paired t-test. BH across cell types within each.",
    "The GLM is the edgeR-style model behind Milo and propeller and handles counting noise; the two transformed "
    "tests are simple and assumption-light. Requiring agreement guards against one model's quirks.",
    "scCODA (Bayesian, used in runs 1–2: with 14 patient covariates its sparsity prior shrinks patient effects); "
    "sccomp / propeller (R); raw-proportion tests (ignore compositionality).",
    "Phipson et al. 2022 Bioinformatics (propeller); Büttner et al. 2021 Nat Commun (scCODA); Aitchison 1982 "
    "(log-ratios); Dann et al. 2022 Nat Biotechnol (Milo).",
)
decisions("composition")

c = q("SELECT * FROM composition_paired ORDER BY glm_log2fc")
units = q("SELECT * FROM composition_units")
if c.empty:
    st.info("No paired composition results in this snapshot.")
    footer()
    st.stop()
nd, ns = int(c.n_donors.iloc[0]), int(c.n_studies.iloc[0])
inference_badge("primary")
st.caption(
    n_note(nd, ns) + " · immune compartment · units with ≥ 200 immune cells, kept only with their partner."
)

# ------------------------------------------------------------------------------------------------ fig 1: forest
st.subheader("Which immune cell types change?")
col = {"higher in tumour": UP_COLOR, "lower in tumour": DOWN_COLOR}
fig = go.Figure(
    go.Scatter(
        x=c.glm_log2fc,
        y=c.cell_type,
        mode="markers",
        marker=dict(
            size=11,
            color=[col.get(v, GREY) for v in c.consistent],
            symbol=["circle" if v in col else "circle-open" for v in c.consistent],
        ),
        error_x=dict(
            type="data",
            symmetric=False,
            array=c.glm_ci_high - c.glm_log2fc,
            arrayminus=c.glm_log2fc - c.glm_ci_low,
            color="#8a8984",
        ),
        customdata=np.stack(
            [
                c.glm_padj,
                c.propeller_padj,
                c.clr_padj,
                c.n_donors_higher_in_tumour,
                c.k_studies_same_direction,
                c.k_studies,
            ],
            axis=1,
        ),
        hovertemplate="%{y}<br>log2 change %{x:.2f}<br>padj GLM %{customdata[0]:.2g} · propeller %{customdata[1]:.2g}"
        " · CLR %{customdata[2]:.2g}<br>higher in tumour in %{customdata[3]} patients"
        "<br>same direction in %{customdata[4]}/%{customdata[5]} studies<extra></extra>",
    )
)
fig.add_vline(x=0, line_color=GREY)
st.plotly_chart(
    style(fig, 90 + 30 * len(c), x="log2 change in share of immune cells, tumour vs adjacent (95 % CI)")
)
_dn = ", ".join(c[c.consistent == "lower in tumour"].sort_values("glm_log2fc").cell_type)
_up = ", ".join(c[c.consistent == "higher in tumour"].sort_values("glm_log2fc", ascending=False).cell_type)
_imm = units[(units.compartment == "Immune") & units.tissue_type.isin(["tumour", "adjacent_normal"])]
_tot = _imm.groupby(["donor_id", "tissue_type"]).n_cells.sum().unstack()
_pairs = _tot.dropna()
_short = _pairs[(_pairs < 200).any(axis=1)]
_drop_note = (
    "; "
    + ", ".join(
        f"{d.split('_')[-1]} ({int(r.min())} immune cells in one tissue)" for d, r in _short.iterrows()
    )
    + " left out (< 200 immune cells in a tissue)"
    if len(_short)
    else ""
)
figure_legend(
    1,
    (
        f"Tumours hold fewer {_dn or '—'} and more {_up or '—'} than adjacent lung of the same patients."
        if _dn or _up
        else "No immune cell type changes consistently between tumour and adjacent lung."
    ),
    "estimate and 95 % confidence interval from the paired negative-binomial model for each immune cell type.",
    "right = larger share of immune cells in tumour. Filled red/blue = all three models agree (BH padj < 0.05, "
    "same direction); open grey = not consistent across models.",
    n_note(nd, ns)
    + "; "
    + f"{int(units[units.compartment == 'Immune'].n_cells.sum()):,} immune cells in paired units"
    + _drop_note
    + ".",
    "NB GLM count ~ patient + tissue with offset log(immune total); cross-checked by arcsin-√ and CLR paired "
    "t-tests; hover for all three adjusted p-values and per-study direction.",
    "only filled markers are claims; hover to see whether the change holds in every study.",
)

# ------------------------------------------------------------------------------------------------ fig 2: slopes
st.subheader("Patient by patient")
pick = st.selectbox(
    "Cell type",
    c.cell_type.tolist(),
    index=int(np.argmax(np.abs(c.glm_log2fc.to_numpy()) * (c.consistent != "not consistent across models"))),
)
im = units[(units.compartment == "Immune") & units.tissue_type.isin(["tumour", "adjacent_normal"])]
tot = im.groupby(["unit", "donor_id", "study", "tissue_type"]).n_cells.sum().rename("tot").reset_index()
tot = tot[tot.tot >= 200]
both = tot.groupby("donor_id").tissue_type.nunique()
tot = tot[tot.donor_id.isin(both[both == 2].index)]
x = im[im.cell_type == pick].groupby("unit").n_cells.sum()
tot["share"] = 100 * tot.unit.map(x).fillna(0) / tot.tot
fig = go.Figure()
for dn, dd in tot.groupby("donor_id"):
    dd = dd.set_index("tissue_type").reindex(["adjacent_normal", "tumour"])
    up = dd.share.iloc[1] > dd.share.iloc[0]
    fig.add_trace(
        go.Scatter(
            x=["Adjacent normal", "LUSC tumour"],
            y=dd.share,
            mode="lines+markers",
            line=dict(color=UP_COLOR if up else DOWN_COLOR, width=1.5),
            marker=dict(size=7),
            name=short_donor(dn),
            showlegend=False,
            hovertemplate=f"{short_donor(dn)}<br>{dd.study.iloc[0]}<br>%{{y:.1f}}%<extra></extra>",
        )
    )
st.plotly_chart(style(fig, 380, x="tissue", y=f"{pick}: % of immune cells"))
row = c.set_index("cell_type").loc[pick]
figure_legend(
    2,
    f"{pick}: higher in tumour in {int(row.n_donors_higher_in_tumour)} of {nd} patients.",
    "one line per patient; red = higher in tumour, blue = lower.",
    "a consistent shift has most lines sloping the same way; a single steep line can dominate a mean but not "
    "this view.",
    n_note(len(tot.donor_id.unique()), tot.study.nunique())
    + f"; higher in tumour in {int(row.n_donors_higher_in_tumour)}.",
    f"GLM log2 change {row.glm_log2fc:+.2f} [{row.glm_ci_low:+.2f}, {row.glm_ci_high:+.2f}], padj "
    f"{row.glm_padj:.2g}; same direction in {int(row.k_studies_same_direction)} of {int(row.k_studies)} studies.",
    "the patient, not the cell, is the unit of this comparison.",
)
show_table(
    c[
        [
            "cell_type",
            "mean_share_adjacent",
            "mean_share_tumour",
            "glm_log2fc",
            "glm_ci_low",
            "glm_ci_high",
            "glm_padj",
            "glm_method",
            "propeller_padj",
            "clr_padj",
            "n_donors_higher_in_tumour",
            "k_studies_same_direction",
            "k_studies",
            "consistent",
        ]
    ].round(4),
    "Paired immune composition, all models. glm_method shows where the NB fit fell back to quasi-Poisson.",
)

# ------------------------------------------------------------------------------------------------ descriptive: compartments
st.subheader("The whole sample, descriptively: what each study captured")
inference_badge("descriptive")
comp = units.groupby(["study", "tissue_type", "compartment"]).n_cells.sum().reset_index()
comp["pct"] = 100 * comp.n_cells / comp.groupby(["study", "tissue_type"]).n_cells.transform("sum")
comp["label"] = comp.study + " · " + comp.tissue_type.str.replace("_", " ")
order = comp.groupby("label").n_cells.sum().sort_values().index
from lib.ui import COMPARTMENT_COLOR  # noqa: E402

fig = go.Figure()
for cp, color in COMPARTMENT_COLOR.items():
    dd = comp[comp.compartment == cp].set_index("label").reindex(order)
    fig.add_trace(go.Bar(y=order, x=dd.pct, orientation="h", name=cp, marker_color=color))
fig.update_layout(barmode="stack")
st.plotly_chart(style(fig, 120 + 20 * len(order), x="% of the study's cells in that tissue"))
figure_legend(
    3,
    "Studies captured very different compartments.",
    "share of each compartment among all cells, per study and tissue.",
    "a bar dominated by one colour is a study that enriched that compartment before sequencing (e.g. CD45+ "
    "sorting = mostly immune).",
    f"{units.donor_id.nunique()} donors, {units.study.nunique()} studies.",
    "cell counts per compartment from the harmonised labels.",
    "this is why composition is tested within the immune compartment and within patients, and why "
    "healthy-vs-tumour proportions are not compared.",
)

# ------------------------------------------------------------------------------------------------ Milo
st.subheader("Neighbourhood-level changes (Milo)")
methods_box(
    "Milo",
    "Builds overlapping neighbourhoods of similar cells on the integrated map and tests, for each, whether "
    "tumour or adjacent tissue contributes more cells (negative-binomial GLM, ~ patient + tissue).",
    "Finds shifts of cell STATES inside a labelled type that a whole-type test averages away.",
    "cluster-level tests (depend on clustering resolution); DA-seq.",
    "Dann et al. 2022 Nat Biotechnol.",
)
mi = q("SELECT nhood_annotation a, logFC, SpatialFDR FROM milo_nhoods")
info = q("SELECT any_value(n_donors) d, any_value(n_samples) u, any_value(design) des FROM milo_nhoods").iloc[
    0
]
# run 3 review: with 0 significant neighbourhoods a count bar chart is empty; show every neighbourhood's
# log2 fold change per cell type (strip plot), significant ones coloured, with the cell-type-level paired GLM
# estimate overlaid, so agreement between the two levels is visible.
sig = mi.SpatialFDR < 0.1
med = mi.groupby("a").logFC.median().sort_values()
keep = mi.groupby("a").size()
med = med[keep.reindex(med.index) >= 10]
cp = q("SELECT cell_type, glm_log2fc FROM composition_paired").set_index("cell_type").glm_log2fc
rng = np.random.default_rng(0)
order = {a: i for i, a in enumerate(med.index)}
m2 = mi[mi.a.isin(order)].copy()
m2["yy"] = m2.a.map(order) + rng.uniform(-0.3, 0.3, len(m2))
fig = go.Figure()
for mask, name, color in [
    (~sig[m2.index], "neighbourhood, SpatialFDR ≥ 0.1", NS_COLOR),
    (sig[m2.index] & (m2.logFC > 0), "more in tumour (SpatialFDR < 0.1)", UP_COLOR),
    (sig[m2.index] & (m2.logFC < 0), "less in tumour (SpatialFDR < 0.1)", DOWN_COLOR),
]:
    d = m2[mask]
    if len(d) or name.startswith("neighbourhood"):
        fig.add_trace(
            go.Scattergl(
                x=d.logFC,
                y=d.yy,
                mode="markers",
                name=f"{name} ({len(d):,})",
                marker=dict(color=color, size=4, opacity=0.45),
                hoverinfo="skip",
            )
        )
fig.add_trace(
    go.Scatter(
        x=med.values,
        y=list(range(len(med))),
        mode="markers",
        name="median neighbourhood log2FC",
        marker=dict(symbol="line-ns", size=16, line=dict(width=3, color=INK)),
    )
)
cc = [(i, cp[a]) for a, i in order.items() if a in cp.index]
if cc:
    fig.add_trace(
        go.Scatter(
            x=[v for _, v in cc],
            y=[i for i, _ in cc],
            mode="markers",
            name="whole-type paired GLM (Figure 1)",
            marker=dict(symbol="diamond", size=10, color="white", line=dict(width=2, color=INK)),
        )
    )
fig.add_vline(x=0, line=dict(color=INK_2, width=1, dash="dot"))
fig.update_yaxes(tickmode="array", tickvals=list(range(len(med))), ticktext=list(med.index))
st.plotly_chart(
    style(fig, 140 + 26 * len(med), x="log2 fold change, tumour vs adjacent (neighbourhood or whole type)")
)
both = [a for a in med.index if a in cp.index]
rho = (  # Spearman = Pearson on ranks (no scipy in the hosted app)
    float(med[both].rank().corr(cp[both].rank())) if len(both) >= 5 else float("nan")
)
n_sig = int(sig.sum())
figure_legend(
    4,
    (
        f"No neighbourhood passes FDR, but neighbourhood shifts point the same way as the whole-type test "
        f"(Spearman ρ = {rho:.2f} across {len(both)} immune types)."
        if n_sig == 0
        else f"{n_sig:,} neighbourhoods change; their direction follows the whole-type test "
        f"(Spearman ρ = {rho:.2f} across {len(both)} immune types)."
    ),
    "each dot is one Milo neighbourhood (a small group of similar cells), placed on the row of its majority "
    "cell type; bars = median per type; diamonds = the whole-type paired estimate from Figure 1 (immune "
    "types only).",
    "a cloud spread across 0 inside one type means states moving in opposite directions; a cloud shifted to "
    "one side agrees with a whole-type change.",
    f"{int(info.d)} patients, {int(info.u)} patient × tissue units; {len(mi):,} neighbourhoods "
    "(types with ≥ 10 neighbourhoods shown).",
    f"Milo, design {info.des}; SpatialFDR (an FDR weighted for overlapping neighbourhoods) < 0.1.",
    "with this many neighbourhoods and 15 patients, per-neighbourhood power is low; read Milo as a map of "
    "where inside a type a shift sits, not as an independent test (issue 18).",
)
st.caption(
    "Healthy donors are not compared here: every healthy donor comes from a different atlas and study than the "
    "tumours (see Cohort), so a difference in proportions could not be told apart from the protocol."
)
footer()
