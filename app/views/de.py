"""Differential expression. Primary = paired tumour vs adjacent within patient (donor × tissue pseudobulk),
with study-grounded robustness (joint per-study model, LOSO, thresholded test). Atlas-confounded comparisons are in
a separate, labelled section. Every number is read from the DB (DEVELOPMENT_LOG §3.6, §3.12)."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import config, has_table, q
from lib.explain import guide
from lib.ui import (
    COMPARTMENT_COLOR,
    DOWN_COLOR,
    GREY,
    INK_2,
    NS_COLOR,
    UP_COLOR,
    caveats,
    cnv_gate_passed,
    decisions,
    figure_legend,
    footer,
    inference_badge,
    methods_box,
    n_note,
    page_header,
    part_tag,
    short_donor,
    show_table,
    style,
)

ROBUST_SQL = (
    "count(*) FILTER (WHERE robust AND padj_thr_global < {alpha} AND NOT coalesce(ambient_suspect, false))"
)
PAIRED = "paired_tumour_vs_adjacent"
dcfg = config().get("de", {})
alpha = dcfg.get("alpha", 0.05)
thr = dcfg.get("lfc_threshold", 0.585)
min_pairs = dcfg.get("min_paired_donors", 6)
loso_a = dcfg.get("loso_alpha", 0.1)
loso_pct = f"{dcfg.get('loso_major_fraction', 0.25):.0%}"

page_header(
    "Differential expression & pathways",
    "Within each cell type, which genes and pathways change in LUSC tumour tissue compared with adjacent "
    "normal lung of the same patients — and do the changes hold study by study?",
    part="evidence",
    part_text="Tumour vs adjacent normal tissue of the same patients, immune cells; the claims of this site come from "
    "here. The comparisons with healthy donors further down this page are Context (descriptive).",
)
guide("de")
caveats("de")

# ------------------------------------------------------------------------------------------------ data
cols = set(q("SELECT * FROM de_pseudobulk LIMIT 0").columns)
has_rob = "robust" in cols
summ = q(
    f"""SELECT "group", any_value(tier) tier, any_value(inference) inference, any_value(n_donors) n_donors,
               any_value(n_studies) n_studies, any_value(studies) studies,
               count(*) n_tested,
               count(*) FILTER (WHERE padj < {alpha}) n_sig,
               count(*) FILTER (WHERE padj_thr < {alpha}) n_thr,
               {ROBUST_SQL.format(alpha=alpha) if has_rob else "0"} n_robust
        FROM de_pseudobulk WHERE contrast = ? GROUP BY 1 ORDER BY n_robust DESC, n_thr DESC""",
    (PAIRED,),
)
primary = summ[summ.tier == "primary"]
_k = (
    q(
        """SELECT "group", max(k_studies) k FROM de_pseudobulk WHERE contrast = ? AND tier = 'primary'
           GROUP BY 1""",
        (PAIRED,),
    )
    if has_rob
    else None
)
k_lo, k_hi = (int(_k.k.min()), int(_k.k.max())) if _k is not None and _k.k.notna().any() else ("?", "?")

# ------------------------------------------------------------------------------------------------ design
st.subheader("The design in one paragraph")
st.markdown(
    "For every patient who gave both a tumour and an adjacent-normal sample, the raw counts of all cells of one "
    "type in one tissue are summed into one **pseudobulk** profile (a patient's tumour regions are pooled, so "
    "each patient counts once per tissue). A negative-binomial model with a term per patient (**~ patient + "
    "tissue**) then compares tumour with adjacent tissue *within* each patient. A cell type is a **primary** "
    f"result only with at least **{min_pairs} patients** from at least **2 studies contributing 2 or more** "
    "each; otherwise it is shown as descriptive."
)
c1, c2 = st.columns(2, gap="medium")
with c1, st.container(border=True, height="stretch"):
    st.markdown("**What pairing cancels**")
    st.markdown(
        "Everything that is the same for both samples of a patient: genotype, age, sex, smoking history, "
        "study, lab and atlas. The sequencing chemistry is also the same within each pair (checked)."
    )
with c2, st.container(border=True, height="stretch"):
    st.markdown("**What it does not cancel**")
    st.markdown(
        "- how the two tissues were handled: digestion, time to processing, blood content, each sample's "
        "ambient RNA (contamination flags cover part of this)\n"
        "- tumour cells inside an 'adjacent' sample (this shrinks real differences toward zero)\n"
        "- changes shared by both tissues, such as smoking damage or a field effect: the comparison cannot "
        "see them"
    )
st.caption(
    "So the result is **tumour relative to adjacent tissue of the same lung, not relative to healthy lung** "
    "(DEVELOPMENT_LOG §3.13, decision 1)."
)
methods_box(
    "pseudobulk + PyDESeq2 (paired)",
    "Counts summed per patient × tissue × cell type; DESeq2's negative-binomial GLM with shrunken dispersions "
    "(Python port, Muzellec et al. 2023).",
    "Treating cells as replicates inflates significance by orders of magnitude (pseudoreplication); pseudobulk "
    "methods are the best-calibrated in benchmarks (Squair et al. 2021; Crowell et al. 2020).",
    "cell-level Wilcoxon/MAST (pseudoreplication); mixed models on cells (NEBULA, slower, similar calibration); "
    "limma-voom/edgeR on pseudobulk (equivalent).",
    "Squair et al. 2021 Nat Commun; Crowell et al. 2020 Nat Commun; Muzellec et al. 2023 Bioinformatics.",
)
methods_box(
    "study-grounded robustness",
    "(1) A thresholded test: is the change larger than 1.5-fold (|log2FC| > "
    f"{thr:.2f}), not merely non-zero? (2) One joint model ~ patient + study:tissue gives each study its own "
    "tumour effect with shared dispersions → fixed-effect meta-analysis with Cochran's Q / I². (3) Leave-one-"
    "study-out refits.",
    "All paired patients come from a few studies, and some hold most of a cell type's pairs. A gene is "
    "**robust** only if it passes the thresholded test, no study with ≥ 2 patients points the other way with a "
    f"95 % CI excluding 0, and it still points the same way (padj < {loso_a}) when each study holding ≥ "
    f"{loso_pct} of that cell type's pairs is left out in turn.",
    "a random effect for how much the tumour effect varies between studies (study differences in baseline "
    f"already cancel through the patient term): with {k_lo}–{k_hi} studies of ≥ 2 patients per cell type, that "
    "variance cannot be estimated reliably (with 2 studies, Hartung–Knapp intervals use a t value on 1 degree "
    "of freedom, 12.7), so results apply to these studies; separate per-study DESeq2 fits (dispersion cannot "
    "be estimated from 1–3 patients).",
    "Nguyen et al. 2023 Nat Commun 14:1570 (meta-analysis vs batch correction for DE); Love et al. 2014 "
    "(thresholded Wald test).",
)
decisions("de")

# ------------------------------------------------------------------------------------------------ fig 1: eligibility
st.subheader("Which cell types can be compared within patients?")
fig = go.Figure()
tier_col = {"primary": "#1c5cab", "descriptive": "#9ec5f4"}
s1 = summ.sort_values("n_donors")
fig.add_trace(
    go.Bar(
        y=s1.group,
        x=s1.n_donors,
        orientation="h",
        marker_color=[tier_col.get(t, GREY) for t in s1.tier],
        text=[
            f"{int(d)} patients · {int(s)} studies" for d, s in zip(s1.n_donors, s1.n_studies, strict=True)
        ],
        textposition="outside",
        customdata=s1.studies,
        hovertemplate="%{y}<br>%{text}<br>%{customdata}<extra></extra>",
    )
)
fig.add_vline(x=min_pairs, line_dash="dash", line_color=INK_2, annotation_text=f"primary needs ≥ {min_pairs}")
fig.update_layout(showlegend=False, xaxis_range=[0, max(s1.n_donors.max() + 5, min_pairs + 3)])
st.plotly_chart(style(fig, max(260, 26 * len(s1) + 80), x="paired patients (tumour + adjacent)", y=None))
figure_legend(
    1,
    "Only immune cell types reach the bar for paired inference.",
    "number of patients contributing both tissues for each cell type (≥ 10 cells per pseudobulk); dark bars = "
    "primary (also ≥ 2 studies with ≥ 2 patients each), light bars = descriptive.",
    "a longer bar means more independent patients; the dashed line is the minimum for a primary result.",
    n_note(int(summ.n_donors.max()), int(summ.n_studies.max())) + " (largest cell type).",
    "eligibility rule fixed in config before the run (senior review).",
    "adjacent-normal tissue held too few epithelial, stromal and endothelial cells, so tumour-epithelium biology "
    "cannot be tested within patient (see the Context section at the end of this page).",
)

# ------------------------------------------------------------------------------------------------ fig 2: counts
st.subheader("How many genes change, and how many hold across studies?")
if not primary.empty:
    p2 = primary.sort_values("n_robust")
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=p2.group,
            x=p2.n_robust,
            orientation="h",
            name="robust (> 1.5-fold at global FDR, no study contradicts, holds when each major study is left out)",
            marker_color="#1c5cab",
        )
    )
    fig.add_trace(
        go.Bar(
            y=p2.group,
            x=p2.n_thr - p2.n_robust,
            orientation="h",
            name="> 1.5-fold, not replicated in every check",
            marker_color="#86b6ef",
        )
    )
    fig.add_trace(
        go.Bar(
            y=p2.group,
            x=(p2.n_sig - p2.n_thr).clip(lower=0),
            orientation="h",
            name="significant, but ≤ 1.5-fold",
            marker_color="#d9d8d3",
        )
    )
    fig.update_layout(barmode="stack")
    st.plotly_chart(style(fig, max(260, 30 * len(p2) + 100), x="number of genes (padj < 0.05)"))
    figure_legend(
        2,
        "Most significant changes are small; a minority are large and hold up in every study check.",
        "per primary cell type, genes significant in the paired model, split by whether the change also exceeds "
        "1.5-fold (thresholded test) and whether it is robust: > 1.5-fold at a global FDR of 5 % across all "
        "primary tests, not contradicted by any study, still pointing the same way (padj < "
        f"{loso_a}) when each study holding ≥ {loso_pct} of the cell type's pairs is left out in turn, and not a "
        "contamination suspect.",
        "read the dark segment as the conservative set to quote; the grey segment is real but small.",
        f"{len(p2)} primary cell types; " + n_note(int(p2.n_donors.max()), int(p2.n_studies.max())) + ".",
        "PyDESeq2 Wald test, BH per cell type (padj); global BH over all primary tests is in the table "
        "(padj_global).",
        "use 'robust' genes for biological claims.",
    )

# ------------------------------------------------------------------------------------------------ volcano
st.subheader("Gene-level view of one cell type")
types = summ.group.tolist()
default = primary.group.iloc[0] if not primary.empty else types[0]
c1, c2 = st.columns([1, 3])
with c1:
    group = st.selectbox("Cell type", types, index=types.index(default))
    row = summ.set_index("group").loc[group]
    inference_badge(row.tier)
    st.caption(n_note(int(row.n_donors), int(row.n_studies)) + f"  \nStudies: {row.studies}")
extra = (
    ", robust, ambient_suspect, ambient_source, ambient_ratio, k_concordant, k_studies, i2" if has_rob else ""
)
d = q(
    f"""SELECT gene, log2FoldChange, lfcSE, ci_low, ci_high, padj, padj_thr, padj_thr_global, baseMean{extra}
        FROM de_pseudobulk WHERE contrast = ? AND "group" = ? AND padj IS NOT NULL""",
    (PAIRED, group),
)
d["nlp"] = -np.log10(d.padj.clip(lower=1e-300))
rob = (
    (d.robust.fillna(False).astype(bool) & (d.padj_thr_global.fillna(1) < alpha))
    if has_rob
    else d.padj_thr < alpha
)
amb = d.ambient_suspect.fillna(False).astype(bool) if has_rob else False
sig = d.padj < alpha
fig = go.Figure()
layers = [
    (~sig, "not significant", NS_COLOR, 4, "circle"),
    (sig & ~rob, "significant, not robust", "#8a8984", 5, "circle"),
    (rob & (d.log2FoldChange < 0) & ~amb, "robust, lower in tumour", DOWN_COLOR, 8, "circle"),
    (rob & (d.log2FoldChange > 0) & ~amb, "robust, higher in tumour", UP_COLOR, 8, "circle"),
    (rob & amb, "robust but contamination suspect", UP_COLOR, 9, "circle-open"),
]
for mask, name, color, size, sym in layers:
    s = d[mask]
    fig.add_trace(
        go.Scattergl(
            x=s.log2FoldChange,
            y=s.nlp,
            mode="markers",
            name=f"{name} ({len(s):,})",
            text=s.gene,
            marker=dict(color=color, size=size, symbol=sym, opacity=0.85),
            hovertemplate="%{text}<br>log2FC %{x:.2f}<br>−log10 padj %{y:.1f}<extra></extra>",
        )
    )
lab = d[rob & ~amb].assign(a=lambda x: x.log2FoldChange.abs()).nlargest(12, "a")
fig.add_trace(
    go.Scatter(
        x=lab.log2FoldChange,
        y=lab.nlp,
        mode="text",
        text=lab.gene,
        textposition="top center",
        showlegend=False,
        hoverinfo="skip",
        textfont=dict(size=11),
    )
)
for v in (-thr, thr):
    fig.add_vline(x=v, line_dash="dot", line_color=GREY)
with c2:
    st.plotly_chart(
        style(
            fig,
            500,
            x="log2 fold change, tumour ÷ adjacent normal (same patients)",
            y="−log10 adjusted p-value",
        )
    )
figure_legend(
    3,
    f"Paired tumour-vs-adjacent changes in {group}.",
    "each dot is a gene; dotted lines mark ±1.5-fold.",
    "right = higher in tumour, up = stronger evidence. Coloured dots are robust; open circles are robust genes "
    "that are ≥ 10× more abundant in another cell type, so ambient RNA or contamination could explain them.",
    n_note(int(row.n_donors), int(row.n_studies)) + f"; {len(d):,} genes tested.",
    "PyDESeq2 Wald test, BH per cell type; labels = the largest robust, non-suspect changes.",
    "quote the labelled genes with their confidence intervals (table below), not the dot cloud.",
)
if has_rob:
    tab = (
        d[sig]
        .sort_values("padj")[
            [
                "gene",
                "log2FoldChange",
                "ci_low",
                "ci_high",
                "padj",
                "padj_thr",
                "robust",
                "k_concordant",
                "k_studies",
                "i2",
                "ambient_suspect",
                "ambient_source",
                "ambient_ratio",
                "baseMean",
            ]
        ]
        .rename(
            columns={
                "log2FoldChange": "log2FC",
                "ci_low": "95% CI low",
                "ci_high": "95% CI high",
                "padj_thr": "padj (>1.5-fold test)",
                "k_concordant": "studies agreeing",
                "k_studies": "studies with ≥2 patients",
                "i2": "I² (heterogeneity)",
            }
        )
    )
    show_table(tab.round(3), f"{group}: significant genes, most significant first.")

# ------------------------------------------------------------------------------------------------ forest
if has_table("de_per_study"):
    st.subheader("Study by study: does a gene's change replicate?")
    genes_here = d[rob].sort_values("padj").gene.tolist() or d.sort_values("padj").gene.head(50).tolist()
    gsel = st.selectbox("Gene (robust genes of this cell type listed first)", genes_here[:200])
    ps = q(
        """SELECT study, n_donors, log2FoldChange l, lfcSE se FROM de_per_study
              WHERE contrast = ? AND "group" = ? AND gene = ? ORDER BY n_donors DESC""",
        (PAIRED, group, gsel),
    )
    pooled = d[d.gene == gsel].iloc[0]
    if not ps.empty:
        ps["lo"], ps["hi"] = ps.l - 1.96 * ps.se, ps.l + 1.96 * ps.se
        ylab = [f"{r.study} ({int(r.n_donors)} pt{'s' if r.n_donors != 1 else ''})" for r in ps.itertuples()]
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=ps.l,
                y=ylab,
                mode="markers",
                name="per-study estimate (joint model)",
                marker=dict(size=9 + 2 * ps.n_donors.clip(upper=6), color=INK_2, symbol="square"),
                error_x=dict(
                    type="data", symmetric=False, array=ps.hi - ps.l, arrayminus=ps.l - ps.lo, color=INK_2
                ),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[pooled.log2FoldChange],
                y=["All paired patients (pooled)"],
                mode="markers",
                name="pooled estimate",
                marker=dict(
                    size=14, symbol="diamond", color=UP_COLOR if pooled.log2FoldChange > 0 else DOWN_COLOR
                ),
                error_x=dict(
                    type="data",
                    symmetric=False,
                    array=[pooled.ci_high - pooled.log2FoldChange],
                    arrayminus=[pooled.log2FoldChange - pooled.ci_low],
                ),
            )
        )
        fig.add_vline(x=0, line_color=GREY)
        for v in (-thr, thr):
            fig.add_vline(x=v, line_dash="dot", line_color=GREY)
        fig.update_yaxes(autorange="reversed")
        st.plotly_chart(
            style(fig, 120 + 45 * (len(ps) + 1), x=f"{gsel} log2 fold change, tumour ÷ adjacent (95 % CI)")
        )
        figure_legend(
            4,
            f"{gsel} in {group}, study by study.",
            "squares = each study's own tumour effect from one joint model (size ∝ patients), diamond = all "
            "paired patients; bars = 95 % confidence intervals; dotted lines = ±1.5-fold.",
            "a replicated change has every study on the same side of zero; studies with one patient have wide "
            "intervals and never count as replication or contradiction.",
            n_note(int(row.n_donors), int(row.n_studies)) + ".",
            "PyDESeq2 with explicit design ~ patient + study:tissue (shared dispersions); fixed-effect meta-analysis "
            f"(I² = {float(pooled.get('i2', np.nan)):.2f} for this gene).",
            "this is the study-grounded view: a result driven by one study is visible here.",
        )

# ------------------------------------------------------------------------------------------------ power
if has_table("de_power"):
    st.subheader("What could these data detect?")
    pw = q("SELECT * FROM de_power WHERE gene_class = 'median-expressed gene' ORDER BY n_pairs DESC")
    fig = go.Figure(
        go.Bar(
            x=pw.group,
            y=2**pw.mde_log2fc,
            marker_color="#86b6ef",
            text=[f"{int(n)} pts" for n in pw.n_pairs],
            textposition="outside",
        )
    )
    fig.add_hline(y=1.5, line_dash="dot", line_color=GREY, annotation_text="1.5-fold")
    st.plotly_chart(style(fig, 320, x="cell type", y="smallest fold change detectable (80 % power)"))
    figure_legend(
        5,
        "With these patient numbers, only fairly large changes are reliably detectable.",
        "approximate minimum detectable fold change for a gene of median expression in each cell type.",
        "taller bars = less power; an absent gene in Figure 3 is not evidence of no change below this size.",
        "patients per cell type on the bars.",
        "Wald approximation for a paired negative-binomial design, SE ≈ √(2(1/μ + φ)/n)/ln 2, two-sided α = 0.05, "
        "power 0.8, μ and φ = median normalised count and dispersion of the cell type's genes.",
        "a non-significant result means 'not detected', not 'unchanged'. These are detection limits for the plain "
        "test against 0 (padj); the robust set additionally requires > 1.5-fold at a global FDR, so its "
        "detection limit is roughly 0.585 log2 higher.",
    )

# ------------------------------------------------------------------------------------------------ checks (§3.13)
fig_no = 6  # figures from here on are numbered as they render (a panel with no data takes no number)
if has_table("de_paired_null") or has_table("adjacent_makeup"):
    st.subheader("Two checks of the paired comparison itself")
    st.markdown(
        "The study-by-study views above ask whether a change *replicates*. These two ask whether the paired "
        "test can be trusted at all: does it produce false hits when there is nothing to find, and what is "
        "actually in the 'normal' samples it uses as controls?"
    )
if has_table("de_paired_null"):
    flag = float(dcfg.get("sign_flip_flag_frac", 0.03))
    pn = q('SELECT "group", seed, n_pairs, n_flipped, n_genes, frac_p_lt_001, flagged FROM de_paired_null')
    agg = (
        pn.groupby("group")
        .agg(mean=("frac_p_lt_001", "mean"), n_pairs=("n_pairs", "max"), flagged=("flagged", "max"))
        .sort_values("mean")
    )
    n_flag = int(agg.flagged.astype(bool).sum())
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            y=agg.index,
            x=100 * agg["mean"],
            orientation="h",
            name="mean of 3 label swaps",
            marker_color=["#a3501f" if f else "#86b6ef" for f in agg.flagged.astype(bool)],
            text=[f"{100 * v:.2f} %" for v in agg["mean"]],
            textposition="outside",
            customdata=agg.n_pairs,
            hovertemplate="%{y}<br>mean %{x:.2f} % of genes at p < 0.01<br>%{customdata} patients used in the "
            "swap<extra></extra>",
        )
    )
    fig.add_trace(
        go.Scatter(
            y=pn.group,
            x=100 * pn.frac_p_lt_001,
            mode="markers",
            name="each label swap",
            marker=dict(color=INK_2, size=7, symbol="line-ns-open", line=dict(width=2, color=INK_2)),
            hovertemplate="%{y}<br>%{x:.2f} %<extra></extra>",
        )
    )
    fig.add_vline(x=1, line_dash="dot", line_color=GREY, annotation_text="1 % expected")
    fig.add_vline(
        x=100 * flag,
        line_dash="dash",
        line_color="#a3501f",
        annotation_text=f"flag above {100 * flag:.0f} %",
        annotation_position="top left",
    )
    fig.update_layout(xaxis_range=[0, max(100 * flag * 1.25, 100 * pn.frac_p_lt_001.max() * 1.15)])
    fig.update_yaxes(categoryorder="array", categoryarray=list(agg.index))
    st.plotly_chart(
        style(fig, max(280, 30 * len(agg) + 110), x="% of genes with p < 0.01 after the label swap")
    )
    figure_legend(
        fig_no,
        (
            "With tumour and adjacent labels swapped in half the patients, no cell type produces too many hits."
            if n_flag == 0
            else f"{n_flag} cell type(s) produce too many hits after the label swap: their p-values are too small."
        ),
        "for each primary cell type, the paired model was refitted after swapping the tumour and adjacent "
        "labels in a random half of the patients (3 different random halves). Bars = mean share of genes "
        "with p < 0.01; ticks = each of the 3 swaps.",
        "after the swap there is no consistent tumour effect left, so about 1 % of genes or fewer should reach p < 0.01 "
        f"by chance. A bar past the dashed line ({100 * flag:.0f} %, fixed before run 4) would mean the paired "
        "test finds differences that are not there.",
        f"{len(agg)} primary cell types, {int(agg.n_pairs.min())}–{int(agg.n_pairs.max())} patients used in the "
        "swap per cell type (with an odd number of paired patients, one random patient is left out so the two "
        "halves are equal, so this can be one fewer than in Figure 1).",
        "PyDESeq2 ~ patient + tissue, the same model as the main test; sign-flip (label-swap) null.",
        f"mean {100 * pn.frac_p_lt_001.mean():.2f} % across all cell types: the paired test is not too liberal, "
        "if anything somewhat conservative. The check itself leans conservative: real tumour effects do not "
        "disappear in the swap, they become extra noise (up in half the patients, down in the other half), "
        "which makes small p-values rarer. So a pass is weaker evidence than a fail would be.",
    )
    fig_no += 1
    show_table(
        pn.assign(pct_p_lt_001=lambda x: (100 * x.frac_p_lt_001).round(2)).drop(columns="frac_p_lt_001"),
        "Paired sign-flip null: one row per cell type and random label swap (table de_paired_null).",
    )
if has_table("adjacent_makeup"):
    st.markdown("**What the controls contain: each patient's adjacent-normal sample**")
    am = q("SELECT * FROM adjacent_makeup ORDER BY study, donor_id")
    am["patient"] = am.donor_id.map(short_donor)  # short names already carry the study
    fig = go.Figure()
    for comp in ("Immune", "Epithelial", "Endothelial", "Stroma", "Unknown"):
        col_ = f"share_{comp.lower()}"
        if col_ not in am or not (am[col_] > 0).any():
            continue
        fig.add_trace(
            go.Bar(
                y=am.patient,
                x=100 * am[col_],
                orientation="h",
                name=comp,
                marker_color=COMPARTMENT_COLOR[comp],
                customdata=am.n_cells,
                hovertemplate="%{y}<br>" + comp + ": %{x:.0f} % of %{customdata:,} cells<extra></extra>",
            )
        )
    fig.update_layout(barmode="stack")
    fig.update_yaxes(autorange="reversed", categoryorder="array", categoryarray=list(am.patient))
    st.plotly_chart(style(fig, 90 + 28 * len(am), x="% of the adjacent sample's cells"))
    called = am[am.n_epithelial_cnv_tested > 0]
    hi = called.sort_values("cnv_call_rate", ascending=False).head(1)
    # which labels do the CNV-called adjacent cells carry? (review: mostly normal airway/alveolar types)
    lab = q("""SELECT harmonised_label, count(*) n FROM cells
               WHERE condition = 'adjacent_normal' AND cnv_status LIKE 'CNV-high%' GROUP BY 1""")
    normal_types = ("Multiciliated", "AT1", "AT2", "Secretory")
    n_normal = int(lab.loc[lab.harmonised_label.isin(normal_types), "n"].sum())
    n_called = int(lab.n.sum())
    figure_legend(
        fig_no,
        "The 'normal' controls are mostly immune cells, and how much else they contain depends on the study.",
        "for each paired patient, the share of cells in the adjacent-normal sample by compartment.",
        "the controls reflect each study's protocol (several enriched immune or endothelial cells before "
        "sequencing), which is why only immune cell types can be compared within patients. The table below "
        "adds what kind of epithelium the sample holds and how much of it looks tumour-like.",
        f"{len(am)} paired patients, {am.study.nunique()} studies, "
        f"{int(am.n_cells.min()):,}–{int(am.n_cells.max()):,} cells per adjacent sample (after QC).",
        "counts of the harmonised cell labels; descriptive, changes no result (DEVELOPMENT_LOG §3.13, check 5b).",
        (
            f"where tested, up to {100 * float(hi.cnv_call_rate.iloc[0]):.0f} % of adjacent epithelial cells "
            f"carry a copy-number call ({hi.patient.iloc[0]}). But {n_normal} of the {n_called} called adjacent "
            "cells are labelled as normal epithelium (ciliated, secretory, AT1 or AT2), which points more to "
            "false-positive calls than to tumour cells in the controls. If a control did hold tumour, it would "
            "likely also hold tumour-infiltrating immune cells, which would shrink differences toward zero; the "
            "'robust' rule requires each result to hold when each major study is left out in turn."
            if len(hi)
            else "no adjacent epithelium had enough cells for a copy-number test."
        ),
    )
    fig_no += 1
    tab = am[
        [
            "patient",
            "n_cells",
            "n_epithelial_labelled",
            "share_airway_of_epithelium",
            "share_blood_myeloid_of_immune",
            "n_epithelial_cnv_tested",
            "cnv_call_rate",
        ]
    ].copy()
    for c in ("share_airway_of_epithelium", "share_blood_myeloid_of_immune", "cnv_call_rate"):
        tab[c] = (100 * tab[c]).round(0)
    st.dataframe(
        tab.rename(
            columns={
                "n_cells": "cells",
                "n_epithelial_labelled": "epithelial cells labelled airway or alveolar",
                "share_airway_of_epithelium": "% airway (basal, ciliated, secretory)",
                "share_blood_myeloid_of_immune": "% blood-type myeloid of immune",
                "n_epithelial_cnv_tested": "epithelial cells CNV-tested",
                "cnv_call_rate": "% CNV-called",
            }
        ),
        hide_index=True,
    )
    st.caption(
        "Airway epithelium is where LUSC arises; alveolar epithelium (AT1, AT2) is normal peripheral lung. "
        "Blood-type myeloid cells (classical and non-classical monocytes, neutrophils) also reach tissue from "
        "blood, so a high share can mean blood content. '% airway' is computed over cells labelled airway or "
        "alveolar; 'CNV-tested' counts all epithelial cells, including unlabelled ones, so the two counts differ. "
        "Compare patients within a study: the shares mostly reflect each study's protocol. Percentages based "
        "on a handful of cells (see the counts) are "
        "unreliable. Empty = not tested or no cells. The malignant-cell calls failed their pre-registered "
        "false-positive check (Malignant cells page), so read the CNV column as a hint, not a diagnosis."
    )

# ------------------------------------------------------------------------------------------------ pathways
st.subheader("Pathway and transcription-factor activity (paired contrast)")
methods_box(
    "decoupler with PROGENy, CollecTRI and MSigDB Hallmark",
    "Scores each pathway or transcription factor from the DE statistic of its target genes (univariate linear "
    "model for PROGENy/CollecTRI; preranked enrichment for Hallmark).",
    "Footprint methods (genes that RESPOND to a pathway) are more reliable than gene sets of pathway members; "
    "scoring the paired DE statistic keeps the paired design.",
    "GSEA with cell-level permutations (pseudoreplication); over-representation on a significant-gene list "
    "(threshold-dependent).",
    "Badia-i-Mompel et al. 2022; Schubert et al. 2018 (PROGENy); Müller-Dott et al. 2023 (CollecTRI); "
    "Liberzon et al. 2015 (Hallmark).",
)
pg = q('SELECT "group", source, score, padj FROM pathways_progeny WHERE contrast = ?', (PAIRED,))
pg = pg[pg.group.isin(primary.group)]
if not pg.empty:
    mat = pg.pivot(index="group", columns="source", values="score")
    star = pg.pivot(index="group", columns="source", values="padj").reindex_like(mat)
    txt = np.where(star < 0.05, "•", "")
    lim = float(np.nanmax(np.abs(mat.to_numpy())))
    fig = go.Figure(
        go.Heatmap(
            z=mat.to_numpy(),
            x=mat.columns,
            y=mat.index,
            text=txt,
            texttemplate="%{text}",
            colorscale=[[0, DOWN_COLOR], [0.5, "#f0efec"], [1, UP_COLOR]],
            zmin=-lim,
            zmax=lim,
            colorbar=dict(title="activity<br>score"),
        )
    )
    st.plotly_chart(style(fig, 80 + 34 * len(mat), x="PROGENy pathway", y="cell type"))
    figure_legend(
        fig_no,
        "Pathway footprints in tumour vs adjacent tissue.",
        "PROGENy activity score per primary cell type and pathway; • = BH padj < 0.05.",
        "red = pathway more active in tumour, blue = less; colour depth = score.",
        f"{len(mat)} primary cell types; same patients as Figure 2.",
        "decoupler ULM on the paired DESeq2 statistic; BH across pathways within each cell type.",
        "read together with Figure 3: a pathway shift should be visible in its target genes.",
    )
    fig_no += 1
c1, c2 = st.columns(2)
for col, tbl, title in [
    (c1, "pathways_hallmark", "Hallmark gene sets"),
    (c2, "pathways_collectri", "Transcription factors (CollecTRI)"),
]:
    with col:
        st.markdown(f"**{title} — {group}**")
        t = q(
            f"""SELECT source, score, padj FROM {tbl} WHERE contrast = ? AND "group" = ? AND padj < 0.05
                  ORDER BY abs(score) DESC LIMIT 15""",
            (PAIRED, group),
        )
        if t.empty:
            st.caption(f"No significant entries for {group}.")
            continue
        fig = go.Figure(
            go.Bar(
                y=t.source,
                x=t.score,
                orientation="h",
                marker_color=[UP_COLOR if v > 0 else DOWN_COLOR for v in t.score],
            )
        )
        fig.update_layout(yaxis=dict(autorange="reversed"))
        st.plotly_chart(style(fig, 460, x="activity score (tumour vs adjacent)"))
        n_up, n_dn = int((t.score > 0).sum()), int((t.score < 0).sum())
        top1 = t.iloc[0]
        figure_legend(
            fig_no,
            f"{group}: strongest {title.split(' (')[0].lower()} shift is {top1.source} "
            f"({'higher' if top1.score > 0 else 'lower'} in tumour); {n_up} of the top {len(t)} are higher, "
            f"{n_dn} lower.",
            "the (up to) 15 strongest significant scores.",
            "red = higher in tumour, blue = lower.",
            n_note(int(row.n_donors), int(row.n_studies)) + ".",
            "decoupler on the paired DE statistic, BH padj < 0.05.",
            "hypotheses for mechanism.",
        )
        fig_no += 1

# ------------------------------------------------------------------------------------------------ QC of pairing
if has_table("de_sex_check"):
    sx = q("SELECT * FROM de_sex_check")
    if "undetermined" not in set(sx.inferred_sex):  # DB built before the tie fix
        sx.loc[(sx.XIST == 0) & (sx.RPS4Y1 == 0), "inferred_sex"] = "undetermined"
    paired = sx.groupby("donor_id").tissue_type.nunique()
    paired = paired[paired > 1].index
    sp_ = sx[sx.donor_id.isin(paired) & (sx.inferred_sex != "undetermined")]
    ok = sp_.groupby("donor_id").inferred_sex.nunique().eq(1)
    per_donor = (
        sx[sx.inferred_sex != "undetermined"]
        .groupby("donor_id")
        .agg(inferred=("inferred_sex", "first"), recorded=("recorded_sex", "first"))
    )
    known = per_donor[~per_donor.recorded.isin(["unknown", "nan", "None"])]
    disc = known[known.inferred != known.recorded]
    n_und = int(sx.groupby("donor_id").inferred_sex.apply(lambda v: (v == "undetermined").all()).sum())
    st.caption(
        f"Sample-swap check (XIST vs RPS4Y1 per patient × tissue): inferred sex agrees between tumour and adjacent "
        f"tissue for **{int(ok.sum())} of {len(ok)}** paired patients"
        + (f" (inconsistent: {', '.join(ok[~ok].index)})" if (~ok).any() else "")
        + f". Inferred vs recorded sex, all donors: {len(known) - len(disc)} agree, {len(disc)} discordant"
        + (
            f" ({', '.join(d.split('_')[-1] for d in disc.index)}; healthy donors only)"
            if len(disc) and not sx[sx.donor_id.isin(disc.index)].tissue_type.ne("healthy").any()
            else (f" ({', '.join(d.split('_')[-1] for d in disc.index)})" if len(disc) else "")
        )
        + f", {len(per_donor) - len(known)} with no recorded sex; {n_und} undetermined (both markers absent)."
    )

# ------------------------------------------------------------------------------------------------ gene lookup
st.subheader("Look up a gene across cell types")
gene = st.text_input("Gene symbol", "SPP1").strip().upper()
g = q(
    """SELECT "group", log2FoldChange l, ci_low, ci_high, padj, tier FROM de_pseudobulk
         WHERE gene = ? AND contrast = ? ORDER BY l""",
    (gene, PAIRED),
)
if g.empty:
    st.caption("Gene not tested in the paired comparison (not expressed enough in any eligible cell type).")
else:
    fig = go.Figure(
        go.Scatter(
            x=g.l,
            y=g.group,
            mode="markers",
            marker=dict(
                size=11,
                color=[
                    (UP_COLOR if v > 0 else DOWN_COLOR) if p < alpha else GREY
                    for v, p in zip(g.l, g.padj, strict=True)
                ],
                symbol=["circle" if t == "primary" else "circle-open" for t in g.tier],
            ),
            error_x=dict(
                type="data",
                symmetric=False,
                array=g.ci_high - g.l,
                arrayminus=g.l - g.ci_low,
                color="#8a8984",
            ),
            hovertemplate="%{y}<br>log2FC %{x:.2f}<extra></extra>",
        )
    )
    fig.add_vline(x=0, line_color=GREY)
    st.plotly_chart(
        style(fig, max(260, 26 * len(g) + 100), x=f"{gene}: log2 fold change, tumour ÷ adjacent (95 % CI)")
    )
    figure_legend(
        fig_no,
        f"{gene} across cell types (paired).",
        "estimate and 95 % CI per cell type; filled = primary tier, open = descriptive; coloured = padj < 0.05.",
        "right = higher in tumour.",
        "patients per cell type: see Figure 1.",
        "paired PyDESeq2 per cell type.",
        "cell-type specificity of a change is visible directly.",
    )
    fig_no += 1
# ------------------------------------------------------------------------------------------------ exploratory
st.divider()
st.subheader("Context: comparisons with healthy donors (different studies, descriptive)")
part_tag(
    "context",
    "Everything below compares patients with healthy donors from other studies. It is descriptive and is "
    "not part of the evidence above.",
)
st.markdown(
    "Healthy donors all come from HLCA studies and tumours from LuCA studies. A model with a random study "
    "effect could be fitted, but healthy and tumour studies also differ systematically in atlas, cell sorting "
    "and chemistry, so it would credit those study differences to the tumour. These comparisons are shown because they are the "
    "only way to ask about tumour **epithelium** (adjacent tissue has too few epithelial cells to pair). Read "
    "them through two lenses: **(a)** how big are pure study differences, and **(b)** does the paired "
    "comparison support a gene?"
)
if has_table("de_cross_study_null"):
    nul = q('SELECT "group", frac_p_lt_001, n_case, n_control, assay FROM de_cross_study_null')
    ex = q(
        """SELECT contrast, "group", avg(CASE WHEN pvalue < 0.01 THEN 1.0 ELSE 0 END) f
              FROM de_pseudobulk WHERE contrast IN ('tumour_vs_healthy', 'adjacent_vs_healthy', ?)
              GROUP BY 1, 2""",
        (PAIRED,),
    )
    common = sorted(set(nul.group) & set(ex[ex.contrast == "tumour_vs_healthy"].group))
    if common:
        fig = go.Figure()
        wsn = q("SELECT * FROM de_within_study_null") if has_table("de_within_study_null") else None
        bars = []
        if wsn is not None and len(wsn):
            bars.append(
                (
                    "healthy vs healthy, random split within one study (calibration)",
                    wsn.groupby("group").frac_p_lt_001.mean(),
                    INK_2,
                )
            )
        for name, dd, color in bars + [
            (
                "healthy vs healthy, two sets of studies (technical null)",
                nul.set_index("group").frac_p_lt_001,
                "#b9b8b2",
            ),
            (
                "adjacent normal (LuCA) vs healthy (HLCA)",
                ex[ex.contrast == "adjacent_vs_healthy"].set_index("group").f,
                "#86b6ef",
            ),
            ("tumour vs healthy", ex[ex.contrast == "tumour_vs_healthy"].set_index("group").f, "#eb6834"),
            (
                "tumour vs adjacent, same patients (primary)",
                ex[ex.contrast == PAIRED].set_index("group").f,
                "#1c5cab",
            ),
        ]:
            fig.add_trace(
                go.Bar(x=common, y=[100 * dd.get(g, np.nan) for g in common], name=name, marker_color=color)
            )
        fig.add_hline(
            y=1, line_dash="dot", line_color=INK_2, annotation_text="1 % expected if nothing differs"
        )
        st.plotly_chart(style(fig, 380, x="cell type", y="% of genes with p < 0.01"))
        figure_legend(
            fig_no,
            "Study differences alone make a large share of genes 'significant'; the test itself is roughly calibrated.",
            "for each cell type, the share of tested genes with p < 0.01 in each comparison.",
            "dark grey = healthy donors of ONE study split at random (nothing differs: should sit near the 1 % "
            "line; it does on average, though not in every cell type, so the unpaired test is roughly calibrated); light grey = healthy donors of two sets of "
            "studies of the same chemistry, lung tissue only (no tumour and, since run 4, no lavage samples; "
            "only study and cohort differ). Everything the light grey "
            "bars show is study effect, and it is as large as the tumour-vs-healthy signal; the sizes and "
            "covariates of the comparisons differ, so compare orders of magnitude, not exact heights.",
            f"technical null: up to {nul.n_case.max()} vs {nul.n_control.max()} donors per cell type, assay "
            f"{nul.assay.iloc[0]}"
            + (
                f"; within-study split: {wsn.study.iloc[0]}, 3 random splits (mean shown)"
                if wsn is not None
                else ""
            )
            + ".",
            "unpaired PyDESeq2; study cannot enter the tumour-vs-healthy models as a fixed effect (it is aliased "
            "with condition).",
            "only genes also supported by the paired comparison are reported from the Context contrasts.",
        )
        fig_no += 1
xs = (
    q(
        """SELECT contrast, "group", count(*) FILTER (WHERE padj < 0.05) n_sig,
                 count(*) FILTER (WHERE padj < 0.05 AND supported_by_paired) n_sup,
                 any_value(inference) inference, any_value(n_case) n_case, any_value(n_control) n_control,
                 any_value(n_donors) n_donors, any_value(n_studies) n_studies
          FROM de_pseudobulk WHERE contrast <> ? GROUP BY 1, 2 ORDER BY 1, 3 DESC""",
        (PAIRED,),
    )
    if "supported_by_paired" in cols
    else None
)
if xs is not None and not xs.empty:
    show_table(
        xs.rename(columns={"n_sig": "significant (padj < 0.05)", "n_sup": "…also supported by paired"}),
        "Context (descriptive, hypothesis-generating) contrasts. 'Supported by paired' = same sign and > 1.5-fold "
        "significant in the paired tumour-vs-adjacent comparison.",
    )
st.markdown("**Tumour epithelium vs healthy basal cells** (the only epithelial comparison available)")
inference_badge("hypothesis" if cnv_gate_passed() else "provisional")
for name, label in [
    ("malignant_vs_healthy_basal", "CNV-called malignant cells"),
    ("cnvlow_tumour_vs_healthy_basal", "tumour epithelium without a CNV call"),
]:
    t = q(
        """SELECT gene, log2FoldChange, ci_low, ci_high, padj, n_case, n_control, n_donors, n_studies
             FROM de_pseudobulk WHERE contrast = ? AND padj < 0.05 AND gene NOT LIKE 'MT%'
             ORDER BY abs(log2FoldChange) DESC LIMIT 25""",
        (name,),
    )
    if t.empty:
        continue
    r0 = t.iloc[0]
    st.caption(
        f"{label} vs healthy basal: {int(r0.n_case)} vs {int(r0.n_control)} pseudobulks, "
        f"{int(r0.n_donors)} donors from {int(r0.n_studies)} studies. Largest changes (mitochondrial "
        "genes and pseudogenes hidden — a chemistry signature):"
    )
    st.dataframe(
        t[["gene", "log2FoldChange", "ci_low", "ci_high", "padj"]].round(3), hide_index=True, height=240
    )

footer()
