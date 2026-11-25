"""Malignant cells: per patient × tissue two-axis CNV calling against a held-out same-sample null, window
mask, pre-registered acceptance gate, hallmark arms, per-patient counts, author-label agreement (DEVELOPMENT_LOG
§3.5, §3.12)."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import has_table, q
from lib.explain import guide
from lib.ui import (
    DOWN_COLOR,
    GREY,
    INK_2,
    MALIGNANT_COLOR,
    TISSUE_COLOR,
    UP_COLOR,
    caveats,
    decisions,
    figure_legend,
    footer,
    methods_box,
    metric_row,
    n_note,
    page_header,
    show_table,
    style,
)

page_header(
    "Malignant cells",
    "Which epithelial cells are tumour cells? Tumour cells carry copy-number changes (gains and losses of "
    "chromosome segments) that can be inferred from expression; normal cells do not.",
    part="context",
    part_text="Descriptive. Tumour cells come mostly from tumour-only patients; the calls failed their "
    "pre-registered false-positive check. Not evidence.",
)
guide("malignant")
caveats("malignant")

cs = q("SELECT * FROM cnv_summary").set_index("cnv_group")
units = q("SELECT * FROM cnv_per_sample") if has_table("cnv_per_sample") else None
mal = cs.loc["malignant (CNV-high)"] if "malignant (CNV-high)" in cs.index else None

# ------------------------------------------------------------------------------------------------ method
st.subheader("How a cell is called malignant")
st.markdown(
    "Genes are ordered along each chromosome and their expression is averaged over windows of neighbouring "
    "genes. Relative to cells that are copy-number normal (immune, endothelial and stromal cells), a gained "
    "segment shows as a run of higher expression and a lost segment as a run of lower expression. Each patient's "
    "tumour and adjacent tissue are analysed **separately, against that patient's own normal cells**, so the "
    "sequencing depth, chemistry and atlas of the reference match the cells being judged."
)
methods_box(
    "inferCNV (infercnvpy) with a two-axis call",
    "Per patient × tissue: (a) CNV signal = mean squared CNV value across chromosome windows; (b) how well a "
    "cell's CNV profile correlates with the consensus profile of its own tumour (built from the other half of "
    "the cells — cross-fitted, so a cell never helps build the profile it is compared with). A fifth of the "
    "reference cells (≥ 100) are held out: their values show what 'no CNV' looks like. A cell is malignant when "
    "both axes exceed the 99th percentile of those held-out cells.",
    "Copy-number changes are a label-free, tumour-specific signal; requiring a clonal pattern (axis b) separates "
    "tumour cells from noisy normal cells (Tirosh 2016; Neftel 2019). Runs 1–2 thresholded infercnvpy's "
    "cnv_score, which is the mean of a cell's cluster, against healthy cells of the other atlas.",
    "CopyKAT, SCEVAN (R; expression-only); Numbat (uses allele counts — the most accurate, but needs read-level "
    "data these atlases do not provide); author labels (differ between atlases).",
    "Patel et al. 2014 Science; Tirosh et al. 2016 Science; Neftel et al. 2019 Cell; infercnvpy (Sturm et al.).",
)
methods_box(
    "window mask for expression programmes",
    "inferCNV is also run on normal epithelium (healthy basal cells + adjacent-normal epithelium); chromosome "
    "windows where it shows a systematic shift are masked before scoring.",
    "Genes of one programme can sit together on a chromosome (keratins on 12q13 and 17q21, the epidermal "
    "differentiation complex on 1q21), so a squamous cell can look 'CNV-high' without any copy-number change.",
    "no mask (inflates squamous/epithelial false positives); a fixed list of loci (incomplete).",
    "§3.12 senior review; Tickle et al. (inferCNV documentation on expression artefacts).",
)
decisions("malignant")

# ------------------------------------------------------------------------------------------------ gate
st.subheader("Pre-registered acceptance check")
st.markdown(
    "Before run 3 started, four criteria were written into the configuration. If any fails, every result that "
    "uses malignant calls is flagged."
)
if "acceptance_passed" in cs.columns:
    checks = [
        (
            "adjacent_call_rate_ok",
            "≤ 5 % of adjacent-normal epithelium is called (false-positive check)",
            f"{float(cs['adjacent_call_rate'].dropna().iloc[0]):.1%}"
            if cs["adjacent_call_rate"].notna().any()
            else "n/a",
        ),
        (
            "min_donors_ok",
            "calls in ≥ 10 patients",
            f"{int(mal.n_donors)} patients" if mal is not None else "n/a",
        ),
        (
            "chr3q_gain_ok",
            "called cells show chr3q gain on average (SOX2/TP63 locus)",
            f"{mal.mean_cnv_chr3q:+.4f}" if mal is not None else "n/a",
        ),
        (
            "chr3p_loss_ok",
            "called cells show chr3p loss on average",
            f"{mal.mean_cnv_chr3p:+.4f}" if mal is not None else "n/a",
        ),
    ]
    rows = [f"| {'✅' if bool(cs[k].iloc[0]) else '❌'} | {txt} | {val} |" for k, txt, val in checks]
    st.markdown("| | Criterion | Observed |\n|---|---|---|\n" + "\n".join(rows))
    if bool(cs["acceptance_passed"].iloc[0]):
        st.success("All four criteria pass: malignant calls are used downstream as planned.")
    else:
        st.error(
            "At least one pre-registered criterion fails: treat every malignant-cell result as provisional."
        )
    if not bool(cs["adjacent_call_rate_ok"].iloc[0]) and units is not None:
        adj = units[(units.tissue_type == "adjacent_normal") & units.resolved.astype(bool)]
        if len(adj):
            by_study = (
                adj.groupby("study")
                .agg(
                    patients=("unit", "size"), epithelial=("n_epithelial", "sum"), called=("n_called", "sum")
                )
                .assign(call_rate=lambda d: 100 * d.called / d.epithelial)
                .sort_values("call_rate", ascending=False)
            )
            st.markdown(
                "**Where the false-positive check fails (post hoc, not a pre-registered test).** Adjacent-normal "
                "epithelium called CNV-high, per study. The criterion stays failed whatever this table shows; it "
                "is shown so the failure can be located. Calls in 'normal' tissue can be real (tumour cells or "
                "field-cancerised epithelium in the adjacent sample) or artefact; the data cannot tell which."
            )
            st.dataframe(
                by_study.reset_index().round({"call_rate": 1}),
                hide_index=True,
                column_config={"call_rate": st.column_config.NumberColumn("call rate (%)", format="%.1f")},
            )

if mal is not None:
    metric_row(
        [
            ("Malignant cells", f"{int(mal.n_cells):,}", None),
            ("Patients with calls", f"{int(mal.n_donors)}", None),
            (
                "Patient × tissue samples scored",
                f"{int(units.resolved.sum())} of {len(units)}" if units is not None else "n/a",
                "patient × tissue units with ≥ 50 epithelial and enough reference cells",
            ),
        ]
    )

# ------------------------------------------------------------------------------------------------ fig 1: two axes
st.subheader("The two axes in one patient's tumour")
if units is not None and units.resolved.any():
    res = units[units.resolved & (units.tissue_type == "tumour")].sort_values("n_called", ascending=False)
    pick = st.selectbox("Patient × tissue unit (most calls first)", res.unit.tolist())
    r = res.set_index("unit").loc[pick]
    dn, tis = pick.rsplit("|", 1)
    cells = q(
        """SELECT cnv_score, cnv_corr, cnv_status FROM cells WHERE donor_id = ? AND tissue_type = ?
                 AND cnv_score IS NOT NULL""",
        (dn, tis),
    )
    fig = go.Figure()
    for status, color, name in [
        ("not called", GREY, "not called"),
        ("CNV-high", UP_COLOR, "malignant (both axes high)"),
    ]:
        d = cells[cells.cnv_status.astype(str) == status]
        fig.add_trace(
            go.Scattergl(
                x=d.cnv_score,
                y=d.cnv_corr,
                mode="markers",
                name=f"{name} ({len(d):,})",
                marker=dict(color=color, size=5, opacity=0.7),
            )
        )
    fig.add_vline(x=r.threshold_score, line_dash="dash", line_color=INK_2, annotation_text="99th pct of null")
    fig.add_hline(y=r.threshold_corr, line_dash="dash", line_color=INK_2)
    fig.update_xaxes(type="log")
    st.plotly_chart(
        style(
            fig,
            420,
            x="CNV signal (mean squared CNV value, log scale)",
            y="correlation with this tumour's CNV profile",
        )
    )
    figure_legend(
        1,
        f"In this patient's tumour, {r.call_rate:.0%} of epithelial cells exceed both same-patient null thresholds "
        f"({dn}, {r.study}; the list starts with the patient with most calls).",
        "each dot is an epithelial cell of this patient's tumour; dashed lines = 99th percentiles of the "
        "held-out normal reference cells of the same patient.",
        "malignant cells sit top-right: a strong CNV signal that also matches the tumour's clonal pattern.",
        f"{int(r.n_epithelial):,} epithelial cells; {int(r.n_reference):,} reference + {int(r.n_heldout_null):,} "
        f"held-out cells ({r.reference_source}); {int(r.n_samples)} sample(s) pooled.",
        "infercnvpy, window 250 genes, step 10; masked windows excluded; cross-fitted consensus.",
        f"{int(r.n_called):,} cells called ({r.call_rate:.0%}).",
    )

# ------------------------------------------------------------------------------------------------ fig 2: per patient
st.subheader("Malignant cells per patient")
pdn = q("""SELECT study, donor_id, sum("sum") malignant, sum(size) cells FROM cnv_malignant_per_donor
           WHERE tissue_type = 'tumour' GROUP BY 1, 2 ORDER BY study, malignant DESC""")
pdn["label"] = pdn.donor_id.str.replace("LuCA:", "").str.replace("HLCA:", "")
fig = go.Figure(
    go.Bar(
        x=pdn.label,
        y=pdn.malignant,
        marker_color=MALIGNANT_COLOR,
        customdata=np.stack([pdn.study, pdn.cells], axis=1),
        hovertemplate="%{x}<br>%{customdata[0]}<br>%{y:,} malignant of %{customdata[1]:,} cells<extra></extra>",
    )
)
for st_name, dd in pdn.groupby("study", sort=False):
    fig.add_annotation(
        x=dd.label.iloc[len(dd) // 2],
        y=1.02,
        yref="paper",
        text=st_name,
        showarrow=False,
        font=dict(size=10, color=INK_2),
    )
fig.update_xaxes(tickangle=-60, tickfont=dict(size=9))
st.plotly_chart(style(fig, 420, x="tumour patient (grouped by study)", y="malignant cells (linear scale)"))
zero = int((pdn.malignant == 0).sum())
top3 = pdn.nlargest(3, "malignant").malignant.sum() / max(pdn.malignant.sum(), 1)
figure_legend(
    2,
    "How evenly malignant calls are spread across patients.",
    "number of called malignant cells in each tumour patient, grouped by study; patients with none are shown as "
    "empty positions.",
    "a result about malignant cells is only as broad as this distribution: a few tall bars mean a few patients "
    "dominate.",
    f"{len(pdn)} tumour patients, {pdn.study.nunique()} studies; {zero} with no calls (most have too few "
    "epithelial cells — several studies enriched immune cells).",
    "per patient × tissue two-axis calls.",
    f"the three largest patients hold {top3:.0%} of all calls.",
)

# ------------------------------------------------------------------------------------------------ fig 3: arms
st.subheader("Do the calls show LUSC's hallmark chromosome changes?")
arms = q("""SELECT donor_id, CASE WHEN malignant THEN 'malignant' WHEN tissue_type = 'adjacent_normal'
                 THEN 'adjacent-normal epithelium' ELSE 'tumour epithelium, not called' END grp,
                 avg(cnv_chr3q) q3, avg(cnv_chr3p) p3, avg(cnv_chr5p) p5, count(*) n
          FROM cells WHERE cnv_score IS NOT NULL GROUP BY 1, 2 HAVING count(*) >= 20""")
fig = go.Figure()
gcol = {
    "adjacent-normal epithelium": TISSUE_COLOR["adjacent_normal"],
    "tumour epithelium, not called": GREY,
    "malignant": UP_COLOR,
}
for g, color in gcol.items():
    d = arms[arms.grp == g]
    for col, lab in [("q3", "chr3q (gain)"), ("p3", "chr3p (loss)"), ("p5", "chr5p (gain)")]:
        fig.add_trace(
            go.Box(
                x=[lab] * len(d),
                y=d[col],
                name=g,
                marker_color=color,
                legendgroup=g,
                showlegend=col == "q3",
                boxpoints="all",
                jitter=0.3,
                pointpos=0,
            )
        )
fig.update_layout(boxmode="group")
fig.add_hline(y=0, line_color=GREY)
st.plotly_chart(
    style(
        fig,
        420,
        x="chromosome arm (expected LUSC change)",
        y="mean inferCNV value over the arm (> 0 gain, < 0 loss; relative units)",
    )
)
figure_legend(
    3,
    "Called cells show the expected LUSC 3q gain and 3p loss more strongly than uncalled tumour epithelium.",
    "each dot is one patient's mean CNV signal over a chromosome arm, separately for called malignant cells, "
    "tumour epithelium not called, and adjacent-normal epithelium (≥ 20 cells per patient and group).",
    "LUSC typically gains 3q (SOX2, TP63, PIK3CA) and 5p and loses 3p: malignant dots should sit above zero for "
    "3q/5p and below for 3p, while adjacent-normal epithelium stays near zero.",
    f"{arms.donor_id.nunique()} patients.",
    "arm windows from window-centre gene positions (GRCh38 centromeres); not used for calling, so this is an "
    "independent check.",
    "consistent with, but not proof of, tumour identity: the shifts are small and the false-positive check failed.",
)

# ------------------------------------------------------------------------------------------------ mask + authors
c1, c2 = st.columns(2)
with c1:
    if has_table("cnv_window_mask"):
        wm = q(
            "SELECT chromosome, sum(CASE WHEN masked THEN 1 ELSE 0 END) m, count(*) n FROM cnv_window_mask GROUP BY 1"
        )
        wm["order"] = wm.chromosome.str.replace("chr", "").replace({"X": "23", "Y": "24"}).astype(int)
        wm = wm.sort_values("order")
        fig = go.Figure(go.Bar(x=wm.chromosome, y=100 * wm.m / wm.n, marker_color="#8a8984"))
        st.plotly_chart(style(fig, 300, x="chromosome", y="% of windows masked"))
        figure_legend(
            4,
            "Windows where normal epithelium deviates.",
            "share of each chromosome's windows masked.",
            "tall bars = expression programmes that mimic CNVs in normal epithelium.",
            f"{int(wm.m.sum())} of {int(wm.n.sum())} windows.",
            "|mean| > median + 3 MAD across windows.",
            "masking these protects against squamous false positives.",
        )
with c2:
    va = q("SELECT * FROM cnv_vs_author_malignant")
    if not va.empty and {"True", "False"} <= set(va.columns):
        t = va.set_index(va.author_malignant.astype(str))
        tp, fn = int(t.loc["True", "True"]), int(t.loc["True", "False"])
        fp, tn = int(t.loc["False", "True"]), int(t.loc["False", "False"])
        st.markdown("**Agreement with the LuCA authors' labels** (tumour tissue, LuCA cells)")
        st.markdown(
            "| | called malignant | not called |\n|---|---|---|\n"
            f"| author: malignant | {tp:,} | {fn:,} |\n| author: other | {fp:,} | {tn:,} |"
        )
        st.caption(
            f"{tp / max(tp + fp, 1):.0%} of our calls are author-malignant; {tp / max(tp + fn, 1):.0%} of "
            "author-malignant cells are called. The authors' labels are themselves partly inferCNV-based, so they are "
            "a reference point, not ground truth."
        )
if units is not None:
    show_table(units, "Per patient × tissue unit: cells, reference source, thresholds, calls.")
_ = (DOWN_COLOR, n_note)
footer()
