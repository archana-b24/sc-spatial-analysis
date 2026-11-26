"""Quality control: what can go wrong inside a droplet, which checks catch it, why the thresholds are set
per sample, how the two doublet detectors work, and how much each sample lost. Every number is computed from the
DB; thresholds come from the config stored in `meta`."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import config, has_table, q
from lib.explain import guide
from lib.ui import (
    ATLAS_COLOR,
    COMPARTMENT_COLOR,
    CONDITION_COLOR,
    CONDITION_LABEL,
    CONDITION_ORDER,
    GREY,
    INK_2,
    caveats,
    decisions,
    figure_legend,
    footer,
    page_header,
    series_color,
    show_table,
    style,
)

page_header(
    "Quality control",
    "Which measurements are not one intact, living cell, and how can they be removed fairly when the samples come "
    "from 20 studies (before QC) sequenced at very different depths?",
    part="built",
)
guide("qc")
caveats("qc")

cfg = config().get("qc", {})
nmad_c, nmad_g = cfg.get("mad_log_counts", "?"), cfg.get("mad_log_genes", "?")
nmad_t, nmad_m = cfg.get("mad_pct_top20", "?"), cfg.get("mad_pct_mt", "?")
max_mt, min_genes = cfg.get("max_pct_mt", "?"), cfg.get("min_genes", "?")
min_cells = cfg.get("min_cells_per_sample", "?")
min_cells_gene = cfg.get("min_cells_per_gene", "?")

qc = q("""SELECT qc.*, s.condition, s.assay, s.atlas, s.study
          FROM qc_per_sample qc
          JOIN (SELECT sample, any_value(condition) condition, any_value(assay) assay,
                       any_value(atlas) atlas, any_value(study) study
                FROM cells GROUP BY 1) s USING (sample)""")
qc["pct_removed"] = 100 * (1 - qc.n_kept / qc.n_cells)
qc["pct_scrublet"] = 100 * qc.n_doublet_scrublet / qc.n_cells
qc["pct_mt_flag"] = 100 * qc.n_mt_outlier / qc.n_cells
mito_floor = cfg.get("mad_floor_pct_mt", 0)
selected = int(q("SELECT sum(n_cells) n FROM cohort_donors").n.iloc[0])
donors_sel = int(q("SELECT count(DISTINCT donor_id) n FROM cohort_donors").n.iloc[0])
donors_post = int(q("SELECT count(DISTINCT donor_id) n FROM cells").n.iloc[0])
cells_in, kept = int(qc.n_cells.sum()), int(qc.n_kept.sum())
small_cells = selected - cells_in
n_out, n_low, n_mt, n_dbl = (
    int(qc[c].sum()) for c in ("n_outlier", "n_low_genes", "n_mt_outlier", "n_doublet")
)
n_scr = int(qc.n_doublet_scrublet.sum())
removed = cells_in - kept

# ---------------------------------------------------------------- QC in one picture
st.subheader("Quality control in one picture")
st.markdown(
    f"**{selected:,} cells** were selected in stage 1. Samples too small to judge on their own were set aside "
    f"first, then every remaining cell was tested against its own sample. **{kept:,} cells "
    f"({kept / cells_in:.1%} of those tested)** passed and go on to every later stage."
)
st.graphviz_chart(
    f"""digraph {{
    rankdir=LR; nodesep=0.25; ranksep=0.4;
    node [shape=box, style="rounded", fontname="Helvetica", fontsize=10, margin="0.14,0.07"];
    edge [color="#888", fontname="Helvetica", fontsize=9, fontcolor="#888"];
    sel [label="Selected in stage 1\\n{selected:,} cells\\n{donors_sel} donors"];
    tested [label="Tested\\n{cells_in:,} cells\\n{len(qc)} samples"];
    flags [label="Flagged (a cell can fail more than one check)\\n{n_out:,} unusual size or complexity\\n{n_low:,} fewer than {min_genes} genes\\n{n_dbl:,} doublets (Scrublet)\\n{n_mt:,} high mitochondrial share", color="#c23b3a", penwidth=2];
    kept [label="Kept\\n{kept:,} cells\\n{donors_post} donors", color="#1c5cab", penwidth=2];
    sel -> tested [label=" −{small_cells:,} cells: samples with\\n fewer than {min_cells} cells ({donors_sel - donors_post} donors lost)"];
    tested -> flags [label=" −{removed:,} cells"];
    tested -> kept [label=" {kept / cells_in:.1%} pass"];
    }}""",
    width="stretch",
)

# ---------------------------------------------------------------- what goes wrong in a droplet
st.subheader("What can go wrong when a cell is captured")
st.markdown(
    "In droplet single-cell sequencing (10x Genomics and similar), each cell is meant to land in its own tiny "
    "droplet with a uniquely barcoded bead, so every RNA molecule read back can be traced to one cell. Three "
    "things break that promise, and each leaves a recognisable fingerprint in the data."
)
FAILS = [
    (
        "Damaged or dying cell",
        "The cell membrane breaks during tissue dissociation and the RNA in the cytoplasm leaks out. RNA "
        "protected inside mitochondria stays behind.",
        "Few molecules and few genes, and an unusually high share of reads from mitochondrial genes.",
        "Would look like a separate, low-quality 'cell type' and distort gene comparisons.",
    ),
    (
        "Empty droplet or debris",
        "A droplet captured no cell, only free-floating RNA from the liquid around the cells (ambient RNA) or a "
        "fragment of one.",
        "Very few genes, with most reads coming from a handful of highly abundant transcripts.",
        "Would add noise and blur the boundaries between real cell types.",
    ),
    (
        "Doublet",
        "Two cells fell into the same droplet and share one barcode, so they are read as one cell.",
        "About twice the usual number of molecules, and marker genes of two different cell types at once "
        "(for example a T-cell and a macrophage).",
        "Would create fake 'hybrid' cell types and fake intermediate states between real ones.",
    ),
]
cols = st.columns(3, gap="medium")
for col, (title, what, sign, harm) in zip(cols, FAILS, strict=True):
    with col, st.container(border=True, height="stretch"):
        st.markdown(f"**{title}**")
        st.markdown(f"*What happened.* {what}")
        st.markdown(f"*Fingerprint in the data.* {sign}")
        st.caption(f"If left in: {harm}")

# ---------------------------------------------------------------- the checks
st.subheader("The checks, one by one")
st.markdown(
    "Each cell gets a few numbers that describe it, and each number is compared with the rest of its own "
    "sample (see the next section for why). The table below lists every check, what it catches and how many "
    "cells it flagged."
)
checks = [
    (
        "Total molecules (UMIs) per cell",
        "How much RNA was captured. A UMI is a unique molecular identifier, so this counts molecules, not reads.",
        "Too low: damaged cell or empty droplet. Too high: likely a doublet.",
        f"more than {nmad_c} MADs from the sample median (log scale)",
    ),
    (
        "Genes detected per cell",
        "How many different genes have at least one molecule.",
        "Too low: damaged cell or empty droplet. Too high: likely a doublet.",
        f"more than {nmad_g} MADs from the sample median (log scale)",
    ),
    (
        "Share of molecules in the top 20 genes",
        "How much of the cell's RNA comes from its 20 most abundant genes.",
        "Too high: a low-complexity profile, typical of debris and empty droplets.",
        f"more than {nmad_t} MADs from the sample median",
    ),
]
rows = [f"| {a} | {b} | {c} | {d} |" for a, b, c, d in checks]
st.markdown(
    "| Check | What it measures | What a bad value means | Flagged when |\n|---|---|---|---|\n"
    + "\n".join(rows)
    + f"\n| | | | **{n_out:,} cells failed at least one of the three checks above** |"
    + f"\n| Minimum genes | Fixed floor, the same for every sample | Too little information to identify the cell"
    f" | fewer than {min_genes} genes: **{n_low:,} cells** |"
    + f"\n| Mitochondrial share | Percentage of molecules from mitochondrial genes, computed on each atlas' full"
    f" gene list before merging | High only: damaged or dying cell | more than {nmad_m} MADs (MAD at least"
    f" {mito_floor} percentage point) above the median of the cell's compartment in its sample, or above"
    f" {max_mt}%: **{n_mt:,} cells** |"
    + f"\n| Doublets | Scrublet, explained below | Two cells read as one | Scrublet call: **{n_dbl:,} cells** |"
)
st.caption(
    f"Genes seen in fewer than {min_cells_gene} cells after QC are also dropped, since they carry no usable "
    "signal. Samples with fewer than "
    f"{min_cells} cells are removed before any of these checks, because a median and spread estimated from a "
    "handful of cells are unreliable and the doublet detectors need enough cells to learn from."
)

# ---------------------------------------------------------------- why per sample
st.subheader("Why the cut-offs are set per sample and per compartment")
st.markdown(
    "The 20 studies (before QC) used different technologies and sequenced to very different depths. A cell from a deeply "
    "sequenced 10x v3 sample may have ten times as many molecules as an equally healthy cell from a shallow "
    "Drop-seq sample. A single fixed cut-off would therefore either let poor cells through in deep samples or "
    "throw away good cells from whole shallow studies. Instead, each sample is judged against itself, and "
    "within the sample each cell is compared with cells of the same broad class (epithelial, immune, "
    "endothelial, stromal): neutrophils naturally contain little RNA and epithelial and tumour cells naturally "
    "carry more mitochondrial RNA, so a sample-wide rule would remove them first (Yates et al. 2025, "
    "*Genome Biology*)."
)
left, right = st.columns([2, 3], gap="large")
with left, st.container(border=True):
    st.markdown("**How a per-sample cut-off works (MAD)**")
    st.markdown(
        "1. Take one sample and one measurement, for example molecules per cell (on a log scale, because these "
        "counts are skewed).\n"
        "2. Find the **median**: the typical cell of that sample.\n"
        "3. Measure the **median absolute deviation (MAD)**: the typical distance of a cell from that median. "
        "It is like a standard deviation but is not pulled around by the very outliers we are trying to find.\n"
        f"4. Flag cells more than **{nmad_c} MADs** away from the median, in either direction (for the "
        "mitochondrial share, only above it: a low share is not a sign of damage).\n\n"
        f"{nmad_c} MADs is deliberately permissive: only clear outliers are removed, so real but unusual cells "
        "(large tumour cells, small lymphocytes) are kept. This follows the single-cell best-practice guide "
        "(Heumos et al. 2023, *Nature Reviews Genetics*)."
    )
with right:
    samples = qc.sort_values("median_counts")["sample"].tolist()
    default = qc.iloc[(qc.median_counts - qc.median_counts.median()).abs().argsort().iloc[0]]["sample"]
    pick = st.selectbox(
        "Example sample", samples, index=samples.index(default), help="Sorted from shallowest to deepest."
    )
    x = np.log10(q("SELECT total_counts FROM cells WHERE sample = ?", (pick,)).total_counts.to_numpy())
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    s = qc.set_index("sample").loc[pick]
    fig = go.Figure(go.Histogram(x=10**x, nbinsx=60, marker_color=INK_2, name="cells that passed"))
    fig.update_xaxes(type="log", title="molecules (UMIs) per cell, log scale")
    fig.update_yaxes(title="number of cells")
    for v, lab in [
        (med, "median"),
        (med - nmad_c * mad, f"−{nmad_c} MAD"),
        (med + nmad_c * mad, f"+{nmad_c} MAD"),
    ]:
        fig.add_vline(
            x=10**v,
            line_dash="solid" if lab == "median" else "dash",
            line_color=GREY if lab == "median" else "#c23b3a",
            annotation_text=lab,
            annotation_position="top",
        )
    fig.update_layout(
        showlegend=False, title=dict(text=f"{pick} · {s.assay} · {s.study}", font=dict(size=13))
    )
    st.plotly_chart(style(fig, 330))
    figure_legend(
        1,
        f"In {pick}, the size cut-offs sit {nmad_c} MADs either side of this sample's own median.",
        "molecules per cell (log scale) of the cells that passed QC in one sample; lines = median and the "
        "± MAD cut-offs.",
        f"pick other samples: the cut-offs move with each sample's depth. {int(s.n_cells):,} cells tested, "
        f"{int(s.n_kept):,} kept.",
        "one sample (choose above).",
        "only kept cells are stored in this database, so the lines are recomputed from them and sit very close "
        "to the values used in the pipeline.",
        "per-sample cut-offs are needed because depth differs by technology (Figure 2).",
    )

st.markdown("**The depth spread across all samples.** Each box is one sample, sorted by its typical depth.")
qs = q("""SELECT sample, any_value(assay) assay, quantile_cont(total_counts, [0.05, 0.25, 0.5, 0.75, 0.95]) qq
          FROM cells GROUP BY 1""")
qs[["p5", "p25", "p50", "p75", "p95"]] = list(qs.qq)
qs = qs.sort_values("p50").reset_index(drop=True)
assays = qs.groupby("assay").p50.median().sort_values().index.tolist()
fig = go.Figure()
for i, a in enumerate(assays):
    d = qs[qs.assay == a]
    fig.add_trace(
        go.Box(
            x=d.index,
            q1=d.p25,
            median=d.p50,
            q3=d.p75,
            lowerfence=d.p5,
            upperfence=d.p95,
            name=f"{a} ({len(d)})",
            marker_color=series_color(i),
            line_width=1,
            hovertext=d["sample"],
            hoverinfo="text+y",
        )
    )
fig.update_xaxes(title="samples, sorted by median molecules per cell", showticklabels=False)
fig.update_yaxes(type="log", title="molecules (UMIs) per cell, log scale")
fig.update_layout(boxgap=0.1, legend_title_text="Technology (samples)")
st.plotly_chart(style(fig, 440))
lo, hi = qs.p50.iloc[0], qs.p50.iloc[-1]
figure_legend(
    2,
    f"Typical sequencing depth differs {hi / lo:.0f}-fold between samples, and the order follows technology.",
    "one box per sample: middle half of cells (box) and 5th–95th percentiles (whiskers) of molecules per cell; "
    "colour = technology.",
    f"medians run from {lo:,.0f} to {hi:,.0f} molecules per cell; a single global threshold would cut into the "
    "left-hand samples while barely touching the right-hand ones.",
    f"{len(qs)} samples.",
    "kept cells; quantiles per sample.",
    "this is why every QC cut-off is set within sample (and compartment).",
)

# ---------------------------------------------------------------- mitochondrial share
st.subheader("The mitochondrial filter")
st.markdown(
    "When a cell's membrane breaks, the RNA in its cytoplasm leaks out but the RNA inside its mitochondria stays, "
    "so a damaged cell shows an unusually **high share of mitochondrial reads**. One atlas (LuCA) was published "
    "without mitochondrial genes, so the share has to be computed **before** the two atlases are merged: from "
    "HLCA's full raw counts, and from LuCA's own per-cell value calculated by its authors on their full gene "
    "list. (In runs 1 and 2 it was computed after merging, when no mitochondrial genes were left, so this filter "
    "silently did nothing; known issue 9.)"
)
mt = q("""SELECT atlas, pct_counts_mt FROM cells USING SAMPLE 60000 (reservoir, 7)""")
src = (
    qc.groupby("mt_source").agg(samples=("sample", "size"), flagged=("n_mt_outlier", "sum")).reset_index()
    if "mt_source" in qc
    else None
)
c1, c2 = st.columns([3, 2], gap="large")
with c1:
    fig = go.Figure()
    for a, col in ATLAS_COLOR.items():
        d = mt[mt.atlas == a]
        fig.add_trace(
            go.Histogram(
                x=d.pct_counts_mt,
                name=f"{a} ({len(d):,} cells shown)",
                marker_color=col,
                opacity=0.7,
                xbins=dict(start=0, end=30, size=0.5),
            )
        )
    fig.add_vline(x=max_mt, line_dash="dash", line_color="#c23b3a", annotation_text=f"ceiling {max_mt}%")
    fig.update_layout(barmode="overlay")
    st.plotly_chart(
        style(fig, 330, x="mitochondrial share of the cell's molecules (%)", y="cells (random sample)")
    )
with c2:
    n_rule = qc["mt_rule"].astype(str).value_counts() if "mt_rule" in qc else None
    zero = int(qc["mt_rule"].astype(str).str.contains("MAD 0").sum()) if "mt_rule" in qc else 0
    st.metric(
        "Cells flagged for a high mitochondrial share", f"{n_mt:,}", help="Before overlap with other checks"
    )
    st.metric("Samples with no mitochondrial signal (ceiling only)", f"{zero} of {len(qc)}")
    if src is not None:
        st.dataframe(src.rename(columns={"mt_source": "where the share comes from"}), hide_index=True)
na = q(
    """SELECT study, count(*) cells, count(DISTINCT sample) samples FROM cells
       WHERE pct_counts_mt IS NULL OR mt_source LIKE 'unavailable%' GROUP BY 1 ORDER BY 2 DESC"""
)
if len(na):
    st.warning(
        f"**No mitochondrial filter for {int(na.cells.sum()):,} kept cells "
        f"({na.cells.sum() / kept:.0%}) in {int(na.samples.sum())} samples from {len(na)} studies** "
        f"({', '.join(na.study)}): these LuCA studies publish an all-zero mitochondrial share, so their cells are "
        "judged by the size, gene-count and doublet checks only (known issue 9, partly resolved)."
    )
figure_legend(
    3,
    "Most kept cells have a low mitochondrial share; the filter removes a thin upper tail.",
    f"distribution of the mitochondrial share in a random sample of kept cells, by atlas; dashed line = the "
    f"{max_mt}% ceiling.",
    "most cells sit near zero; the right tail is where damaged cells would be. LuCA's values are capped at 30 % "
    "by its authors, so both atlases had been filtered once already.",
    f"{len(mt):,} cells shown of {kept:,} kept; {len(qc)} samples.",
    f"per sample and compartment, upper tail only: more than {nmad_m} MADs (at least {mito_floor} pp) above the "
    f"median, or above {max_mt}%; samples whose depositors removed mitochondrial reads (MAD 0) get the ceiling "
    "only.",
    "the filter is lenient on purpose — it removes clear outliers without discarding metabolically active "
    "tumour cells.",
)

# ---------------------------------------------------------------- doublets
st.subheader("Doublets: Scrublet, with a guard against its failure mode")
st.markdown(
    "A doublet cannot be recognised from one number. **Scrublet** (Wolock, Lopez & Klein 2019, *Cell Systems*) "
    "makes artificial doublets on purpose by adding together the profiles of two random cells of the sample, "
    "places real and artificial cells on one map, and gives every real cell a score: the share of its nearest "
    "neighbours that are artificial doublets. A threshold on that score is chosen per sample."
)
d1, d2, d3 = st.columns(3, gap="medium")
n_fallback = (
    int(qc["doublet_rule"].astype(str).str.startswith("fallback").sum()) if "doublet_rule" in qc else 0
)
with d1, st.container(border=True, height="stretch"):
    st.markdown("**Why Scrublet alone**")
    st.markdown(
        "Runs 1 and 2 also used SOLO (a neural-network classifier, Bernstein et al. 2020) and removed a cell only "
        "if both agreed. SOLO flagged about 17 % of cells — far above a plausible rate in data the atlas authors "
        "had already cleaned — because its classifier ignores how rare doublets are. In practice removal "
        "followed Scrublet anyway, so SOLO is dropped."
    )
with d2, st.container(border=True, height="stretch"):
    st.markdown("**The guard**")
    st.markdown(
        "Scrublet's automatic threshold can fail when a sample's score histogram has no clear second peak, and "
        "then calls most of the sample (one sample in run 2: 97 %). If more than 25 % of a sample is called, the "
        "threshold is treated as failed and only the expected doublet rate (0.8 % per 1,000 cells) with the "
        "highest scores is called instead."
    )
    st.caption(f"Samples where the guard was used: {n_fallback} of {len(qc)}.")
with d3, st.container(border=True, height="stretch"):
    st.markdown("**An independent marker check**")
    st.markdown(
        "A cell expressing an immune marker (PTPRC) together with epithelial (EPCAM, KRT8/18) or endothelial "
        "(PECAM1) markers is a likely heterotypic doublet. This rate is reported per sample as a sanity check; "
        "it is not used to remove cells (ambient RNA can also produce it)."
    )
    if "frac_heterotypic_markers" in qc:
        st.caption(f"Median per sample: {100 * qc.frac_heterotypic_markers.median():.2f}% of cells.")
st.markdown(f"**{n_dbl:,} cells ({n_dbl / cells_in:.2%})** were removed as doublets.")
fig = go.Figure()
for cond in ["healthy", "adjacent_normal", "tumour"]:
    d = qc[qc.condition.replace({"tumour_region_unknown": "tumour"}) == cond]
    fig.add_trace(
        go.Box(
            y=d.pct_scrublet,
            name=f"{CONDITION_LABEL[cond]} ({len(d)})",
            marker_color=CONDITION_COLOR[cond],
            boxpoints="all",
            jitter=0.4,
            pointpos=0,
            text=d["sample"] + " · " + d.study,
            hovertemplate="%{text}<br>%{y:.2f}% doublets<extra></extra>",
        )
    )
fig.update_layout(showlegend=False)
st.plotly_chart(style(fig, 320, x="tissue (number of samples)", y="% of the sample's cells called doublets"))
figure_legend(
    4,
    "Doublet rate per sample.",
    "each dot is one sample: the share of its cells Scrublet called doublets.",
    "rates are low because both atlases were cleaned by their authors; a sample far above the rest would point "
    "to a threshold failure (handled by the guard).",
    f"{len(qc)} samples.",
    "Scrublet per sample (scanpy implementation), automatic threshold with the failure guard described above.",
    "doublet removal is a small correction here, as expected for re-processed atlas data.",
)

decisions("qc")

# ---------------------------------------------------------------- result per sample
st.subheader("Was QC fair to every group?")
st.markdown(
    "In the original 2024 version, healthy samples were filtered with the tumour threshold and doublets were "
    "only looked for in tumours, so the groups were not cleaned the same way. Here every sample goes through "
    "identical checks. The chart shows the share of cells each sample lost, by sample group."
)
fig = go.Figure()
for cond in CONDITION_ORDER:
    d = qc[qc.condition == cond]
    if d.empty:
        continue
    fig.add_trace(
        go.Box(
            y=d.pct_removed,
            name=f"{CONDITION_LABEL[cond]} ({len(d)})",
            marker_color=CONDITION_COLOR[cond],
            boxpoints="all",
            jitter=0.4,
            pointpos=0,
            text=d["sample"] + " · " + d.study,
            hovertemplate="%{text}<br>%{y:.1f}% removed<extra></extra>",
        )
    )
fig.update_yaxes(title="% of the sample's cells removed by QC", rangemode="tozero")
fig.update_xaxes(title="sample group (number of samples)")
fig.update_layout(showlegend=False)
st.plotly_chart(style(fig, 380))
by = qc.groupby("condition").apply(
    lambda d: 100 * (1 - d.n_kept.sum() / d.n_cells.sum()), include_groups=False
)
st.caption(
    "Each dot is one sample. Overall removal: "
    + "; ".join(f"{CONDITION_LABEL[c]} {by[c]:.1f}%" for c in CONDITION_ORDER if c in by)
    + ". Most samples lose only a few percent, because both atlases had been cleaned once already by their "
    "authors. Hover over a dot to see which sample and study it is."
)

if has_table("qc_per_compartment"):
    pc = q("SELECT * FROM qc_per_compartment")
    pc["pct"] = 100 * (pc.n_outlier + pc.n_mt_outlier + pc.n_low_genes) / pc.n_cells
    fig = go.Figure(
        go.Bar(
            x=pc.compartment,
            y=pc.pct,
            marker_color=[COMPARTMENT_COLOR.get(c, GREY) for c in pc.compartment],
            text=[f"{v:.1f}%" for v in pc.pct],
            textposition="outside",
        )
    )
    st.plotly_chart(style(fig, 300, x="compartment", y="% of cells flagged (any non-doublet check)"))
    figure_legend(
        5,
        "QC flags by compartment.",
        "share of each compartment's cells flagged by the size/complexity, minimum-gene or mitochondrial checks "
        "(a cell can fail several; doublets not included).",
        "similar bars mean the compartment-aware rule treats cell classes evenly; a much taller bar would mean "
        "one class is systematically lost.",
        f"{int(pc.n_cells.sum()):,} cells tested.",
        "MAD rules within sample × compartment.",
        "immune cells are flagged about twice as often as endothelial or stromal cells; no class loses more than "
        "about 10 %.",
    )

table = (
    qc[
        [
            c
            for c in [
                "sample",
                "study",
                "atlas",
                "condition",
                "assay",
                "n_cells",
                "n_kept",
                "pct_removed",
                "n_outlier",
                "n_low_genes",
                "n_mt_outlier",
                "mt_rule",
                "n_doublet_scrublet",
                "doublet_rule",
                "n_doublet",
                "median_counts",
                "median_genes",
                "median_pct_mt",
            ]
            if c in qc
        ]
    ]
    .rename(
        columns={
            "n_cells": "cells tested",
            "n_kept": "cells kept",
            "pct_removed": "% removed",
            "n_outlier": "size/complexity outliers",
            "n_low_genes": f"< {min_genes} genes",
            "n_mt_outlier": "mitochondrial outliers",
            "n_doublet_scrublet": "Scrublet doublets",
            "n_doublet": "doublets removed",
            "mt_rule": "mitochondrial rule used",
            "doublet_rule": "doublet threshold",
            "median_pct_mt": "median mito %",
            "median_counts": "median UMIs",
            "median_genes": "median genes",
        }
    )
    .sort_values("% removed", ascending=False)
)
show_table(table.round(1), "Per-sample QC, most-filtered first.")
footer()
