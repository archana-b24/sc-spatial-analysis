"""Cell–cell communication: differential LR from the paired DE (gains AND losses), paired tensor programmes
tested within patient, spatial LR co-expression on the LUSC section, per-sample browser with denominators
(DEVELOPMENT_LOG §3.9, §3.12)."""

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from lib.db import has_table, q
from lib.explain import guide
from lib.ui import (
    CONDITION_LABEL,
    DOWN_COLOR,
    GREY,
    TISSUE_COLOR,
    UP_COLOR,
    caveats,
    decisions,
    figure_legend,
    footer,
    inference_badge,
    methods_box,
    page_header,
    short_donor,
    show_table,
    style,
)

_tt = (
    q("SELECT min(padj) p FROM ccc_tensor_condition_test") if has_table("ccc_tensor_condition_test") else None
)
_tensor = (
    "the within-patient test of communication programmes was not run"
    if _tt is None or _tt.p.isna().all()
    else "no communication programme differed between tumour and adjacent tissue in the within-patient test"
    if float(_tt.p.iloc[0]) >= 0.05
    else "the within-patient test found communication programmes that differ between tumour and adjacent tissue"
)
page_header(
    "Cell–cell communication",
    "Which ligand–receptor signals between cell types change in LUSC tumours, within the same patients?",
    part="context",
    part_text=f"Descriptive. Ligand–receptor scores are predictions from mRNA, and {_tensor}. The "
    "differential scores are derived from the paired DE (Evidence), so they add no independent evidence.",
)
guide("ccc")
caveats("ccc")
st.info(
    "All communication here is **predicted from mRNA** of a ligand in a sender cell type and its receptor in a "
    "receiver cell type. It is not measured protein binding; treat single pairs as hypotheses.",
    icon="ℹ️",
)
methods_box(
    "LIANA+",
    "Runs several ligand–receptor scoring methods (CellPhoneDB, CellChat, NATMI, Connectome, log2FC, "
    "SingleCellSignalR) on a consensus resource and aggregates their ranks per sample.",
    "Individual methods agree poorly; a consensus of ranks is more stable (Dimitrov et al. 2022).",
    "one method alone (CellChat or CellPhoneDB); NicheNet (ligand-to-target, different question).",
    "Dimitrov et al. 2022 Nat Commun; Dimitrov et al. 2024 Nat Cell Biol (LIANA+).",
)
decisions("ccc")

# ------------------------------------------------------------------------------------------------ differential
st.subheader("Which signals change between tumour and adjacent tissue?")
inference_badge("primary")
st.markdown(
    "Each ligand–receptor pair is scored from the **paired** pseudobulk DE (Stage 5): the *interaction "
    "statistic* is the mean of the ligand's DE statistic in the sender and the receptor's in the receiver, so it "
    "is high when both go up in tumour and strongly negative when both go down. Only primary-eligible immune "
    "types take part, because only they have paired DE."
)
d = q("""SELECT source, target, ligand_complex, receptor_complex, interaction_stat, interaction_padj,
                ligand_stat, receptor_stat FROM ccc_differential""")
n_sig = int((d.interaction_padj < 0.05).sum())
# run 3 review: a pair is only as strong as its weaker partner's evidence; tag pairs involving a
# descriptive-tier cell type (too few patients/studies) and let the reader hide them
tiers = (
    q(
        """SELECT "group", any_value(tier) tier FROM de_pseudobulk WHERE contrast = 'paired_tumour_vs_adjacent'
       GROUP BY 1"""
    )
    .set_index("group")
    .tier
)
d["both_primary"] = d.source.map(tiers).eq("primary") & d.target.map(tiers).eq("primary")
only_primary = st.toggle("Only pairs whose sender and receiver are both primary-tier cell types", value=True)
if only_primary:
    d = d[d.both_primary]
k = st.slider("Pairs shown per direction", 5, 25, 12, 1)
top = d.nlargest(k, "interaction_stat")
bot = d.nsmallest(k, "interaction_stat")
show = (
    __import__("pandas")
    .concat([top, bot])
    .assign(
        label=lambda x: (
            x.ligand_complex
            + " → "
            + x.receptor_complex
            + "  ("
            + x.source
            + " → "
            + x.target
            + ")"
            + np.where(x.both_primary, "", " · descriptive")
        )
    )
    .sort_values("interaction_stat")
)
fig = go.Figure(
    go.Bar(
        y=show.label,
        x=show.interaction_stat,
        orientation="h",
        marker_color=[UP_COLOR if v > 0 else DOWN_COLOR for v in show.interaction_stat],
        customdata=np.stack([show.ligand_stat, show.receptor_stat, show.interaction_padj], axis=1),
        hovertemplate="%{y}<br>interaction %{x:.2f}<br>ligand %{customdata[0]:.2f} · receptor "
        "%{customdata[1]:.2f}<br>padj %{customdata[2]:.2g}<extra></extra>",
    )
)
fig.add_vline(x=0, line_color=GREY)
st.plotly_chart(
    style(fig, 140 + 22 * len(show), x="interaction statistic (paired DE of ligand and receptor)")
)
figure_legend(
    1,
    "Strongest predicted gains and losses of signalling in tumour.",
    f"the {k} pairs with the highest and the {k} with the lowest interaction statistic"
    + (
        " among pairs whose sender and receiver are both primary-tier types."
        if only_primary
        else "; '· descriptive' marks pairs involving a cell type with too few paired patients or studies."
    ),
    "red = ligand and receptor both higher in tumour; blue = both lower. Hover for the two DE statistics.",
    f"{len(d):,} sender–receiver–pair combinations scored; {n_sig:,} with padj < 0.05; patients as in Stage 5.",
    "LIANA+ df_to_lr on the paired PyDESeq2 statistics (Wald stat; p-values combined over ligand and receptor).",
    "these are the communication hypotheses best supported within patients.",
)
show_table(d.sort_values("interaction_stat", ascending=False).round(3), "All scored interactions.")

# ------------------------------------------------------------------------------------------------ tensor
st.subheader("Communication programmes across patients (Tensor-cell2cell)")
methods_box(
    "Tensor-cell2cell",
    "Stacks each patient-tissue's communication scores into a tensor (patient-tissue × ligand–receptor pair × "
    "sender × receiver) and decomposes it into a few factors; each factor is a programme of pairs used by "
    "particular senders and receivers, with a loading per patient-tissue.",
    "Finds coordinated programmes instead of testing thousands of pairs one by one. Built here on one slice "
    "per patient and tissue for patients with both, with cell types present in both tissues — so a factor cannot "
    "separate tumour from adjacent merely because a cell type exists in only one.",
    "per-pair tests only; NMF on pooled samples (ignores patients).",
    "Armingol et al. 2022 Nat Commun.",
)
if has_table("ccc_tensor_condition_test"):
    t = q("SELECT * FROM ccc_tensor_condition_test ORDER BY factor")
    ctx = q("SELECT * FROM ccc_tensor_contexts")
    factors = [c for c in ctx.columns if c.startswith("Factor")]
    sel = st.selectbox("Factor", factors)
    c1, c2 = st.columns([2, 3])
    with c1:
        long = ctx.assign(
            donor=ctx.element.str.rsplit("|", n=1).str[0], tissue=ctx.element.str.rsplit("|", n=1).str[1]
        )
        fig = go.Figure()
        for dn, dd in long.groupby("donor"):
            dd = dd.set_index("tissue").reindex(["adjacent_normal", "tumour"])
            fig.add_trace(
                go.Scatter(
                    x=["Adjacent", "Tumour"],
                    y=dd[sel],
                    mode="lines+markers",
                    showlegend=False,
                    line=dict(color=UP_COLOR if dd[sel].iloc[1] > dd[sel].iloc[0] else DOWN_COLOR, width=1.4),
                    name=short_donor(dn),
                    hovertemplate=f"{short_donor(dn)}<br>%{{y:.3f}}<extra></extra>",
                )
            )
        st.plotly_chart(style(fig, 360, x="tissue", y=f"{sel} loading"))
        r = t.set_index("factor").loc[sel] if sel in set(t.factor) else None
        if r is not None:
            st.caption(
                f"Paired Wilcoxon: {int(r.n_donors_higher_in_tumour)} of {int(r.n_donors)} patients higher in "
                f"tumour; p = {r.pvalue:.2g}, BH padj = {r.padj:.2g}."
            )
    with c2:
        snd = q("SELECT * FROM ccc_tensor_sender_cells").set_index("element")[sel].sort_values()
        rcv = q("SELECT * FROM ccc_tensor_receiver_cells").set_index("element")[sel].sort_values()
        lr = (
            q("SELECT * FROM ccc_tensor_ligand_receptor_pairs")
            .set_index("element")[sel]
            .nlargest(12)
            .sort_values()
        )
        fig = go.Figure()
        fig.add_trace(go.Bar(y=snd.index, x=snd, orientation="h", name="sender", marker_color="#1c5cab"))
        fig.add_trace(go.Bar(y=rcv.index, x=rcv, orientation="h", name="receiver", marker_color="#86b6ef"))
        fig.update_layout(barmode="group")
        st.plotly_chart(style(fig, 360, x=f"{sel} loading", y="cell type"))
    fig = go.Figure(go.Bar(y=lr.index, x=lr, orientation="h", marker_color="#8a8984"))
    st.plotly_chart(style(fig, 340, x=f"{sel} loading", y="ligand → receptor pair"))
    figure_legend(
        2,
        f"{sel}: who talks to whom, with which pairs, and is it different in tumour?",
        "left: the factor's loading in each patient's adjacent and tumour slice (one line per patient); right: "
        "sender and receiver cell-type loadings; bottom: the 12 highest-loading ligand–receptor pairs.",
        "a programme that differs by tissue has most lines sloping the same way and a small padj.",
        f"{int(t.n_donors.max()) if len(t) else 0} patients (both tissues, every kept cell type present).",
        "tensor rank by elbow of reconstruction error; paired Wilcoxon signed-rank on loadings, BH across factors.",
        "a programme without a tissue difference describes communication shared by tumour and adjacent lung.",
    )
    show_table(t.round(4), "Paired test of each factor's loadings, tumour vs adjacent.")

# ------------------------------------------------------------------------------------------------ spatial
st.subheader("Do ligand and receptor co-occur in the tissue?")
inference_badge("descriptive")
b = q("SELECT interaction, morans FROM ccc_spatial_bivariate ORDER BY morans DESC LIMIT 20").sort_values(
    "morans"
)
fig = go.Figure(go.Bar(y=b.interaction, x=b.morans, orientation="h", marker_color="#8a8984"))
st.plotly_chart(
    style(fig, 520, x="global bivariate Moran's R (spatial co-expression)", y="ligand → receptor")
)
figure_legend(
    3,
    "Ligand–receptor pairs expressed together in neighbouring spots.",
    "the 20 pairs with the highest spatial co-expression across the LUSC Visium section.",
    "R near 1 = where the ligand is high, its receptor is high nearby.",
    "one tissue section (one patient).",
    "LIANA+ bivariate, Gaussian kernel (bandwidth 200 full-resolution image pixels, from Space Ranger coordinates), cosine local score, global Moran's R, 100 permutations.",
    "spatial support for a pair; it does not show which cell types express them.",
)

# ------------------------------------------------------------------------------------------------ browser
st.subheader("Per-sample browser")
st.caption(
    "LIANA+ results per sample, keeping rows with an aggregate magnitude rank ≤ 0.01 (LIANA's consensus rank of "
    "expression strength; lower = stronger; not a percentile). Counts are shown as a fraction of the samples of "
    "each tissue in which both cell types were present."
)
labels = q("SELECT DISTINCT source FROM ccc_by_sample ORDER BY 1").source.tolist()
c1, c2 = st.columns(2)
src = c1.selectbox(
    "Sender",
    labels,
    index=labels.index("Macrophage (monocyte-derived)") if "Macrophage (monocyte-derived)" in labels else 0,
)
tgt = c2.selectbox("Receiver", labels, index=labels.index("CD8 T") if "CD8 T" in labels else 0)
r = q(
    """WITH sc AS (SELECT sample, any_value(tissue_type) tis FROM cells GROUP BY 1),
              present AS (SELECT sc.tis, count(DISTINCT b.sample) n FROM ccc_by_sample b JOIN sc USING (sample)
                          WHERE b.source = ? AND b.target = ? GROUP BY 1)
         SELECT b.ligand_complex || ' → ' || b.receptor_complex AS lr, sc.tis tissue, count(*) n_samples,
                any_value(p.n) n_with_both, median(b.magnitude_rank) median_rank
         FROM ccc_by_sample b JOIN sc USING (sample) JOIN present p ON p.tis = sc.tis
         WHERE b.source = ? AND b.target = ? AND sc.tis IN ('tumour', 'adjacent_normal')
         GROUP BY 1, 2""",
    (src, tgt, src, tgt),
)
if not r.empty:
    r["fraction"] = r.n_samples / r.n_with_both
    w = r.pivot_table(index="lr", columns="tissue", values="fraction", fill_value=0)
    w = w.reindex(columns=["adjacent_normal", "tumour"], fill_value=0)
    w = w.loc[(w.max(axis=1)).sort_values(ascending=False).index[:20]].iloc[::-1]
    fig = go.Figure()
    for tis in ["adjacent_normal", "tumour"]:
        fig.add_trace(
            go.Bar(
                y=w.index,
                x=100 * w[tis],
                orientation="h",
                name=CONDITION_LABEL[tis],
                marker_color=TISSUE_COLOR[tis],
            )
        )
    fig.update_layout(barmode="group")
    st.plotly_chart(
        style(fig, 160 + 24 * len(w), x="% of samples with both cell types where the pair is top-ranked")
    )
    figure_legend(
        4,
        f"{src} → {tgt}: pairs strongly expressed per sample.",
        "for each pair, the share of samples (with both cell types present) where it passes the magnitude-rank "
        "cut, in adjacent-normal and tumour tissue.",
        "fractions, not counts, because tumour samples outnumber adjacent ones.",
        "samples with both types: "
        + ", ".join(
            f"{CONDITION_LABEL.get(t, t)} {int(n)}"
            for t, n in r.groupby("tissue").n_with_both.first().items()
        )
        + ".",
        "LIANA+ rank_aggregate per sample; descriptive (no test).",
        "use Figure 1 for the paired, tested comparison.",
    )
footer()
