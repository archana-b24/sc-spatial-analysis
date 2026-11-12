"""Reader guide on every stage page: What you're looking at / How to read it / What it means / Caveats.

"What" and "How" are curated text (docs/dashboard_content.yaml → `page_guide` table). "What it means" is
COMPUTED here from the database each time the page renders, so an interpretation can never quote a number
the data no longer support (e.g. after the §3.10 fixes are rerun). Each finding function returns short
markdown sentences; numbers are formatted from query results, never typed in.
"""

from __future__ import annotations

import streamlit as st

from .db import q

SEV = {"high": "🔴", "medium": "🟠", "low": "🟡"}


def _one(sql: str, params: tuple | None = None):
    return q(sql, params).iloc[0]


# ------------------------------------------------------------------------------- per-stage findings
def _cohort() -> list[str]:
    pre = _one("SELECT count(DISTINCT donor_id) d, sum(n_cells) c FROM cohort_donors")
    post = _one("SELECT count(DISTINCT donor_id) d, count(*) c, count(DISTINCT study) s FROM cells")
    lusc = _one("""SELECT count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy') n,
                          count(DISTINCT donor_id) FILTER (WHERE condition IN ('tumour', 'adjacent_normal')) ta,
                          count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy' AND atlas = 'HLCA') hlca
                   FROM cohort_donors""")
    paired = _one("""SELECT count(*) n FROM (SELECT donor_id FROM cohort_donors GROUP BY donor_id
                     HAVING bool_or(condition='tumour') AND bool_or(condition='adjacent_normal'))""").n
    return [
        f"HLCA alone holds only **{int(lusc.hlca)}** LUSC donors, which is why LuCA was added: together they "
        f"give **{int(lusc.n)}** LUSC patients ({int(lusc.ta)} with a tumour or adjacent-normal sample; the rest "
        f"only region-unknown tumour tissue), **{int(paired)}** of them with paired tumour + adjacent tissue.",
        f"The selected cohort has **{int(pre.d)}** donors and **{int(pre.c):,}** cells; after QC the atlas has "
        f"**{int(post.d)}** donors and **{int(post.c):,}** cells from **{int(post.s)}** studies. The "
        f"difference of {int(pre.d - post.d)} donors is samples too small for per-sample QC (< 200 cells), "
        "dropped before QC.",
    ]


def _qc() -> list[str]:
    r = _one("""SELECT sum(n_cells) i, sum(n_kept) k, sum(n_doublet) d, sum(n_doublet_scrublet) s,
                       sum(n_doublet_solo) so FROM qc_per_sample""")
    return [
        f"**{r.k / r.i:.1%}** of cells passed QC ({int(r.k):,} of {int(r.i):,}).",
        f"Only **{int(r.d):,}** cells were removed as doublets: SOLO flagged {int(r.so):,} and Scrublet "
        f"{int(r.s):,}, and removal needs both, so doublet removal is effectively Scrublet-only.",
    ]


def _integration() -> list[str]:
    bm = q("SELECT embedding, Total, knn_label_purity, chosen FROM integration_benchmark ORDER BY Total DESC")
    bm["m"] = bm.embedding.str.replace("X_", "", regex=False)
    ch = bm[bm.chosen].iloc[0]
    second = bm[~bm.chosen].iloc[0]
    last = bm.iloc[-1]
    mx = _one("""SELECT count(*) n, count(*) FILTER (WHERE frac_other_atlas_neighbours < expected_if_random) b
                 FROM atlas_mixing_per_label""")
    return [
        f"**{ch.m}** scored highest (scIB total **{ch.Total:.3f}**, next {second.m} {second.Total:.3f}) and keeps "
        f"cell types apart (kNN label purity {ch.knn_label_purity:.3f}).",
        f"Zero-shot foundation-model embedding ({last.m}) scored lowest ({last.Total:.3f}): without "
        "fine-tuning it did not remove study effects.",
        f"In **{int(mx.b)} of {int(mx.n)}** cell types the two atlases mix less than random, partly because the "
        "same label means healthy cells in HLCA and tumour-environment cells in LuCA.",
    ]


def _malignant() -> list[str]:
    cs = q("SELECT cnv_group, n_cells, mean_chr3q FROM cnv_summary").set_index("cnv_group")
    va = q("SELECT * FROM cnv_vs_author_malignant")
    t = va[va.author_malignant.astype(str) == "True"].iloc[0]
    called = int(va["True"].sum())
    out = [
        f"**{int(cs.loc['malignant', 'n_cells']):,}** cells are called malignant. Of these calls, "
        f"**{int(t['True']) / called:.1%}** are also author-labelled malignant (high precision), but only "
        f"**{int(t['True']) / (int(t['True']) + int(t['False'])):.0%}** of author-malignant cells pass the "
        "threshold (low sensitivity): the call is conservative.",
    ]
    if "healthy epithelium" in cs.index:
        ratio = cs.loc["malignant", "mean_chr3q"] / cs.loc["healthy epithelium", "mean_chr3q"]
        out.append(
            f"Mean chr3q CNV score is **{ratio:.1f}×** higher in malignant cells than in healthy epithelium, "
            "the expected LUSC signal."
        )
    d = _one("""SELECT count(*) n, count(*) FILTER (WHERE "sum" >= 10) k FROM cnv_malignant_per_donor
                WHERE condition = 'tumour'""")
    out.append(
        f"Only **{int(d.k)} of {int(d.n)}** tumour donors have ≥ 10 malignant cells, so malignant-cell "
        "results rest on few patients."
    )
    return out


def _de() -> list[str]:
    s = q("""SELECT contrast, "group", count(*) FILTER (WHERE padj < 0.05) n FROM de_pseudobulk
             GROUP BY 1, 2 ORDER BY 1, 3 DESC""")
    p = s[s.contrast == "paired_tumour_vs_adjacent"]
    out = []
    if not p.empty:
        top = p.iloc[0]
        nxt = p[p.group != "Unknown"].iloc[0]
        out.append(
            f"In the paired comparison, **{p.group.nunique()}** cell types were testable. The most changed is "
            f"`{top.group}` ({int(top.n):,} significant genes)"
            + (
                f" — but `Unknown` mixes un-called tumour cells with unmapped cells; the most changed named type "
                f"is **{nxt.group}** ({int(nxt.n):,} genes)."
                if top.group == "Unknown"
                else "."
            )
        )
    g = q("""SELECT gene, log2FoldChange, padj FROM de_pseudobulk WHERE contrast = 'paired_tumour_vs_adjacent'
             AND "group" = 'Macrophage' AND padj < 0.05 AND log2FoldChange > 0 ORDER BY padj LIMIT 3""")
    if not g.empty:
        genes = ", ".join(f"*{r.gene}* (log2FC {r.log2FoldChange:+.2f})" for r in g.itertuples())
        out.append(f"Most significant genes up in tumour macrophages: {genes}.")
    spp1 = q("""SELECT log2FoldChange, padj FROM de_pseudobulk WHERE gene = 'SPP1' AND "group" = 'Macrophage'
                AND contrast = 'paired_tumour_vs_adjacent'""")
    if not spp1.empty:
        r = spp1.iloc[0]
        out.append(
            f"*SPP1* (osteopontin) in macrophages: log2FC **{r.log2FoldChange:+.2f}** (about "
            f"{2**r.log2FoldChange:.0f}× higher in tumour), padj {r.padj:.1e}."
        )
    return out


def _composition() -> list[str]:
    sc = q("""SELECT "Cell Type" ct, "log2-fold change" lfc FROM composition_sccoda
              WHERE Covariate LIKE '%T.tumour]%' AND "Inclusion probability" >= 0.8 ORDER BY lfc""")
    out = []
    if not sc.empty:
        down = sc[sc.lfc < 0]
        up = sc[sc.lfc > 0]
        if len(down):
            out.append(
                f"Credibly **lower** in tumour than adjacent tissue ({len(down)} types): "
                + ", ".join(f"{r.ct} ({r.lfc:+.2f})" for r in down.itertuples())
                + "."
            )
        if len(up):
            out.append(
                "Credibly **higher**: " + ", ".join(f"{r.ct} ({r.lfc:+.2f})" for r in up.itertuples()) + "."
            )
        out.append(
            "These are log2 changes *relative to the reference type `Unknown`*; if `Unknown` itself grows "
            "in tumours (it holds un-called tumour cells), every other type looks lower. Read directions "
            "with care until known issue 2 is fixed."
        )
    mi = _one("SELECT count(*) n, count(*) FILTER (WHERE SpatialFDR < 0.1) s FROM milo_nhoods")
    out.append(
        f"Milo: {int(mi.s)} of {int(mi.n):,} neighbourhoods significant, but it tested the wrong contrast "
        "(known issue 1), so this is not a null result."
    )
    return out


def _trajectory() -> list[str]:
    t = q('SELECT lineage_state, "50%" med, std, count FROM trajectory_summary').set_index("lineage_state")
    if not {"Basal", "Malignant"} <= set(t.index):
        return []
    b, m = t.loc["Basal"], t.loc["Malignant"]
    return [
        f"Median pseudotime is **{m.med:.2f}** for malignant cells vs **{b.med:.2f}** for basal cells "
        f"({int(m['count']):,} vs {int(b['count']):,} cells): malignant cells sit further from healthy basal.",
        f"The spread is narrower for malignant cells (SD {m['std']:.2f} vs {b['std']:.2f}) and the distributions "
        "overlap widely, so pseudotime separates the groups only partly.",
    ]


def _spatial() -> list[str]:
    s = _one("SELECT count(*) n, count(DISTINCT niche) k FROM visium_spots")
    re_ = q('SELECT niche, "Rare epithelial" r FROM spatial_niches ORDER BY r DESC').iloc[0]
    svg = q("SELECT gene, I FROM spatial_svg ORDER BY I DESC LIMIT 3")
    return [
        f"**{int(s.n):,}** spots group into **{int(s.k)}** niches.",
        f"Niche {int(re_.niche)} is **{re_.r:.0%}** 'Rare epithelial' by deconvolution — an artefact of a tiny "
        "reference type absorbing the tumour signal (known issue 3), not biology.",
        "Most spatially clustered genes: "
        + ", ".join(f"*{r.gene}* (I = {r.I:.2f})" for r in svg.itertuples())
        + ".",
    ]


def _ccc() -> list[str]:
    top = q("""SELECT source, target, ligand_complex l, receptor_complex r, interaction_stat s
               FROM ccc_differential ORDER BY interaction_stat DESC LIMIT 3""")
    out = [
        "Strongest predicted gains in tumour: "
        + "; ".join(f"{r.l}→{r.r} ({r.source}→{r.target}, stat {r.s:.1f})" for r in top.itertuples())
        + "."
    ]
    out.append(
        "Tensor factors do not separate tumour from adjacent samples (known issue 4): read them as "
        "sample-level programmes, likely tracking donor and study."
    )
    return out


FINDINGS = {
    "cohort": _cohort,
    "qc": _qc,
    "integration": _integration,
    "malignant": _malignant,
    "de": _de,
    "composition": _composition,
    "trajectory": _trajectory,
    "spatial": _spatial,
    "ccc": _ccc,
}


def findings(stage: str) -> list[str]:
    fn = FINDINGS.get(stage)
    if fn is None:
        return []
    try:
        return fn()
    except Exception as e:  # a guide must never break the page; say what failed instead
        return [f"_(interpretation unavailable: {type(e).__name__})_"]


def guide(stage: str, expanded: bool = True) -> None:
    """The four-part reader guide. Rendered under the page title on every stage page."""
    g = q("SELECT * FROM page_guide WHERE stage = ?", (stage,))
    issues = q("SELECT id, severity, status, title FROM known_issues WHERE stage = ? ORDER BY id", (stage,))
    with st.expander("📖 Guide to this page", expanded=expanded):
        c1, c2 = st.columns(2, gap="large")
        with c1:
            if not g.empty:
                st.markdown(f"**What you're looking at.** {g.iloc[0].what}")
                st.markdown(f"**How to read it.** {g.iloc[0].how}")
        with c2:
            f = findings(stage)
            if f:
                st.markdown("**What it means** *(computed from the data on this page)*")
                st.markdown("\n".join(f"- {s}" for s in f))
            if not issues.empty:
                st.markdown("**Caveats**")
                st.markdown(
                    "\n".join(
                        f"- {SEV.get(r.severity, '⚠️')} Issue {r.id} ({r.status}): {r.title}"
                        for r in issues.itertuples()
                    )
                )
                st.page_link("views/issues.py", label="All known issues", icon="⚠️")
        st.caption("Unfamiliar term? See the glossary.")
