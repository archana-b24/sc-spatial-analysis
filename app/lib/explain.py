"""Reader guide on every stage page: What you're looking at / How to read it / What it means / Caveats.

"What" and "How" are curated text (docs/dashboard_content.yaml → `page_guide` table). "What it means" is
COMPUTED here from the database each time the page renders, so an interpretation can never quote a number
the data no longer support (e.g. after the §3.10 fixes are rerun). Each finding function returns short
markdown sentences; numbers are formatted from query results, never typed in.
"""

from __future__ import annotations

import streamlit as st

from .db import config, has_table, q
from .ui import fixed_in, stage_issues

SEV = {"high": "🔴", "medium": "🟠", "low": "🟡"}


def _one(sql: str, params: tuple | None = None):
    return q(sql, params).iloc[0]


# ------------------------------------------------------------------------------- per-stage findings
def _cohort() -> list[str]:
    pre = _one("""SELECT count(DISTINCT donor_id) d, sum(n_cells) c,
                         count(DISTINCT donor_id) FILTER (WHERE condition = 'healthy') healthy,
                         count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy') lusc,
                         count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy' AND atlas = 'HLCA') lusc_hlca,
                         count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy' AND atlas = 'LuCA') lusc_luca
                  FROM cohort_donors""")
    post = _one("SELECT count(DISTINCT donor_id) d, count(*) c, count(DISTINCT study) s FROM cells")
    lost = _one("""SELECT count(DISTINCT donor_id) FILTER (WHERE condition = 'healthy') h,
                          count(DISTINCT donor_id) FILTER (WHERE condition <> 'healthy') l
                   FROM cohort_donors WHERE donor_id NOT IN (SELECT DISTINCT donor_id FROM cells)""")
    paired = _one("""SELECT count(*) n FROM (SELECT donor_id FROM cohort_donors GROUP BY donor_id
                     HAVING bool_or(condition='tumour') AND bool_or(condition='adjacent_normal'))""").n
    return [
        f"**{int(pre.d)}** donors entered the study: **{int(pre.healthy)}** healthy donors (all from HLCA) and "
        f"**{int(pre.lusc)}** LUSC patients (**{int(pre.lusc_hlca)}** from HLCA, **{int(pre.lusc_luca)}** from "
        f"LuCA). HLCA alone would give only {int(pre.lusc_hlca)} LUSC patients, which is why LuCA was added.",
        f"**{int(paired)}** of the {int(pre.lusc)} LUSC patients gave both a tumour and an adjacent-normal sample "
        "(paired); they carry the main tumour-versus-normal comparison.",
        f"Quality control (next stage) removed **{int(pre.d - post.d)}** donors ({int(lost.h)} healthy, "
        f"{int(lost.l)} LUSC) whose samples had too few cells, leaving **{int(post.d)}** donors and "
        f"**{int(post.c):,}** cells from **{int(post.s)}** studies.",
    ]


def _qc() -> list[str]:
    r = _one("""SELECT sum(n_cells) i, sum(n_kept) k, sum(n_outlier) o, sum(n_low_genes) lg, sum(n_doublet) d,
                       sum(n_mt_outlier) mt, count(*) n, min(median_counts) lo, max(median_counts) hi
                FROM qc_per_sample""")
    return [
        f"**{r.k / r.i:.1%}** of tested cells passed ({int(r.k):,} of {int(r.i):,}, {int(r.n)} samples).",
        f"{int(r.o):,} cells were unusual for their sample and compartment in size or complexity, {int(r.lg):,} "
        f"had too few genes, **{int(r.mt):,}** had a high mitochondrial share (this filter was inactive in runs 1–2) "
        f"and **{int(r.d):,}** were doublets.",
        f"Typical depth differs **{r.hi / r.lo:.0f}-fold** between samples, which is why cut-offs are per sample.",
    ]


def _integration() -> list[str]:
    bm = q("SELECT embedding, Total, knn_label_purity, chosen FROM integration_benchmark ORDER BY Total DESC")
    bm["m"] = bm.embedding.str.replace("X_", "", regex=False)
    ch = bm[bm.chosen].iloc[0]
    second = bm[~bm.chosen].iloc[0]
    out = [
        f"**{ch.m}** scored highest (scIB total **{ch.Total:.3f}**, next {second.m} {second.Total:.3f}) and keeps "
        f"cell types apart (kNN label purity {ch.knn_label_purity:.3f}). scANVI was scored on cells whose labels "
        "were hidden from it during training.",
    ]
    cols = set(q("SELECT * FROM atlas_mixing_per_label LIMIT 0").columns)
    if "scope" in cols:
        mx = q("""SELECT scope, count(*) n, count(*) FILTER (WHERE frac_other_atlas_neighbours < expected_if_random) b
                  FROM atlas_mixing_per_label GROUP BY 1""").set_index("scope")
        for sc_, r in mx.iterrows():
            out.append(f"Atlas mixing below random in **{int(r.b)} of {int(r.n)}** cell types ({sc_}).")
    return out


def _malignant() -> list[str]:
    cs = q("SELECT * FROM cnv_summary").set_index("cnv_group")
    out = []
    mal = cs.loc["malignant (CNV-high)"] if "malignant (CNV-high)" in cs.index else None
    if mal is not None:
        out.append(
            f"**{int(mal.n_cells):,}** epithelial cells from **{int(mal.n_donors)}** patients are called malignant "
            "(both CNV axes above their own sample's null)."
            + (
                " **Provisional:** the pre-registered false-positive check failed (see Malignant cells)."
                if "acceptance_passed" in cs.columns and not bool(cs["acceptance_passed"].iloc[0])
                else ""
            )
        )
        if "adjacent_call_rate" in cs.columns and cs["adjacent_call_rate"].notna().any():
            ok = bool(cs["adjacent_call_rate_ok"].iloc[0]) if "adjacent_call_rate_ok" in cs.columns else True
            out.append(
                f"False-positive check **{'passed' if ok else 'failed'}**: "
                f"**{cs['adjacent_call_rate'].dropna().iloc[0]:.1%}** of adjacent-normal epithelial cells pass the "
                "same rule (pre-registered limit ≤ 5 %; each axis' threshold lets through 1 % of held-out normal cells by "
                "construction, and a call needs both axes)."
            )
        nc = (
            cs.loc["tumour-tissue epithelium, not called"]
            if "tumour-tissue epithelium, not called" in cs.index
            else None
        )
        if nc is not None:
            out.append(
                f"LUSC hallmark arms, malignant vs not-called tumour epithelium: 3q {mal.mean_cnv_chr3q:+.3f} vs "
                f"{nc.mean_cnv_chr3q:+.3f}, 3p {mal.mean_cnv_chr3p:+.3f} vs {nc.mean_cnv_chr3p:+.3f}, 5p "
                f"{mal.mean_cnv_chr5p:+.3f} vs {nc.mean_cnv_chr5p:+.3f} (gain > 0, loss < 0)."
            )
    va = q("SELECT * FROM cnv_vs_author_malignant")
    if not va.empty and "True" in va.columns:
        t = va[va.author_malignant.astype(str) == "True"]
        if len(t):
            t = t.iloc[0]
            called = int(va["True"].sum())
            out.append(
                f"Agreement with the LuCA authors' labels (not ground truth): {int(t['True']) / max(called, 1):.1%} "
                f"of calls are author-malignant; {int(t['True']) / max(int(t['True']) + int(t['False']), 1):.0%} of "
                "author-malignant cells are called."
            )
    d = _one("""SELECT count(*) n, count(*) FILTER (WHERE "sum" >= 10) k FROM cnv_malignant_per_donor
                WHERE tissue_type = 'tumour'""")
    out.append(f"**{int(d.k)} of {int(d.n)}** tumour patients have ≥ 10 malignant cells.")
    return out


PAIRED = "paired_tumour_vs_adjacent"


def _de() -> list[str]:
    p = q(
        """SELECT "group", any_value(n_donors) d, any_value(n_studies) s, any_value(tier) tier,
                    count(*) FILTER (WHERE padj < 0.05) sig, count(*) FILTER (WHERE padj_thr < 0.05) thr,
                    count(*) FILTER (WHERE robust AND padj_thr_global < 0.05
                                       AND NOT coalesce(ambient_suspect, false)) rob
             FROM de_pseudobulk WHERE contrast = ? GROUP BY 1 ORDER BY rob DESC, sig DESC""",
        (PAIRED,),
    )
    out = []
    if p.empty:
        return out
    prim = p[p.tier == "primary"]
    out.append(
        f"**{len(prim)}** cell types meet the primary bar (≥ 6 paired patients, ≥ 2 studies with ≥ 2 each); "
        f"{len(p) - len(prim)} more are shown as descriptive. All are immune cell types: adjacent tissue holds "
        "too few epithelial, stromal or endothelial cells to pair."
    )
    if len(prim):
        top = prim.iloc[0]
        out.append(
            f"Most robust changes: **{top.group}** with **{int(top.rob):,}** robust genes of "
            f"{int(top.thr):,} significant at > 1.5-fold ({int(top.d)} patients, {int(top.s)} studies). Robust = "
            "> 1.5-fold at a global FDR of 5 %, not contradicted by any study, holds when each study with ≥ 25 % "
            "of the pairs is left out in turn, not a contamination suspect."
        )
        g = q(
            """SELECT gene, log2FoldChange l, ci_low, ci_high FROM de_pseudobulk WHERE contrast = ? AND "group" = ?
                 AND robust AND padj_thr_global < 0.05 AND NOT coalesce(ambient_suspect, false)
               ORDER BY abs(log2FoldChange) DESC LIMIT 4""",
            (PAIRED, top.group),
        )
        if len(g):
            out.append(
                f"Largest robust changes in {top.group}: "
                + ", ".join(
                    f"*{r.gene}* {r.l:+.1f} [{r.ci_low:+.1f}, {r.ci_high:+.1f}] (≈ {2 ** abs(r.l):.0f}-fold "
                    f"{'up' if r.l > 0 else 'down'})"
                    for r in g.itertuples()
                )
                + " (log2FC [95 % CI])."
            )
    if has_table("de_paired_null"):
        pn = q("""SELECT "group", avg(frac_p_lt_001) f, bool_or(flagged) flagged FROM de_paired_null
                  GROUP BY 1""")
        if len(pn):
            flag = float(config().get("de", {}).get("sign_flip_flag_frac", 0.03))
            out.append(
                f"Calibration of the paired test (tumour and adjacent labels swapped in a random half of the "
                f"patients): **{pn.f.mean():.1%}** of genes reach p < 0.01 on average (about 1 % expected; a cell "
                f"type above {flag:.0%} would be flagged); **{int(pn.flagged.sum())} of {len(pn)}** cell types flagged."
            )
    if has_table("de_cross_study_null"):
        n = q("SELECT avg(frac_p_lt_001) f FROM de_cross_study_null").f.iloc[0]
        t = q("""SELECT avg(CASE WHEN pvalue < 0.01 THEN 1.0 ELSE 0 END) f FROM de_pseudobulk
                 WHERE contrast = 'tumour_vs_healthy'""").f.iloc[0]
        w = (
            float(q("SELECT avg(frac_p_lt_001) f FROM de_within_study_null").f.iloc[0])
            if has_table("de_within_study_null")
            else None
        )
        if n is not None and t is not None:
            out.append(
                f"Context: healthy vs healthy, two sets of studies of one chemistry (no tumour; study and cohort "
                f"differences only): **{n:.0%}** of genes have "
                f"p < 0.01; tumour vs healthy: **{t:.0%}**. "
                + "Study and cohort differences alone produce hits on the same scale as the tumour comparison "
                "(group sizes and covariates differ, so compare orders of magnitude, not the exact numbers); "
                + (
                    f"random halves of one study give {w:.1%} on average (more in some cell types), so the unpaired "
                    "test itself is roughly calibrated. "
                    if w is not None
                    else ""
                )
                + "Cross-study "
                "comparisons cannot be read gene by gene."
            )
    return out


def _composition() -> list[str]:
    c = q("SELECT * FROM composition_paired ORDER BY glm_log2fc")
    out = []
    if not c.empty:
        d = int(c.n_donors.iloc[0])
        down = c[c.consistent == "lower in tumour"]
        up = c[c.consistent == "higher in tumour"]
        npair = (
            int(_one("SELECT any_value(n_donors) d FROM milo_nhoods").d) if has_table("milo_nhoods") else d
        )
        extra = f" of {npair}; the rest have < 200 immune cells in one tissue" if npair > d else ""
        down = down.reindex(down.glm_log2fc.abs().sort_values(ascending=False).index)
        up = up.reindex(up.glm_log2fc.abs().sort_values(ascending=False).index)
        ctx = f"{d} paired patients{extra}, three models agree; log2 change in share of immune cells"
        if len(up):
            out.append(
                "**Higher** in tumour than adjacent lung of the same patients: "
                + ", ".join(f"{r.cell_type} ({r.glm_log2fc:+.2f})" for r in up.itertuples())
                + f" ({ctx})."
            )
        if len(down):
            out.append(
                "**Lower** in tumour: "
                + ", ".join(f"{r.cell_type} ({r.glm_log2fc:+.2f})" for r in down.itertuples())
                + f" ({ctx})."
            )
        if not len(up) and not len(down):
            out.append(
                f"No immune cell type changes consistently across the three models ({d} paired patients)."
            )
    mi = _one("""SELECT count(*) n, count(*) FILTER (WHERE SpatialFDR < 0.1) s,
                        count(*) FILTER (WHERE SpatialFDR < 0.1 AND logFC > 0) up, any_value(n_donors) d
                 FROM milo_nhoods""")
    out.append(
        f"Milo (paired, {int(mi.d)} patients): **{int(mi.s):,}** of {int(mi.n):,} neighbourhoods change "
        f"({int(mi.up):,} up, {int(mi.s - mi.up):,} down in tumour)."
    )
    return out


def _trajectory() -> list[str]:
    t = q('SELECT lineage_state, "50%" med, std, count FROM trajectory_summary').set_index("lineage_state")
    if not {"Basal", "Malignant"} <= set(t.index):
        return []
    b, m = t.loc["Basal"], t.loc["Malignant"]
    out = [
        f"Median pseudotime **{m.med:.2f}** for malignant cells vs **{b.med:.2f}** for basal cells "
        f"({int(m['count']):,} vs {int(b['count']):,} cells) — a shift of {m.med - b.med:.2f} on a 0–1 scale that "
        "coincides with the atlas difference between root and endpoint. Malignant calls are provisional.",
    ]
    if has_table("trajectory_per_donor"):
        pdn = q("""SELECT lineage_state, count(DISTINCT donor_id) n FROM trajectory_per_donor WHERE n_cells >= 20
                   GROUP BY 1""").set_index("lineage_state")
        if "Malignant" in pdn.index:
            out.append(
                f"Malignant cells come from **{int(pdn.loc['Malignant', 'n'])}** patients with ≥ 20 cells."
            )
    info = q("SELECT value FROM meta WHERE key = 'info_trajectory_info'")
    if len(info):
        import json

        j = json.loads(info.value.iloc[0])
        if j.get("sensitivity_spearman") is not None:
            out.append(
                f"Recomputed on the label-free scVI map, pseudotime correlates at Spearman ρ = "
                f"**{j['sensitivity_spearman']:.2f}** with the main result"
                + (
                    "."
                    if j["sensitivity_spearman"] > 0.75
                    else ": the ordering depends substantially on the map, so read it as a rough picture."
                )
            )
    out.append(
        "Root (HLCA) and malignant cells (LuCA) come from different atlases: read this as descriptive."
    )
    return out


def _spatial() -> list[str]:
    s = _one("SELECT count(*) n, count(DISTINCT niche) k FROM visium_spots")
    ni = q("SELECT * FROM spatial_niches").set_index("niche")
    svg = q("SELECT gene, I FROM spatial_svg ORDER BY I DESC LIMIT 3")
    comp = "; ".join(
        f"niche {int(n)}: "
        + ", ".join(f"{ct} {v:.0%}" for ct, v in row.sort_values(ascending=False).head(2).items())
        for n, row in ni.iterrows()
    )
    out = [
        f"**{int(s.n):,}** spots group into **{int(s.k)}** niches. Top two cell types per niche "
        f"(mean share of spot): {comp}.",
        "Most spatially clustered genes: "
        + ", ".join(f"*{r.gene}* (I = {r.I:.2f})" for r in svg.itertuples())
        + ".",
    ]
    if has_table("c2l_reference_labels"):
        off = q("SELECT label FROM c2l_reference_labels WHERE action <> 'kept'").label.tolist()
        if off:
            out.append(
                "Left out of the deconvolution reference (too few cells or a mixed bucket): "
                + ", ".join(off)
                + ". One section only, so this is descriptive."
            )
    return out


def _ccc() -> list[str]:
    top = q("""WITH t AS (SELECT "group" g FROM de_pseudobulk
                          WHERE contrast = 'paired_tumour_vs_adjacent' AND tier = 'primary' GROUP BY 1)
               SELECT source, target, ligand_complex l, receptor_complex r, interaction_stat s
               FROM ccc_differential WHERE source IN (SELECT g FROM t) AND target IN (SELECT g FROM t)
               ORDER BY interaction_stat DESC LIMIT 3""")
    out = [
        "Strongest predicted gains in tumour among primary-tier cell types (paired DE of both partners): "
        + "; ".join(f"{r.l}→{r.r} ({r.source}→{r.target}, score {r.s:.1f})" for r in top.itertuples())
        + ". Score = mean Wald z of ligand and receptor in the paired DE (> 0: both higher in tumour)."
    ]
    if has_table("ccc_tensor_condition_test"):
        t = q("SELECT * FROM ccc_tensor_condition_test ORDER BY pvalue")
        if len(t):
            sig = t[t.padj < 0.05]
            out.append(
                f"Communication programmes tested within patient ({int(t.n_donors.max())} patients): "
                + (
                    f"**{len(sig)} of {len(t)}** differ between tumour and adjacent tissue (BH padj < 0.05)."
                    if len(sig)
                    else f"none of {len(t)} differs between tumour and adjacent tissue after BH correction."
                )
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


def guide(stage: str, expanded: bool = True) -> None:  # noqa: ARG001 — `expanded` kept for callers
    """The reader guide under the title of every stage page, always visible: what the stage does, why it is
    needed and how it is done (curated), then key results (computed from the data) and open caveats."""
    g = q("SELECT * FROM page_guide WHERE stage = ?", (stage,))
    issues = stage_issues(stage, "id, severity, status, title")
    with st.container(border=True):
        if not g.empty:
            r = g.iloc[0]
            w, y, h = st.columns(3, gap="medium")
            w.markdown(f"**What this stage does**  \n{r.stage_what}")
            y.markdown(f"**Why it is needed**  \n{r.stage_why}")
            h.markdown(f"**How it is done**  \n{r.stage_how}")
        c1, c2 = st.columns([3, 2], gap="large")
        with c1:
            f = findings(stage)
            if f:
                st.markdown("**Key results** *(computed from the data each time the page loads)*")
                st.markdown("\n".join(f"- {s}" for s in f))
        with c2:
            if not issues.empty:
                done = issues["status"].astype(str).str.startswith("resolved")
                if (~done).any():
                    st.markdown("**Caveats**")
                    st.markdown(
                        "\n".join(
                            f"- {SEV.get(r.severity, '⚠️')} Issue {r.id}: {r.title}"
                            for r in issues[~done].itertuples()
                        )
                    )
                if done.any():
                    st.caption(
                        "; ".join(f"{fixed_in(r.status)}: issue {r.id}" for r in issues[done].itertuples())
                        + "."
                    )
        if not g.empty:
            st.caption(f"**This page shows:** {g.iloc[0].what} **Reading the charts:** {g.iloc[0].how}")
        st.caption("Unfamiliar term? See the Glossary page.")
