"""Figures + a Markdown summary assembled ONLY from pipeline outputs (no numbers typed by hand).

Every figure is regenerated from data/processed + results/tables, so the report can never disagree with
the analysis. Missing optional inputs produce a visible "not available" line rather than a crash.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from .config import resolve  # noqa: E402

log = logging.getLogger(__name__)
plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 200, "savefig.bbox": "tight", "font.size": 9})


def _read(cfg: dict, name: str) -> pd.DataFrame | None:
    p = resolve(cfg, "tables", name)
    return pd.read_csv(p) if p.exists() else None


def _fig(cfg: dict, name: str) -> Path:
    p = resolve(cfg, "figures", name)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def fig_sc_overview(cfg: dict) -> list[str]:
    import anndata as ad
    import scanpy as sc

    p = Path(cfg["_repo_root"]) / "data/processed/sc_annotated.h5ad"
    if not p.exists():
        return []
    a = ad.read_h5ad(p, backed="r")
    obs = a.obs.copy()
    umap = a.obsm["X_umap"][:]
    a.file.close()
    small = ad.AnnData(obs=obs)
    small.obsm["X_umap"] = umap
    keys = [
        k
        for k in ["atlas", "condition", "harmonised_label", "de_group", "malignant", "cnv_score"]
        if k in obs
    ]
    sc.pl.umap(small, color=keys, ncols=3, frameon=False, show=False, size=2)
    out = _fig(cfg, "01_umap_overview.png")
    plt.savefig(out)
    plt.close("all")
    return [f"![UMAP overview](../figures/{out.name})"]


def fig_integration(cfg: dict) -> list[str]:
    s = _read(cfg, "integration_benchmark.csv")
    if s is None:
        return ["_Integration benchmark not available._"]
    cols = [c for c in ["Bio conservation", "Batch correction", "Total", "knn_label_purity"] if c in s]
    fig, ax = plt.subplots(figsize=(6, 2.8))
    s.set_index("embedding")[cols].plot.barh(ax=ax)
    ax.set_xlabel("score")
    ax.set_title("Integration benchmark (scib-metrics, batch = study)")
    out = _fig(cfg, "02_integration_benchmark.png")
    fig.savefig(out)
    plt.close(fig)
    chosen = s.loc[s.chosen, "embedding"].iloc[0] if "chosen" in s and s.chosen.any() else "n/a"
    return [
        f"Chosen embedding: **{chosen}**",
        "",
        s.round(3).to_markdown(index=False),
        "",
        f"![benchmark](../figures/{out.name})",
    ]


def fig_volcano(cfg: dict) -> list[str]:
    d = _read(cfg, "de_pseudobulk.csv")
    if d is None or d.empty:
        return ["_DE results not available._"]
    groups = d.group.value_counts().index[:6]
    fig, axes = plt.subplots(1, len(groups), figsize=(3 * len(groups), 3), squeeze=False)
    for ax, g in zip(axes[0], groups, strict=False):
        x = d[d.group == g]
        sig = x.padj < cfg["de"]["alpha"]
        ax.scatter(
            x.log2FoldChange,
            -x.pvalue.clip(lower=1e-300).map(__import__("math").log10),
            s=2,
            c=sig.map({True: "#c0392b", False: "#bbbbbb"}),
        )
        ax.set_title(f"{g}\n{int(sig.sum())} genes padj<{cfg['de']['alpha']}", fontsize=8)
        ax.set_xlabel("log2FC")
    axes[0][0].set_ylabel("-log10 p")
    out = _fig(cfg, "03_volcano_pseudobulk.png")
    fig.savefig(out)
    plt.close(fig)
    n = (
        d.groupby("group")
        .apply(lambda x: int((x.padj < cfg["de"]["alpha"]).sum()))
        .rename("n_sig")
        .reset_index()
    )
    return [n.to_markdown(index=False), "", f"![volcano](../figures/{out.name})"]


def fig_pathways(cfg: dict) -> list[str]:
    p = _read(cfg, "pathways_progeny.csv")
    if p is None:
        return ["_Pathway activities not available._"]
    m = p.pivot_table(index="group", columns="source", values="score")
    fig, ax = plt.subplots(figsize=(0.45 * m.shape[1] + 2, 0.35 * m.shape[0] + 1.5))
    lim = float(abs(m).max().max())
    im = ax.imshow(m, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(m.shape[1]), m.columns, rotation=90)
    ax.set_yticks(range(m.shape[0]), m.index)
    fig.colorbar(im, ax=ax, label="PROGENy activity (ULM t)")
    ax.set_title("Pathway activity: tumour vs adjacent normal")
    out = _fig(cfg, "04_progeny_heatmap.png")
    fig.savefig(out)
    plt.close(fig)
    return [f"![PROGENy](../figures/{out.name})"]


def fig_spatial(cfg: dict) -> list[str]:
    import anndata as ad
    import scanpy as sc

    p = Path(cfg["_repo_root"]) / "data/processed/visium.h5ad"
    if not p.exists():
        return ["_Spatial results not available._"]
    v = ad.read_h5ad(p)
    cts = list(v.uns.get("c2l_cell_types", []))
    keys = [k for k in ["niche", "dominant_type"] if k in v.obs] + cts[:6]
    sc.pl.spatial(v, color=keys, ncols=4, size=1.4, frameon=False, show=False, img_key="hires")
    out = _fig(cfg, "05_spatial_niches_abundance.png")
    plt.savefig(out)
    plt.close("all")
    lines = [f"![spatial](../figures/{out.name})"]
    n = _read(cfg, "spatial_niches.csv")
    if n is not None:
        lines += [
            "",
            "Mean cell-type proportion per niche:",
            "",
            n.set_index("niche").round(2).T.to_markdown(),
        ]
    return lines


def build_report(cfg: dict) -> None:
    root = Path(cfg["_repo_root"])
    sec: list[str] = [
        "# lusc-atlas v2 — results summary",
        "",
        "_Auto-generated by `lusc-atlas report` from pipeline outputs. Do not edit by hand._",
        "",
    ]

    cohort = _read(cfg, "cohort_donors.csv")
    if cohort is not None:
        sec += [
            "## 1. Cohort",
            "",
            cohort.groupby(["atlas", "condition"])
            .agg(donors=("donor_id", "nunique"), cells=("n_cells", "sum"))
            .reset_index()
            .to_markdown(index=False),
            "",
        ]
    xw = _read(cfg, "merge_study_crosswalk.csv")
    if xw is not None:
        sec += ["### Study crosswalk (deduplication)", "", xw.to_markdown(index=False), ""]
    qc = _read(cfg, "qc_per_sample.csv")
    if qc is not None:
        sec += [
            "## 2. QC",
            "",
            f"{int(qc.n_cells.sum())} cells in → {int(qc.n_kept.sum())} kept across {len(qc)} samples "
            f"({int(qc.n_doublet.sum())} consensus doublets removed).",
            "",
        ]
    sec += ["## 3. Integration", ""] + fig_integration(cfg) + [""]
    mix = _read(cfg, "atlas_mixing_per_label.csv")
    if mix is not None:
        sec += [
            "Atlas mixing within harmonised cell types (higher = better integrated; compare with random expectation):",
            "",
            mix.round(3).to_markdown(index=False),
            "",
        ]
    sec += fig_sc_overview(cfg) + [""]
    cnv = _read(cfg, "cnv_summary.csv")
    if cnv is not None:
        sec += ["## 4. Malignant cells (inferCNV)", "", cnv.round(4).to_markdown(index=False), ""]
    sec += (
        ["## 5. Differential expression (pseudobulk DESeq2)", ""]
        + fig_volcano(cfg)
        + [""]
        + fig_pathways(cfg)
        + [""]
    )
    comp = _read(cfg, "composition_sccoda.csv")
    if comp is not None:
        sec += ["## 6. Composition (scCODA)", "", comp.round(3).to_markdown(index=False), ""]
    tr = _read(cfg, "trajectory_summary.csv")
    if tr is not None:
        sec += [
            "## 7. Trajectory (Palantir pseudotime by state)",
            "",
            tr.round(3).to_markdown(index=False),
            "",
            "_Interpretation: transcriptional distance from healthy basal cells; malignant cells are patient-clonal, so this is not lineage tracing._",
            "",
        ]
    sec += ["## 8. Spatial (Visium, n = 1 section — descriptive)", ""] + fig_spatial(cfg) + [""]
    dccc = _read(cfg, "ccc_differential.csv")
    if dccc is not None:
        sec += [
            "## 9. Differential cell–cell communication (top 20)",
            "",
            dccc.head(20)[
                [c for c in ("source", "target", "ligand_complex", "receptor_complex", "ligand_stat",
                             "receptor_stat", "interaction_stat", "interaction_padj") if c in dccc.columns]
            ].round(3).to_markdown(index=False),
            "",
            f"_Ligand/receptor DESeq2 Wald stats from contrast `{cfg['ccc']['de_contrast']}`; "
            "interaction_stat = their mean (LIANA+ `df_to_lr`)._",
            "",
        ]
    sccc = _read(cfg, "ccc_spatial_bivariate.csv")
    if sccc is not None:
        sec += [
            "## 10. Spatially co-expressed ligand–receptor pairs "
            "(Visium, global bivariate Moran's R, top 15)",
            "",
            sccc.head(15)[
                [c for c in ("interaction", "morans", "morans_pvals", "ligand_props", "receptor_props")
                 if c in sccc.columns]
            ].round(3).to_markdown(index=False),
            "",
        ]

    prov_dir = resolve(cfg, "reports", "provenance")
    if prov_dir.exists():
        rows = [json.loads(p.read_text()) for p in sorted(prov_dir.glob("*.json"))]
        sec += [
            "## Provenance",
            "",
            pd.DataFrame(
                [{k: r[k] for k in ("step", "finished", "runtime_s", "config_sha256")} for r in rows]
            ).to_markdown(index=False),
            "",
        ]

    out = resolve(cfg, "reports", "summary.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(sec))
    log.info("report written: %s", out.relative_to(root))
