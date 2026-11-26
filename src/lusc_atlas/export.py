"""Build the dashboard's serving snapshot: pipeline outputs → one read-only DuckDB file.

Design (docs/ROADMAP.md §6):
* Two layers of "source of truth", each with one job (wording reconciled 2026-11-12, docs/design/README.md):
    - computational source of truth = `results/tables/*.csv` and the `.h5ad`/`.zarr` objects, reproducible
      from code + config. Analyses read and write these.
    - serving source of truth = this DuckDB file. The dashboard and the MCP server
      read ONLY this file at runtime (no csv/h5ad/yaml), enforced by tests/test_dashboard.py.
  The DuckDB is DERIVED: `export` rebuilds it from scratch on every run, so it cannot drift from the
  pipeline that produced it, and it is never edited by hand.
* Every table is small enough for a free Hugging Face Space (target < 500 MB file). Large objects are
  reduced on purpose, and each reduction is recorded in the `tables_info` table so the app can say what
  it is showing:
    - cells: all cells' metadata + UMAP, but cell-level expression only for a stratified subsample and a
      marker gene panel (`dashboard.subsample_cells`, `dashboard.marker_genes`);
    - expression: pseudobulk mean log-expression + fraction expressing per gene × cell type × condition
      for every gene;
    - per-sample CCC: only rows passing LIANA's magnitude-rank cut (`dashboard.ccc_magnitude_rank_max`).
* Curated "why" content (decisions, known issues, page guides, glossary) comes from
  `docs/dashboard_content.yaml`, which cites DEVELOPMENT_LOG sections, so every panel in the app points back
  to the written record. It is loaded here, at build time; the app never reads the YAML.
* `SCHEMA_CONTRACT` lists the tables/columns the app and agent rely on; it is written into the DB
  (`schema_contract`) and checked by tests, so a pipeline change that renames a column fails loudly.

The app only ever opens the file read-only (app/lib/db.py).
"""

from __future__ import annotations

import io
import json
import logging
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

SCHEMA_VERSION = (
    5  # run 3 (§3.12): paired/study-grounded tables, per-sample CNV, paired composition, run label
)
# v4: start page: + page_guide.stage_what/why/how (what the stage is, why, how it is done)
# v3: run 2, §3.11: + cells.de_group, c2l_reference_labels, milo design cols; niches open
# v2 (2026-11-12): + page_guide, glossary, step_stamps, schema_contract, meta provenance

# CSVs loaded 1:1 (table name = file stem). Anything not listed here is not shipped.
TABLE_FILES = [
    "cohort_donors",
    "merge_study_crosswalk",
    "merge_exclusions",
    "merge_gene_space",
    "qc_per_sample",
    "integration_benchmark",
    "atlas_mixing_per_label",
    "scanvi_label_confusion",
    "scarches_hlca_label_confusion",
    "cnv_summary",
    "cnv_malignant_per_donor",
    "cnv_vs_author_malignant",
    "de_pseudobulk",
    "pathways_progeny",
    "pathways_collectri",
    "pathways_hallmark",
    "composition_paired",
    "composition_units",
    "milo_nhoods",
    "trajectory_summary",
    "visium_qc",
    "c2l_abundance",
    "spatial_niches",
    "spatial_colocalisation",
    "spatial_nhood_enrichment",
    "spatial_svg",
    "ccc_differential",
    "ccc_spatial_bivariate",
    "ccc_tensor_contexts",
    "ccc_tensor_ligand-receptor_pairs",
    "ccc_tensor_sender_cells",
    "ccc_tensor_receiver_cells",
]
# Shipped when present (run-2 fixes §3.11; run-3 overhaul §3.12); absent in older runs.
OPTIONAL_TABLE_FILES = [
    "c2l_reference_labels",
    "c2l_proportions",
    "qc_per_compartment",
    "cnv_per_sample",
    "cnv_window_mask",
    "de_per_study",
    "de_power",
    "de_sex_check",
    "de_cross_study_null",
    "de_within_study_null",
    "de_paired_null",
    "adjacent_makeup",
    "ccc_tensor_condition_test",
    "trajectory_per_donor",
    "visium_identity",
]

# Tables and columns the app (app/) and the agent (src/lusc_atlas/agent/) depend on. Changing one of these is
# a breaking change: bump SCHEMA_VERSION and update the consumers.
SCHEMA_CONTRACT: dict[str, list[str]] = {
    "cells": [
        "cell_idx",
        "cell_id",
        "atlas",
        "study",
        "sample",
        "donor_id",
        "condition",
        "harmonised_label",
        "de_group",
        "compartment",
        "malignant",
        "cnv_score",
        "umap_1",
        "umap_2",
        "in_subsample",
    ],
    "cohort_donors": ["atlas", "donor_id", "study", "condition", "n_cells"],
    "qc_per_sample": [
        "sample",
        "n_cells",
        "n_kept",
        "n_outlier",
        "n_doublet",
        "n_doublet_scrublet",
        "n_doublet_solo",
        "median_counts",
        "median_genes",
    ],
    "integration_benchmark": [
        "embedding",
        "Total",
        "Batch correction",
        "Bio conservation",
        "knn_label_purity",
        "benchmark_only",
        "chosen",
    ],
    "atlas_mixing_per_label": ["label", "frac_other_atlas_neighbours", "expected_if_random"],
    "cnv_summary": ["cnv_group", "n_cells", "mean_cnv_score", "mean_cnv_chr3q", "acceptance_passed"],
    "cnv_vs_author_malignant": ["author_malignant", "False", "True"],
    "cnv_malignant_per_donor": ["donor_id", "study", "tissue_type", "sum", "size"],
    "de_pseudobulk": [
        "gene",
        "log2FoldChange",
        "padj",
        "stat",
        "baseMean",
        "contrast",
        "group",
        "n_case",
        "n_control",
        "lfcSE",
        "inference",
        "n_donors",
        "n_studies",
    ],
    "pathways_progeny": ["source", "score", "padj", "contrast", "group"],
    "composition_paired": [
        "compartment",
        "cell_type",
        "n_donors",
        "mean_share_adjacent",
        "mean_share_tumour",
        "glm_log2fc",
        "glm_ci_low",
        "glm_ci_high",
        "glm_padj",
        "propeller_padj",
        "clr_padj",
        "consistent",
    ],
    "composition_units": ["unit", "donor_id", "study", "tissue_type", "compartment", "cell_type", "n_cells"],
    "milo_nhoods": ["logFC", "SpatialFDR", "nhood_annotation", "design", "contrast", "n_samples", "n_donors"],
    "c2l_reference_labels": ["label", "n_cells", "action", "target"],
    "trajectory_summary": ["lineage_state", "count", "50%"],
    "lineage_cells": ["lineage_state", "palantir_pseudotime", "umap_1", "umap_2"],
    "visium_spots": ["spot", "niche", "total_counts", "x_fullres", "y_fullres"],
    "spatial_niches": ["niche"],  # + one column per cell2location reference label (set varies by run)
    "spatial_svg": ["gene", "I"],
    "ccc_differential": [
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "interaction_stat",
        "interaction_padj",
    ],
    "ccc_by_sample": ["sample", "source", "target", "ligand_complex", "receptor_complex", "magnitude_rank"],
    "expr_pseudobulk": ["gene", "cell_type", "condition", "mean_lognorm", "frac_expressing", "n_cells"],
    "cell_expr": ["cell_idx", "gene", "value"],
    "genes": ["gene", "in_panel"],
    "decisions": ["stage", "question", "decision", "why", "evidence", "log_ref"],
    "known_issues": ["id", "stage", "severity", "title", "detail", "status", "log_ref"],
    "page_guide": [
        "stage",
        "order",
        "title",
        "one_line",
        "what",
        "how",
        "stage_what",
        "stage_why",
        "stage_how",
    ],
    "glossary": ["term", "definition"],
    "step_runs": ["step", "finished", "runtime_s", "config_sha256", "packages"],
    "step_stamps": ["step", "stamped_at", "stamp"],
    "tables_info": ["table_name", "n_rows", "source", "note"],
    "meta": ["key", "value"],
}

# Always shown on UMAP gene colouring, whatever the DE results say (compartment / LUSC markers).
CANONICAL_MARKERS = [
    "EPCAM",
    "KRT5",
    "KRT14",
    "KRT6A",
    "TP63",
    "SOX2",
    "KRT17",
    "MKI67",
    "TOP2A",
    "CDKN2A",
    "PIK3CA",
    "SFTPC",
    "SFTPB",
    "AGER",
    "FOXJ1",
    "SCGB1A1",
    "MUC5B",
    "PTPRC",
    "CD3E",
    "CD4",
    "CD8A",
    "NKG7",
    "GNLY",
    "MS4A1",
    "JCHAIN",
    "CD68",
    "MARCO",
    "C1QA",
    "LYZ",
    "CD14",
    "FCGR3A",
    "S100A8",
    "CSF3R",
    "TPSAB1",
    "CLEC4C",
    "CLEC9A",
    "PECAM1",
    "VWF",
    "CA4",
    "PROX1",
    "COL1A1",
    "DCN",
    "ACTA2",
    "MYH11",
    "MSLN",
    "PDCD1",
    "CTLA4",
    "TIGIT",
    "HAVCR2",
    "LAG3",
    "FOXP3",
    "CD274",
    "NECTIN2",
    "VEGFA",
    "PLAU",
    "CD44",
    "SPP1",
    "FN1",
]


# ------------------------------------------------------------------------------------------- helpers
def _git(repo: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=repo, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _read_h5ad_parts(path: Path, obsm: list[str], layer: str | None):
    """Read obs, var names, selected obsm and one layer without loading X or the other layers
    (the annotated object is ~3 GB on disk with three matrices; we need one)."""
    import h5py

    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem

    with h5py.File(path, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        emb = {k: np.asarray(read_elem(f["obsm"][k])) for k in obsm if k in f["obsm"]}
        mat = read_elem(f["layers"][layer]) if layer else None
        uns_spatial = None
        if "uns" in f and "spatial" in f["uns"]:
            uns_spatial = read_elem(f["uns"]["spatial"])
    return obs, var, emb, mat, uns_spatial


def _stratified_subsample(obs: pd.DataFrame, n: int, keys: list[str], seed: int) -> np.ndarray:
    """Positional indices of ≤ n cells, proportional per stratum but with every stratum represented
    (rare types such as pDC or malignant cells would vanish from a uniform sample)."""
    rng = np.random.default_rng(seed)
    if n >= len(obs):
        return np.arange(len(obs))
    groups = obs.groupby(keys, observed=True).indices
    floor = max(1, n // (20 * len(groups)))
    frac = n / len(obs)
    picked = []
    for idx in groups.values():
        k = min(len(idx), max(floor, int(round(len(idx) * frac))))
        picked.append(rng.choice(idx, size=k, replace=False))
    out = np.sort(np.concatenate(picked))
    if len(out) > n:
        out = np.sort(rng.choice(out, size=n, replace=False))
    return out


def _pseudobulk(mat, codes: np.ndarray, n_groups: int):
    """Mean expression and fraction of cells > 0 per group, via one sparse indicator product."""
    import scipy.sparse as sp

    mat = sp.csr_matrix(mat)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(n_groups, len(codes)))
    sizes = np.asarray(ind.sum(axis=1)).ravel()
    sums = np.asarray((ind @ mat).todense())
    nz = mat.copy()
    nz.data = np.ones_like(nz.data)
    counts = np.asarray((ind @ nz).todense())
    with np.errstate(divide="ignore", invalid="ignore"):
        return sums / sizes[:, None], counts / sizes[:, None], sizes


def _long_sparse(mat, row_ids: np.ndarray, col_ids: np.ndarray, row_name: str, col_name: str) -> pd.DataFrame:
    import scipy.sparse as sp

    coo = sp.coo_matrix(mat)
    return pd.DataFrame(
        {row_name: row_ids[coo.row], col_name: col_ids[coo.col], "value": coo.data.astype(np.float32)}
    )


def _png(arr: np.ndarray) -> bytes:
    from PIL import Image

    if arr.dtype != np.uint8:
        arr = (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG", optimize=True)
    return buf.getvalue()


# ------------------------------------------------------------------------------------------- builder
class _DB:
    """Thin wrapper that records what each table holds and how it was reduced."""

    def __init__(self, path: Path):
        import duckdb

        self.con = duckdb.connect(str(path))
        self.info: list[dict] = []

    def put(self, name: str, df: pd.DataFrame, source: str, note: str = "") -> None:
        df = df.copy()
        df.columns = [str(c) for c in df.columns]
        self.con.register("_df", df)
        self.con.execute(f'CREATE OR REPLACE TABLE "{name}" AS SELECT * FROM _df')
        self.con.unregister("_df")
        self.info.append({"table_name": name, "n_rows": len(df), "source": source, "note": note})
        log.info("  %-34s %10d rows  ← %s", name, len(df), source)

    def close(self) -> None:
        self.put("tables_info", pd.DataFrame(self.info), "export")
        self.con.execute("CHECKPOINT")
        self.con.close()


def build_dashboard_db(cfg: dict, out: Path) -> None:
    import yaml

    repo = Path(cfg["_repo_root"])
    dcfg = cfg.get("dashboard", {})
    seed = int(cfg["project"]["seed"])
    tables = repo / cfg["paths"]["tables"]
    reports = repo / cfg["paths"]["reports"]
    t0 = time.time()

    tmp = out.with_name(f".{out.name}.partial")
    tmp.unlink(missing_ok=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    db = _DB(tmp)

    # ---- 1. result tables, 1:1
    for stem in TABLE_FILES + [f for f in OPTIONAL_TABLE_FILES if (tables / f"{f}.csv").exists()]:
        path = tables / f"{stem}.csv"
        if stem not in TABLE_FILES and path.stat().st_size < 3:  # an optional analysis that produced no rows
            log.warning(
                "optional table %s is empty — not exported (page shows its 'not available' note)", stem
            )
            continue
        db.put(stem.replace("-", "_"), pd.read_csv(path, low_memory=False), f"results/tables/{stem}.csv")

    # run 3 final review: a second contamination route the per-gene ratio misses. If a gene's source type
    # itself expands in tumour (paired composition: consistently higher), a modest marker ratio (≥ 3) is
    # enough to explain an apparent rise in another type (e.g. IL2RA from expanding Tregs into macrophages).
    tabs = set(db.con.execute("SELECT table_name FROM information_schema.tables").df().table_name)
    if {"de_pseudobulk", "composition_paired"} <= tabs:
        dcols = set(db.con.execute("SELECT * FROM de_pseudobulk LIMIT 0").df().columns)
        if {"ambient_ratio", "ambient_source", "ambient_suspect"} <= dcols:
            min_ratio = float(cfg.get("de", {}).get("ambient_expanding_source_ratio", 3.0))
            db.con.execute(
                f"""ALTER TABLE de_pseudobulk ADD COLUMN ambient_expanding_source BOOLEAN;
                    UPDATE de_pseudobulk SET ambient_expanding_source =
                      log2FoldChange > 0 AND coalesce(ambient_ratio, 0) >= {min_ratio}
                      AND ambient_source IN (SELECT cell_type FROM composition_paired
                                             WHERE consistent = 'higher in tumour');
                    UPDATE de_pseudobulk SET ambient_suspect = TRUE WHERE ambient_expanding_source;"""
            )

    # per-sample LIANA: keep the informative tail only (483 MB CSV → a few % of rows)
    rank_max = float(dcfg.get("ccc_magnitude_rank_max", 0.01))
    keep_cols = [
        "sample",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "lr_means",
        "cellphone_pvals",
        "lrscore",
        "specificity_rank",
        "magnitude_rank",
    ]
    parts = []
    for chunk in pd.read_csv(tables / "ccc_by_sample.csv", usecols=keep_cols, chunksize=500_000):
        parts.append(chunk[chunk["magnitude_rank"] <= rank_max])
    db.put(
        "ccc_by_sample",
        pd.concat(parts, ignore_index=True),
        "results/tables/ccc_by_sample.csv",
        f"filtered to magnitude_rank ≤ {rank_max} (LIANA consensus rank)",
    )

    # ---- 2. single-cell object: metadata, UMAP, pseudobulk, marker panel on a subsample
    sc_path = repo / "data/processed/sc_annotated.h5ad"
    obs, var, emb, lognorm, _ = _read_h5ad_parts(sc_path, ["X_umap"], "lognorm")
    genes = np.asarray(var.index.astype(str))
    cell_cols = [
        "atlas",
        "study",
        "sample",
        "donor_id",
        "assay",
        "condition",
        "harmonised_label",
        "de_group",
        "analysis_group",
        "tissue_type",
        "compartment",
        "fine_label",
        "scanvi_prediction",
        "malignant",
        "cnv_score",
        "cnv_corr",
        "cnv_status",
        "cnv_chr3p",
        "cnv_chr3q",
        "cnv_chr5p",
        "scanvi_heldout",
        "mt_source",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "doublet_score",
        "leiden",
        "age_years",
        "sex",
    ]
    cells = obs[[c for c in cell_cols if c in obs.columns]].copy()
    for c in cells.columns:
        if isinstance(cells[c].dtype, pd.CategoricalDtype) or cells[c].dtype == object:
            cells[c] = cells[c].astype(str)
    cells.insert(0, "cell_id", obs.index.astype(str))
    cells.insert(0, "cell_idx", np.arange(len(cells), dtype=np.int32))
    cells["umap_1"] = emb["X_umap"][:, 0].astype(np.float32)
    cells["umap_2"] = emb["X_umap"][:, 1].astype(np.float32)
    n_sub = int(dcfg.get("subsample_cells", 50_000))
    sub = _stratified_subsample(obs, n_sub, ["harmonised_label", "condition"], seed)
    cells["in_subsample"] = False
    cells.loc[cells.index[sub], "in_subsample"] = True
    db.put(
        "cells", cells.reset_index(drop=True), str(sc_path.relative_to(repo)), f"subsample = {len(sub)} cells"
    )

    # pseudobulk per harmonised label × condition (all genes)
    grp = obs["harmonised_label"].astype(str) + "\t" + obs["condition"].astype(str)
    codes, uniq = pd.factorize(grp, sort=True)
    mean, frac, sizes = _pseudobulk(lognorm, codes, len(uniq))
    lab = [u.split("\t") for u in uniq]
    gi, gj = np.nonzero(frac > 0)
    expr = pd.DataFrame(
        {
            "gene": genes[gj],
            "cell_type": [lab[i][0] for i in gi],
            "condition": [lab[i][1] for i in gi],
            "mean_lognorm": mean[gi, gj].astype(np.float32),
            "frac_expressing": frac[gi, gj].astype(np.float32),
            "n_cells": sizes[gi].astype(np.int32),
        }
    )
    db.put(
        "expr_pseudobulk", expr, "sc_annotated.h5ad layer lognorm", "rows with frac_expressing = 0 dropped"
    )

    # marker panel: canonical markers + top DE genes per group + top HVGs
    de = pd.read_csv(tables / "de_pseudobulk.csv", low_memory=False)
    n_top = int(dcfg.get("marker_top_de_per_group", 15))
    top_de = (
        de.dropna(subset=["padj"])
        .assign(abs_stat=lambda d: d["stat"].abs())
        .sort_values("abs_stat", ascending=False)
        .groupby(["contrast", "group"])
        .head(n_top)["gene"]
    )
    hvg = var.sort_values("highly_variable_rank").index[: int(dcfg.get("marker_top_hvg", 500))]
    panel = pd.Index(CANONICAL_MARKERS).union(pd.Index(top_de)).union(pd.Index(hvg).astype(str))
    panel = panel[panel.isin(genes)]
    gidx = np.flatnonzero(np.isin(genes, panel))
    import scipy.sparse as sp

    sub_mat = sp.csr_matrix(lognorm)[sub][:, gidx]
    db.put(
        "cell_expr",
        _long_sparse(sub_mat, sub.astype(np.int32), genes[gidx], "cell_idx", "gene"),
        "sc_annotated.h5ad layer lognorm",
        f"{len(gidx)} marker genes × {len(sub)} subsampled cells, non-zeros only",
    )
    db.put("genes", pd.DataFrame({"gene": genes, "in_panel": np.isin(genes, panel)}), "sc_annotated.h5ad var")
    del lognorm, sub_mat

    # ---- 3. trajectory object
    lin_path = repo / "data/processed/lineage.h5ad"
    lobs, _, lemb, _, _ = _read_h5ad_parts(lin_path, ["X_umap"], None)
    lcols = ["lineage_state", "palantir_pseudotime", "palantir_entropy", "ct_pseudotime", "term_states_fwd"]
    lin = lobs[[c for c in lcols + ["condition", "atlas", "donor_id"] if c in lobs.columns]].copy()
    for c in lin.columns:
        if isinstance(lin[c].dtype, pd.CategoricalDtype) or lin[c].dtype == object:
            lin[c] = lin[c].astype(str)
    lin.insert(0, "cell_id", lobs.index.astype(str))
    lin["umap_1"], lin["umap_2"] = lemb["X_umap"][:, 0], lemb["X_umap"][:, 1]
    db.put("lineage_cells", lin.reset_index(drop=True), str(lin_path.relative_to(repo)))

    # ---- 4. Visium: spots, H&E image, SVG / LR gene expression
    vis_path = repo / "data/processed/visium.h5ad"
    vobs, vvar, vemb, _, uns_sp = _read_h5ad_parts(vis_path, ["spatial"], None)
    spots = vobs.copy()
    for c in spots.columns:
        if isinstance(spots[c].dtype, pd.CategoricalDtype) or spots[c].dtype == object:
            spots[c] = spots[c].astype(str)
    spots.insert(0, "spot", vobs.index.astype(str))
    spots["x_fullres"], spots["y_fullres"] = vemb["spatial"][:, 0], vemb["spatial"][:, 1]
    db.put("visium_spots", spots.reset_index(drop=True), str(vis_path.relative_to(repo)))

    lib_id, lib = next(iter(uns_sp.items()))
    sf = lib["scalefactors"]
    img_rows = [
        {
            "library_id": lib_id,
            "resolution": res,
            "png": _png(lib["images"][res]),
            "scalef": float(sf[f"tissue_{res}_scalef"]),
            "spot_diameter_fullres": float(sf["spot_diameter_fullres"]),
        }
        for res in ("hires", "lowres")
        if res in lib["images"]
    ]
    db.put("visium_images", pd.DataFrame(img_rows), "visium.h5ad uns/spatial")

    import h5py

    try:
        from anndata.io import read_elem
    except ImportError:
        from anndata.experimental import read_elem
    svg = pd.read_csv(tables / "spatial_svg.csv").sort_values("I", ascending=False)
    lr = pd.read_csv(tables / "ccc_spatial_bivariate.csv")
    vgenes_want = (
        pd.Index(svg["gene"].head(int(dcfg.get("visium_top_svg", 200))))
        .union(pd.Index(lr["ligand"]).union(pd.Index(lr["receptor"])))
        .union(pd.Index(CANONICAL_MARKERS))
    )
    vgenes = np.asarray(vvar.index.astype(str))
    vidx = np.flatnonzero(np.isin(vgenes, vgenes_want))
    with h5py.File(vis_path, "r") as f:
        vX = sp.csr_matrix(read_elem(f["X"]))[:, vidx]  # X = log-normalised (DEVELOPMENT_LOG §3.8)
    db.put(
        "spot_expr",
        _long_sparse(vX, np.asarray(spots["spot"]), vgenes[vidx], "spot", "gene"),
        "visium.h5ad X (log-normalised)",
        f"{len(vidx)} genes: top SVGs by Moran's I + spatial LR genes + markers",
    )

    # ---- 5. JSON reports, provenance, summary
    prov = []
    for p in sorted((reports / "provenance").glob("*.json")):
        r = json.loads(p.read_text())
        prov.append(
            {
                "step": r.get("step", p.stem),
                "finished": r.get("finished"),
                "runtime_s": r.get("runtime_s"),
                "python": r.get("python"),
                "config_sha256": r.get("config_sha256"),
                "packages": json.dumps(r.get("packages", {}), sort_keys=True),
            }
        )
    db.put("step_runs", pd.DataFrame(prov), "results/reports/provenance/*.json")

    info = {}
    for p in [
        reports / "integration_info.json",
        reports / "fm_embed_info.json",
        tables / "trajectory_info.json",
    ]:
        if p.exists():
            info[p.stem] = p.read_text()
    content = yaml.safe_load((repo / "docs/dashboard_content.yaml").read_text())
    db.put("decisions", pd.DataFrame(content["decisions"]), "docs/dashboard_content.yaml")
    db.put("known_issues", pd.DataFrame(content["known_issues"]), "docs/dashboard_content.yaml")
    db.put("page_guide", pd.DataFrame(content["page_guide"]), "docs/dashboard_content.yaml")
    db.put("glossary", pd.DataFrame(content["glossary"]), "docs/dashboard_content.yaml")

    # completion stamps written by scripts/run_all.py (logs/run_all/.done/<step>): which steps had finished
    # cleanly when this snapshot was built
    stamp_dir = repo / "logs" / "run_all" / ".done"
    stamps = (
        [
            {
                "step": p.name,
                "stamped_at": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(p.stat().st_mtime)),
                "stamp": p.read_text().strip(),
            }
            for p in sorted(stamp_dir.glob("*"))
            if p.is_file()
        ]
        if stamp_dir.exists()
        else []
    )
    db.put(
        "step_stamps",
        pd.DataFrame(stamps, columns=["step", "stamped_at", "stamp"]),
        "logs/run_all/.done/*",
        "empty when the pipeline was run by Snakemake (it tracks completion itself)",
    )
    db.put(
        "schema_contract",
        pd.DataFrame(
            [(t, c) for t, cols in SCHEMA_CONTRACT.items() for c in cols],
            columns=["table_name", "column_name"],
        ),
        "src/lusc_atlas/export.py SCHEMA_CONTRACT",
    )

    import hashlib

    cfg_text = (repo / "config/config.yaml").read_text()
    meta = {
        "schema_version": str(SCHEMA_VERSION),
        # which pipeline run produced the numbers (shown in the app footer; §3.12)
        "run_id": str(dcfg.get("run_id", "")),
        "run_label": str(dcfg.get("run_label", "")),
        "role": "serving source of truth for dashboard + agent; derived from pipeline outputs by `export`",
        "config_sha256": hashlib.sha256(cfg_text.encode()).hexdigest(),
        "content_sha256": hashlib.sha256((repo / "docs/dashboard_content.yaml").read_bytes()).hexdigest(),
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_sha": _git(repo, "rev-parse", "HEAD"),
        "git_branch": _git(repo, "rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": str(bool(_git(repo, "status", "--porcelain", "--untracked-files=no"))),
        "config_yaml": cfg_text,
        "summary_md": (reports / "summary.md").read_text(),
        **{f"info_{k}": v for k, v in info.items()},
    }
    db.put("meta", pd.DataFrame({"key": list(meta), "value": list(meta.values())}), "export")
    db.close()
    os.replace(tmp, out)
    log.info("dashboard DB → %s (%.0f MB, %.0fs)", out, out.stat().st_size / 1e6, time.time() - t0)
