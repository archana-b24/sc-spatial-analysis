"""Malignant-cell identification by inferred copy-number variation (DEVELOPMENT_LOG §3.5).

v1 had no malignant/normal distinction (S7): every epithelial cell from a LUSC donor was treated as
"tumour". v2 uses **infercnvpy** (Python re-implementation of inferCNV; Tirosh 2016, Patel 2014):

1. Reference = non-epithelial compartments (immune, endothelial, stroma) — these are diploid in tumours.
2. Smoothed expression along chromosomes → per-cell CNV profile → per-cell CNV score.
3. Threshold = high quantile of the CNV score in epithelium from **normal donors** (an empirical null).
4. Malignant = epithelial ∧ from a LUSC donor ∧ CNV score above that null.
5. Sanity check: LUSC hallmark chr3q gain (SOX2/TP63 locus) should be enriched in malignant cells.
"""

from __future__ import annotations

import gzip
import logging
import re
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def gene_positions_from_gtf(url: str, cache: Path) -> pd.DataFrame:
    """Gene → (chromosome, start, end) from a pinned GENCODE GTF (cached on disk)."""
    cache.parent.mkdir(parents=True, exist_ok=True)
    if not cache.exists():
        log.info("downloading %s", url)
        urllib.request.urlretrieve(url, cache)
    rows = []
    with gzip.open(cache, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8])
            gid = re.search(r'gene_id "([^".]+)', f[8])
            rows.append(
                (name.group(1) if name else None, gid.group(1) if gid else None, f[0], int(f[3]), int(f[4]))
            )
    df = pd.DataFrame(rows, columns=["gene_name", "gene_id", "chromosome", "start", "end"])
    return df[df.chromosome.str.match(r"^chr([0-9]+|X|Y)$")]


def annotate_positions(adata: ad.AnnData, pos: pd.DataFrame) -> None:
    by_id = pos.drop_duplicates("gene_id").set_index("gene_id")
    by_name = pos.drop_duplicates("gene_name").set_index("gene_name")
    key = adata.var["ensembl_id"].str.split(".").str[0] if "ensembl_id" in adata.var else None
    for col in ["chromosome", "start", "end"]:
        vals = key.map(by_id[col]) if key is not None else pd.Series(np.nan, index=adata.var_names)
        vals = vals.where(vals.notna(), adata.var_names.to_series().map(by_name[col]))
        adata.var[col] = vals.values


def run_infercnv(adata: ad.AnnData, ccfg: dict, n_jobs: int = 16, seed: int = 0) -> pd.DataFrame:
    """inferCNV on epithelium + a stratified diploid reference; results written back to ``adata.obs``.

    Adds ``cnv_score``, ``cnv_chr3q``, ``malignant`` (bool) to ``adata.obs`` (NaN/False for cells that
    were not part of the CNV run, i.e. the reference cells beyond the stratified sample).
    """
    import infercnvpy as cnv

    comp = adata.obs[ccfg["compartment_key"]].astype(str)
    epi = comp == ccfg["epithelial_compartment"]
    rng = np.random.default_rng(seed)
    ref_idx = []
    for c in ccfg["reference_compartments"]:
        idx = np.where(comp == c)[0]
        ref_idx.append(
            idx
            if len(idx) <= ccfg["max_ref_cells_per_compartment"]
            else rng.choice(idx, ccfg["max_ref_cells_per_compartment"], replace=False)
        )
    run_idx = np.sort(np.concatenate([np.where(epi)[0], *ref_idx]))
    sub = adata[run_idx].copy()
    sub.obs["cnv_reference"] = np.where(
        sub.obs[ccfg["compartment_key"]].astype(str).isin(ccfg["reference_compartments"]),
        "reference",
        "query",
    )
    cnv.tl.infercnv(
        sub,
        reference_key="cnv_reference",
        reference_cat="reference",
        window_size=ccfg["window_size"],
        step=ccfg["step"],
        n_jobs=n_jobs,
        layer="lognorm",
    )
    cnv.tl.pca(sub)
    cnv.pp.neighbors(sub)
    cnv.tl.leiden(sub, resolution=ccfg["leiden_resolution"])
    cnv.tl.cnv_score(sub)
    sub.obs["cnv_chr3q"] = chr3q_signal(sub, ccfg["window_size"], ccfg["step"])

    sub_epi = sub.obs[ccfg["compartment_key"]].astype(str) == ccfg["epithelial_compartment"]
    null = sub.obs.loc[sub_epi & (sub.obs["condition"] == "healthy"), "cnv_score"]
    if len(null) < 100:
        raise ValueError(
            f"only {len(null)} healthy epithelial cells for the CNV null — cannot set a threshold"
        )
    thr = float(np.quantile(null, ccfg["score_quantile_normal"]))
    in_tumour = sub.obs["condition"].isin(ccfg["tumour_conditions"])
    sub.obs["malignant"] = (sub_epi & in_tumour & (sub.obs["cnv_score"] > thr)).values
    log.info(
        "CNV threshold %.4f (q%.2f of %d healthy epithelial cells) → %d malignant cells",
        thr,
        ccfg["score_quantile_normal"],
        len(null),
        int(sub.obs["malignant"].sum()),
    )

    for col in ["cnv_score", "cnv_chr3q"]:
        adata.obs[col] = np.nan
        adata.obs.loc[sub.obs_names, col] = sub.obs[col].values
    adata.obs["malignant"] = False
    adata.obs.loc[sub.obs_names, "malignant"] = sub.obs["malignant"].values
    adata.obs["cnv_leiden"] = "not_run"
    adata.obs.loc[sub.obs_names, "cnv_leiden"] = sub.obs["cnv_leiden"].astype(str).values
    adata.uns["cnv_threshold"] = thr
    adata.uns["cnv_chr_pos"] = {k: int(v) for k, v in sub.uns["cnv"]["chr_pos"].items()}
    adata.obsm["X_cnv"] = _scatter_rows(sub.obsm["X_cnv"], run_idx, adata.n_obs)

    grp = np.select(
        [
            sub.obs["malignant"],
            sub_epi & (sub.obs["condition"] == "healthy"),
            sub_epi & (sub.obs["condition"] == "adjacent_normal"),
            sub_epi,
        ],
        [
            "malignant",
            "healthy epithelium",
            "adjacent-normal epithelium",
            "tumour-sample epithelium, CNV-low",
        ],
        default="reference (immune/endothelial/stroma)",
    )
    summary = (
        sub.obs.assign(cnv_group=grp)
        .groupby("cnv_group")
        .agg(
            n_cells=("cnv_score", "size"),
            mean_cnv_score=("cnv_score", "mean"),
            mean_chr3q=("cnv_chr3q", "mean"),
        )
        .reset_index()
    )
    summary["threshold"] = thr
    return summary


def author_label_concordance(adata: ad.AnnData, author_label: str) -> pd.DataFrame:
    """Independent check: our CNV calls vs the LuCA authors' 'malignant cell' annotation (LuCA cells only)."""
    o = adata.obs[(adata.obs["atlas"] == "LuCA") & adata.obs["cnv_score"].notna()]
    if o.empty:
        return pd.DataFrame()
    author = o["cell_type"].astype(str) == author_label
    tab = pd.crosstab(author.rename("author_malignant"), o["malignant"].rename("cnv_malignant"))
    return tab.reset_index()


def _scatter_rows(X, rows: np.ndarray, n: int):
    """Place sub-matrix rows back into an n-row sparse matrix (non-run cells stay empty)."""
    import scipy.sparse as sp

    X = sp.csr_matrix(X)
    out = sp.lil_matrix((n, X.shape[1]), dtype=X.dtype)
    out[rows] = X
    return out.tocsr()


CHR3_CENTROMERE_END_GRCH38 = 93_700_000  # chr3 q-arm starts after the centromere (GRCh38)


def chr3q_signal(adata: ad.AnnData, window_size: int, step: int) -> np.ndarray:
    """Mean CNV signal over chr3 q-arm windows (SOX2/TP63 amplicon region in LUSC).

    infercnvpy slides a window of `window_size` genes in steps of `step` over each chromosome's
    position-sorted genes; window j is centred on gene j*step + window_size//2. Windows whose centre
    gene lies beyond the centromere are q-arm windows.
    """
    chr_pos = adata.uns["cnv"]["chr_pos"]
    chroms = sorted(chr_pos, key=chr_pos.get)
    i3 = chroms.index("chr3")
    w0 = chr_pos["chr3"]
    w1 = chr_pos[chroms[i3 + 1]] if i3 + 1 < len(chroms) else adata.obsm["X_cnv"].shape[1]
    genes = adata.var.loc[adata.var["chromosome"] == "chr3"].sort_values("start")["start"].to_numpy()
    centres = np.minimum(np.arange(w1 - w0) * step + window_size // 2, len(genes) - 1)
    q = np.where(genes[centres] > CHR3_CENTROMERE_END_GRCH38)[0] + w0
    return np.asarray(adata.obsm["X_cnv"][:, q].mean(axis=1)).ravel()
