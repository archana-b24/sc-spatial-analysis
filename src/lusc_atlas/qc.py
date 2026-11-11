"""Quality control and doublet removal (DEVELOPMENT_LOG §3.3).

v1 → v2:
* v1 used one global 98th-percentile gene cutoff and the *tumour* threshold on control cells (bug S2).
  v2 uses median-absolute-deviation (MAD) outlier detection **per sample** (Heumos 2023,
  sc-best-practices), so each sample is judged against its own distribution — necessary when HLCA
  combines 10x 3'/5', nuclei and Seq-Well data with very different depth profiles.
* v1 called doublets on tumour cells only (S6). v2 runs Scrublet and scVI-SOLO on **every** sample and
  removes a cell only when both agree — the HLCA input is already QC'd upstream, so we favour precision.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

from .utils import robust_hvg

log = logging.getLogger(__name__)


def _mad_outlier(x: pd.Series, nmads: float) -> pd.Series:
    med = np.median(x)
    mad = median_abs_deviation(x)
    return (x < med - nmads * mad) | (x > med + nmads * mad)


def compute_qc_metrics(adata: ad.AnnData) -> None:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True
    )


def flag_outliers(adata: ad.AnnData, qcfg: dict, sample_key: str = "sample") -> pd.DataFrame:
    """Per-sample MAD flags. Returns a per-sample summary table (also used in the QC report)."""
    o = adata.obs
    flags = pd.Series(False, index=o.index)
    mt_flags = pd.Series(False, index=o.index)
    for _, idx in o.groupby(sample_key, observed=True).groups.items():
        s = o.loc[idx]
        flags.loc[idx] = (
            _mad_outlier(s["log1p_total_counts"], qcfg["mad_log_counts"])
            | _mad_outlier(s["log1p_n_genes_by_counts"], qcfg["mad_log_genes"])
            | _mad_outlier(s["pct_counts_in_top_20_genes"], qcfg["mad_pct_top20"])
        )
        mt_flags.loc[idx] = _mad_outlier(s["pct_counts_mt"], qcfg["mad_pct_mt"]) | (
            s["pct_counts_mt"] > qcfg["max_pct_mt"]
        )
    adata.obs["qc_outlier"] = flags.values
    adata.obs["qc_mt_outlier"] = mt_flags.values
    adata.obs["qc_low_genes"] = (adata.obs["n_genes_by_counts"] < qcfg["min_genes"]).values
    summary = (
        adata.obs.groupby(sample_key, observed=True)
        .agg(
            n_cells=("qc_outlier", "size"),
            n_outlier=("qc_outlier", "sum"),
            n_mt_outlier=("qc_mt_outlier", "sum"),
            n_low_genes=("qc_low_genes", "sum"),
            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),
            median_pct_mt=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    return summary


def call_doublets(adata: ad.AnnData, qcfg: dict, sample_key: str = "sample", seed: int = 0) -> None:
    """Scrublet (per sample) + scVI-SOLO (per sample) → `doublet` = both methods agree."""
    dcfg = qcfg["doublets"]
    adata.obs["doublet_scrublet"] = False
    adata.obs["doublet_solo"] = False
    if dcfg["scrublet"]:
        # scrublet reads .X, which still holds raw counts at this point (normalisation happens later)
        sc.pp.scrublet(adata, batch_key=sample_key, random_state=seed)
        adata.obs["doublet_scrublet"] = adata.obs["predicted_doublet"].astype(bool)
    if dcfg["solo"]:
        import scvi

        tmp = adata.copy()
        tmp.X = tmp.layers["counts"]
        robust_hvg(tmp, 2000, layer=None, batch_key="study", subset=True)
        scvi.model.SCVI.setup_anndata(tmp, batch_key=sample_key)
        vae = scvi.model.SCVI(tmp)
        vae.train(max_epochs=dcfg["solo_max_epochs"], early_stopping=True)
        calls = []
        for b in tmp.obs[sample_key].cat.categories:
            solo = scvi.external.SOLO.from_scvi_model(vae, restrict_to_batch=b)
            solo.train(max_epochs=dcfg["solo_max_epochs"], early_stopping=True)
            pred = solo.predict(soft=False)
            calls.append(pred)
        pred = pd.concat(calls)
        adata.obs["doublet_solo"] = (pred.reindex(adata.obs_names) == "doublet").fillna(False).values
    if dcfg["scrublet"] and dcfg["solo"]:
        adata.obs["doublet"] = adata.obs["doublet_scrublet"] & adata.obs["doublet_solo"]
    else:
        adata.obs["doublet"] = adata.obs["doublet_scrublet"] | adata.obs["doublet_solo"]


def run_qc(
    adata: ad.AnnData, qcfg: dict, seed: int = 0, sample_key: str = "sample"
) -> tuple[ad.AnnData, pd.DataFrame]:
    """Full QC step. Returns filtered AnnData (raw counts in .X and layers['counts']) + per-sample table."""
    compute_qc_metrics(adata)
    # Drop samples too small for per-sample statistics *before* doublet models are trained
    counts = adata.obs[sample_key].value_counts()
    small = counts[counts < qcfg["min_cells_per_sample"]].index
    if len(small):
        log.info("dropping %d samples with < %d cells", len(small), qcfg["min_cells_per_sample"])
        adata = adata[~adata.obs[sample_key].isin(small)].copy()
    adata.obs[sample_key] = adata.obs[sample_key].cat.remove_unused_categories()
    summary = flag_outliers(adata, qcfg, sample_key)
    call_doublets(adata, qcfg, sample_key, seed)
    dsum = (
        adata.obs.groupby(sample_key, observed=True)
        .agg(
            n_doublet_scrublet=("doublet_scrublet", "sum"),
            n_doublet_solo=("doublet_solo", "sum"),
            n_doublet=("doublet", "sum"),
        )
        .reset_index()
    )
    summary = summary.merge(dsum, on=sample_key)
    keep = ~(
        adata.obs["qc_outlier"]
        | adata.obs["qc_mt_outlier"]
        | adata.obs["qc_low_genes"]
        | adata.obs["doublet"]
    )
    adata = adata[keep.values].copy()
    sc.pp.filter_genes(adata, min_cells=qcfg["min_cells_per_gene"])
    summary["n_kept"] = (
        summary[sample_key].astype(str).map(adata.obs[sample_key].astype(str).value_counts()).fillna(0).astype(int)
    )
    log.info("QC kept %d / %d cells", adata.n_obs, int(summary["n_cells"].sum()))
    return adata, summary
