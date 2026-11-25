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


def _mad_outlier(x: pd.Series, nmads: float, upper_only: bool = False) -> pd.Series:
    med = np.median(x)
    mad = median_abs_deviation(x)
    hi = x > med + nmads * mad
    return hi if upper_only else (x < med - nmads * mad) | hi


def compute_qc_metrics(adata: ad.AnnData) -> None:
    """Per-cell QC metrics on the merged gene space.

    ``pct_counts_mt`` is NOT computed here when ``select`` already provides it (§3.10 #9): the merged gene
    space has no MT- genes (LuCA ships none), so recomputing would overwrite the real value with 0 %.
    """
    have_mt = "pct_counts_mt" in adata.obs
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", regex=True)
    if not have_mt and not adata.var["mt"].any():
        raise ValueError(
            "no mitochondrial fraction: select must provide obs['pct_counts_mt'] (MT- genes are absent "
            "from the merged gene space) — see DEVELOPMENT_LOG §3.10 #9"
        )
    qc_vars = ["ribo", "hb"] + ([] if have_mt else ["mt"])
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=[20], log1p=True, inplace=True)


def mito_flags(
    pct: pd.Series, nmads: float, ceiling: float | None, mad_floor: float = 0.0
) -> tuple[pd.Series, str]:
    """Upper-tail mito outliers for ONE group of cells (§3.12).

    Only high mito % indicates damaged cells (cytoplasmic mRNA leaks out, mitochondrial transcripts stay),
    so the
    MAD rule is one-sided. The MAD is floored at ``mad_floor`` percentage points, so a tight distribution
    does not
    turn ordinary cells into outliers. A group with no mito signal (MAD 0: depositors removed MT reads, or
    the value
    is unavailable) gets the hard ceiling only. Returns (flags, rule applied).
    """
    v = pd.to_numeric(pct, errors="coerce")
    if v.isna().all():
        return pd.Series(False, index=pct.index), "none (no mito data)"
    v = v.fillna(0.0)
    hard = v > ceiling if ceiling is not None else pd.Series(False, index=v.index)
    mad = median_abs_deviation(v)
    if mad == 0:
        return hard, "ceiling only (MAD 0)"
    upper = np.median(v) + nmads * max(mad, mad_floor)
    return (v > upper) | hard, "upper MAD + ceiling"


def flag_outliers(
    adata: ad.AnnData, qcfg: dict, sample_key: str = "sample", stratum_key: str | None = "compartment"
) -> pd.DataFrame:
    """Per-sample MAD flags. Returns a per-sample summary table (also used in the QC report).

    Mito outliers are judged within sample × ``stratum_key`` (broad compartment) when that group has at least
    ``mito_min_group`` cells, else within the sample: epithelial and tumour cells legitimately carry more
    mitochondrial RNA than lymphocytes, and a sample-wide MAD removes them preferentially (Yates 2025).
    """
    o = adata.obs
    flags = pd.Series(False, index=o.index)
    mt_flags = pd.Series(False, index=o.index)
    rules = {}
    min_group = qcfg.get("mito_min_group", 50)
    n_nan = int(pd.to_numeric(o["pct_counts_mt"], errors="coerce").isna().sum())
    if n_nan:
        log.info("pct_counts_mt missing for %d cells: judged by the other rules only", n_nan)
    for smp, idx in o.groupby(sample_key, observed=True).groups.items():
        s = o.loc[idx]
        # §3.12: every MAD rule is judged within sample × compartment (small groups pooled per sample). A
        # sample-wide lower cut on counts removes cell types with intrinsically low RNA content first
        # (neutrophils), and a sample-wide mito cut removes metabolically active epithelium first.
        if stratum_key and stratum_key in s:
            strata = s[stratum_key].astype(str)
            small = strata.map(strata.value_counts()) < min_group
            strata = strata.where(~small, "_pooled_small")
        else:
            strata = pd.Series("_all", index=s.index)
        used = set()
        for _, gidx in s.groupby(strata, observed=True).groups.items():
            g = s.loc[gidx]
            flags.loc[gidx] = (
                _mad_outlier(g["log1p_total_counts"], qcfg["mad_log_counts"])
                | _mad_outlier(g["log1p_n_genes_by_counts"], qcfg["mad_log_genes"])
                | _mad_outlier(g["pct_counts_in_top_20_genes"], qcfg["mad_pct_top20"])
            ).values
            f, rule = mito_flags(
                s.loc[gidx, "pct_counts_mt"],
                qcfg["mad_pct_mt"],
                qcfg.get("max_pct_mt"),
                qcfg.get("mad_floor_pct_mt", 0.0),
            )
            mt_flags.loc[gidx] = f.values
            used.add(rule)
        rules[smp] = "; ".join(sorted(used))
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
            p95_pct_mt=(
                "pct_counts_mt",
                lambda x: float(np.nanquantile(x, 0.95)) if x.notna().any() else np.nan,
            ),
        )
        .reset_index()
    )
    summary["mt_rule"] = summary[sample_key].map(rules).values
    if stratum_key and stratum_key in adata.obs:
        adata.uns["qc_per_compartment"] = (
            adata.obs.groupby(stratum_key, observed=True)
            .agg(
                n_cells=("qc_outlier", "size"),
                n_outlier=("qc_outlier", "sum"),
                n_mt_outlier=("qc_mt_outlier", "sum"),
                n_low_genes=("qc_low_genes", "sum"),
            )
            .reset_index()
        )
    if "mt_source" in adata.obs:
        src = adata.obs.groupby(sample_key, observed=True)["mt_source"].agg(lambda x: x.astype(str).iat[0])
        summary["mt_source"] = summary[sample_key].map(src).values
    return summary


def guard_scrublet(
    obs: pd.DataFrame, sample_key: str, max_called: float, rate_per_1000: float = 0.008, rate_cap: float = 0.1
) -> tuple[pd.Series, dict[str, str]]:
    """Per-sample sanity check on Scrublet's automatic threshold (§3.12).

    Scrublet picks a threshold per sample from the bimodal simulated-doublet score histogram. When that
    histogram is not bimodal the threshold can land inside the singlet mode and call most of the sample (run
    2: HLCA:PX9,
    972 / 1,001 cells). A sample with more than ``max_called`` of its cells called is treated as a threshold
    failure; instead, its top ``expected rate`` fraction by doublet score is called, with the 10x rule of
    thumb of ``rate_per_1000`` (≈ 0.8 %) per 1,000 cells in the sample, capped at ``rate_cap``. (The cells
    in a
    sample here are a subsample of the original capture, so this is a lower bound on the true rate:
    conservative.)
    """
    calls = obs["predicted_doublet"].astype(bool).copy()
    rules: dict[str, str] = {}
    for smp, idx in obs.groupby(sample_key, observed=True).groups.items():
        frac = float(calls.loc[idx].mean())
        if frac > max_called:
            rate = min(rate_per_1000 * len(idx) / 1000, rate_cap)
            score = obs.loc[idx, "doublet_score"]
            calls.loc[idx] = (score > score.quantile(1 - rate)).values
            rules[str(smp)] = f"fallback: top {rate:.1%} by score (auto threshold called {frac:.0%})"
        else:
            rules[str(smp)] = "scrublet auto threshold"
    return calls, rules


def heterotypic_rate(adata: ad.AnnData, sample_key: str = "sample", min_umi: int = 2) -> pd.Series:
    """Fraction of cells per sample co-expressing an immune marker (PTPRC) with an epithelial
    (EPCAM/KRT8/KRT18) or endothelial (PECAM1) marker, each ≥ ``min_umi`` UMIs: a marker-based
    heterotypic-doublet check that
    does not depend on Scrublet (§3.12). Reported only; not used to remove cells."""
    import scipy.sparse as sp

    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    X = sp.csr_matrix(X)

    def has(genes):
        cols = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
        if not cols:
            return np.zeros(adata.n_obs, bool)
        return np.asarray((X[:, cols] >= min_umi).sum(axis=1)).ravel() > 0

    co = has(["PTPRC"]) & (has(["EPCAM", "KRT8", "KRT18"]) | has(["PECAM1"]))
    return pd.Series(co, index=adata.obs_names).groupby(adata.obs[sample_key].values, observed=True).mean()


def call_doublets(adata: ad.AnnData, qcfg: dict, sample_key: str = "sample", seed: int = 0) -> None:
    """Scrublet per sample (+ optional scVI-SOLO; off since run 3, §3.12: miscalibrated on atlas data).

    With both enabled, `doublet` = both agree; with one, `doublet` = that method's call."""
    dcfg = qcfg["doublets"]
    adata.obs["doublet_scrublet"] = False
    adata.obs["doublet_solo"] = False
    if dcfg["scrublet"]:
        # scrublet reads .X, which still holds raw counts at this point (normalisation happens later)
        sc.pp.scrublet(adata, batch_key=sample_key, random_state=seed)
        calls, rules = guard_scrublet(
            adata.obs,
            sample_key,
            dcfg.get("max_called_fraction", 1.0),
            dcfg.get("expected_rate_per_1000", 0.008),
            dcfg.get("expected_rate_cap", 0.1),
        )

        adata.obs["doublet_scrublet"] = calls.values
        adata.uns["doublet_rule"] = rules
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
    het = heterotypic_rate(adata, sample_key)
    summary["frac_heterotypic_markers"] = summary[sample_key].map(het).values
    rules = adata.uns.get("doublet_rule", {})
    summary["doublet_rule"] = summary[sample_key].astype(str).map(rules).fillna("n/a").values
    keep = ~(
        adata.obs["qc_outlier"]
        | adata.obs["qc_mt_outlier"]
        | adata.obs["qc_low_genes"]
        | adata.obs["doublet"]
    )
    adata = adata[keep.values].copy()
    sc.pp.filter_genes(adata, min_cells=qcfg["min_cells_per_gene"])
    summary["n_kept"] = (
        summary[sample_key]
        .astype(str)
        .map(adata.obs[sample_key].astype(str).value_counts())
        .fillna(0)
        .astype(int)
    )
    log.info("QC kept %d / %d cells", adata.n_obs, int(summary["n_cells"].sum()))
    return adata, summary
