"""Malignant-cell identification by inferred copy-number variation (DEVELOPMENT_LOG §3.5, §3.12).

v1 had no malignant/normal distinction (S7). v2 runs 1–2 ran one joint infercnvpy job on HLCA-healthy + LuCA
tumour epithelium and thresholded ``cnv_score`` at the 95th percentile of healthy epithelium. Two flaws
(§3.12): ``cnv_score`` is the *mean of the cell's CNV cluster* (so the "percentile" was over ~29 cluster
values and
"Malignant" was exactly 7 clusters, 87 % of calls from 3 patients), and the null came from the other atlas.

Run 3 — per-sample calling with a matched null (Tirosh 2016; Neftel 2019 two-axis rule):

1. **Tissue**: only tumour-patient samples are run (tumour, region-unknown tumour, adjacent normal); healthy
   donors are not used for *calling*, so atlas/chemistry cannot enter the threshold (they do enter the window
   mask, point 3).
2. **Reference** = the same sample's immune/endothelial/stromal cells (diploid), pooled within the study
when a
   sample has too few; a random 20 % of them are held out as **pseudo-query cells**: their CNV statistics, at
   the sample's own depth and chemistry, are the empirical null.
3. **Window mask**: inferCNV also picks up expression programmes whose genes cluster on a chromosome (keratin
   clusters 12q13/17q21, the epidermal differentiation complex 1q21, SPRR/S100). A separate run on *normal*
   squamous-like and alveolar epithelium (healthy basal/suprabasal + adjacent-normal epithelium) finds windows
   where normal epithelium deviates systematically; those windows are masked before scoring.
4. **Per-cell axes**: (a) CNV signal = mean squared CNV value over unmasked windows; (b) correlation of the
   cell's CNV profile with its sample's consensus tumour profile, **cross-fitted**: the consensus for one
   random half of the sample's epithelium is built from the top 10 % (by signal) of the *other* half, so a
   cell never
   contributes to the profile it is compared with.
5. **Call**: malignant = epithelial ∧ tumour tissue ∧ both axes above the sample's 99th percentile of the
   held-out null. Samples with < ``min_epithelial`` epithelial cells are "unresolved".
6. **Validation** (not used for calling): call rate in adjacent-normal epithelium (should be ≈ 1 %), LUSC
   arm-level hallmarks (3q gain, 3p loss, 5p gain) in called vs not-called cells, per-donor counts, and
   agreement with the LuCA authors' labels (themselves inferCNV-informed, so not ground truth).
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


TUMOUR_TISSUE = ("tumour", "tumour_region_unknown")
# GRCh38 centromere bounds (UCSC cytoBand acen), for arm-level validation signals
CENTROMERE = {"chr3": (87_800_000, 94_000_000), "chr5": (46_100_000, 51_400_000)}


def tissue_of(condition: pd.Series) -> pd.Series:
    """tumour (incl. region-unknown tumour) / adjacent_normal / healthy — §3.12: a tumour sample whose
    region was
    not recorded is still tumour tissue (it makes Thienpont's second patient pairable)."""
    c = condition.astype(str)
    return c.where(~c.isin(TUMOUR_TISSUE), "tumour")


def _infercnv(sub: ad.AnnData, ref_mask: np.ndarray, ccfg: dict, n_jobs: int) -> tuple[np.ndarray, dict]:
    import infercnvpy as cnv

    sub.obs["cnv_reference"] = np.where(ref_mask, "reference", "query")
    cnv.tl.infercnv(
        sub,
        reference_key="cnv_reference",
        reference_cat="reference",
        window_size=ccfg["window_size"],
        step=ccfg["step"],
        n_jobs=n_jobs,
        layer="lognorm",
    )
    X = sub.obsm["X_cnv"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return X.astype(np.float32), dict(sub.uns["cnv"]["chr_pos"])


def window_mask(
    adata: ad.AnnData, ccfg: dict, n_jobs: int, seed: int
) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """Windows where NORMAL epithelium shows a systematic CNV-like shift (expression programmes, not copy
    number).

    Normal epithelium = healthy-donor squamous-lineage cells (``mask_normal_labels``, e.g. Basal) +
    adjacent-normal epithelium; reference = immune/endothelial/stromal cells from the same samples. A window
    is masked when
    the absolute mean signal of normal epithelium is a robust outlier (> median + ``mask_nmads`` × MAD
    across windows).
    """
    from scipy.stats import median_abs_deviation

    rng = np.random.default_rng(seed)
    comp = adata.obs[ccfg["compartment_key"]].astype(str)
    lab = adata.obs[ccfg.get("label_key", "harmonised_label")].astype(str)
    tissue = adata.obs["tissue_type"].astype(str)
    epi = comp == ccfg["epithelial_compartment"]
    normal = epi & (
        (tissue == "adjacent_normal") | ((tissue == "healthy") & lab.isin(ccfg["mask_normal_labels"]))
    )
    nidx = np.where(normal)[0]
    if len(nidx) > ccfg["mask_max_cells"]:
        nidx = rng.choice(nidx, ccfg["mask_max_cells"], replace=False)
    samples = set(adata.obs["sample"].iloc[nidx].astype(str))
    refpool = np.where(
        comp.isin(ccfg["reference_compartments"]) & adata.obs["sample"].astype(str).isin(samples)
    )[0]
    if len(refpool) > ccfg["mask_max_cells"]:
        refpool = rng.choice(refpool, ccfg["mask_max_cells"], replace=False)
    run = np.sort(np.concatenate([nidx, refpool]))
    sub = adata[run].copy()
    X, chr_pos = _infercnv(sub, np.isin(run, refpool), ccfg, n_jobs)
    m = X[~np.isin(run, refpool)].mean(0)
    a = np.abs(m)
    cut = np.median(a) + ccfg["mask_nmads"] * median_abs_deviation(a)
    mask = a > cut
    chroms = sorted(chr_pos, key=chr_pos.get)
    starts = [chr_pos[c] for c in chroms]
    win_chr = [chroms[np.searchsorted(starts, w, side="right") - 1] for w in range(len(m))]
    tab = pd.DataFrame({"window": np.arange(len(m)), "chromosome": win_chr, "normal_mean": m, "masked": mask})
    log.info(
        "CNV window mask: %d / %d windows masked (normal epithelium %d cells, |mean| > %.4f); by chromosome: "
        "%s",
        mask.sum(),
        len(mask),
        len(nidx),
        cut,
        tab[tab.masked].chromosome.value_counts().head(8).to_dict(),
    )
    return mask, tab, chr_pos


def _arm_windows(adata_var: pd.DataFrame, chr_pos: dict, n_windows: int, ccfg: dict) -> dict[str, np.ndarray]:
    """Window indices of LUSC hallmark arms (3p, 3q, 5p), from window centre-gene positions."""
    chroms = sorted(chr_pos, key=chr_pos.get)
    out = {}
    for chrom, arm in (("chr3", "p"), ("chr3", "q"), ("chr5", "p")):
        i = chroms.index(chrom)
        w0 = chr_pos[chrom]
        w1 = chr_pos[chroms[i + 1]] if i + 1 < len(chroms) else n_windows
        genes = adata_var.loc[adata_var["chromosome"] == chrom].sort_values("start")["start"].to_numpy()
        centres = np.minimum(np.arange(w1 - w0) * ccfg["step"] + ccfg["window_size"] // 2, len(genes) - 1)
        lo, hi = CENTROMERE[chrom]
        sel = genes[centres] < lo if arm == "p" else genes[centres] > hi
        out[f"{chrom}{arm}"] = np.where(sel)[0] + w0
    return out


def _cross_fit_corr(Xe: np.ndarray, score: np.ndarray, rng, top: float) -> np.ndarray:
    """Correlation of each epithelial cell's profile with a consensus built from the OTHER half's top
    cells."""
    n = len(score)
    half = rng.permutation(n) < n // 2
    corr = np.full(n, np.nan, dtype=np.float32)
    for a in (half, ~half):
        other = ~a
        if other.sum() < 5:
            continue
        k = max(3, int(np.ceil(top * other.sum())))
        o_idx = np.where(other)[0]
        cons = Xe[o_idx[np.argsort(score[o_idx])[-k:]]].mean(0)
        corr[a] = _rowcorr(Xe[a], cons)
    return corr


def _rowcorr(X: np.ndarray, v: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(1, keepdims=True)
    vc = v - v.mean()
    den = np.sqrt((Xc**2).sum(1) * (vc**2).sum())
    return np.divide(Xc @ vc, den, out=np.zeros(len(X), dtype=np.float64), where=den > 0).astype(np.float32)


def run_infercnv(adata: ad.AnnData, ccfg: dict, n_jobs: int = 16, seed: int = 0) -> dict[str, pd.DataFrame]:
    """Per-sample two-axis malignant calling (module docstring). Writes obs columns ``tissue``, ``cnv_score``,
    ``cnv_corr``, ``cnv_chr3p``, ``cnv_chr3q``, ``cnv_chr5p``, ``cnv_status``, ``malignant`` and returns
    tables: per_sample, summary, window_mask."""
    rng = np.random.default_rng(seed)
    adata.obs["tissue_type"] = tissue_of(adata.obs["condition"]).values
    mask, mask_tab, chr_pos = window_mask(adata, ccfg, n_jobs, seed)
    keep_w = ~mask
    arms = _arm_windows(adata.var, chr_pos, len(mask), ccfg)

    obs = adata.obs
    comp = obs[ccfg["compartment_key"]].astype(str)
    epi = (comp == ccfg["epithelial_compartment"]).to_numpy()
    ref_all = comp.isin(ccfg["reference_compartments"]).to_numpy()
    tissue = obs["tissue_type"].astype(str).to_numpy()
    # §3.12 (senior review B4): resolve per patient × tissue (a patient's tumour regions pooled), not per
    # sample — many immune-enriched samples hold too few epithelial cells on their own
    sample = (obs["donor_id"].astype(str) + "|" + obs["tissue_type"].astype(str)).to_numpy()
    study = obs["study"].astype(str).to_numpy()
    run_tissues = {"tumour", "adjacent_normal"}

    n = adata.n_obs
    cols = {
        k: np.full(n, np.nan, dtype=np.float32)
        for k in ["cnv_score", "cnv_corr", "cnv_chr3p", "cnv_chr3q", "cnv_chr5p"]
    }
    status = np.array(["not run"] * n, dtype=object)
    malignant = np.zeros(n, bool)
    null_rows, per_sample = [], []
    for smp in sorted(set(sample[np.isin(tissue, list(run_tissues))])):
        in_s = sample == smp
        e_idx = np.where(in_s & epi)[0]
        tis = tissue[in_s][0]
        st = study[in_s][0]
        if len(e_idx) < ccfg["min_epithelial"]:
            status[e_idx] = f"unresolved (< {ccfg['min_epithelial']} epithelial cells)"
            per_sample.append(
                {"unit": smp, "study": st, "tissue_type": tis, "n_epithelial": len(e_idx), "resolved": False}
            )
            continue
        r_idx = np.where(in_s & ref_all)[0]
        ref_source = "same patient and tissue"
        if len(r_idx) < ccfg["min_reference"]:
            r_idx = np.where((study == st) & ref_all & np.isin(tissue, list(run_tissues)))[0]
            ref_source = "same study (unit has too few)"
        if len(r_idx) > ccfg["max_ref_cells_per_sample"]:
            r_idx = rng.choice(r_idx, ccfg["max_ref_cells_per_sample"], replace=False)
        if len(r_idx) < ccfg["min_reference"]:
            status[e_idx] = "unresolved (too few reference cells)"
            per_sample.append(
                {
                    "unit": smp,
                    "study": st,
                    "tissue_type": tis,
                    "n_epithelial": len(e_idx),
                    "n_reference": len(r_idx),
                    "resolved": False,
                }
            )
            continue
        # held-out null of at least `min_heldout` cells (a q99 from ~20 cells is close to their maximum)
        frac = max(ccfg["heldout_fraction"], ccfg.get("min_heldout", 100) / len(r_idx))
        held = np.zeros(len(r_idx), bool)
        n_held = min(len(r_idx) - 50, int(np.ceil(frac * len(r_idx))))
        held[rng.choice(len(r_idx), n_held, replace=False)] = True
        run = np.concatenate([e_idx, r_idx])
        sub = adata[run].copy()
        is_ref = np.r_[np.zeros(len(e_idx), bool), ~held]
        X, _ = _infercnv(sub, is_ref, ccfg, n_jobs)
        if X.shape[1] != len(mask):
            raise ValueError(f"{smp}: CNV window count {X.shape[1]} ≠ mask {len(mask)}")
        Xk = X[:, keep_w]
        score = (Xk**2).mean(1)
        Xe, se = Xk[: len(e_idx)], score[: len(e_idx)]
        Xh, sh = Xk[len(e_idx) :][held], score[len(e_idx) :][held]
        corr_e = _cross_fit_corr(Xe, se, rng, ccfg["consensus_top_fraction"])
        k = max(3, int(np.ceil(ccfg["consensus_top_fraction"] * len(se))))
        cons_all = Xe[np.argsort(se)[-k:]].mean(0)
        corr_h = _rowcorr(Xh, cons_all)
        q = ccfg["null_quantile"]
        thr_s, thr_c = float(np.quantile(sh, q)), float(np.quantile(corr_h, q))
        call = (se > thr_s) & (corr_e > thr_c)
        cols["cnv_score"][e_idx] = se
        cols["cnv_corr"][e_idx] = corr_e
        for arm, w in arms.items():
            cols[f"cnv_{arm}"][e_idx] = X[: len(e_idx), w].mean(1)
        is_tum = tis == "tumour"
        malignant[e_idx] = call & is_tum
        status[e_idx] = np.where(
            call, "CNV-high" if is_tum else "CNV-high (adjacent: false-positive check)", "not called"
        )
        null_rows.append(pd.DataFrame({"cnv_score": sh, "cnv_corr": corr_h, "unit": smp}))
        per_sample.append(
            {
                "unit": smp,
                "study": st,
                "donor_id": str(obs["donor_id"].iloc[e_idx[0]]),
                "n_samples": int(obs["sample"].iloc[e_idx].nunique()),
                "tissue_type": tis,
                "n_epithelial": len(e_idx),
                "n_reference": int((~held).sum()),
                "n_heldout_null": int(held.sum()),
                "reference_source": ref_source,
                "threshold_score": thr_s,
                "threshold_corr": thr_c,
                "n_called": int(call.sum()),
                "call_rate": float(call.mean()),
                "resolved": True,
            }
        )
        log.info(
            "CNV %s (%s, %s): %d epi, %d ref → %d called (%.1f %%)",
            smp,
            st,
            tis,
            len(e_idx),
            int((~held).sum()),
            int(call.sum()),
            100 * call.mean(),
        )
    for k, v in cols.items():
        adata.obs[k] = v
    adata.obs["cnv_status"] = pd.Categorical(status)
    adata.obs["malignant"] = malignant
    adata.uns["cnv_chr_pos"] = {k: int(v) for k, v in chr_pos.items()}

    null = pd.concat(null_rows, ignore_index=True) if null_rows else pd.DataFrame()
    o = adata.obs
    grp = np.select(
        [
            o["malignant"].to_numpy(),
            epi & (tissue == "tumour") & (status == "not called"),
            epi
            & (tissue == "adjacent_normal")
            & np.isin(status, ["not called", "CNV-high (adjacent: false-positive check)"]),
        ],
        ["malignant (CNV-high)", "tumour-tissue epithelium, not called", "adjacent-normal epithelium"],
        default="",
    )
    rows = []
    for g in ["malignant (CNV-high)", "tumour-tissue epithelium, not called", "adjacent-normal epithelium"]:
        m = grp == g
        rows.append(
            {
                "cnv_group": g,
                "n_cells": int(m.sum()),
                "n_donors": int(o.loc[m, "donor_id"].nunique()),
                **{
                    f"mean_{c}": float(np.nanmean(o.loc[m, c])) if m.any() else np.nan
                    for c in ["cnv_score", "cnv_corr", "cnv_chr3p", "cnv_chr3q", "cnv_chr5p"]
                },
            }
        )
    if len(null):
        rows.append(
            {
                "cnv_group": "held-out reference cells (null)",
                "n_cells": len(null),
                "mean_cnv_score": float(null.cnv_score.mean()),
                "mean_cnv_corr": float(null.cnv_corr.mean()),
            }
        )
    summary = pd.DataFrame(rows)
    adj = (
        (tissue == "adjacent_normal")
        & epi
        & (status != "not run")
        & ~pd.Series(status).str.startswith("unresolved").to_numpy()
    )
    adj_rate = (
        float((status[adj] == "CNV-high (adjacent: false-positive check)").mean()) if adj.any() else np.nan
    )
    summary.attrs["adjacent_call_rate"] = adj_rate
    # pre-registered acceptance gate (config cnv.acceptance, fixed before run 3): if it fails, every page that
    # uses malignant calls says so
    gate = ccfg.get("acceptance", {})
    mal = summary.set_index("cnv_group").loc["malignant (CNV-high)"] if len(summary) else None
    checks = {
        "adjacent_call_rate_ok": bool(adj_rate <= gate.get("max_adjacent_call_rate", 0.05)),
        "chr3q_gain_ok": mal is not None and bool(mal["mean_cnv_chr3q"] > 0),
        "chr3p_loss_ok": mal is not None and bool(mal["mean_cnv_chr3p"] < 0),
        "min_donors_ok": mal is not None and bool(mal["n_donors"] >= gate.get("min_donors_with_calls", 10)),
    }
    for k, v in checks.items():
        summary[k] = bool(v)
    summary["acceptance_passed"] = all(checks.values())
    log.info("CNV acceptance gate: %s → %s", checks, "PASSED" if all(checks.values()) else "FAILED")
    return {"per_sample": pd.DataFrame(per_sample), "summary": summary, "window_mask": mask_tab}


def author_label_concordance(adata: ad.AnnData, author_label: str) -> pd.DataFrame:
    """Independent check: our CNV calls vs the LuCA authors' 'malignant cell' annotation (LuCA cells only)."""
    o = adata.obs[(adata.obs["atlas"] == "LuCA") & adata.obs["cnv_score"].notna()]
    if "tissue_type" in o:
        o = o[o["tissue_type"].astype(str) == "tumour"]
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
