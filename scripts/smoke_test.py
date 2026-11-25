#!/usr/bin/env python
"""End-to-end smoke test on a small real subsample with tiny epoch counts.

Purpose: catch API/shape errors in minutes before a multi-hour full run. Writes nothing into data/ or
results/ — everything goes to $TMPDIR/lusc_smoke. Not a scientific result.

    pixi run python scripts/smoke_test.py [--n 6000] [--with-fm]
    pixi run -e c2l python scripts/smoke_test.py --c2l      # §3.10 #3: reference labels → c2l (tiny epochs)
    pixi run python scripts/smoke_test.py --run3            # §3.12: CNV → DE robustness → composition → CCC

`--run3` / `--c2l` start from a subsample of the existing `data/processed/sc_annotated.h5ad` (whole
tumour/adjacent samples of paired donors + healthy samples) because the upstream steps did not change.
"""

from __future__ import annotations

import argparse
import copy
import logging
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np

from lusc_atlas.config import load_config, set_seed, setup_logging

log = logging.getLogger("smoke")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12000)
    ap.add_argument("--with-fm", action="store_true", help="also run Geneformer on 64 cells")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--c2l", action="store_true", help="cell2location reference prep + tiny fit (c2l env)")
    g.add_argument(
        "--run3", action="store_true", help="§3.12 overhaul: CNV → DE robustness → composition → CCC"
    )
    a = ap.parse_args()
    setup_logging()
    cfg = copy.deepcopy(load_config())
    set_seed(0)
    repo = Path(cfg["_repo_root"])
    out = Path(os.environ.get("TMPDIR", "/tmp")) / "lusc_smoke"
    out.mkdir(parents=True, exist_ok=True)
    if a.c2l or a.run3:
        sub = annotated_subsample(repo, a.n)
        (smoke_c2l if a.c2l else smoke_run3)(cfg, repo, sub)
        return

    # tiny settings
    cfg["qc"]["doublets"]["solo_max_epochs"] = 3
    ic = cfg["integration"]
    ic.update(max_epochs=3, scanvi_max_epochs=2, benchmark_n_cells=1500, n_hvg=1000)
    ic["transfer"]["surgery_epochs"] = 2
    ic["transfer"]["converted_dir"] = str(out / "hlca_ref")
    cfg["cnv"]["max_ref_cells_per_compartment"] = 500
    cfg["trajectory"]["num_waypoints"] = 200

    t = time.time()
    full = ad.read_h5ad(repo / "data/interim/01_cohort.h5ad", backed="r")
    # whole samples (as in real data), ~n/4 cells per condition,
    # so per-sample QC/doublets see real sample sizes
    rng = np.random.default_rng(0)
    idx = []
    by_sample = full.obs.groupby(["condition", "sample"], observed=True).indices
    for cond in full.obs["condition"].cat.categories:
        keys = [k for k in by_sample if k[0] == cond and len(by_sample[k]) >= 250]
        rng.shuffle(keys)
        got = 0
        for k in keys:
            if got >= a.n // 4:
                break
            idx.append(by_sample[k])
            got += len(by_sample[k])
    adata = full[np.sort(np.concatenate(idx))].to_memory()
    full.file.close()
    log.info("smoke subsample %s in %.0fs", adata.shape, time.time() - t)

    from lusc_atlas.qc import run_qc

    cfg["qc"]["min_cells_per_sample"] = 100  # scrublet PCA needs n_cells > 30 components
    adata, qsum = run_qc(adata, cfg["qc"], 0)
    log.info("QC ok: %s; doublets %d", adata.shape, int(qsum.n_doublet.sum()))

    from lusc_atlas.fm import embed
    from lusc_atlas.integrate import atlas_mixing_per_label, benchmark_cells, normalise, run_integration

    normalise(adata)
    fm = None
    if a.with_fm:
        bidx = benchmark_cells(adata.obs_names, ic["benchmark_n_cells"], 0)
        sub = adata[bidx[:64]].copy()
        emb, info = embed(sub, ic["foundation_model"], repo)
        # pad to the benchmark subsample for the shape check
        full_emb = np.zeros((len(bidx), emb.shape[1]), np.float32)
        full_emb[:64] = emb
        full_emb[64:] = np.random.default_rng(0).normal(size=(len(bidx) - 64, emb.shape[1]))
        fm = (adata.obs_names[bidx].to_numpy(), full_emb)
        log.info("FM ok: %s", info)
    scores, chosen, _, extra = run_integration(adata, ic, 0, 16, repo=repo, fm_embedding=fm)
    log.info("integration ok → %s\n%s\nextra=%s", chosen, scores.round(3).to_string(), extra)
    log.info("mixing:\n%s", atlas_mixing_per_label(adata, ic["harmonised_label_key"]).head().to_string())

    from lusc_atlas.cnv import (
        annotate_positions,
        author_label_concordance,
        gene_positions_from_gtf,
        run_infercnv,
    )
    from lusc_atlas.de import make_de_groups

    pos = gene_positions_from_gtf(
        cfg["data"]["gene_positions"]["url"], repo / "data/raw/gencode.v44.basic.annotation.gtf.gz"
    )
    annotate_positions(adata, pos)
    cs = run_infercnv(adata, cfg["cnv"], 16, 0)
    make_de_groups(adata, ic["harmonised_label_key"])
    log.info(
        "CNV ok:\n%s\nconcordance:\n%s",
        cs.to_string(),
        author_label_concordance(adata, cfg["cnv"]["author_label_check"]),
    )

    from lusc_atlas.de import malignant_vs_healthy_basal, pathway_activities, pseudobulk, run_contrast

    cfg["de"]["min_samples_per_condition"] = 2
    pdata = pseudobulk(adata, cfg["de"])
    res = [run_contrast(pdata, cfg["de"], k, v, 8) for k, v in cfg["de"]["contrasts"].items()]
    res.append(malignant_vs_healthy_basal(pdata, cfg["de"], "Basal", 8))
    import pandas as pd

    allres = pd.concat([r for r in res if not r.empty], ignore_index=True)
    log.info(
        "DE ok: %d rows, groups %s",
        len(allres),
        allres.groupby("contrast").group.nunique().to_dict() if len(allres) else {},
    )
    if len(allres):
        pw = pathway_activities(allres, cfg["de"]["pathways"])
        log.info("pathways ok: %s", {k: v.shape for k, v in pw.items()})

    from lusc_atlas.trajectory import run_trajectory, subset_lineage

    sub = subset_lineage(adata, cfg["trajectory"], cfg["trajectory"]["label_key"])
    log.info("lineage subset %s %s", sub.shape, sub.obs.lineage_state.value_counts().to_dict())
    if sub.n_obs > 100:
        info = run_trajectory(sub, cfg["trajectory"], 0)
        log.info("trajectory ok: %s", info)
    log.info("SMOKE TEST PASSED in %.0fs", time.time() - t)


def annotated_subsample(repo: Path, n: int) -> ad.AnnData:
    """Whole samples from sc_annotated.h5ad: paired tumour/adjacent donors first, then healthy samples."""
    full = ad.read_h5ad(repo / "data/processed/sc_annotated.h5ad", backed="r")
    o = full.obs
    ta = o[o.condition.isin(["tumour", "adjacent_normal"])]
    both = ta.groupby("donor_id", observed=True)["condition"].nunique()
    paired = both[both == 2].index
    rng = np.random.default_rng(0)
    in_pairs = o.donor_id.isin(paired) & o.condition.isin(["tumour", "adjacent_normal"])
    samples = list(o.loc[in_pairs, "sample"].unique())
    # healthy samples richest in Basal cells first, so the *_vs_healthy_basal contrasts are exercised
    h = o[o.condition == "healthy"]
    healthy = list(
        (h.harmonised_label.astype(str) == "Basal")
        .groupby(h["sample"], observed=True)
        .sum()
        .sort_values(ascending=False)
        .index
    )
    sizes = o["sample"].value_counts()
    keep, got = [], 0
    for smp in samples + healthy:
        if got >= n and smp in healthy and len([k for k in keep if k in healthy]) >= 8:
            break
        keep.append(smp)
        got += int(sizes[smp])
    cap = max(n // len(keep), 50)  # equal cells per sample, keeping every sample
    by = o.groupby("sample", observed=True).indices
    idx = np.sort(
        np.concatenate([by[k] if len(by[k]) <= cap else rng.choice(by[k], cap, replace=False) for k in keep])
    )
    sub = full[idx].to_memory()
    full.file.close()
    log.info("annotated subsample %s from %d samples", sub.shape, sub.obs["sample"].nunique())
    return sub


def smoke_run3(cfg: dict, repo: Path, adata: ad.AnnData) -> None:
    """§3.12 overhaul on a real annotated subsample: per-sample CNV, DE with study-grounded robustness, paired
    composition, Milo on donor × tissue units, trajectory sensitivity, paired tensor (fixed small rank)."""
    import pandas as pd

    from lusc_atlas import de as D
    from lusc_atlas.ccc import by_sample, paired_loading_test, paired_units, tensor_c2c
    from lusc_atlas.cnv import annotate_positions, gene_positions_from_gtf, run_infercnv

    t = time.time()
    for c in ["malignant", "cnv_score", "de_group", "analysis_group"]:
        if c in adata.obs:
            del adata.obs[c]
    pos = gene_positions_from_gtf(
        cfg["data"]["gene_positions"]["url"], repo / "data/raw/gencode.v44.basic.annotation.gtf.gz"
    )
    annotate_positions(adata, pos)
    cc = dict(cfg["cnv"], mask_max_cells=1500, min_epithelial=30, min_reference=50)
    tabs = run_infercnv(adata, cc, 16, 0)
    log.info(
        "CNV ok:\n%s\nper sample:\n%s", tabs["summary"].to_string(), tabs["per_sample"].head(20).to_string()
    )
    c = cfg["cnv"]
    D.make_de_groups(
        adata,
        cfg["integration"]["harmonised_label_key"],
        cnv_low_label=c["cnv_low_label"],
        author_label=c["author_label_check"],
        tumour_conditions=c["tumour_conditions"],
        refine=cfg["de"]["refine_groups"],
    )
    log.info("de_group: %s", adata.obs["de_group"].value_counts().to_dict())
    dcfg = dict(cfg["de"], min_paired_donors=4, min_paired_studies=2, min_descriptive_donors=2)
    D.add_pseudobulk_unit(adata, dcfg["pseudobulk_sample_key"])
    pdata = D.pseudobulk(adata, dcfg)
    spec = dcfg["contrasts"]["paired_tumour_vs_adjacent"]
    r = D.run_contrast(pdata, dcfg, "paired_tumour_vs_adjacent", spec, 8)
    assert len(r), "paired DE produced nothing"
    assert {"padj_thr", "ci_low", "padj_global", "tier", "inference"} <= set(r.columns)
    loso, per, summ = D.paired_robustness(pdata, dcfg, "paired_tumour_vs_adjacent", spec, r, 8)
    r = r.merge(summ, on=["gene", "group"], how="left")
    amb = D.ambient_flags(pdata, r, dcfg)
    log.info(
        "DE paired ok: %d rows, groups %s; robust %d; LOSO %d rows; per-study %d rows; ambient suspects %d",
        len(r),
        sorted(r.group.unique()),
        int(r.robust.fillna(False).sum()),
        len(loso),
        len(per),
        int(amb.ambient_suspect.sum()) if len(amb) else 0,
    )
    log.info(
        "power:\n%s", D.power_table(r, r.groupby("group")["n_donors"].first().to_dict()).head().to_string()
    )
    log.info("sex check:\n%s", D.sex_check(pdata, dcfg["condition_key"]).head().to_string())
    log.info("cross-study null:\n%s", D.cross_study_null(pdata, dcfg, 8).to_string())
    ex = D.run_contrast(pdata, dcfg, "tumour_vs_healthy", dcfg["contrasts"]["tumour_vs_healthy"], 8)
    if len(ex):
        ex["supported_by_paired"] = D.supported_by_paired(ex, r, dcfg["alpha"])
    acfg = cfg["abundance"]
    cnt = D.composition_counts(adata.obs, acfg, acfg["group_key"], "donor_condition")
    comp = pd.concat(
        [
            D.composition_paired(cnt, {**acfg, **t_, "min_compartment_cells": 50})
            for t_ in acfg["paired_tests"]
        ],
        ignore_index=True,
    )
    log.info("composition ok:\n%s", comp[["cell_type", "glm_log2fc", "glm_padj", "consistent"]].to_string())
    milo = D.milo_da(adata, acfg, "donor_condition", "tissue_type", 0)
    log.info(
        "Milo ok: %d nhoods, %d units, %d donors", len(milo), milo.n_samples.iloc[0], milo.n_donors.iloc[0]
    )
    ccfg = dict(cfg["ccc"], tensor_rank=3, tensor_type_presence=0.6, tensor_outer_fraction=0.5)
    adata.X = adata.layers["lognorm"]
    unit = paired_units(adata.obs, "tumour", "adjacent_normal")
    sub = adata[unit.notna().to_numpy()].copy()
    sub.obs["donor_condition"] = unit.dropna().astype(str).astype("category").to_numpy()
    ures = by_sample(sub, ccfg, "donor_condition", 16)
    try:
        fac = tensor_c2c(sub, ures, "donor_condition", ccfg, 0)
        log.info("tensor ok; paired loading test:\n%s", paired_loading_test(fac, "tumour", "adjacent_normal"))
    except ValueError as exc:  # subsample may be too small for a non-degenerate tensor
        log.warning("tensor skipped on the smoke subsample: %s", exc)
    log.info("SMOKE TEST (run3) PASSED in %.0fs", time.time() - t)


def smoke_c2l(cfg: dict, repo: Path, ref: ad.AnnData) -> None:
    from lusc_atlas.de import make_de_groups
    from lusc_atlas.deconvolve import map_spots, prepare_reference_labels, reference_signatures

    t = time.time()
    c = cfg["cnv"]
    make_de_groups(
        ref,
        cfg["integration"]["harmonised_label_key"],
        cnv_low_label=c["cnv_low_label"],
        author_label=c["author_label_check"],
        tumour_conditions=c["tumour_conditions"],
    )
    cc = dict(cfg["spatial"]["c2l"], ref_max_epochs=3, max_epochs=20, max_ref_cells_per_label=300)
    cc["min_ref_cells"] = 50  # subsample scale
    labels, dec = prepare_reference_labels(
        ref.obs[cc["reference_label"]], cc["min_ref_cells"], cc["merge_small_labels"], cc["exclude_labels"]
    )
    log.info("reference decisions:\n%s", dec.to_string())
    ref = ref[labels.notna().values].copy()
    ref.obs["c2l_label"] = labels.dropna().astype(str).astype("category").values
    sig = reference_signatures(ref, {**cc, "reference_label": "c2l_label"}, 0)
    vis = ad.read_h5ad(repo / "data/interim/visium_qc.h5ad")
    vis = vis[np.random.default_rng(0).choice(vis.n_obs, 300, replace=False)].copy()
    ab, q05, props = map_spots(vis, sig, cc, 0)
    assert set(ab.columns) == set(sig.columns) and "Unknown" not in ab.columns
    assert props is not None and np.allclose(props.sum(axis=1), 1, atol=1e-4), "posterior proportions failed"
    log.info("c2l ok: %s; cell types %s", ab.shape, list(ab.columns))
    log.info("SMOKE TEST (c2l) PASSED in %.0fs", time.time() - t)


if __name__ == "__main__":
    main()
