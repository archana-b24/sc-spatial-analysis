#!/usr/bin/env python
"""End-to-end smoke test on a small real subsample with tiny epoch counts.

Purpose: catch API/shape errors in minutes before a multi-hour full run. Writes nothing into data/ or
results/ — everything goes to $TMPDIR/lusc_smoke. Not a scientific result.

    pixi run python scripts/smoke_test.py [--n 6000] [--with-fm]
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
    a = ap.parse_args()
    setup_logging()
    cfg = copy.deepcopy(load_config())
    set_seed(0)
    repo = Path(cfg["_repo_root"])
    out = Path(os.environ.get("TMPDIR", "/tmp")) / "lusc_smoke"
    out.mkdir(parents=True, exist_ok=True)

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


if __name__ == "__main__":
    main()
