"""Pipeline steps: the ONLY place where file I/O meets analysis functions.

Each step reads declared inputs, writes declared outputs, and records provenance (package versions,
config hash, timings) to results/reports/provenance/<step>.json. Snakemake and run_all.py both call
these functions through `lusc_atlas.cli`, so the two runners cannot diverge.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import tarfile
import time
from importlib import metadata
from pathlib import Path

import pandas as pd

from .config import resolve, set_seed

log = logging.getLogger(__name__)

# ------------------------------------------------------------------ declared I/O (shared with Snakefile)
OUTPUTS = {
    "fetch_visium": ["data/raw/visium/.verified"],
    "fetch_atlases": ["data/raw/atlases_provenance.json"],
    "select": ["data/interim/01_cohort.h5ad", "results/tables/cohort_donors.csv"],
    "qc": ["data/interim/02_qc.h5ad", "results/tables/qc_per_sample.csv"],
    "fm_embed": ["data/interim/fm_geneformer_benchmark.npz"],
    "integrate": ["data/interim/03_integrated.h5ad", "results/tables/integration_benchmark.csv"],
    "cnv": ["data/processed/sc_annotated.h5ad", "results/tables/cnv_summary.csv"],
    "de": ["results/tables/de_pseudobulk.csv", "results/tables/pathways_progeny.csv"],
    "abundance": ["results/tables/composition_sccoda.csv"],
    "trajectory": ["data/processed/lineage.h5ad", "results/tables/trajectory_summary.csv"],
    "spatial_prep": ["data/interim/visium_qc.h5ad", "results/tables/visium_qc.csv"],
    "deconvolve": ["results/tables/c2l_abundance.csv"],
    "spatial": ["data/processed/visium.h5ad", "results/tables/spatial_niches.csv"],
    "ccc": ["results/tables/ccc_by_sample.csv"],
    "report": ["results/reports/summary.md"],
    "export": ["results/dashboard/lusc_atlas.duckdb"],
}
INPUTS = {
    "fetch_visium": [],
    "fetch_atlases": [],
    "select": ["data/raw/atlases_provenance.json"],
    "qc": ["data/interim/01_cohort.h5ad"],
    "fm_embed": ["data/interim/02_qc.h5ad"],
    "integrate": ["data/interim/02_qc.h5ad", "data/interim/fm_geneformer_benchmark.npz"],
    "cnv": ["data/interim/03_integrated.h5ad"],
    "de": ["data/processed/sc_annotated.h5ad"],
    "abundance": ["data/processed/sc_annotated.h5ad"],
    "trajectory": ["data/processed/sc_annotated.h5ad"],
    "spatial_prep": ["data/raw/visium/.verified"],
    "deconvolve": ["data/processed/sc_annotated.h5ad", "data/interim/visium_qc.h5ad"],
    "spatial": ["data/interim/visium_qc.h5ad", "results/tables/c2l_abundance.csv"],
    "ccc": [
        "data/processed/sc_annotated.h5ad",
        "results/tables/de_pseudobulk.csv",
        "data/processed/visium.h5ad",
    ],
    "report": [o for k, v in OUTPUTS.items() if k not in ("report", "export") for o in v],
    # dashboard snapshot (ROADMAP_v2.1 §6): reads tables + processed objects, never feeds back into analysis.
    # docs/dashboard_content.yaml is a tracked source file (like config.yaml), so it is not a declared input.
    "export": [
        "results/reports/summary.md",
        "data/processed/sc_annotated.h5ad",
        "data/processed/lineage.h5ad",
        "data/processed/visium.h5ad",
    ],
}
ORDER = list(OUTPUTS)
ENV = {"deconvolve": "c2l", "export": "dashboard"}  # every other step runs in the `default` pixi env


def _p(cfg: dict, rel: str) -> Path:
    return Path(cfg["_repo_root"]) / rel


def _provenance(cfg: dict, step: str, t0: float, extra: dict | None = None) -> None:
    pkgs = {}
    for name in [
        "scanpy",
        "anndata",
        "scvi-tools",
        "infercnvpy",
        "pydeseq2",
        "decoupler",
        "pertpy",
        "cellrank",
        "palantir",
        "squidpy",
        "spatialdata",
        "cellcharter",
        "liana",
        "cell2location",
        "torch",
        "duckdb",
        "numpy",
        "pandas",
        "lusc-atlas",
    ]:
        try:
            pkgs[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            pass
    cfg_clean = {k: v for k, v in cfg.items() if not k.startswith("_")}
    rec = {
        "step": step,
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "runtime_s": round(time.time() - t0, 1),
        "python": platform.python_version(),
        "config_sha256": hashlib.sha256(json.dumps(cfg_clean, sort_keys=True).encode()).hexdigest(),
        "packages": pkgs,
        **(extra or {}),
    }
    out = resolve(cfg, "reports", "provenance", f"{step}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2))


def _atomic_tmp(path: Path) -> Path:
    """Sibling temp path; written then os.replace()d so a killed step never leaves a truncated
    output that resume logic (run_all.py / Snakemake) could mistake for a finished one."""
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.with_name(f".{path.name}.partial{path.suffix}")


def _write_table(df: pd.DataFrame, path: Path) -> None:
    tmp = _atomic_tmp(path)
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def _write_h5ad(adata, path: Path) -> None:
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str).astype("category")
    tmp = _atomic_tmp(path)
    adata.write_h5ad(tmp, compression="gzip")
    os.replace(tmp, path)


# ================================================================= steps
def fetch_visium(cfg: dict) -> None:
    from .fetch import sha256

    vcfg = cfg["data"]["visium"]
    out = _p(cfg, "data/raw/visium")
    out.mkdir(parents=True, exist_ok=True)
    from .fetch import _try_download

    for name, meta in vcfg["files"].items():
        dest = out / name
        if not (dest.exists() and sha256(dest) == meta["sha256"]):
            if not _try_download(vcfg["base_url"] + meta["url_suffix"], dest):
                raise RuntimeError(f"download failed: {name}")
        got = sha256(dest)
        if got != meta["sha256"]:
            raise RuntimeError(f"checksum mismatch for {name}: {got}")
        if meta.get("extract"):
            with tarfile.open(dest) as tf:
                tf.extractall(out, filter="data")
    (out / ".verified").write_text(json.dumps({k: v["sha256"] for k, v in vcfg["files"].items()}, indent=2))


def fetch_atlases(cfg: dict) -> None:
    """Download each atlas' source h5ad from the pinned Census release; record size + sha256."""
    import cellxgene_census as cc

    from .fetch import sha256

    version = cfg["data"]["census"]["version"]
    records = {}
    for name, acfg in cfg["data"]["atlases"].items():
        dest = _p(cfg, acfg["file"])
        dest.parent.mkdir(parents=True, exist_ok=True)
        uri = cc.get_source_h5ad_uri(acfg["dataset_id"], census_version=version)["uri"]
        if not dest.exists():
            cc.download_source_h5ad(acfg["dataset_id"], to_path=str(dest), census_version=version)
        records[name] = {
            "dataset_id": acfg["dataset_id"],
            "census_version": version,
            "source_uri": uri,
            "bytes": dest.stat().st_size,
            "sha256": sha256(dest),
            "path": acfg["file"],
        }
    out = _p(cfg, OUTPUTS["fetch_atlases"][0])
    out.write_text(json.dumps(records, indent=2))


def select(cfg: dict) -> None:
    from .select import build_merged_cohort

    adata, donors, crosswalk, genes, exclusions = build_merged_cohort(cfg)
    _write_table(exclusions, resolve(cfg, "tables", "merge_exclusions.csv"))
    _write_table(donors, _p(cfg, OUTPUTS["select"][1]))
    _write_table(crosswalk, resolve(cfg, "tables", "merge_study_crosswalk.csv"))
    _write_table(genes, resolve(cfg, "tables", "merge_gene_space.csv"))
    _write_h5ad(adata, _p(cfg, OUTPUTS["select"][0]))


def qc(cfg: dict) -> None:
    import anndata as ad

    from .qc import run_qc

    adata = ad.read_h5ad(_p(cfg, INPUTS["qc"][0]))
    adata, summary = run_qc(adata, cfg["qc"], cfg["project"]["seed"])
    _write_table(summary, _p(cfg, OUTPUTS["qc"][1]))
    _write_h5ad(adata, _p(cfg, OUTPUTS["qc"][0]))


def fm_embed(cfg: dict) -> None:
    """Zero-shot Geneformer embedding of the shared benchmark subsample (CPU-heavy; own step)."""
    import anndata as ad
    import numpy as np

    from .fm import embed
    from .integrate import benchmark_cells

    out = _p(cfg, OUTPUTS["fm_embed"][0])
    out.parent.mkdir(parents=True, exist_ok=True)
    icfg = cfg["integration"]
    if not icfg["foundation_model"]["enabled"]:
        np.savez(
            out, obs_names=np.array([], dtype=str), emb=np.zeros((0, 0)), info=json.dumps({"enabled": False})
        )
        return
    a = ad.read_h5ad(_p(cfg, INPUTS["fm_embed"][0]), backed="r")
    idx = benchmark_cells(a.obs_names, icfg["benchmark_n_cells"], cfg["project"]["seed"])
    sub = a[idx].to_memory()
    a.file.close()
    emb, info = embed(sub, icfg["foundation_model"], Path(cfg["_repo_root"]), cfg["project"]["n_jobs"])
    np.savez(out, obs_names=sub.obs_names.to_numpy().astype(str), emb=emb, info=json.dumps(info))
    resolve(cfg, "reports", "fm_embed_info.json").write_text(json.dumps(info, indent=2))


def integrate(cfg: dict) -> None:
    import anndata as ad
    import numpy as np

    from .integrate import atlas_mixing_per_label, normalise, run_integration

    adata = ad.read_h5ad(_p(cfg, INPUTS["integrate"][0]))
    normalise(adata)
    fm = np.load(_p(cfg, INPUTS["integrate"][1]), allow_pickle=False)
    fm_embedding = (fm["obs_names"], fm["emb"]) if fm["emb"].size else None
    scores, chosen, scanvi_model, extra = run_integration(
        adata,
        cfg["integration"],
        cfg["project"]["seed"],
        cfg["project"]["n_jobs"],
        repo=Path(cfg["_repo_root"]),
        fm_embedding=fm_embedding,
    )
    scanvi_model.save(str(_p(cfg, "data/interim/scanvi_model")), overwrite=True)
    _write_table(scores, resolve(cfg, "tables", "integration_benchmark.csv"))
    resolve(cfg, "reports", "integration_info.json").write_text(json.dumps(extra, indent=2, default=str))
    _write_table(
        atlas_mixing_per_label(adata, cfg["integration"]["harmonised_label_key"]),
        resolve(cfg, "tables", "atlas_mixing_per_label.csv"),
    )
    lab = cfg["integration"]["harmonised_label_key"]
    conf = pd.crosstab(adata.obs[lab].astype(str), adata.obs["scanvi_prediction"].astype(str))
    _write_table(
        conf.rename_axis("harmonised_label").reset_index(),
        resolve(cfg, "tables", "scanvi_label_confusion.csv"),
    )
    if "hlca_ref_prediction" in adata.obs:
        conf2 = pd.crosstab(adata.obs[lab].astype(str), adata.obs["hlca_ref_prediction"].astype(str))
        _write_table(
            conf2.rename_axis("harmonised_label").reset_index(),
            resolve(cfg, "tables", "scarches_hlca_label_confusion.csv"),
        )
    _write_h5ad(adata, _p(cfg, OUTPUTS["integrate"][0]))


def cnv(cfg: dict) -> None:
    import anndata as ad

    from .cnv import annotate_positions, author_label_concordance, gene_positions_from_gtf, run_infercnv
    from .de import make_de_groups

    adata = ad.read_h5ad(_p(cfg, INPUTS["cnv"][0]))
    pos = gene_positions_from_gtf(
        cfg["data"]["gene_positions"]["url"], _p(cfg, "data/raw/gencode.v44.basic.annotation.gtf.gz")
    )
    annotate_positions(adata, pos)
    summary = run_infercnv(adata, cfg["cnv"], n_jobs=cfg["project"]["n_jobs"], seed=cfg["project"]["seed"])
    make_de_groups(adata, cfg["integration"]["harmonised_label_key"])
    _write_table(summary, _p(cfg, OUTPUTS["cnv"][1]))
    _write_table(
        author_label_concordance(adata, cfg["cnv"]["author_label_check"]),
        resolve(cfg, "tables", "cnv_vs_author_malignant.csv"),
    )
    per_donor = adata.obs.groupby(["atlas", "donor_id", "condition"], observed=True)["malignant"].agg(
        ["sum", "size"]
    )
    _write_table(per_donor.reset_index(), resolve(cfg, "tables", "cnv_malignant_per_donor.csv"))
    _write_h5ad(adata, _p(cfg, OUTPUTS["cnv"][0]))


def de(cfg: dict) -> None:
    import anndata as ad

    from .de import malignant_vs_healthy_basal, pathway_activities, pseudobulk, run_contrast

    adata = ad.read_h5ad(_p(cfg, INPUTS["de"][0]))
    pdata = pseudobulk(adata, cfg["de"])
    nj = min(8, cfg["project"]["n_jobs"])
    res = [
        run_contrast(pdata, cfg["de"], name, spec, n_jobs=nj) for name, spec in cfg["de"]["contrasts"].items()
    ]
    res.append(malignant_vs_healthy_basal(pdata, cfg["de"], cfg["trajectory"]["root_label"], n_jobs=nj))
    allres = pd.concat([r for r in res if not r.empty], ignore_index=True)
    _write_table(allres, _p(cfg, OUTPUTS["de"][0]))
    for name, tab in pathway_activities(allres, cfg["de"]["pathways"]).items():
        _write_table(tab, resolve(cfg, "tables", f"pathways_{name}.csv"))


def abundance(cfg: dict) -> None:
    import anndata as ad

    from .de import composition_sccoda, milo_da

    adata = ad.read_h5ad(_p(cfg, INPUTS["abundance"][0]))
    seed = cfg["project"]["seed"]
    # composition per sample (tumour and adjacent-normal of one donor are separate samples)
    eff = composition_sccoda(adata, "de_group", "sample", "condition", seed)
    _write_table(eff, _p(cfg, OUTPUTS["abundance"][0]))
    try:
        milo = milo_da(adata, cfg["abundance"], "sample", "condition", seed)
        _write_table(milo, resolve(cfg, "tables", "milo_nhoods.csv"))
    except Exception as exc:  # noqa: BLE001
        log.warning("Milo failed (%s); scCODA result still written", exc)


def trajectory(cfg: dict) -> None:
    import anndata as ad

    from .trajectory import run_trajectory, subset_lineage

    adata = ad.read_h5ad(_p(cfg, INPUTS["trajectory"][0]))
    sub = subset_lineage(adata, cfg["trajectory"], cfg["trajectory"]["label_key"])
    info = run_trajectory(sub, cfg["trajectory"], cfg["project"]["seed"])
    summ = sub.obs.groupby("lineage_state", observed=True)["palantir_pseudotime"].describe().reset_index()
    _write_table(summ, _p(cfg, OUTPUTS["trajectory"][1]))
    (resolve(cfg, "tables", "trajectory_info.json")).write_text(json.dumps(info, indent=2))
    _write_h5ad(sub, _p(cfg, OUTPUTS["trajectory"][0]))


def spatial_prep(cfg: dict) -> None:
    from .spatial import load_visium, qc_and_cluster, write_spatialdata

    raw = _p(cfg, "data/raw/visium")
    vis = load_visium(raw, cfg["data"]["visium"]["sample_id"])
    summ = qc_and_cluster(vis, cfg["spatial"], cfg["project"]["seed"])
    _write_table(summ, _p(cfg, OUTPUTS["spatial_prep"][1]))
    _write_h5ad(vis, _p(cfg, OUTPUTS["spatial_prep"][0]))
    write_spatialdata(
        raw, _p(cfg, "data/processed/visium.spatialdata.zarr"), cfg["data"]["visium"]["sample_id"]
    )


def deconvolve(cfg: dict) -> None:
    import anndata as ad

    from .deconvolve import map_spots, reference_signatures

    seed = cfg["project"]["seed"]
    ref = ad.read_h5ad(_p(cfg, INPUTS["deconvolve"][0]))
    vis = ad.read_h5ad(_p(cfg, INPUTS["deconvolve"][1]))
    sig = reference_signatures(ref, cfg["spatial"]["c2l"], seed)
    sig.rename_axis("gene").reset_index().to_csv(
        resolve(cfg, "tables", "c2l_reference_signatures.csv"), index=False
    )
    ab = map_spots(vis, sig, cfg["spatial"]["c2l"], seed)
    _write_table(ab.rename_axis("spot").reset_index(), _p(cfg, OUTPUTS["deconvolve"][0]))


def spatial(cfg: dict) -> None:
    import anndata as ad

    from .spatial import niches_cellcharter, spatial_stats

    vis = ad.read_h5ad(_p(cfg, INPUTS["spatial"][0]))
    ab = pd.read_csv(_p(cfg, INPUTS["spatial"][1])).set_index("spot").loc[vis.obs_names]
    cell_types = list(ab.columns)
    props = ab.div(ab.sum(axis=1), axis=0)
    for c in cell_types:
        vis.obs[c] = props[c].values
    vis.obsm["c2l_abundance"] = ab.to_numpy()
    vis.uns["c2l_cell_types"] = cell_types
    for name, tab in spatial_stats(vis, cfg["spatial"], cell_types, cfg["project"]["n_jobs"]).items():
        _write_table(tab, resolve(cfg, "tables", f"spatial_{name}.csv"))
    _write_table(
        niches_cellcharter(vis, cfg["spatial"], cell_types, cfg["project"]["seed"]),
        _p(cfg, OUTPUTS["spatial"][1]),
    )
    _write_h5ad(vis, _p(cfg, OUTPUTS["spatial"][0]))


def ccc(cfg: dict) -> None:
    import anndata as ad

    from .ccc import by_sample, differential, spatial_bivariate, tensor_c2c

    adata = ad.read_h5ad(_p(cfg, INPUTS["ccc"][0]))
    adata.X = adata.layers["lognorm"]
    de_res = pd.read_csv(_p(cfg, INPUTS["ccc"][1]))
    res = by_sample(adata, cfg["ccc"], "sample", cfg["project"]["n_jobs"])
    _write_table(res, _p(cfg, OUTPUTS["ccc"][0]))
    # Each sub-analysis runs even if an earlier one fails, but any failure fails the step
    # (no silent skips — v2.0 run lost differential/tensor CCC to a swallowed API error).
    failures: dict[str, BaseException] = {}
    try:
        _write_table(
            differential(adata, de_res[de_res.contrast == cfg["ccc"]["de_contrast"]], cfg["ccc"]),
            resolve(cfg, "tables", "ccc_differential.csv"),
        )
    except Exception as exc:  # noqa: BLE001
        log.exception("differential CCC failed")
        failures["differential"] = exc
    try:
        for name, f in tensor_c2c(adata, res, "sample", cfg["ccc"], cfg["project"]["seed"]).items():
            _write_table(
                f.rename_axis("element").reset_index(),
                resolve(cfg, "tables", f"ccc_tensor_{name.lower().replace(' ', '_')}.csv"),
            )
    except Exception as exc:  # noqa: BLE001
        log.exception("tensor-cell2cell failed")
        failures["tensor_c2c"] = exc
    try:
        vis = ad.read_h5ad(_p(cfg, INPUTS["ccc"][2]))
        _write_table(spatial_bivariate(vis, cfg["ccc"]), resolve(cfg, "tables", "ccc_spatial_bivariate.csv"))
    except Exception as exc:  # noqa: BLE001
        log.exception("spatial bivariate CCC failed")
        failures["spatial_bivariate"] = exc
    if failures:
        raise RuntimeError(f"ccc sub-analyses failed: {sorted(failures)} — see traceback(s) above")


def report(cfg: dict) -> None:
    from .report import build_report

    build_report(cfg)


def export(cfg: dict) -> None:
    from .export import build_dashboard_db

    build_dashboard_db(cfg, _p(cfg, OUTPUTS["export"][0]))


STEPS = {name: globals()[name] for name in ORDER}


def run_step(cfg: dict, name: str) -> None:
    set_seed(cfg["project"]["seed"])
    t0 = time.time()
    log.info("=== step %s ===", name)
    STEPS[name](cfg)
    _provenance(cfg, name, t0)
    log.info("=== step %s done in %.0fs ===", name, time.time() - t0)


__all__ = ["ORDER", "INPUTS", "OUTPUTS", "ENV", "STEPS", "run_step"]
