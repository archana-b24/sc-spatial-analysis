# lusc-atlas — single-cell + spatial analysis of lung squamous cell carcinoma

Reproducible pipeline that integrates two curated lung single-cell atlases (**HLCA** and **LuCA**) with a
**10x Visium** section of lung squamous cell carcinoma (**LUSC**) to characterise malignant cell states,
the tumour microenvironment, and how they are organised in tissue.

> **v2.0.0** is a ground-up rewrite of the 2024 course project (tag [`v1.0.0`](../../tree/v1.0.0)).
> Every change, and the reason for it, is recorded in **[DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md)**.

## What it does

```
 HLCA (healthy reference) ─┐
                           ├─ dedupe + harmonise ─ per-sample QC ─ integration benchmark ─ inferCNV ─┬─ pseudobulk DE + pathways
 LuCA (45 LUSC donors) ────┘   (§3.3 A–E)          + doublets       scVI/scANVI/sysVI/      malignant  ├─ composition (scCODA, Milo)
                                                                     Harmony → choose          calls    ├─ trajectory (Palantir + CellRank 2)
                                                                                                         └─ cell–cell communication (LIANA+)
 Visium FFPE section ─ QC ─ cell2location (reference = annotated sc data) ─ SVGs, neighbourhoods, CellCharter niches, spatial LR
```

| Question | Method (v2) | Replaces (v1) |
|---|---|---|
| Which cells are malignant? | infercnvpy vs healthy-epithelium null; chr3q check | — (not done) |
| What changes in each cell type? | Pseudobulk PyDESeq2, paired tumour vs adjacent-normal | cell-level Wilcoxon (pseudoreplication) |
| Which pathways / TFs? | decoupler 2 — PROGENy, CollecTRI, Hallmark | GSEA with cells as samples |
| Composition shifts? | scCODA (compositional), Milo (neighbourhoods) | — |
| Basal → malignant axis? | Palantir + CellRank 2 (CytoTRACE + pseudotime kernels) | PAGA graphs |
| Where are cells in tissue? | cell2location deconvolution, CellCharter niches | Scanorama + cosine label per spot |
| Who talks to whom, and what changes? | LIANA+ per sample, differential CCC, Tensor-cell2cell, spatial bivariate | LIANA on tumour cells only |

## Quick start

All environments are pinned in `pixi.lock` ([pixi](https://pixi.sh)).

```bash
pixi install --all                    # default (analysis), c2l (cell2location), workflow (Snakemake)
pixi run test                         # fast unit tests

# Option A — Snakemake (incremental, resumable, cluster-ready)
pixi run -e workflow smk-dry          # show the plan
pixi run -e workflow smk              # run everything

# Option B — plain Python (same steps, same outputs)
pixi run python scripts/run_all.py --dry-run
pixi run run-all
pixi run python scripts/run_all.py --only qc integrate --force
```

Both runners are generated from the same step declarations (`src/lusc_atlas/steps.py`), so they cannot
diverge. Each step can also be run directly: `pixi run lusc-atlas <step>` (`pixi run lusc-atlas list`).

**Resources:** CPU-only is fine (developed on 16 cores / 125 GB RAM, no GPU). Raw data ≈ 40 GB.
**Outputs:** `results/tables/*.csv`, `results/figures/*.png`, `results/reports/summary.md`,
`data/processed/*.h5ad`, `data/processed/visium.spatialdata.zarr`, per-step provenance in
`results/reports/provenance/`.

## Data

All inputs are downloaded by the pipeline from pinned, public sources and checksum-verified.

| Data | Source | Version |
|---|---|---|
| Human Lung Cell Atlas, core + extended (Sikkema et al. 2023, *Nat Med*) | CELLxGENE Census | `2025-11-08` |
| Single-cell Lung Cancer Atlas, extended (Salcher et al. 2022, *Cancer Cell*) | CELLxGENE Census | `2025-11-08` |
| Visium CytAssist FFPE human lung cancer (squamous) | 10x Genomics public dataset | 2.0.1 |
| GENCODE v44 basic annotation | EBI | v44 |

Why both atlases: HLCA contains only 3 LUSC donors from a single study, so it cannot support a
tumour-vs-normal comparison alone. LuCA contributes 45 LUSC donors across studies. HLCA remains the
healthy-lung reference. Merging them safely (deduplication, gene space, label crosswalk, nested batches,
over-correction checks) is described in [DEVELOPMENT_LOG §3.3](DEVELOPMENT_LOG.md#33--merging-two-atlases--what-can-go-wrong-and-what-we-do-about-it).

## Results (first complete run, 2026-11-10)

Full auto-generated report: [`results/reports/summary.md`](results/reports/summary.md) (results are not
committed; regenerate with `pixi run lusc-atlas report`). Interpretation and caveats:
[DEVELOPMENT_LOG §3.4–3.10](DEVELOPMENT_LOG.md#34--qc--doublets--per-sample-both-conditions--resolves-s2-s6).

| Stage | Result |
|---|---|
| Cohort | 249,402 cells from 160 donors, 20 studies (HLCA: 119 healthy + 3 LUSC; LuCA: 38 LUSC; 15 LUSC donors have paired adjacent-normal tissue); 0 cross-atlas duplicate studies |
| QC | 231,504 cells kept across 215 samples (per-sample MAD; 672 Scrublet ∧ SOLO doublets) |
| Integration | 7 embeddings benchmarked with scib-metrics; **scANVI** chosen (total 0.677 vs sysVI 0.663, scVI 0.641, Harmony 0.637, unintegrated PCA 0.564; zero-shot Geneformer 0.385) |
| Malignant cells | 3,942 CNV-called (chr3q score 9× healthy epithelium); 99.9 % agree with author labels, but conservative (13 % of author-labelled malignant cells) |
| Pseudobulk DE | paired tumour vs adjacent-normal: 9,665 significant gene × cell-type results across 15 cell types, strongest in macrophages, neutrophils and T cells |
| Trajectory | Palantir/CellRank 2: Malignant cells sit further from healthy basal cells (median pseudotime 0.67 vs 0.55) |
| Spatial | 6,191 Visium spots deconvolved into 27 cell types; 3 CellCharter niches (tumour-epithelial, immune–stromal, fibro-muscular) |
| Cell–cell communication | Tumour-enriched neutrophil VEGFA/PLAU signalling and the TIGIT–NECTIN2 checkpoint axis; Tensor-cell2cell rank 5; spatial ECM → CD44/LRP1 co-expression |

**Known issues** being fixed before `v2.0.0` is tagged ([§3.10](DEVELOPMENT_LOG.md#310--post-run-review--open-issues-2026-11-11)):
Milo tested the wrong contrast. ~25k CNV-low author-labelled tumour cells sit in the `Unknown` label, which
affects scCODA and the `Unknown` rows of DE/CCC. A 66-cell rare-epithelial reference absorbs tumour signal in
cell2location.

## Repository layout

```
config/config.yaml        every parameter, path, version and checksum (single source of truth)
config/label_crosswalk.csv HLCA ↔ LuCA cell-type mapping to a harmonised level
src/lusc_atlas/           analysis package — one module per stage, `steps.py` wires I/O
workflow/Snakefile        Snakemake rules generated from steps.py
scripts/run_all.py        sequential runner (same steps)
tests/                    unit tests on synthetic data (CI)
legacy/v1/                original notebooks (outputs stripped; originals in tag v1.0.0)
docs/v1_report.pdf        original project report
```

## Limitations

- One Visium section: spatial results are descriptive (n = 1).
- Atlas data are processed counts; ambient RNA cannot be removed without raw droplet matrices.
- Malignant calls are CNV-based; samples with few CNVs may be under-called.
- Pseudotime across donors describes transcriptional distance from healthy basal cells, not lineage.

## Authors

Archana Balasubramanian, Shivaramakrishna Srinivasan, Aaryan Sharma, Aathira Sarath Chandran, Hersh Gupta.

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the HLCA, LuCA and the methods listed in the
development log.
