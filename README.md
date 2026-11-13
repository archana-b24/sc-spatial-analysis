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

### Results dashboard (v2.1, in progress)

An interactive Streamlit site with one page per pipeline stage. Each page pairs the decision with the evidence
behind it and shows the known caveats. It reads a DuckDB snapshot that the `export` step derives from pipeline
outputs (design: `docs/ROADMAP_v2.1.md` §6).

```bash
pixi run -e dashboard export      # pipeline outputs → results/dashboard/lusc_atlas.duckdb (~1 min)
pixi run -e dashboard app         # http://localhost:8501
```

v3 preview (prototypes, `docs/v3/README.md`): a *Start here* guide and computed reader guide on every page, plus the
atlas's read-only DuckDB tools exposed as an MCP server. (The in-app chatbot is parked: the project stays free to
host and run, and a chatbot needs a paid LLM API key.)

```bash
pixi run -e dashboard mcp                                                  # same tools as an MCP server
pixi run -e dashboard cloud-dry-run                                        # rehearse the free deploy locally
```

Hosting is free-only: the app on **Streamlit Community Cloud** (deployed from this public repo, entrypoint
`app/streamlit_app.py`, deps `app/requirements.txt`) and the DuckDB as a **public Hugging Face Dataset**
(`scripts/publish_dataset.sh`, sha256-pinned at app start). Every service used, its free tier and what would
trigger a cost: `docs/v3/C_dashboard_hosting.md` → *Cost: $0*.

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

## Results (run 2, 2026-11-13)

Full auto-generated report: [`results/reports/summary.md`](results/reports/summary.md) (results are not
committed; regenerate with `pixi run lusc-atlas report`). Interpretation and caveats:
[DEVELOPMENT_LOG §3.4–3.11](DEVELOPMENT_LOG.md#34--qc--doublets--per-sample-both-conditions--resolves-s2-s6).

| Stage | Result |
|---|---|
| Cohort | 249,402 cells from 160 donors, 20 studies (HLCA: 119 healthy + 3 LUSC; LuCA: 38 LUSC; 15 LUSC donors have paired adjacent-normal tissue); 0 cross-atlas duplicate studies |
| QC | 231,504 cells kept across 215 samples (per-sample MAD; 672 Scrublet ∧ SOLO doublets) |
| Integration | 7 embeddings benchmarked with scib-metrics; **scANVI** chosen (total 0.677 vs sysVI 0.663, scVI 0.641, Harmony 0.637, unintegrated PCA 0.564; zero-shot Geneformer 0.385) |
| Malignant cells | 3,942 CNV-called (chr3q score 9× healthy epithelium); 99.9 % agree with author labels, but conservative (13 % of author-labelled malignant cells). The other 25,205 author-labelled tumour cells form `Tumour epithelial (CNV-low)`; vs healthy basal they show SOX2, TP63 and CDKN2A up |
| Pseudobulk DE | paired tumour vs adjacent-normal: 8,474 significant gene × cell-type results across 15 cell types, strongest in macrophages, neutrophils and T cells |
| Composition | scCODA (reference cDC): tumour has fewer NK, CD8 T, CD4 T, non-classical monocytes, AT1/AT2 and macrophages than adjacent tissue, and more CNV-low tumour epithelium. Milo (paired, 14 donors): 1,671 of 6,225 neighbourhoods shift, with T-cell states moving in both directions |
| Trajectory | Palantir/CellRank 2: Malignant cells sit further from healthy basal cells (median pseudotime 0.67 vs 0.55) |
| Spatial | 6,191 Visium spots deconvolved into 25 cell types; 9 CellCharter niches (5 tumour-epithelial, 3 stromal–myeloid, 1 fibro-muscular); tumour niches are 40–62 % CNV-low tumour epithelium + 12–14 % Malignant |
| Cell–cell communication | Tumour-enriched neutrophil VEGFA/PLAU signalling and the TIGIT–NECTIN2 checkpoint axis; Tensor-cell2cell rank 6 (sample-level programmes, not tumour-specific); spatial ECM → CD44/LRP1 co-expression |

**Known issues** ([§3.10](DEVELOPMENT_LOG.md#310--post-run-review--open-issues-2026-11-11)): run 1's three
high-severity issues (Milo contrast, tumour cells hidden in `Unknown`, a 66-cell rare-epithelial cell2location
reference) are fixed in run 2; before/after in [DEVELOPMENT_LOG §3.11](DEVELOPMENT_LOG.md) and
[docs/v3/FIXES_RUN2.md](docs/v3/FIXES_RUN2.md). Still open and documented: tensor factors don't separate tumour
from adjacent tissue, SOLO and Scrublet disagree, and cross-atlas mixing is below random within cell types.

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
