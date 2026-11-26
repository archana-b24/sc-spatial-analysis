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
 LuCA (LUSC patients) ─────┘   (§3.3 A–E)          (mito, MAD,      scVI/scANVI/sysVI/      per-patient├─ paired composition, Milo
                                                    Scrublet)        Harmony → choose        calls      ├─ trajectory (Palantir + CellRank 2)
                                                                                                        └─ cell–cell communication (LIANA+)
 Visium FFPE LUSC section ─ identity gate ─ QC ─ cell2location ─ SVGs, neighbourhoods, CellCharter niches, spatial LR
```

**Design in one line:** pool both atlases for *representation* (map, labels, references), but draw *inferences*
only within patients — tumour vs adjacent-normal tissue of the same person, replicated across studies — because
healthy donors and tumour patients never come from the same study ([DEVELOPMENT_LOG §3.12 A](DEVELOPMENT_LOG.md)).
The dashboard is split accordingly into **Evidence** (same-patient comparisons: DE & pathways, Composition) and
**Context** (descriptive: malignant cells, trajectory, spatial, cell–cell communication, and every comparison with
healthy donors) ([§3.13](DEVELOPMENT_LOG.md)).

| Question | Method (v2, run 4) | Replaces (v1) |
|---|---|---|
| Which cells are malignant? | infercnvpy per patient × tissue, same-sample held-out null, two axes (signal + correlation), pre-registered acceptance check | — (not done) |
| What changes in each cell type? | Pseudobulk PyDESeq2 on patient × tissue units, paired tumour vs adjacent-normal; thresholded test, per-study model, leave-one-study-out | cell-level Wilcoxon (pseudoreplication) |
| Which pathways / TFs? | decoupler 2 — PROGENy, CollecTRI, Hallmark | GSEA with cells as samples |
| Composition shifts? | paired NB GLM within the immune compartment + propeller + CLR; Milo (paired) | — |
| Basal → malignant axis? | Palantir + CellRank 2 (CytoTRACE + pseudotime kernels); scVI sensitivity | PAGA graphs |
| Where are cells in tissue? | cell2location deconvolution, CellCharter niches | Scanorama + cosine label per spot |
| Who talks to whom, and what changes? | LIANA+ per sample, differential CCC, Tensor-cell2cell on paired slices, spatial bivariate | LIANA on tumour cells only |

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

### Results dashboard

An interactive Streamlit site with one page per pipeline stage. Each page pairs the decision with the evidence
behind it and shows the known caveats. It reads a DuckDB snapshot that the `export` step derives from pipeline
outputs (design: `docs/ROADMAP.md` §6).

```bash
pixi run -e dashboard export      # pipeline outputs → results/dashboard/lusc_atlas.duckdb (~1 min)
pixi run -e dashboard app         # http://localhost:8501
```

Also included (design notes in `docs/design/README.md`): a *Start here* guide and computed reader guide on every page, plus the
atlas's read-only DuckDB tools exposed as an MCP server. (The in-app chatbot is parked: the project stays free to
host and run, and a chatbot needs a paid LLM API key.)

```bash
pixi run -e dashboard mcp                                                  # same tools as an MCP server
pixi run -e dashboard cloud-dry-run                                        # rehearse the free deploy locally
```

Hosting is free-only: the app on **Streamlit Community Cloud** (deployed from this public repo, entrypoint
`app/streamlit_app.py`, deps `app/requirements.txt`) and the DuckDB as a **public Hugging Face Dataset**
(`scripts/publish_dataset.sh`, sha256-pinned at app start). Every service used, its free tier and what would
trigger a cost: `docs/design/C_dashboard_hosting.md` → *Cost: $0*.

## Data

All inputs are downloaded by the pipeline from pinned, public sources and checksum-verified.

| Data | Source | Version |
|---|---|---|
| Human Lung Cell Atlas, core + extended (Sikkema et al. 2023, *Nat Med*) | CELLxGENE Census | `2025-11-08` |
| Single-cell Lung Cancer Atlas, extended (Salcher et al. 2022, *Cancer Cell*) | CELLxGENE Census | `2025-11-08` |
| Visium CytAssist FFPE human lung squamous cell carcinoma (`CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma`) | 10x Genomics public dataset | Space Ranger 2.0.0 |
| GENCODE v44 basic annotation | EBI | v44 |

Why both atlases: HLCA contains only 3 LUSC donors from a single study, so it cannot support a
tumour-vs-normal comparison alone. LuCA contributes most LUSC patients, across several studies. HLCA remains the
healthy-lung reference. Merging them safely (deduplication, gene space, label crosswalk, nested batches,
over-correction checks) is described in [DEVELOPMENT_LOG §3.3](DEVELOPMENT_LOG.md#33--merging-two-atlases--what-can-go-wrong-and-what-we-do-about-it).

## Results (run 4, 2026-11-17)

Full auto-generated report: [`results/reports/summary.md`](results/reports/summary.md) (results are not
committed; regenerate with `pixi run lusc-atlas report`). What changed from run 2 and why:
[docs/design/OVERHAUL_RUN3.md](docs/design/OVERHAUL_RUN3.md) §5; interpretation and caveats:
[DEVELOPMENT_LOG §3.12](DEVELOPMENT_LOG.md); run 4 (healthy = lung tissue from people without cancer; two checks
of the paired test) and its results: [§3.13](DEVELOPMENT_LOG.md). Run 4 changed only Context results, by design:
the paired samples, and so the Evidence, are the same as in run 3.

| Stage | Result |
|---|---|
| Cohort | 128 donors, 15 studies, 219,557 cells selected → 121 donors, 200,848 cells after QC. Healthy: 87 donors selected, 81 after QC (lung tissue from people without cancer, 9 studies). 15 LUSC patients with tumour and adjacent-normal tissue, from 4 studies (the paired design) |
| QC | 200,848 of 217,462 cells kept (per sample × compartment MAD; 5,655 high-mitochondrial cells; 898 Scrublet doublets) |
| Integration | scANVI chosen on held-out labels (scIB total 0.683; sysVI 0.670, scVI 0.648, Harmony 0.608, PCA 0.562) |
| Malignant cells | 16,174 CNV-called cells; 29 patients with ≥ 10, with the expected 3q gain / 3p loss — but the **pre-registered false-positive check failed** (7.9 % of adjacent-normal epithelium called, limit 5 %), so every malignant-cell result is provisional |
| Paired DE | 10 immune cell types tested at the primary tier; 126 robust changes (> 1.5-fold at a global FDR of 5 %, no study contradicting, surviving leave-one-study-out, not a contamination suspect), e.g. *CXCL13*, *ENTPD1* in CD8 T cells, *SPP1*, *MMP12* in macrophages. No epithelial or stromal type has enough paired patients |
| Checks of the paired test | label-swap (sign-flip) null: 0.51 % of genes at p < 0.01 on average, no cell type above the 3 % flag, so the paired test is not too liberal; a per-patient table shows what each adjacent-normal control contains |
| Cross-study null | two sets of *healthy* studies (lung tissue only) differ in 6.4–60.8 % of genes at p < 0.01: comparisons across studies (tumour vs healthy) are Context, hypothesis-generating only |
| Composition | within patients, tumour has more regulatory T, B and plasma cells and fewer alveolar macrophages and non-classical monocytes (three models agree); Milo agrees in direction but no neighbourhood passes FDR |
| Trajectory | malignant cells sit further from healthy basal cells (median pseudotime 0.28 vs 0.19); descriptive, root and endpoint come from different atlases, and the ordering depends substantially on the embedding (scVI sensitivity ρ = 0.59) |
| Spatial | 10x LUSC section (3,816 spots after QC) deconvolved into 26 cell types; 8 CellCharter niches, two of them 58–67 % tumour epithelium and one 30 %; n = 1 patient |
| Cell–cell communication | strongest paired changes involve neutrophils (VEGFA, PLAU; descriptive tier — 5 patients, 1 study); the paired communication tensor (11 patients) shows no programme that differs between tissues |

**Known issues:** listed with severity and status in the dashboard (*Caveats & limitations* on each page) and in
[DEVELOPMENT_LOG §3.10/§3.12](DEVELOPMENT_LOG.md). The ones that limit conclusions: the paired design is
immune-only, four studies carry it, malignant calls failed their false-positive check, and the spatial data are
one patient.

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
- Malignant calls are CNV-based; samples with few CNVs may be under-called, and the calls of runs 3 and 4 failed
  the pre-registered adjacent-tissue false-positive check.
- Only immune cell types can be compared within patients; tumour-epithelium results cross atlases.
- Pseudotime across donors describes transcriptional distance from healthy basal cells, not lineage.

## Authors

Archana Balasubramanian, Shivaramakrishna Srinivasan, Aaryan Sharma, Aathira Sarath Chandran, Hersh Gupta.

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the HLCA, LuCA and the methods listed in the
development log.
