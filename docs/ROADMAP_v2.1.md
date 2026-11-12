# Roadmap — v2.1 candidates

> Changes **considered but deliberately deferred** from v2.0.0, with the reason they were deferred and what
> would make them worth doing. v2.0.0 scope and decisions live in [DEVELOPMENT_LOG.md](../DEVELOPMENT_LOG.md).
> An item moves into a release only when it answers a concrete scientific question with an adequate design.

Status legend: 💡 idea · 🔍 scoping · ✅ accepted for v2.1 · ❌ rejected (with reason)

| # | Item | Status | Effort | Depends on |
|---|---|---|---|---|
| 1 | Histology CNN / pathology foundation-model features ↔ expression | 💡 | M | more Visium sections (ideally) |
| 2 | Additional LUSC Visium / Xenium sections | 💡 | M | public data search |
| 3 | Fine-tuned (not zero-shot) foundation-model classifier for malignant cells | 💡 | M–L | GPU |
| 4 | Ambient-RNA sensitivity analysis | 💡 | S | — |
| 5 | Survival / TCGA-LUSC projection of malignant programmes | 💡 | M | — |
| 6 | **Results dashboard** — DuckDB results store + Streamlit app on Hugging Face Spaces (§6) | ✅ | M–L | v2.0 outputs |
| 7 | v3 groundwork: Claude Skill wrapping the pipeline steps | 💡 | M | v2.0 release |

---

## 1. Histology CNN ↔ spatial expression (deferred from v2.0 on 2026-11-09)

**Idea.** Extract a feature vector for the H&E image patch under each Visium spot with a pretrained image
model, then ask: *does tissue morphology predict the cell2location abundances (malignant, immune, stromal) and
the CellCharter niches?* Agreement would validate the deconvolution; disagreement would flag spots where
expression and morphology tell different stories.

**Candidate encoders.**
- ImageNet ResNet-50 via `squidpy.im.calculate_image_features` — trivially available, weak for pathology.
- Pathology foundation models (UNI / UNI2, CONCH, Virchow, H-optimus) — much better representations, but
  model weights are gated (licence acceptance) and inference wants a GPU.

**Design.** Full-resolution image (`*_tissue_image.tif`, 3.8 GB, 10x 2.0.1 release) → 55 µm patches at spot
centres → frozen encoder → ridge / gradient-boosted regression to abundances. **Spatially blocked
cross-validation** (contiguous tissue blocks held out), otherwise neighbouring spots leak information and
accuracy is grossly inflated.

**Why deferred.** With **one** section, the model is trained and evaluated on the same tissue from the same
patient; even with blocked CV it cannot show generalisation, so any number is exploratory. It also adds a GPU-
and licence-dependent component to an otherwise CPU-reproducible pipeline.

**What would make it worth doing.** ≥ 3 LUSC sections from different patients (item 2) for leave-one-section-out
evaluation, or access to a GPU + a pathology foundation model licence.

**Sketch of deliverables.** `src/lusc_atlas/histology.py` (patching + encoder), step `histology` in its own pixi
env, `results/tables/histology_cv.csv`, figure of predicted vs observed abundance maps.

## 2. More spatial sections
Search 10x / GEO / HTAN / Zenodo for public LUSC Visium, Visium HD or Xenium data. Turns spatial results from
descriptive (n = 1) into comparative, and unlocks item 1.

## 3. Fine-tuned foundation model for malignant-cell calling
v2.0 benchmarks Geneformer **zero-shot** embeddings (integration benchmark). A fine-tuned classifier on CNV-
derived malignant labels — validated against LuCA author labels and on held-out studies — would test whether
foundation models add value over inferCNV. Needs a GPU for reasonable runtime.

## 4. Ambient-RNA sensitivity
Atlas data are processed counts, so CellBender/SoupX cannot run. A sensitivity analysis — re-running the primary
DE after removing genes with high ambient signatures (e.g. SFTPC, SCGB1A1, HB genes in non-expressing types) —
would quantify how much ambient contamination could move the conclusions.

## 5. TCGA-LUSC projection
Score malignant programmes from v2 on TCGA-LUSC bulk RNA-seq (decoupler) and test association with survival and
the TCGA molecular subtypes (classical / basal / secretory / primitive).

## 6. Results dashboard (accepted for v2.1 on 2026-11-10) — 🚧 built and running locally (2026-11-11), not yet deployed

**Goal.** Make the analysis *readable*. Anyone, with or without a programming background, should be able to see what is in
the data, what each step did and **why**, and the evidence behind each decision. This is the interactive form of
DEVELOPMENT_LOG §3 and §5: each page pairs a decision with the data that justified it.

**Decisions (user, 2026-11-10).** Public hosting on **Hugging Face Spaces** (free tier, 16 GB RAM). All inputs are
public atlases, so derived results can be public. Python-only stack: **Streamlit**. R Shiny would mean
reimplementing the stack in R; Dash stays the fallback if we need linked multi-view interactions.

### Data architecture
```
pipeline steps ──► results/tables/*.csv, *.h5ad, *.zarr     (source of truth: reproducible from code + config)
        │
        └─► new step `export` ──► results/dashboard/lusc_atlas.duckdb   (derived, versioned serving snapshot)
                                          │
                                          └─► Streamlit app (read-only)
```
- **Why DuckDB is a serving layer, not the source of truth.** It suits tables (columnar, SQL, one file, no
  server), but the 249k × 17k sparse expression matrix and the SpatialData objects belong in `.h5ad`/`.zarr`.
  Making the DB the primary store would let it drift from the code that produces it. Instead, `export` rebuilds it
  deterministically from pipeline outputs, and every row traces back to a run.
- **Planned schema** (one table per result CSV, plus):

| Table | Grain | Approx. size | Feeds |
|---|---|---|---|
| `runs` | one row per pipeline run: git SHA, config sha256, step runtimes (from `reports/provenance/`) | tiny | footer on every page |
| `cells` | one row per cell: sample, donor, study, atlas, condition, cell type, malignant flag, QC metrics, UMAP x/y | 249k rows, ~50 MB | UMAP explorer, cohort |
| `samples` / `donors` | cohort metadata, pairing | ~200 rows | cohort page |
| `expr_pseudobulk` | gene × cell type × condition: mean log-expression, fraction expressing | ~17k × 27 × 4 | gene lookup / dot plots |
| `visium_spots` | spot x/y, QC, cell2location abundances, niche | 6k rows | spatial page |
| `decisions` | step, question, decision, rationale, evidence link (curated from DEVELOPMENT_LOG) | ~30 rows | "why" panels |
| result tables | DE, pathways, scCODA, Milo, CNV, trajectory, CCC (by sample, differential, tensor, spatial), integration benchmark | from CSVs | per-stage pages |

### Pages (one per pipeline stage, each answering *question → decision → evidence → caveats*)
0. **Overview**: the question, cohort at a glance, the pipeline as a diagram, and links to each stage
1. **Cohort and merge**: HLCA + LuCA, study crosswalk, exclusions, pairing (why HLCA alone was not enough)
2. **QC**: per-sample MAD thresholds, doublet agreement, cells kept/removed per sample
3. **Integration**: the 7-candidate benchmark (scIB metrics), UMAP coloured by atlas, study or cell type, and why the winner was chosen
4. **Malignant cells**: CNV scores, threshold from healthy epithelium, chr3q positive control, agreement with author labels
5. **DE and pathways**: volcano per cell type (paired design), PROGENy / CollecTRI / Hallmark heatmaps, gene lookup
6. **Composition**: scCODA credible effects, Milo neighbourhoods
7. **Trajectory**: pseudotime by state, with the "not lineage tracing" caveat
8. **Spatial**: Visium H&E with deconvolved abundances and niches (n = 1, descriptive)
9. **Cell–cell communication**: differential LR, tensor factors vs condition, spatial LR co-expression, and the Malignant-not-in-tensor limitation
10. **Provenance**: runs, versions, config, environment

### Engineering notes
- 249k-point UMAP: Plotly `scattergl` with a 50k stratified subsample by default, full set on demand (or datashader).
- `st.cache_resource` for the read-only DuckDB connection; the app never writes.
- The DuckDB file (est. < 500 MB) ships in the HF Space via Git LFS, or is downloaded from an HF Dataset repo at startup. The second
  option keeps the app repo small and versions data separately (preferred).
- App code lives in `app/` in this repo. CI builds the DB from a smoke-test run and checks every page renders (`streamlit.testing`).
- No reimplementation of a general single-cell browser. For free-form cell-level exploration, link out to cellxgene
  (processed `.h5ad`); consider embedding Vitessce for the spatial page.

### Build order
1. `export` step → DuckDB + schema test → 2. app skeleton (overview, provenance) → 3. stage pages in pipeline
order → 4. HF Space + HF Dataset deploy → 5. README badge/link.

### Open questions → defaults taken while building (2026-11-11, open to revision)
- **Expression at cell level:** both. Pseudobulk (mean log-expression and fraction expressing per gene × cell type ×
  condition) for **all** genes and cells, so the dot plot covers any gene. Per-cell values only for a marker panel
  (canonical markers + top-15 DE genes per contrast × group + top-500 HVGs ≈ 1.1k genes) on a 50k stratified
  subsample, for UMAP colouring. Every cell type × condition stratum is represented. Size knobs live in
  `config.yaml → dashboard`.
- **`decisions` source:** curated by hand in `docs/dashboard_content.yaml`, with every entry citing its
  DEVELOPMENT_LOG section. Parsing the log would need a machine-readable block in the log first. That is still the
  better single-source design, deferred until the log's §3 settles after the §3.10 fixes.

### As built (2026-11-11)
| Piece | Where | Notes |
|---|---|---|
| `export` step | `src/lusc_atlas/export.py`, registered in `steps.py` (runs in the new pixi env `dashboard`) | 44 tables, **118 MB**, 42 s. Atomic write (`.partial` → `os.replace`). Records git SHA / dirty flag, config, summary.md and per-table row counts plus any reduction (`tables_info`). Per-sample LIANA shipped only at magnitude rank ≤ 0.01 (2.17 M → 248k rows) |
| pixi env `dashboard` | `pixi.toml` (Python 3.12: streamlit, duckdb, plotly, anndata, h5py) | separate env, so the analysis lock is untouched (lock diff = insertions only). `app/requirements.txt` = app-only subset for HF Spaces |
| App | `app/streamlit_app.py`, `app/views/*.py` (13 pages), `app/lib/` | overview, known issues, 9 stage pages, cell & gene explorer, provenance. Each stage page = question → "why" expander → known-issue banners → evidence charts → table view |
| Tests | `tests/test_dashboard.py` | schema checks + every page renders through the real navigation (`streamlit.testing`). Skipped when streamlit or the DB is absent, so default-env CI is unchanged |
| Colour | `app/lib/ui.py` | validated categorical palette in fixed order. Condition and atlas colours are constant across pages. 27 cell types are shown by compartment colour or highlight-one-grey-the-rest, never 27 hues |

Run locally: `pixi run -e dashboard export` (or `run_all.py --only export`), then `pixi run -e dashboard app` → http://localhost:8501.

### Remaining build steps
1. Review pass on wording and figures (user).
2. HF Dataset repo for the DuckDB file plus an HF Space for `app/` (`LUSC_DB_HF_REPO` env var is already supported in `app/lib/db.py`).
3. CI job: dashboard env + `export` on smoke-test outputs + `tests/test_dashboard.py`.
4. README badge/link.
5. Rebuild the snapshot after the §3.10 fixes. The known-issue banners update from `docs/dashboard_content.yaml`.

## 7. v3 groundwork
See DEVELOPMENT_LOG §4. **2026-11-12: three v3 concepts prototyped — agentic "Ask the atlas" (A), NVIDIA Nemotron
assessment (B), dashboard explainers + single-source-of-truth enforcement + hosting bundle (C) — index, verdicts
and open questions in [docs/v3/README.md](v3/README.md).** Hosting note for item 6: Hugging Face now requires a
paid plan for Docker Spaces; options compared in [docs/v3/C_dashboard_hosting.md](v3/C_dashboard_hosting.md).

---

## Changelog of this document
| Date | Change |
|---|---|
| 2026-11-09 | Created; item 1 (histology CNN) deferred from v2.0 by decision |
| 2026-11-10 | Item 6 widened from "results browser" to the results dashboard (DuckDB + Streamlit on HF Spaces); accepted for v2.1, design under review |
| 2026-11-12 | §7 points to docs/v3 (concepts A–C); item 6 hosting re-evaluated (HF Docker Spaces need PRO; Streamlit Community Cloud as free fallback) |
