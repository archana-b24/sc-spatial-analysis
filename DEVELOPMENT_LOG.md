# Development Log — v2.0.0

> Living record of **what** changed between v1 and v2, **why** it changed, and **where** it lives.
> Every entry is written so a reviewer can audit a decision without reading the code.
> Newest decisions are appended inside each section; the [Timeline](#timeline) at the bottom is chronological.

| | |
|---|---|
| **Baseline** | tag `v1.0.0` — original BENG207 course notebooks (Feb–May 2024), preserved byte-for-byte |
| **Working branch** | `v2-dev` |
| **Target** | tag `v2.0.0` — reproducible pipeline + corrected science with current (2025–26) methods |
| **Next** | AI-native layer (Claude Skills / agent interface) — see [§4](#4-looking-ahead-ai-native-layer) |
| **Long-running jobs** | tmux session `lusc-v2` (windows: `main`, `env`, `pipeline`) — `tmux attach -t lusc-v2` |
| **Storage rule** | Everything (envs, caches, temp, data) lives on `/home/abalasubramanian/storage`; `source ../storage_env.sh` first (see §2.3) |

Status legend: ✅ done · 🚧 in progress · ⏳ planned · ❌ dropped (with reason)

---

## 0. Audit of v1 (why v2 exists)

The v1 audit (2026-11-09) found three classes of problems. Each v2 entry below references the issue ID it resolves.

### Reproducibility breaks
| ID | Problem in v1 | Evidence |
|---|---|---|
| R1 | Notebooks do not chain: files are read that no notebook writes | `03` reads `data/sclcCells.h5ad` but `02` writes `sclcSinglets.h5ad` → **doublet removal never reached downstream analysis**; `05`/`06` read `sclcHV.h5ad` (never created); `07` reads `sclcSingletwGenes.h5ad`, `normalCellsOld.h5ad`, uses undefined `sclcCells_1` |
| R2 | Hardcoded / inconsistent paths | `/Users/shivaramakrishna/...`, `hlca/` vs `data/`, `BENG207/...` |
| R3 | Saved errors & stubs | 11 cells with tracebacks; `04_trajectory` is 3 cells referencing undefined objects |
| R4 | No environment specification | Python 3.9 and 3.11 mixed; `%pip install` inside notebooks; README says "Python ≥3.6" |
| R5 | Manual data access | Google Drive link, no checksums, no provenance |
| R6 | Repo bloat | 33 MB notebook outputs, 25 MB `.h5` tracked despite `.gitignore`, `.DS_Store` files, git pack 227 MB |

### Scientific issues
| ID | Problem in v1 | Consequence |
|---|---|---|
| S1 | "SCLC" used for *squamous* cell lung carcinoma | SCLC conventionally = **small**-cell lung cancer. Squamous = **LUSC**. Misleading to any reader |
| S2 | Control QC used tumour threshold (`upper_lim` instead of `upper_lim_ctrl`) | Control cells filtered with wrong cutoff |
| S3 | Cell-level Wilcoxon / t-test tumour-vs-normal DE | Pseudoreplication — cells from the same donor are not independent; p-values badly inflated |
| S4 | GSEA with `permutation_type='phenotype'` using cells as samples | Same pseudoreplication problem in pathway statistics |
| S5 | No batch / study correction across HLCA studies | Tumour and normal come from different studies → condition confounded with batch |
| S6 | Doublets called on tumour only, without batch key | Asymmetric QC between conditions |
| S7 | No malignant-cell identification | "Progression" cannot be studied without separating malignant from normal epithelium |
| S8 | Trajectory = PAGA graphs only | No pseudotime / fate inference |
| S9 | Visium "integration" = Scanorama + cosine label transfer | Visium spots (55 µm) contain multiple cells — requires deconvolution, not 1:1 labels |
| S10 | No spatial statistics | No SVGs, neighbourhood enrichment, niches, spatial CCC |
| S11 | LIANA on tumour only | No differential (tumour vs normal) communication |
| S12 | Deprecated tools | `diffxpy` (unmaintained), legacy `decoupler`/`liana` APIs |

### Engineering gaps
| ID | Gap |
|---|---|
| E1 | No package / functions duplicated across notebooks |
| E2 | No tests, no CI |
| E3 | No LICENSE, CITATION, CHANGELOG, versioning |
| E4 | No workflow manager / single entry point |

---

## 1. Architecture

### 1.1 ✅ Preserve v1, then restructure — *resolves R6, E3*
- **What:** Tagged current `main` as `v1.0.0`; moved all v1 notebooks to `legacy/v1/` (outputs stripped), v1 report to `docs/v1_report.pdf`; untracked `.DS_Store`, Space Ranger CSVs, `.h5`, result CSVs.
- **Why:** v2 rewrites nearly everything. A tag gives an immutable, citable reference for "what v1 did" so we can strip, move and delete freely without losing history. Stripped legacy notebooks keep the code readable while dropping ~50 MB of embedded images.
- **Not done (needs sign-off):** rewriting git history to purge the 227 MB pack. Requires a force-push to `origin` and coordination with co-authors.

### 1.2 ✅ Layout: package + generated runners
```
sc-spatial-analysis/
├── config/config.yaml          # every path, threshold, seed — single source of truth
├── src/lusc_atlas/             # importable package: all analysis logic lives here
├── workflow/                   # Snakemake: Snakefile + rules/*.smk + envs/
├── scripts/run_all.py          # dependency-free sequential runner (same steps as Snakemake)
├── notebooks/                  # thin, narrative notebooks that only call the package and plot
├── tests/                      # unit tests + tiny synthetic fixtures for CI
├── data/{raw,interim,processed}/  # never committed; created by the fetch step
├── results/{tables,figures,reports}/
├── legacy/v1/                  # frozen v1 code for comparison
└── docs/
```
- **Why package + thin orchestrators:** v1 duplicated logic across notebooks, and they drifted (R1). With one function per analysis step, Snakemake, `run_all.py`, notebooks and tests all call the *same* code — two runners cannot disagree.
- **Why both Snakemake and `run_all.py`:** Snakemake gives incremental re-runs, provenance and cluster execution; `run_all.py` is the zero-learning-curve path for someone who just wants `python scripts/run_all.py`. Both read `config/config.yaml`.
- **How they stay identical:** `src/lusc_atlas/steps.py` declares each step's `INPUTS`, `OUTPUTS` and pixi `ENV`. The Snakefile *generates* one rule per step from those declarations; `run_all.py` iterates the same `ORDER` and skips steps whose declared outputs exist. A unit test asserts every input is produced by an earlier step — the exact failure that broke v1 (R1).
- **Provenance:** every step writes `results/reports/provenance/<step>.json` (runtime, Python, package versions, sha256 of the resolved config).

| Step | Env | Purpose |
|---|---|---|
| `fetch_visium` | default | 10x download + sha256 + extract |
| `fetch_atlases` | default | HLCA + LuCA source h5ads from pinned Census; size + sha256 |
| `select` | default | two-atlas dedup, conditions, count check, gene space, labels, caps (§3.3) |
| `qc` | default | per-sample MAD QC + Scrublet ∧ SOLO doublets |
| `integrate` | default | PCA / Harmony / scVI / scANVI / sysVI → scib benchmark + guard → choose |
| `cnv` | default | infercnvpy malignant calls, chr3q control, author-label concordance |
| `de` | default | pseudobulk PyDESeq2 contrasts + decoupler pathways/TFs |
| `abundance` | default | scCODA + Milo |
| `trajectory` | default | Palantir + CellRank 2 on Basal → Malignant |
| `spatial_prep` | default | Visium QC/clustering + SpatialData export |
| `deconvolve` | **c2l** | cell2location reference signatures + spot mapping |
| `spatial` | default | SVGs, neighbourhood enrichment, CellCharter niches |
| `ccc` | default | LIANA+ by sample, differential, Tensor-cell2cell, spatial bivariate |
| `report` | default | figures + `results/reports/summary.md` from outputs only |

### 1.3 ✅ Dashboard hosting: free services only (decided 2026-11-12)
- **Decision (user):** "I have to keep everything free … I don't get sucked into a paid plan." Nothing in the repo
  may create, require or nudge towards a paid plan or paid API.
- **What changed:** the Hugging Face **Docker Space** path (`build_space.py`, `deploy_space.sh`,
  `space_dry_run.sh`, tasks `space-build`/`space-dry-run`) was removed, because creating a Docker or Gradio Space
  now needs HF PRO. The scripts are parked outside the repo (`../.tooling/parked_paid_deploy_2026-11-12/`).
  Hosting = **Streamlit Community Cloud** (free; deployed by the user from the public GitHub repo, entrypoint
  `app/streamlit_app.py`, deps `app/requirements.txt`, Python 3.12) + the DuckDB as a **public HF Dataset**
  (free), downloaded at start-up through `hf_hub_download`, pinned to the upload's commit and verified by sha256.
  `app/lib/db.py` reads the `LUSC_DB_*` settings from env **or** `st.secrets` (Community Cloud's Secrets box).
  `.streamlit/config.toml` moved to the repo root (Community Cloud's working directory).
- **New scripts:** `scripts/publish_dataset.sh` (user runs it with `HF_TOKEN`: public Dataset only, refuses private
  repos and files > 1 GB, never creates a Space, prints the secrets TOML) and `scripts/cloud_dry_run.sh` (clean
  clone of exactly what git would push → fresh Python 3.12 venv from `app/requirements.txt` → secrets TOML → DB via
  `hf_hub_download` from a seeded offline cache → every page rendered → wrong-sha negative control → health check).
- **Options considered:** HF Docker Space (paid, rejected), MotherDuck (free tier exists but adds an account and a
  network DB for a 119 MB read-only file; rejected), bundling the DB in git (119 MB derived binary in history;
  rejected), Community Cloud + HF Dataset (**chosen**).
- **Evidence:** rehearsal passes, peak RSS 444 MB vs Community Cloud's 690 MB–2.7 GB allocation (docs, Feb 2024
  figures). Every service, its free tier and what would trigger a cost: `docs/design/C_dashboard_hosting.md` → *Cost: $0*.

---

## 2. Reproducibility

### 2.1 ✅ Long-running work in tmux
- Session `lusc-v2`: `main` (interactive), `env` (environment builds), `pipeline` (pipeline runs). Survives SSH disconnects; attach with `tmux attach -t lusc-v2`.

### 2.2 ✅ Data acquisition with provenance — *resolves R5*

| Source | Pinned version | Size | How it is verified |
|---|---|---|---|
| HLCA core+extended (Sikkema 2023) | CELLxGENE Census `2025-11-08` (stable/LTS), dataset `9f222629-…` | 21.9 GB | byte size = S3 `Content-Length`; sha256 recorded in `data/raw/hlca/provenance.json` |
| LuCA extended (Salcher 2022) | Census `2025-11-08`, dataset `1e6a6ef9-…` | 17.6 GB | same |
| Visium CytAssist FFPE lung SCC | 10x public release **2.0.1** | 64 MB + 34 MB | sha256 pinned in `config/config.yaml` |
| GENCODE v44 basic GTF | EBI FTP | — | gene positions for inferCNV |

- **Finding:** the v1 Visium matrix (`sclcVisiumFiltered.h5`, 25 MB) is **not** byte-identical to the official 10x release (64 MB; different sha256; v1 also used the older `tissue_positions_list.csv` layout). v2 uses the official release as the canonical, publicly verifiable source; the v1 file remains in tag `v1.0.0` for comparison only.
- 10x's `2.0.0` path returns HTTP 403 to scripted clients; the `2.0.1` path is public.
- **Why pinned Census releases:** "latest" changes monthly; a pinned LTS build means the same cells, same metadata, forever.

### 2.3 ✅ Pinned environments (pixi) — *resolves R4*

- **Incident (2026-11-09):** the first env build (conda, in `~`) filled the 250 GB root volume. Rule since then: **all** envs, package caches, temp files and data live on `/storage`. `storage_env.sh` (one level above the repo) exports `PIXI_HOME`, `PIXI_CACHE_DIR`, `UV_CACHE_DIR`, `PIP_CACHE_DIR`, `CONDA_PKGS_DIRS`, `TMPDIR`, `MPLCONFIGDIR`, `NUMBA_CACHE_DIR`, `XDG_CACHE_HOME`. Every tmux window sources it.
- **Why pixi** (over conda + `environment.yml`): one manifest (`pixi.toml`) + one lockfile (`pixi.lock`) for several *isolated* environments, conda and PyPI resolved together, and exact reproduction with `pixi install --frozen`.

| Env | Python | Purpose | Why separate |
|---|---|---|---|
| `default` | 3.13 | the whole analysis | pertpy 1.5 requires ≥ 3.13 |
| `c2l` | 3.12 | cell2location only | needs scvi-tools 1.3.x + pyro < 1.9.2; main stack uses scvi-tools 1.5 |
| `workflow` | 3.13 | Snakemake only | orchestrator must not share deps with analysis (see conflicts) |

Conflicts found and how they were resolved (each was a real solver failure):

| Conflict | Resolution |
|---|---|
| `pygpcca` (← cellrank) hard-pins `jinja2==3.0.3`; conda-forge `pytorch` and `snakemake` need `jinja2≥3.1` | torch from the PyPI **CPU** wheel index; snakemake in its own env |
| `cellcharter 0.3.9` needs `scipy≥1.18`; `pygam` (← cellrank) needs `scipy<1.17` | cellcharter left floating → solver picks 0.3.7 |
| `snakemake 9.27` needs `packaging<26`; scanpy env had 26.3 | resolved by moving snakemake out (above) |
| `cell2location 0.1.5` has no upper bound on scvi-tools and breaks on 1.5 | `c2l` env pins `scvi-tools>=1.3,<1.4` |
| `pyro-ppl 1.9.2` removed `deep_getattr` used by cell2location | `c2l` env pins `pyro-ppl<1.9.2` |

Verification: `scripts/check_env.py` imports every package in each env (all pass).

Missing optional extra found by the full run (2026-11-10): `pertpy` imports cleanly without JAX, but `pt.tl.Sccoda()` raises `ImportError` at call time — scCODA's NUTS sampler needs the `[jax]` extra (numpyro, flax, optax, ott-jax). `check_env.py` could not catch this because it only imports packages. Fix: `pertpy = { version = ">=1.5", extras = ["jax"] }`; reinstalling left torch, scanpy, numpy, scipy and anndata on the same versions. Lesson for `check_env.py`: also instantiate tools whose dependencies are loaded lazily.

Two more found by the resumed run (2026-11-10):
- **`tabulate` missing**: `report` uses `DataFrame.to_markdown()`, which needs the optional dependency `tabulate`. It was never exercised because the smoke test stops before `report`. Added to `pixi.toml`; nothing else changed.
- **LIANA+ API drift (1.x → 2.x)**: `pixi.toml` already required `liana>=2.1.1`, but `ccc.py` was written against the 1.x API (`li.multi.*`, `li.ut.*`), which 2.x reorganised. Decision: **port to 2.x, not pin to 1.x**. The manifest already chose 2.x, 1.x is no longer developed, and pinning back risks new conflicts with the Python 3.13 stack. Mapping: `li.multi.df_to_lr` → `li.mt.df_to_lr` (now expects the DE table *indexed by gene* with the groupby column, which the old code did not supply); `li.multi.to_tensor_c2c` → `li.ms.to_tensor_c2c` (`inverse_fun` → `inverse_fn`, whose default changed to −log10, so `1 − x` stays explicit); `li.ut.spatial_neighbors` → `li.pp.spatial_neighbors`. Lesson: as with scCODA above, import checks can't catch API drift. Only exercising a step does, so `ccc` and `report` should join the smoke test (v2.1).

### 2.4 ✅ Notebook & repo hygiene — *resolves R6*
`.pre-commit-config.yaml`: **nbstripout** (outputs can no longer be committed), **ruff** lint + format, `check-added-large-files` (2 MB) so data cannot re-enter git. Data/results directories are git-ignored.

### 2.5 ✅ Tests + CI — *resolves E2*
`tests/` runs in ~1 s on synthetic data — no downloads, no training. What is tested is chosen by *what would silently corrupt results*: DAG wiring (R1), per-sample QC (S2), design collinearity and paired blocking (S3/S5), the over-correction guard, chr3q window mapping, duplicate-study detection, origin import, raw-count detection, stratified caps, checksums. `.github/workflows/ci.yml` runs lint, tests, the Snakemake dry-run and the `run_all` plan with `pixi install --frozen`.

### 2.6 ✅ Resumable runs — a step is "done" only if it finished cleanly
**Problem.** `run_all.py` originally skipped a step whenever its declared output files existed. That would have caused two silent failures on resume:
1. *Truncated outputs:* a step killed while writing an `.h5ad` or CSV leaves a file that exists but is incomplete. Resume would skip the step, and every later step would read corrupt data.
2. *Stale outputs:* outputs from an older run (e.g. before a config or code change) are treated as current.

This matters here because the pipeline runs unattended for 20+ hours (`deconvolve` alone took 10.6 h on the full data), so being interrupted and resumed is the normal case.

**Fix (two independent layers):**

| Layer | Where | What it guarantees |
|---|---|---|
| Atomic writes | `steps._write_table`, `steps._write_h5ad` | write to `.<name>.partial<ext>` in the same directory, then `os.replace()`. On one filesystem the rename is atomic, so a declared output is always either complete or absent. |
| Completion stamps | `scripts/run_all.py` | `logs/run_all/.done/<step>` is **deleted when the step starts** and **written only after exit code 0** (it records the timestamp and runtime). `done(step)` = stamp exists **and** all declared outputs exist. A crash, kill or non-zero exit therefore always means re-run. |

- `--adopt STEP …` stamps steps whose outputs were produced before stamps existed. It refuses if any output is missing. It was used once (2026-11-10) for the 12 steps that had completed in the first full run, each checked against its ✓ line in `logs/run_full.log` or a traceback-free step log. It is never needed for new runs.
- `--force` still re-runs regardless, and `--dry-run` shows the plan using the same logic.
- The **Snakemake** path is unchanged: Snakemake already marks a job's outputs as incomplete when the job fails or is killed, and re-runs it. The atomic writes protect that path too.
- Run logs are appended to (`tee -a logs/run_full.log`, with a `=== resume <timestamp>` header), so the history of a multi-session run stays in one file.
- Not covered: outputs written by third-party code straight to their final path (the `.npz` from `fm_embed`, cell2location's tables, the SpatialData `.zarr`). The stamp layer still protects these, since the step has to exit 0 for its stamp to be written.

---

## 3. Correct science

Each entry: **v1 → v2**, **why**, **evidence/citation**.

### 3.1 🚧 Naming: SCLC → LUSC — *resolves S1*
All v2 code, files and figures say **LUSC** (lung squamous cell carcinoma). "SCLC" is reserved in the field for small-cell lung cancer.

### 3.2 ✅ Cohort: HLCA alone cannot answer the question → HLCA + LuCA

**Finding (2026-11-09, from the HLCA source h5ad):** HLCA's entire LUSC content is **3 donors from one study** (`Thienpont_2018`, i.e. Lambrechts et al. 2018), 20,631 cells, of which 84 % are immune and only ~1,950 epithelial (442 basal). No normal donor shares that study. v1's tumour-vs-normal comparison was therefore *n = 3 vs other studies' donors*, with condition perfectly confounded with study — a fact v1 never surfaced.

**Decision (user, 2026-11-09):** combine **HLCA** (healthy-lung reference, 1.3 M normal cells, gold-standard harmonised annotation) with **LuCA extended** (Salcher 2022; 45 LUSC donors / 190 k LUSC cells across several studies, plus normal tissue). This turns an n = 3 anecdote into a multi-study cohort, at the cost of a cross-atlas integration problem that §3.3 treats explicitly.

### 3.3 🚧 Merging two atlases — what can go wrong, and what we do about it

Merging two curated atlases is *not* `anndata.concat`. Each atlas is itself an integration of many studies with its own processing, QC thresholds, gene annotation and cell-type ontology. The checklist below is the design; each item names the failure mode, how we detect it, and how we handle it.

**A. Duplicated cells and donors (the silent killer).**
- *Failure mode:* the same patients appear in both atlases (LuCA re-distributes public datasets that HLCA also contains). Duplicates (i) count one patient twice → pseudoreplication in donor-level statistics, and (ii) form a perfect artificial "batch" between identical cells that any integration method will happily "correct", **inflating integration metrics**.
- *Evidence:* CELLxGENE flags `is_primary_data = False` on **360,842 LuCA cells** (incl. 8,373 LUSC cells) — they are copies of cells primary elsewhere. Donor-ID strings share **zero** overlap between atlases, so naive donor matching would miss every duplicate.
- *Handling:* keep only `is_primary_data == True` cells per atlas, then build an explicit **study crosswalk** (`results/tables/merge_study_crosswalk.csv`) and assert that no study contributes donors from both atlases. Any study present in both is taken from exactly one atlas (the one whose annotation we rely on for that study).

**B. Are the counts actually counts?**
- *Failure mode:* atlases store `X` normalised and raw counts in `.raw`, but some contributing studies only deposited normalised or non-integer data (e.g. ambient-corrected). Count models (scVI, DESeq2, cell2location) silently mis-model non-integers.
- *Handling:* per-study check that the matrix used is non-negative integer counts; studies failing the check are excluded and listed in the crosswalk with the reason.

**C. Gene space.**
- *Failure mode:* symbol-based joins break on renamed genes / GENCODE version drift; per-study gene filtering upstream makes "missing" look like "zero".
- *Handling:* join on **Ensembl IDs (version stripped)**, never symbols; take the **intersection** of genes measured in both atlases for integration; respect CELLxGENE `feature_is_filtered`; report genes lost. (For Visium — a probe-based FFPE panel of ~18 k genes — cell2location uses the sc ∩ Visium intersection.)

**D. Metadata harmonisation.**
- Donor and sample IDs are namespaced by atlas (`HLCA:…`, `LuCA:…`).
- One `condition` vocabulary: `tumour` (LUSC tumour tissue), `adjacent_normal` (non-tumour tissue from LUSC patients, LuCA `origin`), `healthy` (HLCA/LuCA normal donors). The v1 binary tumour/normal hid this distinction, which matters: adjacent-normal lung of a smoker with cancer is *not* healthy lung.
- Tissue restricted to lung primary sites; lymph-node / effusion / metastasis samples excluded.

**E. Cell-type labels are two different ontologies.**
- *Failure mode:* HLCA `ann_finest_level` and LuCA `cell_type` / `cell_type_major` were built independently. Comparing fine labels across atlases without a crosswalk is meaningless; LuCA has tumour-specific states (e.g. tumour cells, TANs) that HLCA lacks.
- *Handling:* an explicit, reviewed **label crosswalk** (`config/label_crosswalk.csv`) maps both schemes to one harmonised coarse level (anchored on Cell Ontology IDs, which both atlases carry). Fine labels are kept per atlas and never compared directly. A semi-supervised **scANVI** model uses the harmonised coarse label (unmapped → `Unknown`) and its predictions are checked against the original labels (confusion matrix per atlas).

**F. Batch structure is nested and crossed.**
- `sample ⊂ donor ⊂ study ⊂ atlas`, with `assay` (10x 3′ v1/v2/v3, 5′, BD Rhapsody, GEXSCOPE, inDrop, Seq-Well…) crossing studies.
- *Choice:* `batch_key = sample` (the finest unit of technical variation), with `study` and `assay` as categorical covariates. **Not** `batch = atlas`: that is the coarsest grouping and would leave within-atlas study effects untouched while risking removal of biology that differs between atlases by design (HLCA ≈ healthy, LuCA ≈ tumour).

**G. Which integration method, and how we decide.**
- Candidates: **scVI** (Lopez 2018), **scANVI** (Xu 2021; label-aware, best bio-conservation in scIB), **sysVI** (Hrovatin et al., VampPrior + cycle-consistency; designed for *substantial* batch effects such as cross-technology/cross-atlas, built into scvi-tools ≥ 1.3), and **Harmony on PCA** as a classical baseline. Unintegrated PCA is the floor.
- *Benchmark:* **scib-metrics** (Luecken 2022) on a stratified subsample — bio-conservation (cell-type ASW, NMI/ARI, isolated labels, cLISI) vs batch removal (iLISI, kBET, graph connectivity, PCR), computed with batch = study.
- *Over-correction guard (the part most pipelines skip):* integration must not erase real biology. Checks: (i) malignant cells must not co-cluster with healthy epithelium; (ii) iLISI computed **within** each harmonised cell type (mixing of the same cell type across atlases is good, mixing of different types is not); (iii) known tumour-only states remain distinct.
- *Decision rule:* the method with the best aggregate scIB score **that passes the over-correction checks** is used downstream; results of all candidates are kept in `results/tables/integration_benchmark.csv`.

**H. Integration is for representation, not for statistics.**
- The integrated latent space is used for neighbours, UMAP, clustering, label transfer, trajectories and Milo. **Differential expression always uses raw counts with batch terms in the model**, never batch-"corrected" expression (corrected values are not counts and their variance is distorted).
- Primary DE contrast is **paired within donor**: LUSC tumour vs adjacent-normal from the *same* patients (`~ donor + condition`), which cancels donor and study effects by design. Tumour vs healthy (`~ study + assay + condition`) is a secondary, explicitly confounded-by-design sensitivity analysis.

**I. Uniform reprocessing.**
Both atlases were QC'd with different thresholds. v2 re-applies one per-sample MAD QC and per-sample doublet detection to every sample from both atlases (§3.4), so cell inclusion criteria are identical across the merged object.

**J. Composition differences between atlases.**
HLCA healthy lung is AT2/AT1/macrophage-rich; tumours are immune- and malignant-rich. Integration methods can mis-align non-overlapping populations; scANVI's label guidance and the over-correction checks (G) are the mitigations. Differential abundance is modelled compositionally (scCODA) rather than by comparing raw proportions.

**K. Validation of the merge.**
- Atlas mixing per harmonised cell type (should be high) vs across types (should be low).
- Marker sanity per compartment in both atlases.
- **Leave-one-study-out** DE: the sign of top tumour-vs-adjacent genes must be reproducible when any single study is removed.

> **Results below come from the first complete run** (finished 2026-11-10 21:13 PST; every number is from `results/reports/summary.md` or `results/tables/`). The post-run review in [§3.10](#310--post-run-review--open-issues-2026-11-11) lists what the numbers revealed that still needs fixing. ⚠️ marks a result that depends on one of those issues.

### 3.4 ✅ QC & doublets — per sample, both conditions — *resolves S2, S6*
v1: global 98th-percentile cutoff, control filtered with the tumour threshold (bug), doublets on tumour only. v2: MAD outliers **per sample** (Heumos 2023), Scrublet + scVI-SOLO **per sample on every sample**, cell removed only when both agree.

**Result:** `select` gave 249,402 cells from 160 donors; 26 samples with < 200 cells (`qc.min_cells_per_sample`; 3,052 cells) were dropped before per-sample QC, which removed 9 donors entirely. 246,350 cells in → **231,504 kept** across **215 samples** (**151 donors**, 20 studies). *(Corrected 2026-11-12: this line previously said 160 donors, the pre-QC cohort count; verified against the post-QC object and `logs/run_all/qc.log`, see docs/design/C_dashboard_hosting.md §4.)* 13,203 MAD outliers; 672 consensus doublets. The two doublet callers disagree strongly: SOLO flags 41,531 cells (17 %), Scrublet 2,082 (0.8 %), so the "both agree" rule effectively removes only what Scrublet finds (see §3.10 item 5).

### 3.5 ✅ Malignant-cell calling (infercnvpy) — *resolves S7*
v1: none. v2: inferCNV against immune/endothelial/stromal reference; threshold = 95th percentile of CNV score in healthy-donor epithelium; malignant = epithelial ∧ tumour sample ∧ above threshold; LUSC chr3q gain checked as a positive control.

**Result:** threshold 0.0129; **3,942 malignant cells**. The chr3q positive control works: mean chr3q score is 0.039 in malignant cells vs 0.004 in healthy and adjacent-normal epithelium.

| Check | Value | Reading |
|---|---|---|
| Precision vs author labels | 3,903 / 3,907 CNV calls with an author label are author-"malignant" (99.9 %) | calls are trustworthy |
| Sensitivity vs author labels | 3,903 / 29,172 author-"malignant" cells called (13 %) | **very conservative** |
| CNV-low tumour-sample epithelium | 27,577 cells, mean chr3q 0.025 (6× healthy) | many under-called malignant cells sit here |
| Tumour donors with ≥ 1 / ≥ 10 malignant calls | 22 / 8 of 40 | malignant-cell analyses rest on few donors |

⚠️ The 95th-percentile null is strict by design (false positives would contaminate the malignant signature), but at 13 % sensitivity the un-called malignant cells fall into the `Unknown` label and leak into other analyses (§3.10 item 2).

### 3.6 ✅ DE, pathways, abundance — *resolves S3, S4, S5, S12*
v1: cell-level Wilcoxon + phenotype-permutation GSEA with cells as samples (pseudoreplication). v2: pseudobulk per sample × cell type → PyDESeq2; decoupler 2.x PROGENy / CollecTRI / Hallmark on the Wald statistic; scCODA + Milo for abundance.

**Result (padj < 0.05):**

| Contrast | Design | Cell-type groups | Significant gene × group rows |
|---|---|---|---|
| `paired_tumour_vs_adjacent` (primary) | `~ donor + condition` | 15 | 9,665 |
| `tumour_vs_healthy` (sensitivity, confounded by design) | `~ study + assay + condition` | 20 | 20,174 |
| `malignant_vs_healthy_basal` | unpaired | 1 | 300 |

In the paired contrast the strongest signal is in myeloid and lymphoid cells: Macrophage 2,038 genes (14 donor pairs), Neutrophil 1,652 (5), CD4 T 924 (11), Monocyte classical 753 (12), CD8 T 706 (12). No epithelial group is tested at all, because none has ≥ 10 cells in ≥ 3 samples per condition (`min_cells_per_pseudobulk`, `min_samples_per_condition`). This matches the earlier power note that adjacent-normal tissue holds few epithelial cells. The epithelial question is answered by `malignant_vs_healthy_basal` instead. ⚠️ `Unknown` (2,181 paired genes) is a mixed bucket and should not be interpreted (§3.10 item 2). *(Run 1.)* **Run 2 (§3.11):** paired 8,474 rows; only `Unknown` changes (2,181 → 990, now mostly unannotated HLCA cells), every other group is identical; `tumour_vs_healthy` 19,452; new `cnvlow_tumour_vs_healthy_basal` 1,040.

**scCODA** (reference cell type chosen automatically = `Unknown`; credible = inclusion probability ≥ 0.8): relative to adjacent-normal, tumour samples have lower CD4 T, CD8 T, NK, Macrophage, Monocyte non-classical, AT1 and AT2 proportions, and no credible change in Malignant. Healthy donors have fewer CD4 T / CD8 T / NK / B / classical monocytes and more Multiciliated cells than adjacent-normal tissue. ⚠️ With `Unknown` (mostly un-called tumour cells) as the reference, these are shifts *relative to the tumour-cell bucket*. That explains why nearly everything goes down in tumour, so the direction must be read with care (§3.10 item 2). *(Run 1.)* **Run 2 (reference fixed to cDC, §3.11):** tumour vs adjacent, credible: NK −1.15, CD8 T −0.80, CD4 T −0.66, Monocyte non-classical −0.43, AT2 −0.37, Macrophage −0.33, AT1 −0.29, EC lymphatic −0.12 (log2 fold change), and `Tumour epithelial (CNV-low)` +2.23. The run-1 directions hold with a stable reference, so they were not an artefact of the `Unknown` reference; the effect sizes are slightly larger.

**Milo:** 17,511 neighbourhoods, **0 significant** (SpatialFDR < 0.1). ⚠️ This is a bug, not a null result: the wrong contrast was tested (§3.10 item 1). *(Run 1.)* **Run 2 (paired, `~donor_id + condition`, 14 donors / 48 samples):** **1,671 of 6,225** neighbourhoods significant (917 enriched in tumour, 754 depleted). By annotation: CD4 T 280 up / 59 down, B 261 / 0, CD8 T 155 / 154, Macrophage 56 / 221, Monocyte classical 45 / 102, NK 26 / 112, CNV-low tumour epithelium 18 / 0. Milo finds what scCODA cannot: *states* within a type move in opposite directions (CD8 T, CD4 T), while NK, macrophages and classical monocytes are mostly depleted, consistent with scCODA.

### 3.7 ✅ Trajectory — *resolves S8*
v1: PAGA graphs, stub notebook. v2: Palantir pseudotime rooted in healthy basal cells + CellRank 2 (CytoTRACE + pseudotime kernels). Reported as transcriptional distance from normal basal (malignant cells are patient-clonal; not lineage tracing).

**Result:** 10,093 cells (6,151 Basal, 3,942 Malignant); root = one HLCA healthy basal cell; CellRank terminal states: Malignant, Basal_1, Basal_2. Median pseudotime is 0.55 for Basal and 0.67 for Malignant, with a narrower spread for Malignant (SD 0.17 vs 0.28). Malignant cells sit further from healthy basal, as expected, but the distributions overlap widely.

### 3.8 ✅ Spatial — *resolves S9, S10*
v1: Scanorama + cosine label transfer (one label per multi-cell spot). v2: cell2location deconvolution (separate env), squidpy SVGs / neighbourhood enrichment, CellCharter niches, SpatialData export. n = 1 section → descriptive only.

**Result:** 6,195 → **6,191 spots** kept (18,061 genes; median 66k UMIs, 10k genes per spot). cell2location ran for 10.6 h on CPU and estimated abundances for **27 harmonised cell types** per spot. CellCharter found **3 niches**: an epithelial/tumour niche (0), a mixed immune–stromal niche (1) and a fibroblast + smooth-muscle niche (2: Fibroblast 22 %, SMC/pericyte 41 %). Spatial bivariate LR is dominated by ECM → CD44/LRP1 pairs (FN1, collagens, CCN1/2, SERPING1; Moran's R 0.74–0.84). Moran's I is significant for 1,998 of the 2,000 genes tested, which is typical for Visium; genes are ranked by I rather than filtered by p.

⚠️ **Rare epithelial** (only 66 healthy HLCA cells: ionocytes, neuroendocrine, brush cells) gets 66 % of niche 0 and 41 % of niche 1, while Malignant gets 11 % and 5 %. This is almost certainly a deconvolution artefact, not biology (§3.10 item 3), so per-type spatial abundances should not be quoted until it is fixed. *(Run 1.)*

✅ **Run 2 (2026-11-12/13, §3.11):** reference labels with < 200 cells are dropped (Rare epithelial 66, Mesothelial 91) and `Unknown` (28,224) is excluded, so the reference has **25 signatures**, including `Tumour epithelial (CNV-low)` (decisions in `results/tables/c2l_reference_labels.csv`). cell2location took 11.2 h. CellCharter's stability-based auto-K now picks **9 niches** (3 in run 1): five epithelial/tumour niches (0, 2, 3, 6, 8: CNV-low tumour epithelium 40–62 % + Malignant 12–14 %, with Multiciliated 11–14 % and Plasma 6–9 %), a smooth-muscle/fibroblast niche (7: SMC/pericyte 42 %, Fibroblast 24 %, like run 1's niche 2), and three mixed stromal–myeloid niches (1, 4, 5: Fibroblast, Neutrophil, classical Monocyte). Rare epithelial no longer absorbs the tumour signal. Caveats: CNV-low vs Malignant within a spot is not separable with confidence (their signatures overlap); Mesothelial cells can no longer be mapped; Multiciliated and Plasma in tumour niches may be real (airway-like differentiation, plasma-cell infiltrate) or partial absorption — unverified on n = 1 section.

### 3.9 ✅ Cell–cell communication — *resolves S11*
v1: LIANA on pooled tumour cells only. v2: LIANA+ per sample in every condition, differential CCC from pseudobulk DE, Tensor-cell2cell programmes, spatial bivariate LR on Visium.

Issues found in the first complete `ccc` run (2026-11-10) and their fixes:

| Issue | Why it mattered | Fix |
|---|---|---|
| Sub-analyses wrapped in `try/except → log.warning` | the step "succeeded" with differential CCC and the tensor silently missing | every sub-analysis still runs, but any failure raises at the end of the step, so it gets no completion stamp (§2.6) |
| Differential CCC filter used a stale group name (`Malignant_vs_normal_epithelium`) | the DE table holds two overlapping cell-type contrasts (`paired_tumour_vs_adjacent`, `tumour_vs_healthy`), giving 27k duplicate gene × group rows, so LR scores would have mixed paired and unpaired designs | explicit `ccc.de_contrast: paired_tumour_vs_adjacent` (the primary paired design, §3.6); `differential()` raises on duplicate gene × group rows |
| Tensor built with `how="outer_cells"` over all 215 samples | that option intersects LR pairs across *all* samples. Only **1 of 3,159** LR pairs (`PKM^CD44`) is present in every sample, so all 8 factors were trivial. The step still exited 0 | `how="outer"`, which keeps cells and LR pairs present in ≥ `tensor_outer_fraction` (1/3) of samples, with missing entries set to 0. A guard raises on a degenerate tensor (< 50 LR pairs, < 3 cell types or < 10 samples) |
| Tensor scope: 133 of 215 samples are healthy donors, and an `Unknown` label sat in 79 % of samples | the factors would mostly capture variation between healthy donors, not tumour vs normal | `ccc.tensor_conditions: [tumour, tumour_region_unknown, adjacent_normal]` (82 samples) and `tensor_exclude_labels: [Unknown]`. Result: 82 samples × 868 LR pairs × 10 senders × 10 receivers (user decision) |
| Elbow rank search hard-coded to `upper_rank=12` | on the rescoped tensor each rank took longer than the last on CPU (rank 1: 11.5 min, rank 2: ~30 min; each rank = 10 random-init fits × ≤ 100 iterations, the cell2cell `regular` preset). Extrapolated, ranks 3–12 would take > 15 h, and the step had run 24 min before rescoping | new config `ccc.tensor_upper_rank: 8` (read in `ccc.py`). The rank is still chosen from the elbow of reconstruction error, so it stays data-driven, just over a shorter range. With 10 sender/receiver types, 8 factors is already about as many as can be interpreted one by one. Fits per rank left at 10: lowering them would mean calling `elbow_rank_selection` directly. **Rejected:** a fixed rank (subjective) and leaving 12 (≈ a day of compute for ranks that are unlikely to be picked). **Check after the run:** if the elbow lands at 8, the true optimum may be higher, so rerun with a higher ceiling on GPU or with fewer fits per rank (user decision, 2026-11-10) |

Known limitation: **Malignant cells are not in the tensor.** Only ~4 % of per-sample LIANA results have ≥ `min_cells` (10) malignant cells, which is far below the 1/3 sample-presence cut. Malignant-cell communication is covered by the per-sample LIANA results and the spatial bivariate analysis instead. Follow-up options: a tumour-only tensor with a lower `outer_fraction`, or a lower `min_cells` for Malignant.

**Result (run finished 2026-11-10 21:12, `ccc` 5.3 h):**
- **Elbow rank = 5**, below the ceiling of 8, so the shorter search did not cut off the optimum and no rerun is needed.
- **Differential CCC** (paired tumour vs adjacent): the top interactions are neutrophil-derived VEGFA/PLAU signalling (→ CD44 on NK, PLAUR on macrophages) and the TIGIT–NECTIN2 checkpoint axis between CD4 T cells and macrophages/monocytes. CEACAM1 → HAVCR2 (TIM-3) on CD8 T cells also ranks in the top 20.
- **Tensor factors** (82 samples × 868 LR pairs × 10 × 10): Factor 1 is macrophage/cDC-sender (antigen presentation: B2M–CD1A/B, GRN–TNFR); Factor 3 is fibroblast/EC-sender to cytotoxic T (LCK–CD8, GZMA–F2R); Factor 5 is fibroblast/monocyte-sender to macrophages (FN1, HGF, MMP1 → CD44). ⚠️ **No factor separates tumour from adjacent-normal samples** (Mann–Whitney on sample loadings, 58 tumour vs 24 adjacent: all p ≥ 0.31). The factors describe communication programmes that vary between samples, not tumour-specific ones. Most likely they track donor and study, since LuCA studies differ in immune enrichment. Reported as such; a study-adjusted test of the loadings is a follow-up item. *(Run 1.)*

**Run 2 (2026-11-13, §3.11):** `ccc` 7.0 h. The tensor now has **82 samples × 951 LR pairs × 11 × 11**: `Tumour epithelial (CNV-low)` is present in enough samples to enter, so tumour-epithelial communication is in the tensor for the first time (`Malignant` still is not). Elbow rank **6** (ceiling 8). The CNV-low group loads most on Factor 5 (sender 0.42, with macrophages 0.55 and cDC 0.40; top LR pairs LAMB3–CD44, CXCL8–CD79A, GNAI2–EGFR, S100A9–CD68). Factor 6 is a chemokine programme from CD8 T / NK / classical monocytes (CCL5/CCL4 → CCR1). Still **no factor separates tumour from adjacent** (52 vs 24 samples, all p ≥ 0.34), so §3.10 #4 stands. Differential CCC is unchanged except that `Unknown`-sender rows (e.g. NECTIN4 → TIGIT in run 1's top 5) drop out: those were tumour cells hiding in `Unknown`.

### 3.10 ✅ Post-run review — open issues (2026-11-11)

The pipeline ran end to end, but reading the outputs (not just the exit codes) showed the issues below. They are listed by how much they change conclusions. None of them is a crash, so all of them would have shipped silently without this review.

| # | Issue | Evidence | Effect | Proposed fix | Status |
|---|---|---|---|---|---|
| 1 | **Milo tested the wrong contrast** | `milo.da_nhoods(design="~condition")` with no `model_contrasts`. pertpy then tests the *last* level alphabetically, `tumour_region_unknown` (2 HLCA donors), against the reference `adjacent_normal` | the 0 significant neighbourhoods are meaningless, not a null result | `model_contrasts="conditiontumour - conditionadjacent_normal"` on tumour + adjacent samples, with donor in the design (paired, like DE). Rerun Milo only (minutes; scCODA unaffected) | ✅ resolved, run 2 (§3.11) |
| 2 | **`Unknown` holds ~25k un-called malignant cells** | 57,335 `Unknown` cells; 29,172 have the LuCA author label `malignant cell`, and only 3,903 of those pass the CNV threshold | `Unknown` is mostly tumour cells mixed with unmapped HLCA cells. It contaminates scCODA (chosen as the reference type), DE (7,947 genes) and CCC (top-20 rows). | split it: label CNV-low author-malignant cells as a separate group (e.g. `Tumour epithelial (CNV-low)`), keep `Unknown` for truly unmapped cells, and fix the scCODA reference to a stable type (e.g. Fibroblast or EC capillary). Reruns cnv → downstream | ✅ resolved, run 2 (§3.11) |
| 3 | **Rare epithelial absorbs tumour signal in cell2location** | 66-cell reference (ionocyte/PNEC/brush) → 66 % of the tumour niche | spatial cell-type abundances are wrong for epithelium | drop or merge reference labels below a minimum size (e.g. < 200 cells) before building signatures. Fixing #2 also gives the reference a proper tumour-epithelium signature. Reruns deconvolve (~10 h CPU) → spatial → ccc spatial part | ✅ resolved, run 2 (§3.8, §3.11) |
| 4 | Tensor factors do not separate conditions | see §3.9 | interpretation only | report as sample-level programmes; follow-up: study-adjusted test | documented; still true in run 2 (all p ≥ 0.34) |
| 5 | SOLO vs Scrublet disagreement | 41,531 vs 2,082 flagged; consensus 672 | doublet removal is effectively Scrublet-only; SOLO's rate looks miscalibrated on atlas data that the original authors already QC'd | check SOLO's score distribution per study. Consensus is the conservative choice, so no rerun for v2.0 | documented |
| 6 | Atlas mixing below random for every cell type | e.g. Fibroblast 0.125 vs 0.493 expected | partly expected (HLCA ≈ healthy, LuCA ≈ tumour, so the same label differs in state), but it limits cross-atlas claims | add a condition-matched mixing metric (adjacent-normal LuCA vs healthy HLCA within type) | documented |
| 7 | `report` provenance row is one run stale | the report reads `provenance/report.json` before writing its own | cosmetic | write own provenance first, or exclude self | ✅ fixed, run 2 (§3.11) |

### 3.11 ✅ Run 2 — fixes for §3.10 #1–#3 (2026-11-12/13)

Decisions taken unattended under defaults agreed with the user ("start on #1–#3"), each config-driven and
reversible. Open questions are collected in [docs/design/FIXES_RUN2.md](docs/design/FIXES_RUN2.md).

**#1 Milo contrast.** `de.milo_design()` keeps only `tumour` + `adjacent_normal` samples of donors who gave both
(14 donors, 48 samples post-QC), designs `~donor_id + condition` and passes the explicit
`model_contrasts="conditiontumour - conditionadjacent_normal"` (config `abundance.milo_contrast`). The kNN graph and
neighbourhoods are rebuilt on that subset (healthy and region-unknown cells would otherwise sit in neighbourhoods
that the test can't use). Checked against the installed pertpy 1.5.0: with a contrast, `da_nhoods` drops the
intercept (`+ 0`), so donor gets full dummy coding and condition one treatment-coded column; the contrast resolves
to that column (unit test `test_milo_contrast_matches_pertpy_coefficients` builds the design with formulaic and runs
pertpy's own `_contrast_vector`). Default solver = pydeseq2 (R-free); 48 samples vs 15 coefficients leaves 33
residual df. **Rejected:** a random donor intercept `(1 | donor_id)` (pertpy's mixed model; slower, and the fixed
donor term matches the primary pseudobulk DE exactly); keeping all conditions with `subset_samples` (neighbourhoods
built partly from cells that are never tested). Milo failures now fail the step (run 1 swallowed them with a
warning; same rule as `ccc`, §3.9).

**#2 `Unknown` split + scCODA reference.**
- `make_de_groups()` (cnv step) gives cells whose LuCA author label is `malignant cell`, that sit in a tumour sample
  (`cnv.tumour_conditions`) and fail the CNV threshold, the new group **`Tumour epithelial (CNV-low)`**
  (`cnv.cnv_low_label`). `Malignant` (CNV-high) is unchanged; `Unknown` keeps only cells without a usable label
  (mostly 25,951 unannotated HLCA cells). The **64** author-malignant cells in *adjacent-normal* samples stay
  `Unknown`: a tumour label in normal tissue would enter the paired DE and CCC as if it had a counterpart. The split
  lives in `de_group` (the analysis grouping), **not** in `harmonised_label`, which steered scANVI: changing it would
  force re-integration (2 h) for no gain, since the crosswalk deliberately maps the author label to `Unknown` so it
  never steers integration (§3.3 E). Expected size from run-1 objects: 25,205 cells.
- **scCODA reference = `cDC`** (`abundance.sccoda_reference`, was `automatic` → `Unknown`). Chosen from per-sample
  proportions in run-1 objects (215 samples, new grouping):

  | Candidate | present in samples (min over conditions) | CV of proportion | mean prop. adjacent / healthy / tumour |
  |---|---|---|---|
  | **cDC** | 87 % (81 %) | 1.21 | 0.024 / 0.019 / 0.021 |
  | Macrophage | 89 % (84 %) | 1.05 | 0.201 / 0.277 / 0.111 |
  | Fibroblast | 71 % (50 %) | 3.51 | 0.008 / 0.028 / 0.024 |
  | EC capillary | 56 % (46 %) | 3.25 | 0.053 / 0.019 / 0.014 |

  The brief suggested Fibroblast or EC capillary; the data rule them out: several LuCA studies are CD45-enriched, so
  stromal/endothelial cells are absent from 29–44 % of samples, and a reference must be present almost everywhere.
  Macrophage is the most ubiquitous but itself changes 2.5× between tumour and healthy (a reference is assumed
  constant). cDC is nearly as ubiquitous with a flat mean across conditions. No type reaches scCODA's own 95 %
  presence rule, which is why "automatic" fell back to the mixed `Unknown`.
- **Downstream treatment of the new group:**

  | Analysis | Treatment | Why |
  |---|---|---|
  | Paired / unpaired DE | not tested (`de.tumour_only_groups`, like `Malignant`) | no counterpart in control tissue |
  | New DE contrast `cnvlow_tumour_vs_healthy_basal` | CNV-low group vs healthy basal, unpaired, same code as `malignant_vs_healthy_basal` (`de.vs_healthy_basal`) | shows whether the CNV-low cells look like tumour epithelium (validates the label) |
  | scCODA, Milo | included as a cell type | it is a real, tumour-only population |
  | LIANA by-sample, tensor | included (tensor still excludes `Unknown`); enters the tensor only if present in ≥ 1/3 of tensor samples | tumour-cell communication was missing from the run-1 tensor |
  | Differential CCC | absent (no DE stats; not in the paired contrast) | same as `Malignant` |
  | Trajectory | **not** included (Basal + CNV-called Malignant only) | the endpoint must stay CNV-confirmed; adding unconfirmed cells would blur it (possible follow-up) |
  | cell2location reference | included | gives the reference a real tumour-epithelium signature (#3) |

**#3 cell2location reference.** `deconvolve.prepare_reference_labels()` runs before signature building
(config `spatial.c2l.min_ref_cells: 200`, `merge_small_labels: {}`, `exclude_labels: [Unknown]`) and writes
`results/tables/c2l_reference_labels.csv`. Counts are taken on the full reference before the per-label subsample.
Expected from run-1 objects: `Rare epithelial` (66 cells) and `Mesothelial` (91) dropped. No merge: the crosswalk
gives no obvious parent (ionocytes/PNEC/brush cells are not secretory or basal; mesothelium is not fibroblast).
`Unknown` (28k cells after the split) is excluded because a signature of a mixed bucket is exactly the kind of
catch-all that absorbed tumour signal. Rejected: lowering `max_ref_cells_per_label` or reweighting (doesn't fix a
66-cell signature), keeping `Rare epithelial` with a prior (no principled prior from one section).

**#7 report provenance** (trivial, fixed): `steps.report` writes a preliminary `provenance/report.json` before
building the report, so the provenance table shows the current run. `run_step` overwrites it with the final runtime.
The report also gained a Milo section (design, contrast, counts by cell type) and prints the scCODA reference.

**Verification before the rerun:** 9 unit tests (`tests/test_fixes_run2.py`); `scripts/smoke_test.py --run2`
(12k real cells from 56 whole samples incl. all 48 paired samples: de_group split → all DE contrasts incl. both
vs-healthy-basal → pathways → scCODA with reference cDC → Milo paired, 48 samples / 14 donors) passed in 132 s;
`pixi run -e c2l python scripts/smoke_test.py --c2l` (reference decisions → RegressionModel → Cell2location, tiny
epochs; `Unknown` absent from the output) passed in 202 s.

**Rerun:** run 1's small outputs archived to `results/_archive_run1_2026-11-10/` (tables, reports, figures, the
run-1 DuckDB, stamps, log; git-ignored). `run_all.py --only cnv de abundance trajectory deconvolve spatial ccc report
export --force --keep-going` in tmux `lusc-v2:pipeline`, started 2026-11-12 11:18 PST. `--only` + `--force` reruns
exactly the steps downstream of the change (upstream `qc`/`integrate`/`spatial_prep` keep their stamps);
`trajectory` is rerun although its definition is unchanged, so no output depends on a stale input.

**Result (2026-11-13 07:07 PST, ALL DONE, no failures, 19.8 h wall):** cnv 10 min · de 4 min · abundance 77 min
· trajectory 4 min · **deconvolve 11.2 h** · spatial 2 min · **ccc 7.0 h** · report 1 min · export 43 s.

| | Run 1 (2026-11-10) | Run 2 (2026-11-13) |
|---|---|---|
| Milo (SpatialFDR < 0.1) | 0 / 17,511, wrong contrast | **1,671 / 6,225** (917 up, 754 down in tumour), paired, 14 donors / 48 samples |
| scCODA reference | automatic → `Unknown` (mixed) | `cDC` (fixed) |
| scCODA tumour vs adjacent (credible log2FC) | NK −0.95, CD8 T −0.82, CD4 T −0.67, Mono NC −0.33, AT2 −0.30, Mac −0.22, AT1 −0.19 | NK −1.15, CD8 T −0.80, CD4 T −0.66, Mono NC −0.43, AT2 −0.37, Mac −0.33, AT1 −0.29, EC lymphatic −0.12; CNV-low tumour +2.23 |
| `Unknown` cells | 57,335 (incl. ~25k tumour cells) | 28,224 (unannotated HLCA + 64 adjacent author-malignant) |
| `Tumour epithelial (CNV-low)` | — (inside `Unknown`) | 25,205 cells |
| `Malignant` (CNV-called) | 3,942 | 3,942 (identical; threshold 0.0129 reproduced) |
| Paired DE significant rows | 9,665 (`Unknown` 2,181) | 8,474 (`Unknown` 990; all other groups identical) |
| Tumour epithelium vs healthy basal | Malignant 300 genes | Malignant 300; **CNV-low 1,040** (40 vs 22 samples) |
| cell2location reference | 27 labels incl. Rare epithelial (66 cells) | 25 labels (Rare epithelial, Mesothelial dropped; `Unknown` excluded) |
| Tumour-niche composition | Rare epithelial 66 %, Malignant 11 % | CNV-low tumour 40–62 %, Malignant 12–14 % |
| CellCharter niches (auto-K) | 3 | 9 (5 epithelial/tumour, 1 SMC–fibroblast, 3 stromal–myeloid) |
| Tensor | 82 × 868 × 10 × 10, rank 5 | 82 × 951 × 11 × 11 (+ CNV-low), rank 6 |
| Tensor factor vs condition | all p ≥ 0.31 | all p ≥ 0.34 (unchanged, §3.10 #4) |
| Differential CCC top 10 | incl. `Unknown` → CD4 T NECTIN4–TIGIT | same rows minus the `Unknown` one |

**Does the CNV-low label hold up?** Against healthy basal cells it shows the LUSC programme: **SOX2** (3q26, the
defining LUSC amplicon) log2FC +2.24, **TP63** +1.40, **CDKN2A** +5.56, and DLX5 is the top non-technical gene in
both the CNV-low and Malignant contrasts; 58 of Malignant's 118 up-genes are also up in CNV-low. These are tumour
cells whose CNV signal is below a deliberately conservative threshold, not mislabelled normal epithelium. Caveat:
the top-ranked CNV-low genes are mitochondrial pseudogenes (MTCO3P18, MTND4P33, …), a study/chemistry artefact
that the `~assay + epi_state` design does not fully absorb; read the list from the marker genes, not the top rows.

**Found by the post-run checks:** CellCharter picked 9 niches, one more than the dashboard's 8-colour palette, so the
Spatial page crashed (`IndexError`, caught by `tests/test_dashboard.py`). Fix: `ui.series_color(i)` (8 validated
colours + 4 extra, then cycle) and a regression test.

**Checks after the run:** `pixi run lint` and `ruff format --check` clean; default tests 24 passed / 2 skipped;
dashboard + agent tests 45 passed; `agent-eval --check` 24 items, 0 grader failures (truths recomputed on the run-2
DB, no LLM); `scripts/cloud_dry_run.sh` passed on the run-2 DB (15/15 pages, DB via secrets + `hf_hub_download`,
peak RSS 440 MB = 16 % of 2.7 GB, wrong sha256 refused, health OK; `.tooling/cloud_dry_run_4.log`, after re-exporting with the updated known-issues text). App restarted
in `lusc-v2:app` on the new DB.

Still open (documented, not fixed here): #4 tensor not condition-specific, #5 SOLO/Scrublet, #6 atlas mixing;
Mesothelial can't be mapped spatially; trajectory still covers Basal + CNV-called Malignant only.

---

## 4. Looking ahead: AI-native layer

🔍 **Concepts prototyped 2026-11-12 in [docs/design/](docs/design/README.md)** (agentic "Ask the atlas" with read-only DuckDB tools, provenance and an eval set; Nemotron assessment; dashboard explainers + hosting bundle). Original note: ⏳ To be scoped after v2.0.0 ships. Direction: each pipeline step is already a pure function with declared inputs/outputs (`lusc_atlas.steps`) — exactly the shape an agent tool or **Claude Skill** needs. Candidates: a `lusc-atlas` Skill that runs/inspects steps and explains results; an MCP server over the processed AnnData/SpatialData objects; agent-assisted annotation review (cf. CellAgent, scAgent, Biomni — all early, human review required).

---

## 5. Project claims → evidence → v2 upgrade

The project's public summary (CV / portfolio) makes specific claims. This section keeps every claim
**traceable to code and output**, records where v1's evidence was weaker than the wording, and shows how
v2 turns each one into something that survives expert questioning. Checked against tag `v1.0.0` on 2026-11-09.

| Claim (v1 wording) | v1 evidence | Verdict | v2 upgrade (where) |
|---|---|---|---|
| Integrated squamous lung cancer profiles with **HLCA** reference data | LUSC cells *came from* HLCA (3 donors, 1 study); no reference mapping | ⚠️ thin | Two-atlas integration **HLCA + LuCA** (45 + 3 LUSC donors), deduplicated, label crosswalk, benchmarked (§3.2–3.3; `select.py`, `integrate.py`) |
| **Scanpy/AnnData** | all notebooks | ✅ | packaged pipeline on scanpy 1.12 / anndata 0.13 |
| **scVI** | scVI used only as SOLO backbone | ✅ (narrow) | scVI **and** scANVI **and** sysVI as integration candidates, chosen by scib-metrics + over-correction guard (`integrate.py`) |
| **Seurat** | only `flavor='seurat_v3'` HVG selection inside Scanpy; no Seurat (R) code | ⚠️ overstated | not added (Python-only stack). Accurate phrasing: "Seurat-v3 HVG selection" — or drop |
| QC + **VAE/SOLO doublet removal** | SOLO on tumour cells only (126 doublets / 20,631); output file never read downstream (R1) | ⚠️ bug | per-sample MAD QC + **Scrublet ∧ scVI-SOLO on every sample of both atlases** (`qc.py`) |
| Cluster and annotate **18,381** single-cell profiles | `03_analysis` 18,381 cells after QC | ✅ number verified | annotation = harmonised HLCA/LuCA labels via Cell Ontology crosswalk + scANVI; cell count regenerated by the pipeline |
| Integrated with **3,858 Visium spots** using **Scanorama** | v1 matrix 3,858 × 18,085 (an earlier 10x release); Scanorama + cosine label transfer | ✅ number / ⚠️ method | official 10x 2.0.1 release (6,195 spots); **cell2location** deconvolution (`deconvolve.py`) — Scanorama retired: a 55 µm spot holds several cells |
| **FDR-adjusted DGE** | `pvals_adj` from cell-level Wilcoxon/t-test | ⚠️ pseudoreplication | **pseudobulk PyDESeq2**, paired tumour vs adjacent-normal within donor (`de.py`) |
| **GSEA** | gseapy phenotype permutation with cells as samples | ⚠️ pseudoreplication | decoupler 2 **Hallmark GSEA + PROGENy + CollecTRI** on DESeq2 statistics (`de.py`) |
| **PAGA trajectory** | PAGA graphs; `04_trajectory` never ran | ⚠️ partial | **Palantir pseudotime + CellRank 2** (CytoTRACE + pseudotime kernels), rooted in healthy basal cells (`trajectory.py`) |
| **LIANA** ligand–receptor | rank_aggregate on pooled tumour cells | ✅ (narrow) | **LIANA+ per sample, all conditions**, differential CCC, Tensor-cell2cell, spatial bivariate LR (`ccc.py`) |
| Disease-associated **immune, epithelial, stromal** states | compartment-level UMAPs/markers | ✅ | + CNV-defined **malignant** state, scCODA/Milo abundance shifts, CellCharter spatial **niches** |
| **Mapped 36 cell types into tissue space** | per-spot cosine scores for 27 type columns (`05_integration`) | ⚠️ count/method | per-spot **cell-type abundances** for 25 harmonised types incl. Malignant and CNV-low tumour epithelium (`results/tables/c2l_abundance.csv`; run 2 fixed the §3.10 #3 artefact) |
| Led a **five-member team** | 5 authors, git history (2 committers) | ✅ | unchanged; authors in `CITATION.cff` |

**Suggested v2 wording.** Numbers are copied from `results/reports/summary.md` / `results/tables/`
(run 2, 2026-11-13), not typed from memory:

> *Spatial transcriptomics of squamous lung cancer (v2, reproducible pipeline).* Led a five-member team;
> rebuilt the analysis as a tested, Snakemake/pixi-pinned Python package. Merged two curated atlases (HLCA +
> LuCA; **231,504** cells after QC, **151** donors, **20** studies) with explicit deduplication, gene-space and
> Cell-Ontology label harmonisation; benchmarked seven embeddings (incl. scVI/scANVI/sysVI/Harmony, scArches and
> a Geneformer foundation-model baseline) with scib-metrics and an over-correction guard (scANVI selected).
> Called **3,942** malignant cells by inferred CNV (chr3q-gain validated), ran donor-paired pseudobulk DESeq2
> with PROGENy/CollecTRI/Hallmark activity, compositional (scCODA) and donor-paired neighbourhood (Milo:
> **1,671** of 6,225 neighbourhoods shifted in tumour) abundance, Palantir/CellRank 2 trajectories, and LIANA+
> differential, tensor (Tensor-cell2cell) and spatial cell–cell communication; deconvolved **6,191** Visium spots
> with cell2location into **25** cell types and **9** CellCharter niches.

Run 2 resolved §3.10 #1–#3, so Milo, the scCODA directions and the per-type spatial abundances can now be quoted.
Still not to be claimed: tumour-specific communication programmes from the tensor (#4: no factor separates
tumour from adjacent).

---

## Timeline

| Date | Change |
|---|---|
| 2026-11-09 | v1 audit; tag `v1.0.0`; branch `v2-dev`; tmux session `lusc-v2`; legacy move; this log created |
| 2026-11-09 | Literature review of 2024–26 methods (stage-by-stage recommendations → §3) |
| 2026-11-09 | Root disk filled by conda build → all tooling moved to `/storage`; pixi envs `default`/`c2l`/`workflow` built and import-verified |
| 2026-11-09 | Visium switched to official 10x 2.0.1 release (v1 file differs) |
| 2026-11-09 | HLCA found to hold only 3 LUSC donors → decision: HLCA + LuCA; cross-atlas merge design (§3.3) |
| 2026-11-09 | Project claims audited against v1 outputs; claim → evidence → v2 traceability (§5) |
| 2026-11-09 | `select` run: 249,402 cells × 17,457 genes; 20 studies kept, 0 cross-atlas duplicates after `is_primary_data`; 15/40 LUSC donors paired; Smart-seq2 excluded (non-UMI); exclusions in `merge_exclusions.csv` |
| 2026-11-09 | Added integration candidates: **scArches transfer learning** onto the published HLCA scANVI model and **zero-shot Geneformer-V2-104M_CLcancer** (benchmark subsample, 20k cells). External pickles loaded only via `safe_load.py` (sha256-pinned, whitelisted types after static `pickletools` scan) — user-approved |
| 2026-11-09 | Smoke test (`scripts/smoke_test.py`, 16k real cells, tiny epochs) found and fixed 5 bugs before the full run: scrublet `layer` kwarg; per-sample seurat_v3 LOESS singular on small samples → HVGs batch-aware per **study**; categorical `fillna`; CytoTRACE needs an imputed layer (kNN smoothing); smoke sampler realism. Passed end-to-end in 323 s |
| 2026-11-09 | Observation: adjacent-normal tissue has only 1,290 epithelial cells (several LuCA studies are immune-enriched) → paired tumour-vs-adjacent DE is well powered for immune/stromal types, thin for epithelium; malignant-vs-healthy-basal contrast covers the epithelial question |
| 2026-11-09 | Full run launched unattended in tmux `lusc-v2:pipeline` (`run_all.py --keep-going`); roadmap started (`docs/ROADMAP.md`) |
| 2026-11-10 | First full run: ✓ qc 2.5 h · fm_embed 2.0 h · integrate 2.0 h · cnv 10 min · de 4 min · trajectory 4 min · **deconvolve 10.6 h** · spatial 2.4 min; ✗ abundance (scCODA needs `pertpy[jax]`, §2.3). Stopped during ccc to fix the env |
| 2026-11-10 | Resume hardening (§2.6): atomic writes + completion stamps + `--adopt`; 12 completed steps adopted; resumed run covers abundance → ccc → report |
| 2026-11-10 | Resumed run: ✓ abundance 2.3 h (scCODA + Milo); ✗ ccc (LIANA 1.x API on 2.1.1) · ✗ report (`tabulate`). Ported `ccc.py` to LIANA 2.x, added `tabulate` (§2.3); subsample-tested; reran → ✓ ccc 24 min · ✓ report 1 min (**ALL DONE**) |
| 2026-11-10 | Output review found a degenerate tensor (1 LR pair) and silent failure / contrast-mixing risks in `ccc` → fixes in §3.9; tensor rescoped to tumour + adjacent samples (82 × 868 × 10 × 10); `ccc` + `report` rerun |
| 2026-11-10 | Rescoped tensor's rank search (`upper_rank=12`) projected > 15 h (rank 1 11.5 min, rank 2 ~30 min) → stopped; ceiling made configurable and set to 8 (§3.9); `ccc` + `report` relaunched 15:52 |
| 2026-11-10 | ✓ ccc 5.3 h (elbow rank 5, below the ceiling of 8) · ✓ report 65 s → **full pipeline complete** (21:13 PST); all 15 steps stamped |
| 2026-11-11 | Post-run output review → results written into §3.4–3.9; 7 open issues in §3.10 (Milo contrast bug; `Unknown` holds ~25k un-called malignant cells; Rare-epithelial cell2location artefact; tensor factors not condition-specific; SOLO/Scrublet disagreement; atlas mixing; stale report provenance); §5 wording filled; README results section |
| 2026-11-11 | Dashboard built on branch `v2.1-dashboard` (ROADMAP §6): new `export` step → DuckDB serving snapshot (118 MB, derived, never source of truth) in a separate pixi env `dashboard`; Streamlit app with one page per stage (question → decision → evidence → known-issue banners); page-render tests. Building it surfaced a count discrepancy: the post-QC object holds **151** donors, not the 160 quoted in §3.4 (160 = `cohort_donors`, pre-QC; 8 small healthy HLCA donors + 1 LuCA tumour donor have no cells in `sc_annotated.h5ad`) — to verify and correct |
| 2026-11-12 | **Zero-cost decision:** project must stay free → "Ask the atlas" chatbot page removed from the dashboard and the deploy bundle (needs a paid LLM key; page parked outside the repo); anthropic/openai dropped from app/requirements.txt; hosting = Streamlit Community Cloud + HF Dataset. Agent library, eval and MCP server kept, unwired (docs/design/README.md decision note). |
| 2026-11-12 | AI-layer concepts (docs/design/README.md), additive and uncommitted: **A** agentic layer `src/lusc_atlas/agent/` (11 read-only DuckDB tools incl. 8 curated, caveat-aware; Claude default + OpenAI-compatible Nemotron provider; tool-use loop with SQL provenance and a post-hoc number-grounding check; MCP server; "Ask the atlas" page; 24-question eval with SQL-computed truths, not yet run live: no keys); **B** Nemotron assessment (open-weights fallback, no GPU here, hosted API prototyping-only); **C** DuckDB enforced as the single serving source of truth (schema contract + runtime-read tests, schema 2: page_guide, glossary, step_stamps, config/content hashes), Start-here page + computed reader guide on every stage page + glossary, HF Docker-Space bundle + `deploy_space.sh` (not run) + local dry deploy (16/16 pages, 435 MB). Finding: HF Docker/Gradio Spaces now need a paid plan; Streamlit SDK deprecated 2025-04-30. §3.4 donor count corrected 160 → 151 post-QC (wording only) |
| 2026-11-12 | **Free-only hosting implemented** (§1.3): HF Docker-Space bundle parked outside the repo; `scripts/publish_dataset.sh` (DuckDB → public HF Dataset, never a Space; user-run) and `scripts/cloud_dry_run.sh` (clean clone, Python 3.12 venv from `app/requirements.txt`, secrets, DB via `hf_hub_download`, all pages, wrong-sha negative control) — passed on the run-1 DB (15/15 pages, 444 MB); `.streamlit/config.toml` moved to repo root; DB location read from env or `st.secrets` |
| 2026-11-12 | §3.10 #1–#3 + #7 fixed in code (§3.11): paired Milo contrast; `Tumour epithelial (CNV-low)` split out of `Unknown` + new vs-healthy-basal contrast; scCODA reference `cDC`; cell2location reference ≥ 200 cells, `Unknown` excluded; report provenance current. 9 unit tests + `smoke_test.py --run2` / `--c2l` passed. Run-1 outputs archived; run 2 launched 11:18 PST |
| 2026-11-13 | **Run 2 ALL DONE** 07:07 PST (19.8 h; deconvolve 11.2 h, ccc 7.0 h): Milo 1,671 / 6,225 significant (run 1: 0); scCODA directions hold with cDC reference; CNV-low group 25,205 cells (SOX2/TP63/CDKN2A up vs healthy basal); tumour niches now CNV-low tumour 40–62 % + Malignant 12–14 % (run 1: Rare epithelial 66 %); 9 niches; tensor 11 × 11 incl. CNV-low, rank 6, still not condition-specific. Dashboard palette crash on 9 niches fixed; all checks + cloud dry run (440 MB) pass. Write-up: `docs/design/FIXES_RUN2.md` |
