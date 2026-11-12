# C · A dashboard that explains the data, one source of truth, a public link

*Concept C of [docs/v3](README.md). Written 2026-11-12. Hosting facts verified that day on the official docs
(URLs inline).*

## Verdict

* **Single source of truth: done, and now enforced.** At runtime the app and the agent read only the DuckDB file.
  A test fails if any app file reads CSV/h5ad/YAML/JSON, and a schema-contract test fails if the DB lacks a table or
  column the app or agent relies on. A Docker-free "dry deploy" proves the bundle runs standalone.
* **Explain and interpret: done.** There is a new *Start here* page (workflow map + one computed headline per stage),
  a reader guide on every stage page (*What you're looking at · How to read it · What it means · Caveats*), a
  glossary, and *Ask the atlas*. Every interpretation sentence is computed from the DB when the page renders, so none
  can go stale after the §3.10 reruns.
* **Public access: ready to publish, with one finding that changes the earlier hosting decision.** Hugging Face now
  requires a **paid plan (PRO)** to create Docker or Gradio Spaces, and it deprecated the Streamlit SDK on
  2025-04-30. The bundle therefore targets a **Docker Space** and also runs unchanged on **Streamlit Community
  Cloud** (free). **Recommendation:** publish the DuckDB as a free HF Dataset either way. Host the app on an HF
  Docker Space if you have (or want) PRO, otherwise on Streamlit Community Cloud. Open question for you (README).

---

## 1. One source of truth (what changed and why)

**Wording, reconciled with the 2026-11-10 decision.** There are two layers, each with one job:

| Layer | What | Who reads it | Truth for |
|---|---|---|---|
| Computational | `data/processed/*.h5ad`, `*.zarr`, `results/tables/*.csv` | pipeline steps | the analysis (reproducible from code + config) |
| Serving | `results/dashboard/lusc_atlas.duckdb` | the dashboard, the Ask-the-atlas agent, the MCP server | everything a reader sees or asks |

The 2026-11-10 decision ("DuckDB is a derived snapshot, not the source of truth; h5ad/zarr stay for matrices")
still holds for the *computational* layer. What is new is that the DuckDB is the **only** serving source: nothing at
runtime reads anything else, and it is rebuilt from scratch by `export`, never edited.

**Audit result.** Before this session the app already read only DuckDB (curated content was loaded into the DB by
`export`, not read from YAML at runtime). What was missing was *enforcement* and *provenance*:

| Added | Where | Why |
|---|---|---|
| `test_app_reads_only_duckdb` | `tests/test_dashboard.py` | Fails if any `app/` file calls `read_csv`/`read_h5ad`/`yaml`/`open(`/`json.load`/`h5py`/`anndata` (two reviewed exceptions carry an explicit marker: decoding PNG bytes from the DB, and hashing the DB file) |
| `SCHEMA_CONTRACT` + `schema_contract` table + `test_schema_contract` | `src/lusc_atlas/export.py`, tests | Lists the tables/columns the app and agent depend on; a pipeline change that renames one now fails loudly instead of breaking a page silently |
| `meta`: `config_sha256`, `content_sha256`, `role`, schema v2 | export | Which config and which curated text produced this snapshot |
| `step_stamps` table | export (from `logs/run_all/.done/*`) | Which pipeline steps had finished cleanly when the snapshot was built |
| `page_guide`, `glossary` tables | export (from `docs/dashboard_content.yaml`) | Reader-guide text lives in the DB like everything else |
| `LUSC_DB_SHA256` pin | `app/lib/db.py` | A deployed app refuses to serve a snapshot whose sha256 differs from the one it was built and tested with |

## 2. Explain and interpret

| Page / element | What it does |
|---|---|
| **Start here** (new default page) | What the site is; how to use it in three steps; the workflow as a diagram (data → cohort → QC → integration → malignant → DE / composition / trajectory → spatial → CCC → DuckDB); one card per stage with its one-line purpose, a page link and its first computed finding; how to query the DB yourself |
| **Reader guide on every stage page** | *What you're looking at* and *How to read it* (curated, plain language); **What it means** (2–3 sentences computed from the DB at render time, e.g. "Only 7 of 38 tumour donors have ≥ 10 malignant cells"); **Caveats** (that stage's known issues with severity and status, linked to the Known-issues page). The existing question → decision → evidence → caveat banners stay |
| **Glossary** | 26 plain-language terms (pseudobulk, padj, scCODA, Moran's I, grounding check…), searchable |
| **Ask the atlas** | Concept A's agent with full provenance; works without a key as a tool runner |
| **Provenance** | Now also shows the config/content hashes, the step stamps and the schema contract |

**Honesty about §3.10.** The computed interpretations name the problem where it applies: `Unknown` heads the paired
DE table but holds un-called tumour cells; scCODA's "decreases" are relative to `Unknown`; Milo's 0 hits are a
wrong-contrast artefact, not a null result; the 66 % "Rare epithelial" niche is a deconvolution artefact. The
dashboard does not fix them (the user deferred §3.10); it makes sure nobody can read a number without the caveat.

**Charts.** No new chart types. The existing palette is the dataviz skill's validated reference palette, unchanged.
The palette validator could not be re-run because Node.js is not installed on this host. Condition colours stay
fixed across pages; cell types use highlight-one-grey-the-rest rather than 27 hues.

## 3. Public access

### 3.1 Options compared (facts retrieved 2026-11-12)

| | **HF Docker Space** + HF Dataset | **Streamlit Community Cloud** + HF Dataset | **MotherDuck** (hosted DuckDB) |
|---|---|---|---|
| What it is | Container on HF; app on port 8501 via `sdk: docker`, `app_port: 8501` | Streamlit's free hosting, deployed from a GitHub repo | Cloud DuckDB service; a DB can be shared read-only (`ATTACH 'md:_share/<db>/<token>'`) |
| Cost | **Needs a paid HF plan to create** Docker/Gradio Spaces (PRO for personal accounts; https://huggingface.co/docs/hub/spaces-overview). Datasets are free | Free | Free "Lite" tier: 10 GB storage, 10 compute-unit hours / month (fees addendum, Feb 2026; sources disagree on user caps) |
| Resources | CPU Basic: 2 vCPU, 16 GB RAM, 50 GB ephemeral disk; sleeps when idle | ~0.08–2 cores, 690 MB–2.7 GB RAM, ≤ 50 GB storage (docs say approximate and subject to change); sleeps after 12 h without traffic (https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app) | n/a (database only, still needs an app host) |
| Fits this app? | Yes: dry run peaked at **435 MB** RSS rendering every page | Yes: 435 MB is below the 2.7 GB ceiling | Would replace the file download with a network DB, but the app still needs hosting. Adds an account and a token, and the serving DB would no longer be a single checksummed file |
| Secrets (LLM keys) | Space *Secrets* → env vars | `st.secrets` / app settings | — |
| Streamlit support | Streamlit SDK **deprecated 2025-04-30**; Docker template instead (https://huggingface.co/docs/hub/en/spaces-changelog) | native | — |
| Verdict | **Best if you have PRO**: one HF account for app + data, Docker = reproducible runtime | **Best free option** | **Not recommended now**: a 119 MB read-only file is exactly what a file download + sha256 does well. Revisit if the DB grows past a few GB or needs live updates |

**Data: Dataset + sha256 pin, not bundled.** The 119 MB DuckDB is published as a public **HF Dataset** and downloaded
by the app at startup, pinned to a dataset revision and verified by sha256. Why: (1) the dataset *is* the "download
and query it yourself" link; (2) data can be re-published without rebuilding the app and the other way round; (3) the
pin keeps the guarantee that the site serves exactly the tested snapshot. `--bundle-db` / `BUNDLE_DB=1` puts it
inside the Space instead (simplest single-repo option).

### 3.2 What is ready (free-only, revised 2026-11-12)

The HF Docker Space path (`build_space.py`, `deploy_space.sh`, `space_dry_run.sh`) was **removed** after the
zero-cost decision: creating a Docker Space needs a paid HF plan. What remains uses free tiers only.

| File | Purpose |
|---|---|
| `scripts/publish_dataset.sh` (`pixi run -e dashboard publish-dataset`) | **You run it, with `HF_TOKEN`.** Creates or updates the **public** HF Dataset (card + DuckDB) and prints the Streamlit secrets to paste. It never creates a Space or sets hardware, refuses files over 1 GB, asks before uploading and never prints the token |
| `scripts/cloud_dry_run.sh` (`pixi run -e dashboard cloud-dry-run`) | Local rehearsal of Community Cloud: copies what git would push to a clean folder, makes a fresh Python 3.12 venv from `app/requirements.txt` only, writes a secrets file, fetches the DB through the same `hf_hub_download` call (offline cache), renders every page and checks the health endpoint |
| `docs/v3/dataset_card.md` | README for the HF Dataset |
| `app/lib/db.py` | Reads `LUSC_DB_HF_REPO` / `LUSC_DB_HF_REVISION` / `LUSC_DB_SHA256` from env vars or Streamlit secrets; refuses a file whose sha256 doesn't match |

### Cost: $0

| Service | Used for | Free tier | What would cost money (avoid) |
|---|---|---|---|
| GitHub (public repo) | code, CI | free for public repos | private repo with heavy Actions use |
| Hugging Face Dataset (public) | the DuckDB file (~119 MB) | public storage is free | Spaces with Docker/Gradio SDK or upgraded hardware, PRO subscription, private storage beyond the quota |
| Streamlit Community Cloud | the app | free; apps sleep after inactivity and wake on the next visit | nothing: there is no paid tier to upgrade into |
| LLM APIs | nothing (chatbot removed) | n/a | any API key on the public app |

### 3.3 How anyone can try it

1. **Once deployed:** open the public `*.streamlit.app` link. No login is needed to browse. If the app was
   asleep, click *wake up* and wait about 30 seconds.
2. **Locally (full repo):**
   ```bash
   git clone https://github.com/archana-b24/sc-spatial-analysis && cd sc-spatial-analysis
   pixi install -e dashboard
   # either build the DB from pipeline outputs …
   pixi run -e dashboard export
   # … or download the published snapshot and point the app at it
   export LUSC_DB=/path/to/lusc_atlas.duckdb
   pixi run -e dashboard app                       # http://localhost:8501
   ```
3. **Just the data:** download `lusc_atlas.duckdb` from the HF Dataset, then
   ```bash
   duckdb -readonly lusc_atlas.duckdb
   SELECT * FROM tables_info;                       -- what every table is and where it came from
   SELECT * FROM known_issues;                      -- read these before quoting numbers
   ```
   or use it from Python (`duckdb.connect("lusc_atlas.duckdb", read_only=True)`), or from an AI assistant through the
   MCP server (`pixi run -e dashboard mcp`; setup in `src/lusc_atlas/agent/mcp_server.py`).

### 3.4 Steps for you to deploy (free)

```bash
export HF_TOKEN=hf_...                       # write token, free account
pixi run -e dashboard publish-dataset        # confirms, uploads, prints three secrets
```
Then on share.streamlit.io: *Create app* → *Deploy a public app from GitHub* → repo
`archana-b24/sc-spatial-analysis`, branch `v2-dev`, main file `app/streamlit_app.py`; *Advanced settings*:
Python 3.12 and paste the three printed lines (`LUSC_DB_HF_REPO`, `LUSC_DB_HF_REVISION`, `LUSC_DB_SHA256`)
into *Secrets*. Community Cloud installs `app/requirements.txt` from the entry-point folder.

## 4. Pending verification resolved: 151 vs 160 donors

**Both numbers are right; they describe different stages.** Verified 2026-11-12 against the DuckDB and
`logs/run_all/qc.log`:

* `select` produced the cohort: **160 donors, 249,402 cells, 20 studies** (`cohort_donors`).
* QC first drops samples with fewer than `qc.min_cells_per_sample: 200` cells; the log says "dropping 26 samples with
  < 200 cells" (3,052 cells). Nine donors had *only* such samples (8 small healthy HLCA donors from Banovich_Kropski,
  Gomperts, Schiller, Kaminski, Teichmann_Meyer, Schultze_unpubl; and LuCA Zilionis_Klein_2019 patient_2).
* So per-sample QC saw 246,350 cells in 215 samples, and the post-QC atlas holds **231,504 cells from 151 donors**
  (114 HLCA + 37 LuCA) in 215 samples and 20 studies.

DEVELOPMENT_LOG §3.4 said "231,504 kept across 215 samples (160 donors…)", which mixed the two stages. It and the §5
resume wording are corrected (wording only, no science changed). The cohort page, the agent's `cohort_summary`
note and eval item q02 (which fails an answer that says "160 donors" for the post-QC atlas) now carry the
distinction.
