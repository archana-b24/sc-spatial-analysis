# C · A dashboard that explains the data, one source of truth, a public link

*Concept C of [docs/v3](README.md). Written 2026-11-12; §3 and "Cost: $0" revised later that day for the
free-only deploy (fixes session, [FIXES_RUN2.md](FIXES_RUN2.md)). Hosting facts verified on the official docs (URLs
inline).*

## Verdict

* **Single source of truth: done, and now enforced.** At runtime the app and the agent read only the DuckDB file.
  A test fails if any app file reads CSV/h5ad/YAML/JSON, and a schema-contract test fails if the DB lacks a table or
  column the app or agent relies on. A Docker-free "dry deploy" proves the bundle runs standalone.
* **Explain and interpret: done.** There is a new *Start here* page (workflow map + one computed headline per stage),
  a reader guide on every stage page (*What you're looking at · How to read it · What it means · Caveats*), a
  glossary, and *Ask the atlas* (since removed: zero-cost rule). Every interpretation sentence is computed from the DB when the page renders, so none
  can go stale after the §3.10 reruns.
* **Public access: free-only, ready for you to publish** (revised 2026-11-12 in the fixes session after your
  "keep everything free" decision). Hugging Face now requires a **paid plan (PRO)** to create Docker or Gradio Spaces
  and deprecated its Streamlit SDK on 2025-04-30, so the HF Space path was **removed** (scripts parked outside the
  repo). The app goes on **Streamlit Community Cloud** (free) and the DuckDB on a **public HF Dataset** (free),
  downloaded at app start and sha256-pinned. A local rehearsal (`cloud-dry-run`) passes. See [Cost: $0](#cost-0)
  for every service, its free tier and what would trigger a cost.

---

## 1. One source of truth (what changed and why)

**Wording, reconciled with the 2026-11-10 decision.** There are two layers, each with one job:

| Layer | What | Who reads it | Truth for |
|---|---|---|---|
| Computational | `data/processed/*.h5ad`, `*.zarr`, `results/tables/*.csv` | pipeline steps | the analysis (reproducible from code + config) |
| Serving | `results/dashboard/lusc_atlas.duckdb` | the dashboard, the MCP server (and the parked agent) | everything a reader sees or asks |

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
| ~~**Ask the atlas**~~ | Removed from the site 2026-11-12 (needs a paid LLM key); page parked outside the repo |
| **Provenance** | Now also shows the config/content hashes, the step stamps and the schema contract |

**Honesty about §3.10.** The computed interpretations name the problem where it applies: `Unknown` heads the paired
DE table but holds un-called tumour cells; scCODA's "decreases" are relative to `Unknown`; Milo's 0 hits are a
wrong-contrast artefact, not a null result; the 66 % "Rare epithelial" niche is a deconvolution artefact. The
dashboard does not fix them (the user deferred §3.10); it makes sure nobody can read a number without the caveat.
*Update 2026-11-13:* run 2 fixed these three ([FIXES_RUN2.md](FIXES_RUN2.md)); the interpretations now read
the run-2 tables, and the Known-issues page shows them as resolved.

**Charts.** No new chart types. The existing palette is the dataviz skill's validated reference palette, unchanged.
The palette validator could not be re-run because Node.js is not installed on this host. Condition colours stay
fixed across pages; cell types use highlight-one-grey-the-rest rather than 27 hues.

## 3. Public access (free-only; revised 2026-11-12, fixes session)

### 3.1 Options compared

Facts re-checked on 2026-11-12 on the official docs (links inline). Community Cloud's resource figures are the
ones its docs publish ("as of February 2024", "may change at any time without notice").

| | **Streamlit Community Cloud** + HF Dataset (**chosen**) | HF Docker Space (**removed**) | MotherDuck |
|---|---|---|---|
| Cost | Free | **Paid**: creating Docker/Gradio Spaces needs HF PRO (https://huggingface.co/docs/hub/spaces-overview); Streamlit SDK deprecated 2025-04-30 (https://huggingface.co/docs/hub/en/spaces-changelog) | Free "Lite" tier with usage caps; paid above them (not re-verified today; not needed) |
| Resources | CPU 0.078–2 cores, **RAM 690 MB–2.7 GB**, storage ≤ 50 GB (https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app) | — | database only; the app still needs a host |
| Idle behaviour | Apps without traffic for **12 h go to sleep**; any visitor can wake it ("Yes, get this app back up!"), first load then takes a minute or two (same page) | — | — |
| Source | GitHub repo; **public repo → public app**; only **one private app** at a time (https://docs.streamlit.io/deploy/streamlit-community-cloud/share-your-app) | — | — |
| Python | Default **3.12**; any version still receiving security updates, chosen under *Advanced settings* (https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/deploy) | — | — |
| Dependencies | Looks in the **entrypoint's directory first, then the repo root**; first file wins (`uv.lock` > `Pipfile` > `environment.yml` > `requirements.txt` > `pyproject.toml`); `requirements.txt` is installed with uv, falling back to pip (https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies) → our `app/requirements.txt` | — | — |
| Config / cwd | Working directory is **always the repo root**; `.streamlit/config.toml` must be at the root (https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/file-organization) → moved from `app/.streamlit/` to `.streamlit/` | — | — |
| Settings | *Advanced settings → Secrets* takes a `secrets.toml` body, read via `st.secrets` → `app/lib/db.py` reads `LUSC_DB_*` from env **or** `st.secrets` | — | — |
| Fits? | Yes: rehearsal peak **444 MB** RSS rendering all 15 pages in one process = 16 % of the 2.7 GB ceiling and below even the 690 MB minimum allocation | — | — |

**Data: public HF Dataset + sha256 pin, not committed to git.** The 119 MB DuckDB is git-ignored (derived output). It
is published as a public **HF Dataset** and downloaded by the app at startup through `huggingface_hub.hf_hub_download`,
pinned to the upload's **commit sha** (immutable) and verified by sha256. Why: (1) the dataset is also the "download
and query it yourself" link; (2) data and code can be updated independently; (3) the pin guarantees the site serves
the exact snapshot that was tested. HF public storage for free accounts is "best-effort" with no fixed quota, and the
docs ask that large uploads be useful to others (https://huggingface.co/docs/hub/storage-limits); 119 MB is far below
anything that page treats as large. Private storage (100 GB free) is not used.

### 3.2 What is ready

| File | Purpose |
|---|---|
| `scripts/publish_dataset.sh` (`pixi run -e dashboard publish-dataset`) | **You run it** (needs `HF_TOKEN`). Creates the *public* dataset if missing (refuses if it exists as private), uploads `docs/v3/dataset_card.md` as README and the DuckDB, then prints the secrets TOML to paste into Community Cloud (repo, upload commit sha, file sha256). **Never creates a Space**, never touches hardware, refuses files > 1 GB (`MAX_MB`), never prints the token. Not run here (rule: no uploads) |
| `scripts/cloud_dry_run.sh` (`pixi run -e dashboard cloud-dry-run`) | **Local rehearsal, no upload, no account.** (1) copies exactly what git would push (tracked + untracked-not-ignored) to a clean directory, so the DB and `.pixi/` are absent; (2) fresh **Python 3.12** venv from `app/requirements.txt` only; (3) writes `.streamlit/secrets.toml` in the clone; (4) seeds an offline HF cache, so the app fetches the DB through the same `hf_hub_download` call as in production, then checks the pin; (5) renders every page with AppTest from the repo-root cwd and reports peak RSS vs 2.7 GB; (5b) **negative control**: a wrong sha256 must stop the app; (6) serves the clone on :8502 and checks `/_stcore/health` |
| `app/lib/db.py` | `LUSC_DB`, `LUSC_DB_HF_REPO`, `LUSC_DB_HF_REVISION`, `LUSC_DB_SHA256` from env or `st.secrets` |
| `app/requirements.txt` | 8 pins, exactly the pixi `dashboard` env versions (streamlit 1.67.0, duckdb 2.0.0, plotly 7.1.0, pandas 2.3.3, numpy 2.5.4, pyarrow 26.0.1, pillow 13.0.0, huggingface-hub 2.4.0). No LLM SDKs |
| `.streamlit/config.toml` | Theme + headless; at the repo root where Community Cloud reads it |
| `docs/v3/dataset_card.md` | README of the HF Dataset |

**Rehearsal result (2026-11-12, run-1 DB):** 15/15 pages OK, DB resolved via secrets + `hf_hub_download` from the
seeded cache, wrong sha256 refused, health OK, **peak RSS 444 MB**. Re-run after the run-2 export: see
[FIXES_RUN2.md](FIXES_RUN2.md). What it cannot reproduce: Community Cloud's container and real limits, the real HF
download, GitHub access.

### 3.3 How anyone can try it

1. **Once deployed:** open the `*.streamlit.app` link. No login is needed to browse. If it has slept (12 h without
   visitors), click the wake-up button and wait a minute or two.
2. **Locally (full repo):**
   ```bash
   git clone https://github.com/archana-b24/sc-spatial-analysis && cd sc-spatial-analysis
   pixi install -e dashboard
   pixi run -e dashboard export          # build the DB from pipeline outputs …
   export LUSC_DB=/path/to/lusc_atlas.duckdb   # … or point at the downloaded snapshot
   pixi run -e dashboard app             # http://localhost:8501
   ```
3. **Just the data:** download `lusc_atlas.duckdb` from the HF Dataset, then
   ```bash
   duckdb -readonly lusc_atlas.duckdb
   SELECT * FROM tables_info;                       -- what every table is and where it came from
   SELECT * FROM known_issues;                      -- read these before quoting numbers
   ```
   or from Python (`duckdb.connect("lusc_atlas.duckdb", read_only=True)`), or from Claude Code / Claude Desktop
   through the MCP server (`pixi run -e dashboard mcp`; setup in `src/lusc_atlas/agent/mcp_server.py`).

### 3.4 Steps for you to deploy (all free)

1. Push the branch with the dashboard to the **public** GitHub repo (your call when; nothing was pushed).
2. Publish the data:
   ```bash
   export HF_TOKEN=hf_...            # free account; write token from https://huggingface.co/settings/tokens
   scripts/publish_dataset.sh        # shows the plan, asks, uploads; prints three secrets lines
   ```
3. On https://share.streamlit.io (sign in with GitHub): *Create app* → repo `archana-b24/sc-spatial-analysis`,
   your branch, main file `app/streamlit_app.py`; *Advanced settings* → Python **3.12** and paste the three printed
   lines into *Secrets*:
   ```toml
   LUSC_DB_HF_REPO = "<you>/lusc-atlas-db"
   LUSC_DB_HF_REVISION = "<commit sha printed by the script>"
   LUSC_DB_SHA256 = "<sha256 printed by the script>"
   ```
4. Deploy. First boot installs 8 packages and downloads 119 MB, then the app checks the sha256.
5. After any re-export: rerun step 2 and update the two secrets lines (revision, sha). Code-only changes redeploy
   automatically on push.

## Cost: $0

Rule (you, 2026-11-12): everything stays free; nothing may create, require or nudge towards a paid plan or API.
Audit of the deploy path as it stands after the fixes session:

| Service | Used for | Free tier we rely on | What would trigger a cost (avoid) |
|---|---|---|---|
| **Streamlit Community Cloud** | hosting the app | free; public app from a public repo; 0.078–2 CPU, 690 MB–2.7 GB RAM, ≤ 50 GB; sleeps after 12 h idle | Nothing in Community Cloud itself is billed. Costs would only come from moving to a different host (Snowflake-hosted Streamlit, a VM, an HF Space with paid hardware). Note: only **one private** app is allowed, so keep the repo public |
| **Hugging Face Hub, Dataset repo** | the 119 MB DuckDB + card | free account, public repo, "best-effort" public storage | Making the repo **private** (counts against a 100 GB private quota; above it only paid plans), buying a storage add-on, uploading GB-scale files repeatedly. `publish_dataset.sh` refuses private repos and files > 1 GB |
| **Hugging Face Spaces** | **not used** | — | Creating a Docker/Gradio Space needs PRO; any paid hardware tier is billed per hour. **No script in the repo creates a Space** (the old ones are parked outside the repo) |
| **GitHub** | public repo, CI (Actions) | public repos: Actions minutes free on standard runners (GitHub's long-standing policy; not re-verified today) | Private repo beyond the free Actions minutes, larger runners, Git LFS beyond the free quota (we don't use LFS) |
| **LLM APIs** (Anthropic, OpenAI, NVIDIA) | **not used** by the site or the pipeline | — | The chatbot page is removed; `anthropic`/`openai` stay only in the local pixi env for the parked agent and are **not** in `app/requirements.txt`. No key is read by the app. `agent-eval` was never run (it would cost ≈ $2–3) |
| **MotherDuck** | **not used** | — | Hosted DuckDB; not needed for a 119 MB read-only file |
| Local compute (this VM) | pipeline, export, dry run | your existing machine | none added; no cloud compute was used |

**How to keep it at $0:** keep the HF account free and the dataset public; never click "Create Space" or pick
hardware on HF; don't add a payment method anywhere for this project; if Community Cloud ever asks for an upgrade,
the app does not need it (440 MB peak on the run-2 DB vs 2.7 GB).

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
