# Fixes session: free-only deploy + pipeline run 2

*2026-11-12 → 2026-11-13, run unattended on branch `v2.1-dashboard`. Everything is **uncommitted** for your
review. Nothing was deployed, uploaded, pushed or created on any external service. The detailed engineering
record is [DEVELOPMENT_LOG §1.3](../../DEVELOPMENT_LOG.md) (hosting) and §3.11 (run 2); this page is the summary.*

## TL;DR

1. **Hosting is now free-only, by construction.** The HF Docker-Space bundle (which needs HF PRO) is parked
   outside the repo. The only path left: you publish the DuckDB snapshot as a **public HF Dataset**
   (`scripts/publish_dataset.sh`, run by you) and deploy the app on **Streamlit Community Cloud** from the public
   GitHub repo. A local rehearsal (`scripts/cloud_dry_run.sh`) reproduces the cloud install and DB download and
   passed on the final DB: 15/15 pages, peak memory 440 MB of the 2.7 GB ceiling.
2. **The three high-severity issues from the run-1 review (§3.10 #1–#3) are fixed, and the full downstream pipeline
   was rerun** (run 2: 19.8 h, no failures). Milo went from a meaningless 0 to **1,671 of 6,225** neighbourhoods
   changing between tumour and adjacent tissue. The 25k tumour cells hidden in `Unknown` now form their own group.
   The Rare-epithelial artefact in the tissue maps is gone. #7 (stale report provenance) is fixed too.
3. **The run-1 conclusions mostly survive.** The scCODA directions hold with a stable reference, the differential
   CCC top hits are unchanged, and the tensor still does not separate tumour from adjacent tissue (#4 stays open).

## 1. Free-only deploy (Part 1 of the brief)

| Change | Why |
|---|---|
| Parked `build_space.py`, `deploy_space.sh`, `space_dry_run.sh`, the old `app/.streamlit/config.toml`, pre-edit doc copies → `../.tooling/parked_paid_deploy_2026-11-12/` (with a README) | HF Docker/Gradio Spaces need a paid plan; keeping a script that creates one is a nudge toward paying. Parked, not deleted (recoverable) |
| `scripts/publish_dataset.sh` (`pixi run -e dashboard publish-dataset`) | Uploads the card + DuckDB to a **public** dataset only; refuses private repos and files > 1 GB; never creates a Space or touches hardware; token only from `HF_TOKEN`, never printed; prints the Streamlit secrets TOML (repo, pinned revision, sha256) |
| `scripts/cloud_dry_run.sh` (`pixi run -e dashboard cloud-dry-run`) | Rehearses Community Cloud offline: clean copy of what git would push, fresh Python 3.12 venv from `app/requirements.txt` only, `secrets.toml`, DB through the real `hf_hub_download` call from a seeded offline cache, every page rendered, peak RSS vs the 2.7 GB ceiling, a wrong-sha negative control, `/_stcore/health` |
| `app/lib/db.py`: `setting()` reads env vars, then `st.secrets` | Community Cloud passes configuration via its Secrets box, not env vars |
| `.streamlit/config.toml` moved to the repo root | Community Cloud runs from the repo root and reads config from there |
| Docs: C_dashboard_hosting §3 rewritten + new **"Cost: $0"** table, README, dataset card, ROADMAP, DEVELOPMENT_LOG §1.3 | Each service's free-tier limits and what would trigger a cost |

**Cost audit result: $0.** Nothing in the repo creates, requires or nudges toward a paid plan or API:
- No script creates a Space or sets hardware.
- The dataset is public by enforcement.
- No LLM API is called by the site or the pipeline: `anthropic`/`openai` are not in `app/requirements.txt`, and the
  chatbot was removed earlier. `agent-eval` was never run live; `--check` needs no LLM.
- No cloud compute was used; everything ran on this VM.
- Old open questions in `docs/v3/README.md` and `B_nemotron_assessment.md` that suggested paid API spend or HF PRO
  are struck through and marked superseded (kept for the record).

Free-tier limits relied on (verified 2026-11-12 on the providers' docs):
- **Community Cloud:** 0.078–2 CPU, 690 MB–2.7 GB RAM, ≤ 50 GB storage; sleeps after 12 h idle, and anyone can wake
  it; a public repo gives a public app (one private app only).
- **HF:** public storage is free (best effort); private storage has a 100 GB quota.

Full table: [C_dashboard_hosting.md → Cost: $0](C_dashboard_hosting.md#cost-0).

**What you do to go live (all free):** `export HF_TOKEN=…; pixi run -e dashboard publish-dataset`, push the repo
to a public GitHub repo, create the app on share.streamlit.io (entrypoint `app/streamlit_app.py`, Python 3.12),
and paste the printed TOML into *Advanced settings → Secrets*. Steps: C §3.4.

## 2. Pipeline fixes (Part 2: §3.10 #1–#3, #7)

| # | Problem in run 1 | Fix (config-driven) | Why this option |
|---|---|---|---|
| 1 | Milo tested `tumour_region_unknown` (2 donors) vs adjacent: the wrong contrast | Only donors with both tumour and adjacent samples (14 donors, 48 samples); design `~donor_id + condition`; explicit contrast `conditiontumour - conditionadjacent_normal`; graph rebuilt on that subset; failures now fail the step | Matches the primary paired pseudobulk DE. A random donor effect (slower) and testing within all conditions (neighbourhoods from untestable cells) were rejected. A unit test checks the contrast against pertpy's own design coding |
| 2 | 25k author-labelled tumour cells below the CNV threshold sat in `Unknown`, which was also scCODA's automatic reference | New analysis group `Tumour epithelial (CNV-low)` (in `de_group`; `harmonised_label` untouched, so no re-integration); scCODA reference fixed to **cDC**; new contrast `cnvlow_tumour_vs_healthy_basal` | cDC is present in 87 % of samples (≥ 81 % per condition) with a flat mean proportion across conditions. **Fibroblast / EC capillary (the brief's suggestion) are absent from 29–44 % of samples**, because several LuCA studies are CD45-enriched |
| 3 | A 66-cell `Rare epithelial` signature absorbed 66 % of the tumour niche in cell2location | Reference labels with < 200 cells are dropped (`spatial.c2l.min_ref_cells`), `Unknown` is excluded, and the decisions go to `c2l_reference_labels.csv` | No sensible parent to merge into; a mixed bucket is exactly the kind of catch-all that absorbs signal |
| 7 | The report showed the previous run's provenance | The report writes its own provenance first | Trivial |

**How the new group is treated downstream:**
- **Excluded** from the paired and unpaired DE, like `Malignant`, since it has no counterpart in control tissue.
- **Included** in scCODA, Milo, LIANA and the tensor.
- **Not used** for the trajectory, which stays on CNV-confirmed cells.

Verification before the 20 h rerun:
- 9 new unit tests (`tests/test_fixes_run2.py`);
- `scripts/smoke_test.py --run2` (12k real cells including all 48 paired samples, every changed stage);
- `--c2l` (reference decisions → both cell2location models, tiny epochs).

## 3. Before / after

| | Run 1 (2026-11-10) | Run 2 (2026-11-13) |
|---|---|---|
| Milo significant neighbourhoods | 0 / 17,511 (wrong contrast) | **1,671 / 6,225**: 917 up, 754 down in tumour (CD4 T 280↑/59↓, B 261↑, CD8 T 155↑/154↓, Macrophage 56↑/221↓, NK 26↑/112↓) |
| scCODA tumour vs adjacent (credible log2FC) | ref `Unknown`: NK −0.95, CD8 −0.82, CD4 −0.67, Mono NC −0.33, AT2 −0.30, Mac −0.22, AT1 −0.19 | ref `cDC`: NK −1.15, CD8 −0.80, CD4 −0.66, Mono NC −0.43, AT2 −0.37, Mac −0.33, AT1 −0.29, EC lymphatic −0.12; **CNV-low tumour +2.23** |
| `Unknown` | 57,335 cells | 28,224 (unannotated HLCA + 64 adjacent-tissue author-malignant cells) |
| `Malignant` (CNV) | 3,942 | 3,942 (threshold reproduced exactly) |
| Paired DE significant rows | 9,665 (`Unknown` 2,181) | 8,474 (`Unknown` 990; every other group identical) |
| Tumour epithelium vs healthy basal | Malignant 300 genes | Malignant 300; CNV-low 1,040 (SOX2 +2.2, TP63 +1.4, CDKN2A +5.6 log2FC) |
| cell2location reference | 27 labels | 25 (Rare epithelial 66, Mesothelial 91 dropped; `Unknown` excluded) |
| Tumour niche | Rare epithelial 66 %, Malignant 11 % | CNV-low tumour 40–62 %, Malignant 12–14 % |
| CellCharter niches | 3 | 9 (auto-K) |
| Tensor | 82 × 868 × 10 × 10, rank 5 | 82 × 951 × 11 × 11 (CNV-low enters), rank 6 |
| Tensor factor vs condition | all p ≥ 0.31 | all p ≥ 0.34 (still not tumour-specific) |
| Differential CCC top 10 | 1 row from `Unknown` (NECTIN4→TIGIT) | the same rows minus that one |

How to read this:
- **Milo vs scCODA.** Milo adds what scCODA cannot see. CD8 and CD4 T-cell *states* move in both directions, so
  some T-cell phenotypes are enriched in tumour even though total T-cell proportions fall.
- **The CNV-low group is tumour.** Its SOX2/TP63 profile (the 3q LUSC programme) shows these are tumour cells with a
  weak CNV signal, not mislabelled normal epithelium.

## 4. Things found and fixed after the run

- **Dashboard crash:** 9 niches are more than the 8-colour palette, so the Spatial page raised `IndexError`. The
  dashboard tests caught it. Fixed with `ui.series_color()` (8 validated colours + 4 extra, then cycle) and a
  regression test.
- **Content refresh:** `docs/dashboard_content.yaml` was updated (issue 3 resolved, rank 6), and the DB was
  re-exported so the Known-issues page shows it.
- **Final checks, all passing:**
  - `pixi run lint`;
  - `ruff format --check`;
  - `pytest -q` (24 passed, 2 skipped);
  - dashboard + agent tests (45 passed);
  - `agent-eval --check` (24 items, 0 failures, truths recomputed on the run-2 DB);
  - `cloud_dry_run.sh` on the final DB.
- The app in tmux `lusc-v2:app` runs the new version at http://localhost:8501.

## 4b. Updating the live site to run 2

The site you deployed during this session runs the **run-1** snapshot from the `sc-deploy` worktree (branch
`v2-dev`, which holds the dashboard and free-deploy scripts but none of the §3.10 fix code). To move it to run 2:

1. Bring this branch's changes onto `v2-dev`. Expect small overlaps in `app/lib/db.py`, `scripts/cloud_dry_run.sh`,
   `README.md`, `docs/v3/README.md` and `docs/v3/C_dashboard_hosting.md`; the versions here are the newer ones.
   Schema v3 adds `cells.de_group`, `c2l_reference_labels` and new Milo/scCODA columns, so the app code and the DB
   must move together.
2. Copy `results/dashboard/lusc_atlas.duckdb` (run 2, 120 MB) into the worktree and run `cloud-dry-run` there.
3. Run `publish-dataset` again. It prints a new revision and sha256; paste both into the app's Secrets. The sha
   pin makes the app refuse a mismatched DB rather than serve it silently.

## 5. Honest limitations (still open)

- **#4:** the tensor factors describe sample-level programmes, not tumour-specific ones. A study-adjusted test of the
  loadings is a v2.1 item.
- **#5 (SOLO vs Scrublet) and #6 (atlas mixing):** documented, not changed.
- **Spot-level CNV-low vs Malignant:** the split within a spot is not reliable (their signatures overlap); quote the
  two together as tumour epithelium.
- **Mesothelial** cells can no longer be mapped spatially (91 reference cells).
- **Tumour niches** still hold Multiciliated (11–14 %) and Plasma (6–9 %) cells. That may be real or partial
  absorption; one section can't tell.
- **CNV-low vs basal DE:** the top-ranked genes are mitochondrial pseudogenes (a study/chemistry artefact); read the
  marker genes.
- **Adjacent tissue:** the 64 author-malignant cells in adjacent-normal samples stay `Unknown`.
- **Trajectory** still covers Basal + CNV-called Malignant only.

## 6. Open questions for you

1. **Licence.** The repo has no LICENSE file. Suggested: **MIT** for the code and **CC BY 4.0** for the published
   DuckDB snapshot. The snapshot is derived from CELLxGENE (HLCA, LuCA: CC BY 4.0) and 10x Genomics data, so check
   the 10x dataset terms before publishing.
2. **Commit split.** Suggested three commits:
   - (a) v2.1 dashboard;
   - (b) v3 concepts (agent library, docs/v3 A–C);
   - (c) this fixes session (free-only deploy + run-2 fixes, tests, docs).

   Nothing is committed. Generated results stay git-ignored.
3. **Concept 4 (Claude Skill).** Not started. It works through Claude Code at no extra cost when you want it.
4. **scCODA reference.** I used `cDC` instead of the suggested Fibroblast / EC capillary (data in §2). Optional
   sensitivity check: rerun scCODA with `Macrophage` as the reference (about 1 h, `abundance` step only).
5. **Trajectory.** Add the CNV-low group as a second endpoint (v2.1 option), or keep it CNV-confirmed only?
6. **Adjacent-normal author-malignant cells (64).** Keep them as `Unknown` (current choice) or drop them?

## Files touched in this session

- **Added:**
  - `scripts/publish_dataset.sh`, `scripts/cloud_dry_run.sh`, `.streamlit/config.toml` (moved);
  - `tests/test_fixes_run2.py`;
  - `docs/v3/FIXES_RUN2.md`.
- **Changed:**
  - `config/config.yaml`;
  - `src/lusc_atlas/{de,deconvolve,steps,report,export}.py`, `src/lusc_atlas/agent/{__init__,tools}.py`;
  - `app/lib/{db,ui,explain}.py`, `app/views/{composition,spatial}.py`;
  - `scripts/smoke_test.py`, `tests/test_dashboard.py`, `tests/agent_eval/questions.yaml`;
  - `pixi.toml` (tasks only; `pixi.lock` unchanged, backup in `../.tooling/`), `app/requirements.txt` (header only);
  - `.gitignore` (`logs/smoke/`);
  - `README.md`, `DEVELOPMENT_LOG.md`, `docs/dashboard_content.yaml`, `docs/ROADMAP_v2.1.md`;
  - `docs/v3/{README,A_agentic_research,B_nemotron_assessment,C_dashboard_hosting,dataset_card}.md`.
- **Parked (outside the repo):** `../.tooling/parked_paid_deploy_2026-11-12/`.
- **Archived (git-ignored):** `results/_archive_run1_2026-11-10/`, holding run 1's tables, reports, figures, DuckDB,
  stamps and log.
