# Design notes: agentic atlas · Nemotron · dashboard that explains the data

*Session of 2026-11-12, run unattended. Everything is **additive and uncommitted** on branch `v2.1-dashboard`, for
your review. Nothing was deployed, uploaded or pushed. The §3.10 pipeline issues were deliberately **not** touched
(deferred by you); the pipeline was not rerun. Only `results/dashboard/` was regenerated, via the `export` step.*

| Doc | Concept | One-line verdict |
|---|---|---|
| [A_agentic_research.md](A_agentic_research.md) | A · make it agentic | Build an **interpreter with provenance**, not an "AI scientist": read-only tools, caveats attached to data, SQL behind every number, mechanical number-grounding. **Built and tested; not yet run against a live model (no keys).** |
| [B_nemotron_assessment.md](B_nemotron_assessment.md) | B · NVIDIA Nemotron | **Yes as a pluggable open-weights second brain; no as default today.** No GPU here, the free hosted API is prototyping-only, and there is no evidence yet that it matches Claude on this task. Provider switch implemented; same eval decides. |
| [C_dashboard_hosting.md](C_dashboard_hosting.md) | C · dashboard + one source of truth + public link | **DuckDB is now the enforced single serving source of truth**; every page explains itself with interpretations computed from the data; a deploy bundle is ready and dry-run tested. **New finding: HF Docker Spaces now need a paid plan** → Streamlit Community Cloud is the free fallback. |

> **Decision 2026-11-12 (after this session): everything stays free, and the chatbot is out for now.**
> The *Ask the atlas* page is removed from the dashboard and the deploy bundle (the page file is parked outside
> the repo at `../.tooling/parked_chatbot/ask.py`), because a public chatbot needs a paid LLM API key. Hosting
> = Streamlit Community Cloud + a free HF Dataset (no HF PRO). The agent library (`src/lusc_atlas/agent/`), its
> tests, the eval set and the MCP server are kept but not wired into the site; the MCP server and the planned
> Claude Skill (concept 4) work through Claude Code with no extra cost. Mentions of *Ask the atlas* below
> describe the prototype as built, not the current site.

> **Update 2026-11-12 (fixes session, [FIXES_RUN2.md](FIXES_RUN2.md)): the paid deploy path is gone.** The HF
> Docker-Space bundle (`build_space.py`, `deploy_space.sh`, `space_dry_run.sh`, tasks `space-build` /
> `space-dry-run`) is parked outside the repo (`../.tooling/parked_paid_deploy_2026-11-12/`). The free-only
> replacements are `scripts/publish_dataset.sh` (`pixi run -e dashboard publish-dataset`: DuckDB → public HF
> Dataset, never a Space) and `scripts/cloud_dry_run.sh` (`pixi run -e dashboard cloud-dry-run`: clean clone, fresh
> Python 3.12 venv from `app/requirements.txt`, secrets TOML, DB fetched through `hf_hub_download`, every page
> rendered, wrong-sha negative control). See [C §3 and "Cost: $0"](C_dashboard_hosting.md#cost-0). Mentions of the
> Space bundle below describe what this session originally built. Open question 2 is answered: Community Cloud.

## Executive summary

**A: agentic application.** I surveyed ~25 agent systems (CellAgent, CellVoyager, SpatialAgent, scChat, Biomni,
GeneAgent, Virtual Lab, AI co-scientist, Robin/PaperQA2, ToolUniverse, BioMANIA, AutoBA, CRISPR-GPT, Kosmos,
POPPER …), the text-to-SQL and MCP literature, the evaluation benchmarks (BixBench 17 %, ScienceAgentBench 32 %,
DiscoveryBench ~25 % end-to-end) and the documented failure modes (hallucinated numbers, LLM-driven wrong
conclusions, p-hacking, prompt injection). The consistent lesson is that the agents that hold up on real data
**restrict the model to typed tools and verify before showing**. So the prototype in `src/lusc_atlas/agent/` gives
an LLM 11 read-only tools over the DuckDB: 3 generic, including guarded SQL (one SELECT, ≤ 200 rows, no file/network
access, config locked), and 8 curated ones (top DE genes, gene look-up, composition, communication, expression,
cohort, contrasts, caveats). Tools **attach the relevant known issues to their rows**. Every answer comes back with
the exact SQL and tables, and a **grounding check** flags any number in the answer that is not in a tool result.
The same tools are exposed as an **MCP server** and on an **Ask the atlas** page, and a **24-question eval** with
SQL-computed ground truth is ready to compare providers. The standout versus the landscape is caveat-aware tools,
mechanical grounding, one source of truth for chart and agent, and identical guardrails over MCP. Phase 2
(re-running pipeline steps) is proposed with human approval gates, after the eval and the §3.10 fixes.

**B: Nemotron.** It is genuinely interesting: open weights under commercial-use licences, released pre- and
post-training data, published recipes, and a hybrid Mamba-Transformer MoE line (Nemotron 3 Nano/Super/Ultra, 3.5
Lightning/Bolt) with a reasoning on/off switch. For *this* project it makes sense as a **reproducible, self-hostable
fallback brain**, not as the default. The tool-capable sizes need datacentre GPUs (Super: B200 / DGX Spark minimum).
The CPU-runnable 4B-class GGUFs are weak at multi-turn tool use by NVIDIA's own numbers (Tau2 29 for Nano-4B). The free
build.nvidia.com API is for prototyping only. There is no independent head-to-head against current Claude models.
The provider switch (`--provider nvidia`, default `nvidia/nemotron-3-super-120b-a12b`, correct reasoning flag per
generation) is implemented and mocked in tests; one ~$3 eval run with both keys will settle it on our own task.
Bio stack: BioNeMo (now bionemo-recipes) does include **Geneformer**, which this repo already benchmarks zero-shot.
Its value would be GPU fine-tuning, and rapids-singlecell would accelerate the pipeline. Both are GPU-only, so both
are parked until GPU time exists.

**C: dashboard.** The app already read only DuckDB at runtime. That is now **enforced** by a test, and a **schema
contract** test pins the tables and columns the app and agent need. Provenance grew too: config/content hashes,
pipeline step stamps, schema 2, and a sha256 pin so a deployed site refuses an unverified snapshot. Wording is
reconciled: h5ad/zarr/CSV are the *computational* source of truth, and DuckDB is the single *serving* source of
truth. The new **Start here** page shows the workflow map and one computed headline per stage. Every stage page now
has a four-part **reader guide** (what / how to read / what it means / caveats) whose interpretations are computed
from the DB at render time and name the §3.10 problems where they apply. There is also a **glossary**. Public
access: the Space bundle builder, a one-command `deploy_space.sh` (needs only `HF_TOKEN`, not run) and a local
**dry deploy** (clean venv, all 16 pages render from the bundle, health check OK, 433 MB peak) are ready. Verified
today: HF deprecated the Streamlit SDK (2025-04-30) and **creating Docker Spaces needs HF PRO**. So the recommendation
is the DuckDB as a free HF Dataset (sha256-pinned), plus the app on an HF Docker Space if you have PRO, otherwise on
Streamlit Community Cloud. The pending **151 vs 160 donors** question is resolved: 160 is the pre-QC cohort, and
QC dropped 26 samples under 200 cells, removing 9 donors, so 151 is post-QC. The log wording is corrected.

## Open questions for you

1. **What was the 4th concept?** Your message said "4 concepts" but listed three; I did those three.
2. ~~**Hosting:** do you have, or want, **HF PRO** (needed for a Docker Space; current price at https://huggingface.co/pricing, not verified here)?
   If not, I recommend Streamlit Community Cloud + HF Dataset (free); steps are in C §3.4. This revises the
   2026-11-10 "HF Spaces" decision only in *where the app runs*; the data plan is unchanged.~~
   **Answered 2026-11-12: free only → Community Cloud + public HF Dataset; no HF PRO.**
3. ~~**LLM keys and spend:** may I run the eval once with your `ANTHROPIC_API_KEY` (≈ $2–3 on Opus 5.5) and an
   `NVIDIA_API_KEY` (free trial) to get the first real numbers? And for the public site, should *Ask the atlas* use
   **your** key (capped at 10 questions per browser session), **visitors' own keys only** (current default when no
   secret is set), or be hidden?~~
   **Superseded 2026-11-12: everything stays free; the chatbot page is removed and no LLM API is called.** The
   eval stays unrun; it would only be run with a free option (e.g. a local open-weights model).
4. **Licence:** the repo has no LICENSE file, so the Space card asserts none. Which licence do you want for code and
   for the published DuckDB snapshot (whose content derives from CELLxGENE/10x data)?
5. ~~**Default model:** I used `claude-opus-5-5` (current default Opus). Prefer Sonnet 5.5 for cost on the public site?~~
   **Moot: no hosted model is used.**
6. **Agent phase 2** (propose pipeline re-runs with a config diff + your approval): in scope, or keep the agent to
   interpretation?
7. **Commit plan:** split into (a) the dashboard (already pending) and (b) the agent concepts, or one commit? Nothing
   is committed.
8. Minor, pre-existing: README says "LuCA (45 LUSC donors)" in the diagram and prose while the cohort uses 38 LuCA
   LUSC donors (45 is likely LuCA's availability before dedup and selection). Worth a wording check during review.

## Every file added or changed in this session

**Added**

| File | What |
|---|---|
| `src/lusc_atlas/agent/__init__.py`, `tools.py`, `llm.py`, `loop.py`, `mcp_server.py`, `evaluation.py`, `__main__.py` | The agent: tools, providers, loop + grounding, MCP server, eval grading, CLI |
| `app/views/start.py`, `app/views/glossary.py`, `app/views/ask.py` | Start here, Glossary, Ask the atlas pages |
| `app/lib/explain.py` | Reader guide + computed interpretations per stage |
| `tests/test_agent.py` | 24 tests: guards, loop with fake LLM, mocked Anthropic/OpenAI adapters, reasoning flags, eval integrity, MCP |
| `tests/agent_eval/questions.yaml` | 24-question eval with `truth_sql` |
| `scripts/run_agent_eval.py` | Eval runner (`--check` without LLM; `--providers anthropic,nvidia`) |
| `scripts/build_space.py`, `scripts/deploy_space.sh`, `scripts/space_dry_run.sh` | Space bundle, one-command deploy (not run), local dry deploy |
| `docs/design/README.md`, `A_agentic_research.md`, `B_nemotron_assessment.md`, `C_dashboard_hosting.md`, `dataset_card.md` | This documentation + the HF Dataset card |

**Changed**

| File | Change |
|---|---|
| `src/lusc_atlas/export.py` | Schema 2: `SCHEMA_CONTRACT` (+ `schema_contract` table), `page_guide`, `glossary`, `step_stamps`, meta `config_sha256` / `content_sha256` / `role`; docstring wording on the two sources of truth |
| `docs/dashboard_content.yaml` | + `page_guide` (9 stages) and `glossary` (26 terms); provenance decision reworded |
| `app/streamlit_app.py` | Navigation (Start here default, Glossary, Ask the atlas); finds the agent package in a plain checkout or bundle |
| `app/lib/db.py` | `LUSC_DB_HF_REVISION`, `LUSC_DB_SHA256` pin (refuses mismatching snapshot) |
| `app/views/{cohort,qc,integration,malignant,de,composition,trajectory,spatial,ccc}.py` | Reader guide inserted under each title (plus `ruff format` reflow) |
| `app/views/provenance.py` | Hashes, step stamps, schema contract; serving-vs-computational wording |
| `app/views/spatial.py` | Marker comment on the in-memory PNG decode (allowed by the runtime-read test) |
| `app/requirements.txt` | Pinned to the exact `dashboard` env versions; + anthropic, openai |
| `tests/test_dashboard.py` | + schema-contract, runtime-reads-only-DuckDB, computed-guides, sha256-pin tests |
| `pixi.toml` / `pixi.lock` | `dashboard` env + `anthropic`, `openai`, `mcp`, `huggingface-hub` (lock diff: **insertions only**; backup at `../.tooling/pixi.lock.pre-v3`); tasks `ask`, `mcp`, `agent-eval`, `space-build`, `space-dry-run` |
| `pyproject.toml` | E501 ignore extended to `app/lib/explain.py` (reader prose, same rule as `app/views`) |
| `.gitignore` | `build/`, `logs/agent_eval/` |
| `DEVELOPMENT_LOG.md` | §3.4 and §5 donor count 160 → 151 post-QC (wording only, with explanation); §4 pointer; timeline row 2026-11-12 |
| `docs/ROADMAP.md` | §7 pointer to docs/design; hosting note; changelog row |
| `README.md` | Preview commands under the dashboard section |
| `results/dashboard/lusc_atlas.duckdb` | Regenerated by `pixi run -e dashboard export` (git-ignored; 48 tables, schema 2, 119 MB) |

Outside the repo, all under /storage: `.tooling/pixi.lock.pre-v3` (backup), `.tooling/space_venv/` (dry-run venv,
deletable), `.tooling/space_home/`, `.tooling/space_dry_run.log`; `build/space/` inside the repo is git-ignored.

## How to try each piece

```bash
cd sc-spatial-analysis && source ../storage_env.sh
pixi run -e dashboard app                          # http://localhost:8501 → "Start here" (running now in tmux lusc-v2:app)
pixi run -e dashboard pytest -q tests/test_agent.py tests/test_dashboard.py   # 45 tests
pixi run -e dashboard agent-eval --check           # eval set + grader self-test, no LLM
export ANTHROPIC_API_KEY=…                          # then:
pixi run -e dashboard ask "How conservative is the malignant-cell call?"
pixi run -e dashboard agent-eval --providers anthropic,nvidia   # needs both keys; ≈ $2–3 on Opus 5.5
pixi run -e dashboard mcp                          # MCP stdio server (Claude Desktop / Claude Code setup in mcp_server.py)
pixi run -e dashboard cloud-dry-run                # rehearse the free Community Cloud deploy locally (no upload)
HF_TOKEN=… scripts/publish_dataset.sh              # YOU run it: DuckDB → public HF Dataset (asks first; never a Space)
```

## Checks run at the end of the session
* `pixi run lint` (ruff check src scripts tests app): **pass**; `ruff format --check`: **pass**.
* `pixi run test` (default env, as CI runs it): **15 passed, 2 skipped** (dashboard/agent tests skip there).
* `pixi run -e dashboard pytest tests/test_dashboard.py tests/test_agent.py`: **45 passed**.
* `agent-eval --check`: 24/24 truth queries run; grader self-test 0 failures.
* `space_dry_run.sh`: 16/16 pages render from the bundle in a clean venv; health OK; peak RSS 433 MB.
* Live app restarted on the new version in `lusc-v2:app`; `curl localhost:8501/_stcore/health` → `ok`.

## Recommended next steps
1. Answer the open questions (4th concept, hosting plan, keys).
2. Run `agent-eval` once with keys; record the table in A §6 and B §5 (replacing the token estimates).
3. Review the wording on *Start here* and the reader guides in the running app; commit (your call on the split).
4. Deploy (Dataset first with `publish_dataset.sh`, then the app on Streamlit Community Cloud; C §3.4).
5. Then return to the §3.10 fixes. Afterwards `export` → the computed interpretations and the eval truths update
   automatically, and the eval doubles as a regression check that the fixes changed what they should.
