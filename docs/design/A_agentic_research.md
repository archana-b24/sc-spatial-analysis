# A · Making the atlas agentic: research, design, prototype

*Concept A of [docs/design](README.md). Written 2026-11-12. Literature checked with web search/fetch that day; every
citation lists its ID/URL, and anything confirmed only from a secondary source is marked **(sec.)**, anything not
confirmed **UNVERIFIED**.*

**Verdict.** The most useful and defensible agent for this project is one that **interprets the atlas, does not
redo the science**: a small tool-using LLM that can only read the DuckDB snapshot through read-only tools, attaches
the known caveats to the rows it reads, returns the exact SQL behind every number, and is checked after the fact so
that no number in an answer is made up. That is built and tested (§5). Re-running or re-parameterising pipeline
steps (the "AI scientist" direction) is a later phase, gated on the eval in §6 and on the §3.10 fixes.

---

## 1. Why this shape, in one paragraph

Published agent benchmarks for biology and data analysis put end-to-end success at **17–35 %** (BixBench 17 % open
answer; ScienceAgentBench 32.4 %; DSBench 34 %; DA-Code 30.5 %; DiscoveryBench ~25 %; table §3). Evaluations of
"AI scientist" systems report coding failures in 42 % of experiments and wrong conclusions in roughly a third of
hypotheses when LLMs sit in the analysis loop (§4). The systems that work best on data **restrict the model to
typed, documented tools** (scChat, SpatialAgent's skill templates, ToolUniverse, BioMANIA) and **verify before
showing** (GeneAgent, CellTypeAgent, CASSIA; Kosmos cites every statement to code). So the agent here is narrow
on purpose: it is an *interpreter with provenance*, sitting on a database that is already the single serving
source of truth (Concept C).

---

## 2. Landscape: LLM agents for single-cell, spatial and biomedical analysis

### 2a. Single-cell and spatial agents (closest to this project)

| System | Citation | What it does | Architecture | Evaluation | What we borrow |
|---|---|---|---|---|---|
| **CellAgent** | Xiao Y … Peng J. bioRxiv 10.1101/2024.05.13.593861; arXiv:2407.09811 (2024); revised version **ICLR 2026** ("…Natural Language-Based Single-Cell Analysis"). https://arxiv.org/abs/2407.09811 | End-to-end scRNA-seq (ICLR version adds spatial) from natural-language requests | Multi-agent: planner, executor, evaluator; code execution; evaluator-triggered re-runs | >50 scRNA-seq datasets vs experts and methods (no numbers in the abstract) | An evaluator step before display → our post-hoc grounding check |
| **scAgent** | Mao Y … Gao Y (Zhejiang). arXiv:2504.04698 (2025). Not the same as *scAgents* (Gerstein lab, ICML 2025 workshop) | Universal cell-type annotation, novel-type detection | Single agent (planning + memory) + scGPT with tissue-specific LoRA plug-ins | 160 cell types / 35 tissues | Tissue-specific tools selected by metadata |
| **CellVoyager** | Alber S … Zou J (Stanford). bioRxiv 10.1101/2025.06.03.657517 (2025); *Nature Methods* 2026 **(sec.)**. https://github.com/zou-group/CellVoyager | Re-analyses published scRNA-seq data to propose analyses the paper did not run | Single agent writing and running code in a Jupyter notebook; self-review | CellBench (50 studies in the preprint): up to +23 % over GPT-4o/o3-mini at predicting the authors' analyses | Every agent action is a reproducible artefact → we log tool calls + SQL per answer |
| **SpatialAgent** | Wang H … Leskovec J, Regev A (Genentech/Stanford). bioRxiv 10.1101/2025.04.03.646459 (2025). https://github.com/Genentech/SpatialAgent | Spatial biology end to end (panel design, annotation, CCI, mapping) | Plan-Act-Conclude loop, 72 tools, **17 skill templates**, model-agnostic | ~2M cells (brain, heart, colitis); abstract claims expert-level | **Skill templates** → our curated semantic tools |
| **STAgent** | Liu lab, Harvard. bioRxiv 10.1101/2025.04.01.646731 (2025); first author UNVERIFIED | Spatial transcriptomics analysis + report | Multimodal agent that looks at its own plots | One case study | Later: let the agent look at the chart it cites |
| **scChat** | Li group, Purdue. bioRxiv 10.1101/2024.10.01.616063 (2024); first author UNVERIFIED | Conversational scRNA-seq co-pilot | Multi-agent + RAG; **every quantitative step is a structured function call** | Annotation accuracy, case studies | **Closest analogue**: the LLM never computes numbers itself |
| **CASSIA** | Xie E … Kendziorski C. bioRxiv 10.1101/2024.12.04.626476; *Nat Commun* 2025 **(sec.)** | Reference-free annotation with rationale and quality score | Annotator / validator / scorer agents | >970 populations | Per-answer quality signal (our grounding report) |
| **CellTypeAgent** | UNC. arXiv:2505.08844 (2025) | Annotation verified against CELLxGENE expression | Generate → verify against data | 9 datasets, 303 types | Verify claims against data before showing |

### 2b. General biomedical agents, tool layers and agent science

| System | Citation | What it does | Architecture | Evaluation | What we borrow |
|---|---|---|---|---|---|
| **Biomni** | Huang K … Leskovec J (Stanford). bioRxiv 10.1101/2025.05.30.656746 (2025); *Science* 2026, 10.1126/science.adz4351 **(sec.)**. https://github.com/snap-stanford/Biomni | General biomedical agent | Environment of ~150 tools / 105 packages / 59 databases; code-as-action; retrieval-based tool selection | LAB-Bench DbQA 74.4 % (human 74.7 %), SeqQA 81.9 % **(sec.)** | Tool retrieval once the tool set grows |
| **BioDiscoveryAgent** | Roohani Y … Leskovec J. arXiv:2405.17631; **ICLR 2025** | Designs CRISPR screen rounds | Agent + critic agent + literature/code tools | +21 % hit detection over Bayesian optimisation (v3) | Feed prior results into the next step |
| **GeneAgent** | Wang Z … Lu Z (NCBI/NLM). *Nature Methods* 2025, 10.1038/s41592-025-02748-6; arXiv:2405.16205 | Gene-set function naming with fewer hallucinations | Generate → **self-verify each claim against databases** → revise | 1,106 gene sets; experts: 92 % of verification decisions correct | Claim-level verification → our number-grounding check |
| **Virtual Lab** | Swanson K, Wu W, Bulaong NL, Pak JE, Zou J. *Nature* 2025, 10.1038/s41586-025-09442-9 | AI "lab" that designed SARS-CoV-2 nanobodies | PI agent + specialist agents + **scientific critic**, human-guided | Wet lab: >90 % of 92 nanobodies expressed; 2 new binders | Critic pass (future: LLM-judge with rubric, §6) |
| **AI co-scientist** | Gottweis J … Natarajan V (Google). arXiv:2502.18864 (2025); *Nature* 2026, 10.1038/s41586-026-10644-y **(sec.: press release)** | Hypothesis generation and ranking | Gemini multi-agent: generation, reflection, Elo tournament ranking, evolution | Expert ratings; wet-lab AML repurposing, liver fibrosis | Pairwise ranking if we ever rank candidate findings |
| **Robin / PaperQA2** (FutureHouse) | Ghareeb AE … Rodriques SG. arXiv:2505.13400; *Nature* 2026, 10.1038/s41586-026-10652-y **(sec.)**. PaperQA2: Skarlinski MD … White AD, arXiv:2409.13740 | Literature → hypothesis → data analysis (Finch sub-agent; no standalone Finch paper found) | Separate literature agents (PaperQA-based) and data agent | Ripasudil for dry AMD (wet lab); LitQA2 precision 85.2 % | **Keep literature and data tools separate** (prompt-injection boundary, §4) |
| **ToolUniverse** | Gao S … Zitnik M (Harvard). arXiv:2509.23426 (2025–26) | Tool ecosystem for "AI scientists" | Typed tool schemas, any LLM; MCP servers | Case studies | Our tools are typed JSON-schema tools, exposed over MCP too |
| **BioMANIA** | Dong Z, Zhong V, Lu YY. bioRxiv 10.1101/2023.10.29.564479 (journal version UNVERIFIED) | Turns a Python library into a chatbot | API retrieval + argument filling + error correction | 12 tools; taxonomy of tool-design issues | Generate question → tool-call pairs as tests |
| **AutoBA** | Zhou J … Gao X (KAUST). *Advanced Science* 2024, 10.1002/advs.202407094 | Runs omics pipelines from a description | Planner, code generator, **automated code repair** | End-to-end 65 % → 87.5 % with error feedback | Tool errors go back to the model (implemented) |
| **CRISPR-GPT** | Qu Y … Cong L. *Nat Biomed Eng* 2025, 10.1038/s41551-025-01463-z | Gene-editing experiment design | Planner + state machine + user checkpoints + safety guardrails | Wet lab by novices (A549, a lung adenocarcinoma line, not LUSC) | Explicit checkpoints before any action with side effects (phase 2) |
| **Kosmos** | Mitchener L et al. (Edison Scientific/FutureHouse). arXiv:2511.02824 (2025) | Autonomous discovery runs | Data + literature agents; **every report statement cited to code or literature** | 79.4 % of statements judged accurate; data-analysis statements ~85 % vs interpretive ~58 % (breakdown from search) | Provenance per statement; hold interpretation to a higher bar than lookup |
| **POPPER** | Huang K, Jin Y, Li R, Li MY, Candès E, Leskovec J. **ICML 2025**, arXiv:2502.09858 | Hypothesis validation by sequential falsification | Agent designs falsification tests; sequential testing with Type-I error control | Matches experts, ~10× faster | Phase-2 pattern for hypothesis questions |

Name clashes worth knowing: *scAgent* ≠ *scAgents*; *GeneAgent* (NCBI) ≠ *GenoAgent/GenoMAS* (Wang lab, arXiv:2507.21035);
CellAgent exists in three versions with two titles; no system called "OmicsAgent" was found (nearest: OmicsNavigator,
bioRxiv 10.1101/2025.07.21.665821, and the agent benchmark arXiv:2508.13201); ether0 (FutureHouse) is a chemistry
reasoning *model*, not an agent.

### 2c. Text-to-SQL and the tool layer

| Item | Citation | Key point for us |
|---|---|---|
| BIRD | Li J et al., NeurIPS 2023 D&B, arXiv:2305.03111 | GPT-4 54.9 % vs humans 93 % execution accuracy on dirty, large DBs **(search summary)** |
| Spider 2.0 | Lei F et al., ICLR 2025, arXiv:2411.07763 | Best agent 21.3 % on enterprise schemas (>1,000 columns) → **don't expose wide raw tables as the main interface** |
| BIRD-Interact | Huo N et al., ICLR 2026, arXiv:2510.05318 | Multi-turn, ambiguous requests; frontier models low (exact numbers UNVERIFIED) |
| DIN-SQL / MAC-SQL / CHESS / ReFoRCE | arXiv:2304.11015; 2312.11242; 2405.16755; 2502.00675 | Schema linking, decomposition, execution-feedback repair, column exploration |
| DuckDB-NSQL-7B | MotherDuck + Numbers Station, 2024, https://huggingface.co/motherduckdb/DuckDB-NSQL-7B-v0.1 | Small DuckDB-specific model; dated (DuckDB 0.9.2) → baseline only |
| MotherDuck MCP server | https://github.com/motherduckdb/mcp-server-motherduck | Read-only by default, row/char caps; README warns read-only alone is not sufficient |
| Securing DuckDB | https://duckdb.org/docs/2.0/operations_manual/securing_duckdb/overview.html | `enable_external_access=false`, no extension autoload/autoinstall, then `lock_configuration=true` (**implemented**) |
| MCP | Anthropic, 25 Nov 2024, https://www.anthropic.com/news/model-context-protocol; donated to the Linux Foundation's Agentic AI Foundation 9 Dec 2025 **(search)**; Python SDK `mcp` 2.3.0 (2 Oct 2026), where `FastMCP` became `MCPServer` | Our MCP server targets `mcp` 2.x |
| Building effective agents | Schluntz E, Zhang B, Anthropic, 19 Dec 2024, https://www.anthropic.com/engineering/building-effective-agents | Workflows vs agents; "poka-yoke" tools that are hard to misuse |
| Writing effective tools for agents | Aizawa K, Anthropic, 11 Sep 2025, https://www.anthropic.com/engineering/writing-tools-for-agents | Few high-impact tools, readable fields, truncation with steering messages, eval-driven |

---

## 3. Evaluation landscape

| Benchmark | Citation | Measures | Headline |
|---|---|---|---|
| BixBench | Mitchener L et al., arXiv:2503.00096 (2025) | Multi-step bioinformatics analysis, open answers | 17 % open-answer accuracy (GPT-4o, Claude 3.5 Sonnet) |
| LAB-Bench | Laurent JM et al., arXiv:2407.10362 (2024) | >2,400 MCQs: literature, DB access, sequences, protocols | Human PhD baselines |
| ScienceAgentBench | Chen Z et al., ICLR 2025, arXiv:2410.05080 | 102 data-driven discovery tasks from 44 papers | Best agent 32.4 % (34.3 % with expert knowledge) |
| DiscoveryBench | Majumder BP et al., ICLR 2025, arXiv:2407.01725 | Hypothesis search and verification on datasets | ~25 % |
| DSBench / DA-Code | arXiv:2409.07703 (ICLR 2025) / arXiv:2410.07331 (EMNLP 2024) | Data-science agent tasks | 34.1 % / 30.5 % |
| BioCoder / BioDSBench | arXiv:2308.16458 (*Bioinformatics* 2024) / arXiv:2410.21591 | Bioinformatics code; biomedical data science | LLMs fail on data understanding and practice |

Lesson: aim for **verifiable partial answers** and build a **domain eval with gold answers computed from the data**
(BixBench / DiscoveryBench style), not autonomous discovery.

## 4. Failure modes and what each guardrail answers

| Failure mode | Evidence | Guardrail in this repo |
|---|---|---|
| Hallucinated numbers | Kosmos: interpretive statements far less accurate than data statements; Spider 2.0 / BIRD SQL errors | Every number must come from a tool result: **post-hoc grounding check** (`loop.check_grounding`) flags any number that does not match a result from this run; shown in the UI and scored in the eval |
| Wrong conclusions / p-hacking | Baumann J, Röttger P, Urman A et al., "LLM Hacking", arXiv:2509.08825 (~31 % wrong conclusions with state-of-the-art models); Asher, Malzahn et al. 2026 working paper on Claude Code/Codex p-hacking (design verified, results UNVERIFIED); Kosch & Feger, arXiv:2504.14571 | The agent **cannot run new statistics**: it can only read the pipeline's pre-specified tests (pseudobulk DESeq2, scCODA, LIANA). The system prompt requires padj and the contrast name; eval item q24 asks it to report a raw-p "hit" and expects a refusal |
| Ignoring known problems | §3.10 issues (Unknown cells, Milo contrast, Rare-epithelial c2l) | **Caveats travel with the data**: tools attach the relevant open issues to their rows; eval items q07, q11, q12, q17 require the caveat |
| Provenance loss | CellVoyager, AutoBA, Kosmos all keep code/citations per step | Every tool result carries the exact SQL and tables; the app shows them under the answer; the CLI prints them |
| Destructive or exfiltrating actions; prompt injection via data | Greshake K et al., arXiv:2302.12173; Willison's "lethal trifecta" (16 Jun 2025) | Read-only DB, single-SELECT parser check, external access off and config locked, no network/write tools, row caps, a query timeout; tool output is labelled as data in the prompt; no literature tool shares the data agent (removes the exfiltration leg) |
| Answering from memory | General LLM behaviour | System prompt + eval scoring: an answer with zero tool calls fails; negative questions (q21, q22: survival, EGFR mutations) must say the atlas does not contain this |
| Error loops / runaway cost | AutoBA repair loop; tool-use agents | Tool errors are returned to the model (bounded repair); hard step cap (10) with a forced final answer; per-session question cap on the public site |

---

## 5. Design for this repository (as built)

### 5.1 What the agent is for
* **Phase 1 (built): interpret and explore the atlas.** "Which genes go up in tumour macrophages?", "How much should I
  trust the composition result?", "Is SPP1 changed, and where?". Read-only, cited, caveat-aware.
* **Phase 2 (not built; proposal): re-run or re-parameterise pipeline steps.** Every step is already a pure function
  with declared inputs/outputs (`lusc_atlas.steps`), so a tool like `propose_rerun(step, config_patch)` could write a
  config diff and a dry-run plan for **human approval** (CRISPR-GPT-style checkpoint), never execute directly.
  Prerequisites: the eval below passing, and the §3.10 fixes, so the agent is not let loose on known-broken outputs.

### 5.2 Architecture

```
question ─▶ ChatSession (Claude | Nemotron/OpenAI-compatible | Fake)
               │ tool calls (JSON schema, strict)
               ▼
           AtlasTools ──▶ DuckDB (read_only, external access off, config locked)
               │  rows + SQL + tables + caveats
               ▼
           loop.ask ──▶ answer + provenance + grounding report ──▶ Streamlit page / CLI / eval
                                                      MCP server ──▶ Claude Desktop / Claude Code
```

| Piece | File | Notes |
|---|---|---|
| Tool layer | `src/lusc_atlas/agent/tools.py` | 3 generic tools (`list_tables`, `describe_table`, guarded `run_sql`) + 8 curated semantic tools (`cohort_summary`, `describe_contrasts`, `top_de_genes`, `gene_de_lookup`, `gene_expression`, `composition_change`, `ccc_top_interactions`, `caveats_and_decisions`). Poka-yoke: enums for contrasts/conditions/stages, errors list the valid values (e.g. tested DE groups) |
| Providers | `src/lusc_atlas/agent/llm.py` | `AnthropicSession` (official SDK, `claude-opus-5-5`, adaptive thinking, explicit effort, strict tools, prompt caching of the tools+system prefix, server-side refusal fallback enabled), `OpenAICompatSession` (build.nvidia.com / NIM / any OpenAI-compatible endpoint; Nemotron reasoning switch per generation), `FakeSession` (tests). Keys from env (or pasted per browser session in the app), never logged |
| Loop | `src/lusc_atlas/agent/loop.py` | Manual tool-use loop (needed to drive two provider APIs identically and record provenance per call); all tool results returned in one message; step cap; `check_grounding` |
| MCP server | `src/lusc_atlas/agent/mcp_server.py` | Same 11 tools over stdio, annotated read-only/idempotent/closed-world; `pixi run -e dashboard mcp` |
| UI | `app/views/ask.py` | Answer, grounding badge, every query with its SQL/rows/caveats; disabled with instructions when no key; "run a tool yourself" works with no LLM at all |
| CLI | `python -m lusc_atlas.agent "…"` (`pixi run -e dashboard ask`) | Prints answer + provenance; `--json` for the full record |

**Why `src/lusc_atlas/agent/` and not `app/agent/`:** four consumers (Streamlit page, MCP server, eval runner, CLI)
import it, and it depends only on duckdb/pandas plus optional SDKs, never on the analysis stack. The Space bundle that
vendored it (`scripts/build_space.py`) was parked on 2026-11-12 under the zero-cost rule.

**Why curated tools *and* SQL:** the curated tools encode the decisions a reader should not have to rediscover
(which contrast is primary, that scCODA effects are relative to `Unknown`, that padj not p is the cut). SQL is
the escape hatch for questions nobody anticipated, kept safe by the guards above.

### 5.3 What makes it stand out
1. **Caveat-aware tools.** Most agents return rows; these return rows plus the open issues that change how those
   rows should be read. The model cannot quote the scCODA result without being handed issues 1 and 2.
2. **Mechanical number grounding.** Not "please don't hallucinate", but a checker that matches each number in the
   answer against this run's tool results and shows any mismatch to the reader (and fails the eval item).
3. **One source of truth, two front doors.** The dashboard and the agent read the same DuckDB, under a schema
   contract test, so an answer and a chart cannot disagree.
4. **Same tools over MCP.** Any MCP client gets identical guardrails; there is no second, laxer path.
5. **Provider-neutral eval.** The same 24 questions, graded deterministically against SQL-computed truths, compare
   Claude and Nemotron on accuracy, grounding, caveat use and tokens (Concept B).

## 6. Evaluation plan

* **Set:** `tests/agent_eval/questions.yaml`: 24 questions across lookup (5), aggregate (4), DE (4), composition (2),
  CCC (2), spatial (2), caveat (1), negative (2), adversarial (2). Ground truth is a `truth_sql` executed on the same
  snapshot at run time, so the eval moves with the data when the pipeline is rerun.
* **Grading (deterministic, no LLM judge):** value (number within tolerance / required entities present), required
  mentions (caveats, refusals, padj), forbidden statements (e.g. "160 donors" for the post-QC atlas), grounding (no
  untraceable numbers), and at least one tool call. Pass = all five.
* **Runner:** `pixi run -e dashboard agent-eval --check` validates every truth query and self-tests the grader with
  a synthetic right and wrong answer (no LLM; 24/24 OK). With keys:
  `agent-eval --providers anthropic,nvidia` writes per-question JSON to `logs/agent_eval/` and prints a comparison
  table (pass rate, value accuracy, grounding, caveat rate, calls/question, tokens).
* **Not yet run live:** no API keys exist on this machine (by design). Expected cost of one full run is in
  [B_nemotron_assessment.md §5](B_nemotron_assessment.md#5-cost-and-hardware).
* **Next:** add a rubric-based LLM judge (position-randomised, never the model under test) for interpretation
  quality, and track pass rate in CI as a non-blocking job once a key is available as a CI secret.

## 7. Limitations (honest list)
* The grounding check is lexical: it proves a number *appears* in a tool result, not that it is the *right* number
  for the claim. The value grader in the eval covers correctness for the 24 questions; the reader still has to look.
* Small integers (≤ 10) are not checked (too many false alarms from "top 5", "3 niches").
* The agent inherits every §3.10 issue in the data. It is told about them; it cannot fix them.
* No live provider run exists yet, so there is no measured accuracy. Treat the prototype as tested plumbing with an
  untested brain until `agent-eval` has run.
