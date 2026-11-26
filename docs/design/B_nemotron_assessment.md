# B · Does NVIDIA Nemotron make sense here?

*Concept B of [docs/design](README.md). Facts retrieved 2026-11-12 from NVIDIA model cards (huggingface.co/nvidia),
research.nvidia.com, build.nvidia.com / docs.api.nvidia.com, NVIDIA GitHub, arXiv and the Anthropic pricing page.
Benchmarks are **NVIDIA's own** unless a row says otherwise. **UNVERIFIED** = not confirmed on an official
source; **ESTIMATE** = my reasoning, no source.*

## Verdict

**Yes, as a pluggable, open-weights second brain for the atlas agent. No, not as the default today, and not for
anything that has to run on this machine.**

* **Where it fits.** Nemotron is unusually open: weights under commercial-use licences, large parts of the
  pre- and post-training data released, and training recipes on GitHub. For a reproducibility-minded project that is
  the strongest argument: an agent whose model you can pin, self-host and audit, with no data leaving your
  infrastructure. Its reasoning on/off switch also makes a clean experiment: does thinking change grounding?
* **Where it doesn't (yet).**
  1. **No GPU here.** Every Nemotron that is good at multi-step tool use needs datacentre GPUs (Super-class: a B200
     or DGX Spark minimum). The CPU-feasible sizes (4B-class, GGUF) are, by NVIDIA's own numbers, weak at multi-turn
     tool use (Nano-4B: Tau2 29.2, 14.9 at Q4).
  2. **The free hosted API is for prototyping only.** NVIDIA's FAQ restricts build.nvidia.com trial access to
     prototyping/research/testing; production use needs NVIDIA AI Enterprise. A public Space backed by a free NVIDIA
     key is therefore not a legitimate production setup.
  3. **No evidence yet that it matches a frontier model on our task.** There is no independent head-to-head of
     Nemotron 3 against current Claude models on tool use; the one independent index found (Artificial Analysis)
     compared open-weights models only. Our own eval is the right way to decide, and it is ready to run.
* **"Overall" (the phenomenon).** Nemotron is NVIDIA's bet that an open model *plus* its data *plus* its recipes is
  worth more than weights alone, and the hybrid Mamba-Transformer MoE designs (Nemotron-H → Nemotron 3) are a
  serious, published architecture line aimed at long context and cheap inference. For science that is genuinely
  interesting: you can, in principle, reproduce or fine-tune the brain of your analysis agent.

**What was implemented:** the agent's provider switch (`--provider nvidia`, `OpenAICompatSession`), defaulting to
`nvidia/nemotron-3-super-120b-a12b` on `https://integrate.api.nvidia.com/v1`, with the correct reasoning switch per
Nemotron generation, and the same 24-question eval runnable against both providers. **Nothing was called live**
(no keys on this machine, by design; unit tests use mocked clients).

---

## 1. The Nemotron families (what exists, verified)

### Current generation: Nemotron 3 / 3.5 (hybrid Mamba-2 + Transformer + MoE)

| Model | Size (total / active) | Context | Reasoning switch | Tool calling | Licence | Released | Source |
|---|---|---|---|---|---|---|---|
| Nemotron 3 Nano 30B-A3B | ~30B / ~3.5B (research page: 31.6B / 3.2B) | up to 1M (256k default) | `chat_template_kwargs.enable_thinking` | yes (vLLM `qwen3_coder` parser) | NVIDIA Nemotron Open Model License | 2025-12-15 | https://huggingface.co/nvidia/NVIDIA-Nemotron-3-Nano-30B-A3B-BF16 ; tech report arXiv:2512.20848; family paper arXiv:2512.20856 |
| Nemotron 3 Nano 4B | ~4B (`nemotron_h`) | up to 262k | system prompt | trained on tool-calling data; not documented as a feature | Nemotron OML | 2026-03 | https://huggingface.co/nvidia/NVIDIA-Nemotron-3-Nano-4B-GGUF (Q4_K_M GGUF 2.84 GB) |
| **Nemotron 3 Super 120B-A12B** | 120B / 12B | up to 1M | `enable_thinking` (+ low-effort mode) | yes (`qwen3_coder`) | Nemotron OML | 2026-03-11 | https://huggingface.co/nvidia/NVIDIA-Nemotron-3-Super-120B-A12B-NVFP4 ; arXiv:2604.12374 |
| Nemotron 3 Ultra 550B-A55B | 550B / 55B | 1M | `enable_thinking` (+ medium effort) | yes | **OpenMDW-1.1** | 2026-06-04 | https://research.nvidia.com/labs/nemotron/Nemotron-3-Ultra/ |
| Nemotron 3.5 Lightning 30B-A3B | 30B / 3B | up to 1M | `enable_thinking` | yes | OpenMDW-1.1 | 2026-08-11 | https://huggingface.co/nvidia/NVIDIA-Nemotron-3.5-Lightning-30B-A3B-BF16 |
| Nemotron 3.5 Bolt E4B | 8B / 4B "effective"; **post-trained from Google Gemma-4-E4B**, not a Nemotron architecture | 64k | `enable_thinking` | yes (native) | OpenMDW-1.1 + Gemma 4 Apache-2.0 | 2026-10 | https://huggingface.co/nvidia/NVIDIA-Nemotron-3.5-Bolt-E4B-4GB (Q4 GGUF 6.18 GB) |

Also in the 3.5 line: Super VL 120B (multimodal) and Bolt 74B-A9B (aggressively quantised GGUF for 24–32 GB GPUs).

### Earlier generations (still served)

| Model | Notes | Reasoning switch | Source |
|---|---|---|---|
| Llama-3.1-Nemotron-Nano-8B-v1, Nano-4B-v1.1 | Llama-derived; BFCL v2 Live ~64–68 | `detailed thinking on/off` | HF cards; arXiv:2505.00949 |
| Llama-3.3-Nemotron-Super-49B-v1.5 | 49B dense (NAS from Llama-3.3-70B); BFCL v3 71.75; single H100-80GB target | `/think`, `/no_think` | https://huggingface.co/nvidia/Llama-3_3-Nemotron-Super-49B-v1_5 |
| Llama-3.1-Nemotron-Ultra-253B-v1 | 253B dense; BFCL v2 Live 74.1; 8×H100 BF16 | `detailed thinking on/off` | HF card; arXiv:2505.00949 |
| Nemotron-H 8B/47B/56B | First hybrid Mamba-Transformer line; base checkpoints carry a **research-only** licence | — | arXiv:2504.03624 |
| NVIDIA-Nemotron-Nano-9B-v2 / 12B-v2 | Hybrid Mamba-2; thinking budget; BFCL v3 66.9 | `/think`, `/no_think` | arXiv:2508.14444 |
| Nemotron-Mini-4B-Instruct | 2024; RAG / function calling | — | https://huggingface.co/nvidia/Nemotron-Mini-4B-Instruct |

**Licences have changed three times:** NVIDIA Open Model License (+ Llama licence) for Llama-derived models →
NVIDIA Nemotron Open Model License (Nemotron 3 Nano/Super) → OpenMDW-1.1 (Nemotron 3 Ultra and 3.5). All allow
commercial use of the *weights*; read the exact card before relying on one.

## 2. Openness: data, recipes, tooling

| Item | Fact | Source |
|---|---|---|
| Pre-training data | Nemotron-CC (v1, v2, v2.1: +2.5T English tokens), Nemotron-CC-Code, Nemotron-Pretraining-Code/Specialized sets | https://research.nvidia.com/labs/nemotron/Nemotron-3 ; arXiv:2412.02595 |
| Post-training data | Llama-Nemotron / Nemotron-Post-Training-Dataset v1–v3, SFT and RL data. v2 is mostly CC-BY-4.0 and **gated** on HF, and its card warns of inherited Qwen/DeepSeek terms | https://huggingface.co/datasets/nvidia/Nemotron-Post-Training-Dataset-v2 |
| Caveat | The Nano 30B card cites ~25T training tokens, but its listed datasets sum to ~11–13T, so **not everything is released** | HF Nano 30B card |
| Recipes | GitHub NVIDIA-NeMo/Nemotron (training recipes); NeMo / Megatron stack | https://research.nvidia.com/labs/nemotron/Nemotron-3-Ultra/ |
| Agent tooling | NeMo Agent Toolkit (formerly AgentIQ / AIQ), Apache-2.0, v1.9.0 (2026-09-10), MCP client and server support | https://github.com/NVIDIA/NeMo-Agent-Toolkit |

## 3. Access paths

| Path | How | Fit for this project |
|---|---|---|
| **Hosted API (build.nvidia.com)** | `https://integrate.api.nvidia.com/v1`, OpenAI-compatible chat completions, `NVIDIA_API_KEY`; model id e.g. `nvidia/nemotron-3-super-120b-a12b` | **Good for the eval comparison.** Trial terms: "prototyping, research, development and testing purposes only"; production needs NVIDIA AI Enterprise (from $4,500 per GPU per year). The old credit system was replaced by per-model rate limits that NVIDIA does not publish (the often-quoted "40 RPM" is UNVERIFIED). Hosted `tools` support was not tested live (UNVERIFIED per model) |
| **Self-hosted NIM** | GPU containers (CUDA only; no CPU NIM found); OpenAI-compatible; tool calling needs `--enable-auto-tool-choice` + the right parser | Not possible on this host (no GPU). Free for R&D on up to 16 GPUs under the Developer Program |
| **vLLM / SGLang** | Open weights; Nemotron 3 parsers documented | Same GPU requirement |
| **Ollama / llama.cpp (GGUF)** | Official tags `nemotron-3-nano` (4b, 30b), `nemotron-3.5-lightning`, `nemotron-3.5-bolt`; OpenAI-compatible local server | **The only path that runs here (CPU).** Point `NVIDIA_BASE_URL` at the local server (`NVIDIA_API_KEY` can be any non-empty string). Quality caveat below |
| **Third-party hosts (e.g. OpenRouter)** | OpenAI-compatible; Nemotron 3 Super ~$0.075 / $0.40 per MTok in/out (OpenRouter models API, 2026-11-12) | Cheapest legitimate hosted option for a public app; set `NVIDIA_BASE_URL` + the host's model id (ids not verified here) |

## 4. NVIDIA's bio stack, and what applies to a single-cell atlas

| Component | What it is (verified) | Relevance here |
|---|---|---|
| **BioNeMo** | `github.com/NVIDIA/bionemo-framework` now redirects to **NVIDIA-BioNeMo/bionemo-recipes** (Apache-2.0, v3.0.0, 2026-06-24): TransformerEngine-accelerated model classes + training recipes (ESM-2, AMPLIFY, Evo2, CodonFM, **Geneformer**). GPU-only | **Geneformer is in BioNeMo** (checkpoints 10M and 106M, trained on CZ CELLxGENE). This repo already benchmarks zero-shot **Geneformer-V2-104M_CLcancer** as an integration candidate (from the Geneformer authors' HF release, loaded CPU-only through `safe_load.py`), where it scored lowest (scIB total 0.385). BioNeMo's value would be **fine-tuning** Geneformer for LUSC, which needs a GPU. https://docs.nvidia.com/bionemo-recipes/2.5/models/geneformer/ |
| BioNeMo NIMs | Hosted/containerised Boltz-2, OpenFold2/3, msa-search, Evo2-40B, ESM2-650M, AlphaFold2, GenMol, DiffDock | Protein/structure/chemistry. Not used by a transcriptomic atlas |
| BioNeMo Agent Toolkit | Agent "skills" for structure prediction, docking, generative chemistry, genomics, protein design, biomarker discovery; Nemotron for reasoning (announced 2026-06-23) | **No single-cell tools named.** Interesting as an MCP-compatible tool source for a future drug-target extension |
| **rapids-singlecell** | scverse project (MIT, v0.18.1): GPU-accelerated, mostly Scanpy-compatible; parts of Squidpy, decoupler, pertpy | The natural NVIDIA speed-up for this pipeline (QC → integration → neighbours ran ~2–2.5 h per step on CPU), but **requires an NVIDIA GPU** |
| Clara Parabricks | GPU genomics (alignment, variant calling), free container | Not relevant (no reads are processed here) |
| scGPT | No NVIDIA recipe or NIM found (UNVERIFIED / likely none) | — |

**Bottom line for bio:** NVIDIA's single-cell-relevant pieces (Geneformer fine-tuning, rapids-singlecell) are
GPU accelerators for the **pipeline**, not for the dashboard or agent. They become relevant only if the project
gets GPU time, and then they are worth a benchmark row in the integration step (fine-tuned Geneformer vs scANVI).

## 5. Cost and hardware

Token use per question is an **ESTIMATE**: system prompt + tool schemas ≈ 1.7k tokens (6,859 characters measured),
tool results 0.3–2k tokens each (measured 1.1–8.1k characters), ~3 model calls per question → **~15k input + ~2k
output tokens per question** (output includes reasoning). Replace with measured usage after the first `agent-eval` run.

| Option | Where it runs | Hardware | Price basis (retrieved 2026-11-12) | ≈ per question | ≈ 24-question eval |
|---|---|---|---|---|---|
| Claude Opus 5.5 (`claude-opus-5-5`, default) | Anthropic API | none local | $4 / $20 per MTok in/out (cache reads $0.20) | ~$0.10 | ~$2.40 |
| Claude Sonnet 5.5 | Anthropic API | none local | $2 / $10 | ~$0.05 | ~$1.20 |
| Claude Haiku 5.5 | Anthropic API | none local | $0.10 / $0.50 (prompts ≤100k) | ~$0.003 | ~$0.06 |
| Nemotron 3 Super via build.nvidia.com | NVIDIA hosted | none local | free trial, unpublished rate limits; **prototyping only** | $0 | $0 |
| Nemotron 3 Super via OpenRouter (third party) | hosted | none local | ~$0.075 / $0.40 | ~$0.002 | ~$0.05 |
| Nemotron 3 Super, self-hosted NIM/vLLM | your GPU | ≥1× B200 or DGX Spark (NVFP4); tested 1–8× H100/H200 | GPU time + AI Enterprise for production ($4,500 / GPU / year) | GPU-hour cost | — |
| Nemotron 3 Nano 30B-A3B, self-hosted | your GPU | 1× H100/A100-80GB | GPU time | — | — |
| Nemotron 3 Nano 4B / 3.5 Bolt E4B, Ollama on **this CPU host** | local | 16 CPU, 125 GB RAM: fits easily (2.8 / 6.2 GB GGUF) | $0 marginal | ESTIMATE: tens of seconds per turn (no NVIDIA CPU throughput figure exists) | slow; weak tool use (Tau2 29.2 / 55.1 per NVIDIA) |

Prompt caching (on by default in the Claude adapter) cuts the repeated ~1.7k-token prefix to the cache-read rate after
the first call, so real Claude costs should land below these estimates.

## 6. Plan (recommended, mostly reversible)

1. **Done now:** pluggable provider; Nemotron default model + reasoning switch; mocked tests; eval runner that
   compares providers on identical, data-derived ground truth.
2. **Not planned (free-only decision, 2026-11-12):** the paid half of this step is off. Free variant: compare
   against a local GGUF Nemotron via `NVIDIA_BASE_URL` only. Original plan, kept for the record — *when you have keys (~$3 total):* `ANTHROPIC_API_KEY=… NVIDIA_API_KEY=… pixi run -e dashboard agent-eval
   --providers anthropic,nvidia`. Run Nemotron twice, reasoning on and off (`reasoning=` in `make_session`; add a CLI
   flag if useful). Decide on measured pass rate, grounding rate and caveat rate, not on benchmarks.
3. *(Moot while the site has no chatbot.)* **If Nemotron is within a few points of Claude:** offer it as the public site's default through a production-legal
   host (OpenRouter or self-hosted), keep Claude as the quality option. If not: keep it as the documented
   open-weights fallback.
4. **If GPU time appears:** fine-tune Geneformer with BioNeMo recipes as an extra integration-benchmark candidate,
   and try rapids-singlecell for the CPU-bound steps. Both are pipeline changes and belong after the §3.10 fixes.

## 7. Uncertainties
* Hosted per-model `tools` support and rate limits on build.nvidia.com were not tested (no key).
* CPU throughput for GGUF Nemotron on this host is an estimate; nothing was downloaded (multi-GB weights, by brief).
* Licences and model ids change quickly (three licence regimes in 18 months); re-check the card before any deployment.
* Benchmarks above are vendor-reported except Artificial Analysis (Nemotron 3 Ultra 47.7, Super 36.0 on its
  Intelligence Index, 2026-06-04, open-weights comparison only).
