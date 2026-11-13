"""Provider-agnostic chat sessions for the tool-use loop.

A ``ChatSession`` keeps the provider's NATIVE message history (so provider-specific blocks such as Claude's
thinking blocks are passed back unchanged) and exposes two calls the loop needs:

    turn = session.send_user(text)
    turn = session.send_tool_results([(call, json_text, is_error), ...])

Each returns a ``Turn``: final text and/or requested tool calls, the stop reason and token usage.

Providers
---------
* ``anthropic`` (default) - Claude via the official Anthropic SDK. Key: ``ANTHROPIC_API_KEY``.
* ``nvidia``    - any OpenAI-compatible endpoint; defaults to NVIDIA's hosted API (build.nvidia.com) for
                  Nemotron. Key: ``NVIDIA_API_KEY``; ``NVIDIA_BASE_URL`` points it at a self-hosted NIM.
* ``fake``      - scripted turns for tests and offline demos; never touches the network.

SDKs are imported lazily: the dashboard, MCP server and tests run without them. Keys are read from the
environment only and are never logged, echoed or stored.
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

# Defaults are overridable per call or via LUSC_AGENT_MODEL / LUSC_AGENT_PROVIDER.
DEFAULT_MODELS = {
    "anthropic": "claude-opus-5-5",
    # Nemotron 3 Super (120B-A12B hybrid Mamba-Transformer MoE, tool calling) on build.nvidia.com;
    # see docs/design/B_nemotron_assessment.md for the alternatives and why this one
    "nvidia": "nvidia/nemotron-3-super-120b-a12b",
    "fake": "fake-scripted",
}
KEY_ENV = {"anthropic": "ANTHROPIC_API_KEY", "nvidia": "NVIDIA_API_KEY"}
NVIDIA_BASE_URL = "https://integrate.api.nvidia.com/v1"


class ProviderUnavailable(RuntimeError):
    """Raised when a provider cannot be used (no key / SDK missing). Message is user-facing."""


@dataclass
class ToolCall:
    id: str
    name: str
    args: dict


@dataclass
class Turn:
    text: str = ""
    tool_calls: list[ToolCall] = field(default_factory=list)
    stop_reason: str = ""
    usage: dict = field(default_factory=dict)


class ChatSession:
    provider = "base"

    def __init__(self, model: str, system: str, tools: list[dict]):
        self.model, self.system, self.tools = model, system, tools

    def send_user(self, text: str) -> Turn:  # pragma: no cover - interface
        raise NotImplementedError

    def send_tool_results(self, results: list[tuple[ToolCall, str, bool]]) -> Turn:  # pragma: no cover
        raise NotImplementedError


def available_providers() -> dict[str, str | None]:
    """provider → None if usable, else the reason it is not (shown in the UI). Never reveals key values."""
    out: dict[str, str | None] = {}
    for p, env in KEY_ENV.items():
        out[p] = None if os.environ.get(env) else f"set {env} to enable"
    return out


# ------------------------------------------------------------------------------------------ Claude
class AnthropicSession(ChatSession):
    """Manual tool-use loop on the Messages API (not the beta tool runner) because the same loop must
    drive non-Anthropic providers and record provenance per call. Adaptive thinking is on; effort is set
    explicitly. Server-side refusal fallbacks are enabled (``fallbacks="default"``)."""

    provider = "anthropic"

    def __init__(
        self,
        model: str,
        system: str,
        tools: list[dict],
        effort: str = "high",
        max_tokens: int = 16000,
        api_key: str | None = None,
    ):
        super().__init__(model, system, tools)
        if not (api_key or os.environ.get("ANTHROPIC_API_KEY")):
            raise ProviderUnavailable("Claude is not configured: set ANTHROPIC_API_KEY to enable.")
        try:
            import anthropic
        except ImportError as e:  # pragma: no cover - env-dependent
            raise ProviderUnavailable("the `anthropic` package is not installed in this environment") from e
        self._anthropic = anthropic
        # an explicit key (pasted by a visitor, kept only in their session) wins; else the SDK's env lookup
        self.client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()
        self.effort, self.max_tokens = effort, max_tokens
        self.messages: list[dict] = []
        # stable prefix (tools → system) so prompt caching hits across turns and questions
        self._tools = [dict(t, strict=True) for t in tools]

    def _call(self) -> Turn:
        a = self._anthropic
        try:
            resp = self.client.beta.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                system=[{"type": "text", "text": self.system, "cache_control": {"type": "ephemeral"}}],
                tools=self._tools,
                messages=self.messages,
                thinking={"type": "adaptive"},
                output_config={"effort": self.effort},
                betas=["server-side-fallback-2026-07-01"],
                fallbacks="default",
            )
        except a.AuthenticationError as e:
            raise ProviderUnavailable("Claude rejected the API key (ANTHROPIC_API_KEY).") from e
        except a.RateLimitError as e:
            raise ProviderUnavailable("Claude rate limit reached; try again shortly.") from e
        except a.APIConnectionError as e:
            raise ProviderUnavailable("could not reach the Claude API (network).") from e
        # keep the full content (thinking + tool_use blocks) in history, unchanged
        self.messages.append({"role": "assistant", "content": resp.content})
        text = "".join(b.text for b in resp.content if b.type == "text")
        calls = [ToolCall(b.id, b.name, dict(b.input)) for b in resp.content if b.type == "tool_use"]
        if resp.stop_reason == "refusal":
            text = text or "The model declined to answer this question."
        if resp.stop_reason == "max_tokens" and not calls:
            text = (text + "\n\n[answer cut off at the output token limit]").strip()
        u = resp.usage
        usage = {"input_tokens": u.input_tokens, "output_tokens": u.output_tokens}
        for k in ("cache_read_input_tokens", "cache_creation_input_tokens"):
            if getattr(u, k, None):
                usage[k] = getattr(u, k)
        return Turn(text=text, tool_calls=calls, stop_reason=resp.stop_reason or "", usage=usage)

    def send_user(self, text: str) -> Turn:
        self.messages.append({"role": "user", "content": text})
        return self._call()

    def send_tool_results(self, results):
        blocks = [
            {
                "type": "tool_result",
                "tool_use_id": c.id,
                "content": out,
                **({"is_error": True} if err else {}),
            }
            for c, out, err in results
        ]
        self.messages.append({"role": "user", "content": blocks})  # all results in ONE user message
        return self._call()


# ------------------------------------------------------------------------- OpenAI-compatible (NIM)
class OpenAICompatSession(ChatSession):
    """Any OpenAI-compatible chat-completions endpoint with function calling. Default target is NVIDIA's
    hosted API for Nemotron; set NVIDIA_BASE_URL for a self-hosted NIM container (same API).

    Nemotron generations switch reasoning differently (verified on the model cards, 2026-11-12):
    Nemotron 3 / 3.5 → ``chat_template_kwargs.enable_thinking``; Super-49B v1.5 and Nano-9B-v2 →
    ``/think`` / ``/no_think`` in the system prompt; Llama-Nemotron v1 → ``detailed thinking on|off``.
    ``reasoning`` ("on" / "off" / None = model default) is translated per model."""

    provider = "nvidia"

    def __init__(
        self,
        model: str,
        system: str,
        tools: list[dict],
        reasoning: str | None = None,
        max_tokens: int = 8192,
        base_url: str | None = None,
        api_key_env: str = "NVIDIA_API_KEY",
        api_key: str | None = None,
    ):
        super().__init__(model, system, tools)
        if not (api_key or os.environ.get(api_key_env)):
            raise ProviderUnavailable(f"Nemotron is not configured: set {api_key_env} to enable.")
        try:
            import openai
        except ImportError as e:  # pragma: no cover
            raise ProviderUnavailable("the `openai` package is not installed in this environment") from e
        self._openai = openai
        self.client = openai.OpenAI(
            base_url=base_url or os.environ.get("NVIDIA_BASE_URL", NVIDIA_BASE_URL),
            api_key=api_key or os.environ[api_key_env],
        )
        self.max_tokens = max_tokens
        prefix, self.extra_body = reasoning_controls(model, reasoning)
        self.messages: list[dict] = [{"role": "system", "content": (prefix + "\n" + system).strip()}]
        self._tools = [
            {
                "type": "function",
                "function": {
                    "name": t["name"],
                    "description": t["description"],
                    "parameters": t["input_schema"],
                },
            }
            for t in tools
        ]

    def _call(self) -> Turn:
        o = self._openai
        try:
            resp = self.client.chat.completions.create(
                model=self.model,
                messages=self.messages,
                tools=self._tools,
                tool_choice="auto",
                max_tokens=self.max_tokens,
                temperature=0.0,
                **({"extra_body": self.extra_body} if self.extra_body else {}),
            )
        except o.AuthenticationError as e:
            raise ProviderUnavailable("the endpoint rejected NVIDIA_API_KEY.") from e
        except o.RateLimitError as e:
            raise ProviderUnavailable("rate limit reached on the NVIDIA endpoint; try again shortly.") from e
        except o.APIConnectionError as e:
            raise ProviderUnavailable(
                "could not reach the NVIDIA endpoint (network / NVIDIA_BASE_URL)."
            ) from e
        choice = resp.choices[0]
        msg = choice.message
        calls = []
        for tc in msg.tool_calls or []:
            try:
                args = json.loads(tc.function.arguments or "{}")
            except json.JSONDecodeError:
                args = {"__invalid_json__": tc.function.arguments}
            calls.append(ToolCall(tc.id, tc.function.name, args))
        self.messages.append(msg.model_dump(exclude_none=True))
        text = _strip_think(msg.content or "")
        usage = {}
        if resp.usage:
            usage = {"input_tokens": resp.usage.prompt_tokens, "output_tokens": resp.usage.completion_tokens}
        return Turn(text=text, tool_calls=calls, stop_reason=choice.finish_reason or "", usage=usage)

    def send_user(self, text: str) -> Turn:
        self.messages.append({"role": "user", "content": text})
        return self._call()

    def send_tool_results(self, results):
        for c, out, _err in results:
            self.messages.append({"role": "tool", "tool_call_id": c.id, "content": out})
        return self._call()


def reasoning_controls(model: str, reasoning: str | None) -> tuple[str, dict]:
    """(system-prompt prefix, extra request body) that switch a Nemotron model's reasoning on/off."""
    if reasoning not in ("on", "off"):
        return "", {}
    on = reasoning == "on"
    m = model.lower()
    if "nemotron-3" in m:  # Nemotron 3 / 3.5: chat-template flag
        return "", {"chat_template_kwargs": {"enable_thinking": on}}
    if "super-49b-v1.5" in m or "super-49b-v1_5" in m or "nano-9b-v2" in m or "nano-12b-v2" in m:
        return ("/think" if on else "/no_think"), {}
    if "nemotron" in m:  # Llama-Nemotron v1 family
        return f"detailed thinking {'on' if on else 'off'}", {}
    return "", {}


def _strip_think(text: str) -> str:
    """Reasoning models may inline <think>…</think>; the user sees only the answer."""
    if "</think>" in text:
        text = text.split("</think>", 1)[1]
    return text.strip()


# ------------------------------------------------------------------------------------------- fake
class FakeSession(ChatSession):
    """Deterministic stand-in. ``script`` is either a list of Turns (played in order) or a callable
    ``(history) -> Turn`` that can look at earlier tool results; ``history`` is a list of
    ("user", text) / ("tool", ToolCall, output, is_error) events."""

    provider = "fake"

    def __init__(self, model: str, system: str, tools: list[dict], script: list[Turn] | Callable = ()):
        super().__init__(model, system, tools)
        self.script = script
        self.history: list[tuple] = []
        self._i = 0

    def _next(self) -> Turn:
        if callable(self.script):
            return self.script(self.history)
        if self._i >= len(self.script):
            return Turn(text="(fake script exhausted)", stop_reason="end_turn")
        t = self.script[self._i]
        self._i += 1
        return t

    def send_user(self, text: str) -> Turn:
        self.history.append(("user", text))
        return self._next()

    def send_tool_results(self, results):
        for c, out, err in results:
            self.history.append(("tool", c, out, err))
        return self._next()


def make_session(
    provider: str | None, system: str, tools: list[dict], model: str | None = None, **kw: Any
) -> ChatSession:
    provider = (provider or os.environ.get("LUSC_AGENT_PROVIDER") or "anthropic").lower()
    model = model or os.environ.get("LUSC_AGENT_MODEL") or DEFAULT_MODELS.get(provider)
    if provider == "anthropic":
        return AnthropicSession(model, system, tools, **kw)
    if provider in ("nvidia", "nemotron", "openai_compat"):
        return OpenAICompatSession(model, system, tools, **kw)
    if provider == "fake":
        return FakeSession(model, system, tools, **kw)
    raise ProviderUnavailable(f"unknown provider '{provider}' (anthropic | nvidia | fake)")
