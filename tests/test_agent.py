"""Agent: tool guardrails, the tool-use loop (fake LLM), provider adapters (mocked SDK clients, never
the network), grounding check, eval-set integrity and the MCP server.

Needs the `dashboard` pixi env and a built DB; skipped otherwise (the default-env CI job is unaffected).
"""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("duckdb")
DB = Path(__file__).resolve().parents[1] / "results" / "dashboard" / "lusc_atlas.duckdb"
pytestmark = pytest.mark.skipif(not DB.exists(), reason="dashboard DB not built")

from lusc_atlas.agent import llm  # noqa: E402
from lusc_atlas.agent.llm import FakeSession, ProviderUnavailable, ToolCall, Turn  # noqa: E402
from lusc_atlas.agent.loop import SYSTEM_PROMPT, ask, check_grounding  # noqa: E402
from lusc_atlas.agent.tools import TOOL_NAMES, TOOL_SPECS, AtlasTools, ToolResult  # noqa: E402


@pytest.fixture(scope="module")
def tools():
    return AtlasTools(DB)


@pytest.fixture(autouse=True)
def no_keys(monkeypatch):
    for k in ("ANTHROPIC_API_KEY", "NVIDIA_API_KEY", "LUSC_AGENT_PROVIDER", "LUSC_AGENT_MODEL"):
        monkeypatch.delenv(k, raising=False)


# ------------------------------------------------------------------------------------------ tools
def test_specs_match_methods(tools):
    for s in TOOL_SPECS:
        assert hasattr(tools, f"t_{s['name']}"), s["name"]
        assert s["input_schema"]["additionalProperties"] is False


@pytest.mark.parametrize(
    "sql, msg",
    [
        ("DROP TABLE cells", "only SELECT"),
        ("DELETE FROM known_issues", "only SELECT"),
        ("SELECT 1; SELECT 2", "exactly one"),
        ("CREATE TABLE x AS SELECT 1", "only SELECT"),
        ("ATTACH 'x.db'", "only SELECT"),
    ],
)
def test_run_sql_rejects_non_select(tools, sql, msg):
    r = tools.call("run_sql", {"sql": sql})
    assert r.error and msg in r.error


def test_run_sql_blocks_file_access(tools):
    r = tools.call("run_sql", {"sql": "SELECT * FROM read_csv('/etc/passwd')"})
    assert r.error and "disabled" in r.error.lower()


def test_run_sql_row_cap_and_provenance(tools):
    r = tools.call("run_sql", {"sql": "SELECT cell_id FROM cells", "max_rows": 5})
    assert r.error is None and r.n_rows == 5 and r.truncated
    assert r.tables == ["cells"] and "LIMIT 6" in r.sql[0]


def test_unknown_tool_and_bad_args(tools):
    assert "unknown tool" in tools.call("rm_rf", {}).error
    assert "bad arguments" in tools.call("top_de_genes", {"gene": "X"}).error


def test_cohort_summary_matches_db(tools):
    r = tools.call("cohort_summary")
    row = r.rows[0]
    n = tools._query("SELECT count(*) c, count(DISTINCT donor_id) d FROM cells").iloc[0]
    assert row["cells"] == n.c and row["donors"] == n.d
    pre = tools._query("SELECT count(DISTINCT donor_id) d FROM cohort_donors").iloc[0]
    assert f"{int(pre.d)} donors" in r.note  # the pre-QC vs post-QC distinction travels with the number


def test_composition_carries_caveats(tools):
    r = tools.call("composition_change", {"condition": "tumour"})
    assert {c["issue"] for c in r.caveats} >= {13, 18}
    assert all("consistent" in row and "glm_ci_low" in row for row in r.rows)


def test_top_de_rejects_untested_group(tools):
    r = tools.call("top_de_genes", {"contrast": "paired_tumour_vs_adjacent", "cell_type": "Basal"})
    assert r.error and "Tested groups" in r.error


def test_model_json_is_bounded(tools):
    r = tools.call("run_sql", {"sql": "SELECT * FROM cells", "max_rows": 200})
    assert len(r.to_model_json(max_chars=5000)) <= 6000


# ------------------------------------------------------------------------------------- grounding
def test_grounding_accepts_rounding_and_percent():
    r = ToolResult(tool="t", args={}, rows=[{"n": 231504, "frac": 0.13379, "lfc": 4.454334}])
    g = check_grounding("231,504 cells; 13.4% pass; log2FC 4.45 in CD8 T and chr3q [cells].", [r])
    assert g["ok"], g


def test_grounding_flags_invented_numbers():
    r = ToolResult(tool="t", args={}, rows=[{"n": 231504}])
    g = check_grounding("About 240,000 cells, 47% of them T cells.", [r])
    assert not g["ok"] and set(g["ungrounded"]) == {"240,000", "47%"}


# ------------------------------------------------------------------------------------------ loop
def _scripted(question_calls):
    """A fake model that calls the given tools, then answers by copying a number from the result."""

    def script(history):
        n_tool = sum(1 for h in history if h[0] == "tool")
        if n_tool < len(question_calls):
            name, args = question_calls[n_tool]
            return Turn(tool_calls=[ToolCall(f"c{n_tool}", name, args)], stop_reason="tool_use")
        out = json.loads(history[-1][2])
        return Turn(text=f"There are {out['rows'][0]['cells']:,} cells [cells].", stop_reason="end_turn")

    return script


def test_loop_with_fake_session(tools):
    s = FakeSession("fake", SYSTEM_PROMPT, TOOL_SPECS, script=_scripted([("cohort_summary", {})]))
    res = ask("How many cells?", tools, session=s)
    n = int(tools._query("SELECT count(*) c FROM cells").c.iloc[0])
    assert res.answer == f"There are {n:,} cells [cells]."
    assert res.steps == 1 and [c.tool for c in res.calls] == ["cohort_summary"]
    assert "cells" in res.tables and res.sql and res.grounding["ok"]


def test_loop_reports_tool_errors_to_model(tools):
    seen = {}

    def script(history):
        if len(history) == 1:
            return Turn(tool_calls=[ToolCall("x", "run_sql", {"sql": "DROP TABLE cells"})])
        seen["err"] = history[-1][3]
        return Turn(text="I cannot modify the database; it is read-only.")

    res = ask("drop it", tools, session=FakeSession("f", "", [], script=script))
    assert seen["err"] is True and res.calls[0].error
    assert tools.call("run_sql", {"sql": "SELECT count(*) n FROM cells"}).rows[0]["n"] > 0


def test_loop_step_cap(tools):
    forever = lambda h: Turn(tool_calls=[ToolCall(str(len(h)), "list_tables", {})])  # noqa: E731
    res = ask("loop", tools, session=FakeSession("f", "", [], script=forever), max_steps=3)
    assert res.steps == 3 and res.stop_reason == "step_limit"


def test_providers_unavailable_without_keys(tools):
    assert all(v for v in llm.available_providers().values())
    for p in ("anthropic", "nvidia"):
        with pytest.raises(ProviderUnavailable, match="_API_KEY"):
            ask("hi", tools, provider=p)
    with pytest.raises(ProviderUnavailable):
        llm.make_session("nope", "", [])


# ---------------------------------------------------------------------- provider adapters (mocked)
def test_anthropic_adapter_message_shapes(monkeypatch):
    pytest.importorskip("anthropic")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-not-a-real-key")
    calls = []
    replies = [
        SimpleNamespace(
            content=[SimpleNamespace(type="tool_use", id="tu1", name="cohort_summary", input={})],
            stop_reason="tool_use",
            usage=SimpleNamespace(input_tokens=10, output_tokens=5),
        ),
        SimpleNamespace(
            content=[SimpleNamespace(type="text", text="done")],
            stop_reason="end_turn",
            usage=SimpleNamespace(input_tokens=20, output_tokens=3),
        ),
    ]
    s = llm.make_session("anthropic", "sys", TOOL_SPECS)

    def fake_create(**kw):
        calls.append(json.loads(json.dumps(kw, default=lambda o: "<obj>")))
        return replies[len(calls) - 1]

    monkeypatch.setattr(s.client.beta.messages, "create", fake_create)
    t1 = s.send_user("q")
    assert t1.tool_calls[0].name == "cohort_summary"
    t2 = s.send_tool_results([(t1.tool_calls[0], '{"ok":1}', False)])
    assert t2.text == "done"
    kw = calls[0]
    assert kw["model"] == "claude-opus-5-5" and kw["thinking"] == {"type": "adaptive"}
    assert all(t["strict"] for t in kw["tools"]) and kw["fallbacks"] == "default"
    last = calls[1]["messages"][-1]
    assert last["role"] == "user" and last["content"][0]["type"] == "tool_result"
    assert last["content"][0]["tool_use_id"] == "tu1"


def test_openai_compat_adapter_message_shapes(monkeypatch):
    pytest.importorskip("openai")
    monkeypatch.setenv("NVIDIA_API_KEY", "test-not-a-real-key")
    calls = []

    def msg(content=None, tool_calls=None):
        m = SimpleNamespace(content=content, tool_calls=tool_calls)
        m.model_dump = lambda exclude_none=True: {"role": "assistant", "content": content or ""}
        return m

    tc = SimpleNamespace(id="call1", function=SimpleNamespace(name="list_tables", arguments="{}"))
    replies = [
        SimpleNamespace(
            choices=[SimpleNamespace(message=msg(tool_calls=[tc]), finish_reason="tool_calls")], usage=None
        ),
        SimpleNamespace(
            choices=[SimpleNamespace(message=msg("<think>x</think>answer"), finish_reason="stop")], usage=None
        ),
    ]
    s = llm.make_session("nvidia", "sys", TOOL_SPECS, reasoning="off")

    def fake_create(**kw):
        calls.append({**kw, "messages": list(kw["messages"])})  # snapshot: the session appends afterwards
        return replies[len(calls) - 1]

    monkeypatch.setattr(s.client.chat.completions, "create", fake_create)
    t1 = s.send_user("q")
    t2 = s.send_tool_results([(t1.tool_calls[0], "{}", False)])
    assert t2.text == "answer"
    assert calls[0]["model"] == "nvidia/nemotron-3-super-120b-a12b"
    assert calls[0]["extra_body"] == {"chat_template_kwargs": {"enable_thinking": False}}
    assert calls[0]["tools"][0]["type"] == "function"
    assert calls[1]["messages"][-1] == {"role": "tool", "tool_call_id": "call1", "content": "{}"}


def test_reasoning_controls_per_generation():
    rc = llm.reasoning_controls
    assert rc("nvidia/nemotron-3-nano-30b-a3b", "on") == (
        "",
        {"chat_template_kwargs": {"enable_thinking": True}},
    )
    assert rc("nvidia/llama-3.3-nemotron-super-49b-v1.5", "off") == ("/no_think", {})
    assert rc("nvidia/llama-3.1-nemotron-ultra-253b-v1", "on") == ("detailed thinking on", {})
    assert rc("anything", None) == ("", {})


# ---------------------------------------------------------------------------------- eval + MCP
def test_eval_set_self_check(tools):
    from lusc_atlas.agent.evaluation import compute_truths, grade_value, load_items

    items = compute_truths(load_items(), tools)
    assert 15 <= len(items) <= 30 and len({i.id for i in items}) == len(items)
    for it in items:
        if it.grade["type"] == "number":
            assert grade_value(it, f"value {float(it.truth):.6g}")
        elif it.grade["type"] == "set":
            assert grade_value(it, " ".join(it.truth)) and not grade_value(it, "nothing")


def test_mcp_server_lists_read_only_tools(tools):
    pytest.importorskip("mcp")
    import asyncio

    from lusc_atlas.agent.mcp_server import build_server

    server = build_server(tools)
    listed = asyncio.run(server.list_tools())
    assert sorted(t.name for t in listed) == sorted(TOOL_NAMES)
    assert all(t.annotations.read_only_hint for t in listed)
    out = asyncio.run(server.call_tool("run_sql", {"sql": "DELETE FROM cells"}))
    assert "only SELECT" in out.content[0].text
