"""Grade agent answers against tests/agent_eval/questions.yaml.

Ground truth comes from `truth_sql` run on the live snapshot; grading is deterministic (no LLM judge), so
two providers are compared on exactly the same yardstick and a run costs only the agent's own tokens.
A question passes only if ALL of these hold:

* value      - the expected number / strings are in the answer (``grade``);
* mentions   - every ``must_mention`` group is satisfied (caveats, refusals, statistical honesty);
* forbidden  - no ``must_not_mention`` string appears;
* grounded   - every number in the answer is traceable to a tool result (loop.check_grounding);
* used_tools - at least one tool call (an answer from memory fails even if it is right).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .loop import _NUM, AgentResult, _parse

EVAL_FILE = Path(__file__).resolve().parents[3] / "tests" / "agent_eval" / "questions.yaml"


@dataclass
class Item:
    id: str
    category: str
    question: str
    grade: dict
    truth_sql: str | None = None
    must_mention: list[list[str]] = field(default_factory=list)
    must_not_mention: list[str] = field(default_factory=list)
    truth: Any = None


def load_items(path: Path = EVAL_FILE) -> list[Item]:
    raw = yaml.safe_load(Path(path).read_text())["questions"]
    return [Item(**q) for q in raw]


def compute_truths(items: list[Item], tools) -> list[Item]:
    """Fill ``item.truth`` from the DB (via the tools' own read-only connection)."""
    for it in items:
        if not it.truth_sql:
            continue
        df = tools._query(it.truth_sql)
        if df.empty:
            raise ValueError(f"{it.id}: truth_sql returned no rows")
        it.truth = df.iloc[0, 0] if it.grade["type"] == "number" else [str(v) for v in df.iloc[:, 0]]
    return items


def _numbers(text: str) -> list[tuple[float, bool]]:
    out = []
    for m in _NUM.finditer(text):
        p = _parse(m.group(0))
        if p:
            out.append((p[0], p[2]))
    return out


def grade_value(it: Item, answer: str) -> bool:
    t = it.grade["type"]
    a = answer.lower()
    if t == "text":
        return True
    if t == "set":
        hits = sum(1 for s in it.truth if re.search(rf"(?<![a-z0-9]){re.escape(s.lower())}(?![a-z0-9])", a))
        return hits >= min(int(it.grade.get("min_hits", len(it.truth))), len(it.truth))
    if t == "number":
        truth = float(it.truth)
        tol = max(float(it.grade.get("tol_abs", 0.0)), abs(truth) * float(it.grade.get("tol_rel", 0.0)))
        for val, pct in _numbers(answer):
            cands = [val]
            if it.grade.get("as_percent") or pct:
                cands.append(val / 100)
            # a number written with fewer decimals than tol allows still counts if it rounds to the truth
            if any(abs(c - truth) <= tol + 1e-9 for c in cands):
                return True
            if tol == 0 and any(abs(c - truth) < 0.5 for c in cands):  # integer counts written as-is
                return True
        return False
    raise ValueError(f"unknown grade type {t}")


def grade(it: Item, res: AgentResult) -> dict:
    a = res.answer.lower()
    mentions = all(any(alt.lower() in a for alt in group) for group in it.must_mention)
    forbidden = not any(s.lower() in a for s in it.must_not_mention)
    value = grade_value(it, res.answer)
    used = len(res.calls) > 0
    grounded = bool(res.grounding.get("ok", False))
    return {
        "id": it.id,
        "category": it.category,
        "pass": value and mentions and forbidden and used and grounded and not res.error,
        "value": value,
        "mentions": mentions,
        "forbidden_ok": forbidden,
        "used_tools": used,
        "grounded": grounded,
        "ungrounded": res.grounding.get("ungrounded", []),
        "n_calls": len(res.calls),
        "n_tool_errors": sum(1 for c in res.calls if c.error),
        "steps": res.steps,
        "usage": res.usage,
        "duration_s": round(res.duration_s, 2),
        "error": res.error,
        "truth": it.truth if not hasattr(it.truth, "item") else it.truth.item(),
    }


def summarise(rows: list[dict]) -> dict:
    n = len(rows) or 1
    s = {
        "n": len(rows),
        "pass_rate": sum(r["pass"] for r in rows) / n,
        "value_acc": sum(r["value"] for r in rows) / n,
        "grounded_rate": sum(r["grounded"] for r in rows) / n,
        "caveat_rate": sum(r["mentions"] for r in rows) / n,
        "mean_calls": sum(r["n_calls"] for r in rows) / n,
        "input_tokens": sum(r["usage"].get("input_tokens", 0) for r in rows),
        "output_tokens": sum(r["usage"].get("output_tokens", 0) for r in rows),
        "by_category": {},
    }
    for r in rows:
        c = s["by_category"].setdefault(r["category"], [0, 0])
        c[0] += r["pass"]
        c[1] += 1
    return s
