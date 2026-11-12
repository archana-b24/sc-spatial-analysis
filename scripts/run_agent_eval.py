#!/usr/bin/env python
"""Run the "Ask the atlas" eval set against one or more LLM providers.

    pixi run -e dashboard agent-eval --check                       # no LLM: truths compute, grader self-test
    pixi run -e dashboard agent-eval --providers anthropic         # needs ANTHROPIC_API_KEY
    pixi run -e dashboard agent-eval --providers anthropic,nvidia  # Claude vs Nemotron, same questions

Live runs make real, billed API calls (≈ 24 questions × a few tool steps). Providers without a key are
skipped with a message, never faked. Results go to logs/agent_eval/<timestamp>_<provider>.json plus a
markdown comparison table (logs/ is git-ignored; copy a report into docs/ if you want to keep it).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from lusc_atlas.agent.evaluation import EVAL_FILE, compute_truths, grade, grade_value, load_items, summarise
from lusc_atlas.agent.llm import ProviderUnavailable, available_providers
from lusc_atlas.agent.loop import ask
from lusc_atlas.agent.tools import AtlasTools

REPO = Path(__file__).resolve().parents[1]


def self_check(items, tools) -> int:
    """Validate the eval set without any LLM: every truth query runs, and the grader accepts a
    synthetic correct answer and rejects a synthetic wrong one for each value-graded item."""
    compute_truths(items, tools)
    bad = 0
    for it in items:
        t = it.grade["type"]
        if t == "number":
            good, wrong = (
                f"The value is {float(it.truth):.6g}.",
                f"The value is {float(it.truth) * 3 + 7:.6g}.",
            )
        elif t == "set":
            good, wrong = "Answer: " + ", ".join(it.truth), "Answer: none of them."
        else:
            good = wrong = None
        ok = True
        if good is not None:
            ok = grade_value(it, good) and not grade_value(it, wrong)
        print(f"{'OK ' if ok else 'BAD'} {it.id:28s} [{it.category}] truth={it.truth!r}"[:150])
        bad += not ok
    print(f"\n{len(items)} items, {bad} grader self-test failures")
    return 1 if bad else 0


def run(items, tools, provider: str, out_dir: Path) -> dict | None:
    compute_truths(items, tools)
    rows, records = [], []
    for it in items:
        try:
            res = ask(it.question, tools, provider=provider)
        except ProviderUnavailable as e:
            print(f"[{provider}] skipped: {e}")
            return None
        g = grade(it, res)
        rows.append(g)
        records.append({"grade": g, "result": res.to_dict()})
        print(
            f"[{provider}] {'PASS' if g['pass'] else 'fail'} {it.id} "
            f"({g['n_calls']} calls, {g['duration_s']}s)"
        )
    s = summarise(rows)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{stamp}_{provider}.json").write_text(
        json.dumps({"provider": provider, "summary": s, "records": records}, indent=2, default=str)
    )
    return s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="validate the eval set; no LLM calls")
    ap.add_argument("--providers", default="anthropic", help="comma list: anthropic,nvidia")
    ap.add_argument("--only", default=None, help="comma list of question ids")
    ap.add_argument("--db", default=None)
    ap.add_argument("--questions", default=str(EVAL_FILE))
    ap.add_argument("--out", default=str(REPO / "logs" / "agent_eval"))
    a = ap.parse_args(argv)

    tools = AtlasTools(a.db)
    items = load_items(Path(a.questions))
    if a.only:
        keep = set(a.only.split(","))
        items = [i for i in items if i.id in keep]
    if a.check:
        return self_check(items, tools)

    avail = available_providers()
    table = []
    for p in a.providers.split(","):
        p = p.strip()
        if avail.get(p, None if p == "fake" else "unknown provider"):
            print(f"[{p}] skipped: {avail.get(p, 'unknown provider')}")
            continue
        s = run(items, tools, p, Path(a.out))
        if s:
            table.append((p, s))
    if not table:
        print("No provider ran. Set ANTHROPIC_API_KEY and/or NVIDIA_API_KEY (see docs/v3/README.md).")
        return 2
    print("\n| provider | pass | value | grounded | caveats | calls/q | in tok | out tok |")
    print("|---|---|---|---|---|---|---|---|")
    for p, s in table:
        print(
            f"| {p} | {s['pass_rate']:.0%} | {s['value_acc']:.0%} | {s['grounded_rate']:.0%} | "
            f"{s['caveat_rate']:.0%} | {s['mean_calls']:.1f} | {s['input_tokens']:,} | "
            f"{s['output_tokens']:,} |"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
