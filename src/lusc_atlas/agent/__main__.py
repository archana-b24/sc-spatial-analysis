"""CLI: ask the atlas one question from the terminal.

pixi run -e dashboard ask "Which genes go up most in tumour macrophages?"
pixi run -e dashboard ask --provider nvidia "..."      # Nemotron via build.nvidia.com / NIM
pixi run -e dashboard ask --json "..."                 # full provenance as JSON
"""

from __future__ import annotations

import argparse
import json
import sys

from .llm import ProviderUnavailable
from .loop import ask, format_provenance_markdown
from .tools import AtlasTools


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m lusc_atlas.agent", description=__doc__)
    ap.add_argument("question")
    ap.add_argument("--provider", default=None, help="anthropic (default) | nvidia")
    ap.add_argument("--model", default=None)
    ap.add_argument("--db", default=None, help="path to lusc_atlas.duckdb (default: $LUSC_DB or repo build)")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    try:
        res = ask(a.question, AtlasTools(a.db), provider=a.provider, model=a.model)
    except ProviderUnavailable as e:
        print(f"Agent unavailable: {e}", file=sys.stderr)
        return 2
    if res.error:
        print(f"Agent unavailable: {res.error}", file=sys.stderr)
        return 2
    if a.json:
        print(json.dumps(res.to_dict(), indent=2, default=str))
    else:
        print(res.answer, "\n\n---\n" + format_provenance_markdown(res))
    return 0


if __name__ == "__main__":
    sys.exit(main())
