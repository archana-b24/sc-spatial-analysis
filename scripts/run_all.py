#!/usr/bin/env python
"""Run the whole pipeline sequentially — the zero-learning-curve alternative to Snakemake.

    pixi run run-all                      # everything, skipping steps whose outputs already exist
    pixi run python scripts/run_all.py --from integrate --force
    pixi run python scripts/run_all.py --only qc
    pixi run python scripts/run_all.py --dry-run

Every step runs as `pixi run -e <env> lusc-atlas <step>` in a fresh process (exactly what Snakemake
does), so memory is released between steps and the cell2location step gets its own environment.
A step is skipped only if its declared outputs (`lusc_atlas.steps.OUTPUTS`, shared with the Snakefile)
exist AND its completion stamp `logs/run_all/.done/<step>` exists — the stamp is written only after a
clean exit, so interrupted or failed steps always re-run on resume (DEVELOPMENT_LOG §2.6).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from lusc_atlas.steps import ENV, INPUTS, ORDER, OUTPUTS  # noqa: E402

STAMPS = REPO / "logs" / "run_all" / ".done"


def stamp(step: str) -> Path:
    return STAMPS / step


def done(step: str) -> bool:
    """A step is complete only if its declared outputs exist AND it finished cleanly (stamp written
    after exit 0). Outputs alone are not trusted: a step killed mid-way, or outputs left over from an
    older run, would otherwise be skipped on resume."""
    return stamp(step).exists() and all((REPO / o).exists() for o in OUTPUTS[step])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--only", choices=ORDER, nargs="+", help="run only these steps")
    g.add_argument("--from", dest="start", choices=ORDER, help="start at this step")
    ap.add_argument("--force", action="store_true", help="re-run even if outputs exist")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and exit")
    ap.add_argument(
        "--keep-going",
        action="store_true",
        help="on failure, continue with steps whose inputs exist (unattended runs)",
    )
    ap.add_argument("--config", default=None)
    ap.add_argument(
        "--adopt",
        choices=ORDER,
        nargs="+",
        metavar="STEP",
        help="mark steps whose outputs exist as done without re-running "
        "(for runs that pre-date completion stamps; check their logs first)",
    )
    a = ap.parse_args()

    if a.adopt:
        STAMPS.mkdir(parents=True, exist_ok=True)
        for s in a.adopt:
            if not all((REPO / o).exists() for o in OUTPUTS[s]):
                print(f"✗ cannot adopt {s}: outputs missing", file=sys.stderr)
                return 1
            stamp(s).write_text("adopted (completed before stamps existed)\n")
            print(f"  adopted {s}")
        return 0

    steps = a.only or (ORDER[ORDER.index(a.start) :] if a.start else ORDER)
    plan = [(s, "skip (outputs exist)" if (done(s) and not a.force) else "run") for s in steps]
    width = max(len(s) for s in steps)
    for s, what in plan:
        print(f"  {s:<{width}}  [{ENV.get(s, 'default'):>7}]  {what}")
    if a.dry_run:
        return 0

    log_dir = REPO / "logs" / "run_all"
    log_dir.mkdir(parents=True, exist_ok=True)
    failed: list[str] = []
    for s, what in plan:
        if what != "run":
            continue
        missing = [i for i in INPUTS[s] if not (REPO / i).exists()]
        if missing:
            print(f"✗ {s}: missing inputs {missing}", file=sys.stderr, flush=True)
            if a.keep_going:
                failed.append(s)
                continue
            return 1
        cmd = [
            "pixi",
            "run",
            "--manifest-path",
            str(REPO / "pixi.toml"),
            "-e",
            ENV.get(s, "default"),
            "lusc-atlas",
            s,
        ]
        if a.config:
            cmd += ["--config", a.config]
        t0 = time.time()
        stamp(s).unlink(missing_ok=True)  # invalidate before touching outputs
        print(f"▶ {s}  (log: logs/run_all/{s}.log)", flush=True)
        with open(log_dir / f"{s}.log", "w") as fh:
            rc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            print(
                f"✗ {s} failed (exit {rc}) after {time.time() - t0:.0f}s — see logs/run_all/{s}.log",
                file=sys.stderr,
                flush=True,
            )
            if a.keep_going:
                failed.append(s)
                continue
            return rc
        STAMPS.mkdir(parents=True, exist_ok=True)
        stamp(s).write_text(f"{time.strftime('%Y-%m-%dT%H:%M:%S')}\t{time.time() - t0:.0f}s\n")
        print(f"✓ {s}  {time.time() - t0:.0f}s", flush=True)
    print("ALL DONE" if not failed else f"FINISHED WITH FAILURES: {failed}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
