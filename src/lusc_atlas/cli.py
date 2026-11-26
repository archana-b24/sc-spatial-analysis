"""Command-line entry point: `lusc-atlas <step> [--config path]`.

Both runners call this:
* Snakemake  → one rule per step, `pixi run -e <env> lusc-atlas <step>`
* run_all.py → loops over steps in order, same command.
"""

from __future__ import annotations

import argparse
import sys

from .config import load_config, setup_logging
from .steps import ORDER, run_step


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="lusc-atlas", description=__doc__)
    ap.add_argument("step", choices=ORDER + ["list"], help="pipeline step to run ('list' prints the order)")
    ap.add_argument("--config", default=None, help="path to config.yaml (default: config/config.yaml)")
    args = ap.parse_args(argv)
    if args.step == "list":
        print("\n".join(ORDER))
        return 0
    setup_logging()
    run_step(load_config(args.config), args.step)
    return 0


if __name__ == "__main__":
    sys.exit(main())
