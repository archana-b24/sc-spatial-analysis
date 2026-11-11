"""Import every analysis dependency and print versions — the gate before trusting an environment."""

import importlib
import sys

mods = sys.argv[1:]
bad = []
for m in mods:
    try:
        mod = importlib.import_module(m)
        print(f"{m:18s} {getattr(mod, '__version__', '?')}")
    except Exception as e:  # noqa: BLE001
        bad.append(m)
        print(f"{m:18s} FAILED: {type(e).__name__}: {str(e)[:200]}")
sys.exit(1 if bad else 0)
