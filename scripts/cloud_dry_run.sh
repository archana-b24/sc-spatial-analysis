#!/usr/bin/env bash
# Local rehearsal of the free deployment (Streamlit Community Cloud + public HF Dataset), with no upload
# and no account:
#
#   scripts/cloud_dry_run.sh                 # or: pixi run -e dashboard cloud-dry-run
#   KEEP_SERVING=1 scripts/cloud_dry_run.sh  # leave the :8502 server running afterwards
#
# It reproduces what Community Cloud does, as far as it can be reproduced offline:
#   1. "clone": copies exactly the files git would push (tracked + untracked-but-not-ignored) to a clean
#      directory — so anything git-ignored (results/, the DuckDB, .pixi/, secrets.toml) is absent;
#   2. "install": a fresh venv on Python $PY_VERSION (Community Cloud's default) from app/requirements.txt
#      only — no pixi env, no `pip install -e .`, no repo src/ on PYTHONPATH;
#   3. "secrets": writes .streamlit/secrets.toml in the clone (Community Cloud's Secrets box) with the
#      dataset repo, a pinned revision and the DB sha256; no LUSC_* environment variables are set;
#   4. "fetch": the app downloads the DB through the same huggingface_hub.hf_hub_download call as in
#      production, served here from a pre-seeded offline HF cache (HF_HUB_OFFLINE=1), then checks the
#      sha256 pin;
#   5. renders every page with Streamlit's AppTest in one process and reports peak RSS against the
#      Community Cloud memory ceiling; then serves the clone on :$PORT and checks /_stcore/health.
# What it cannot reproduce: Community Cloud's container, its real resource limits and the HF download.
set -euo pipefail
cd "$(dirname "$0")/.."
source ../storage_env.sh
REPO=$PWD
DB=$REPO/results/dashboard/lusc_atlas.duckdb
PY_VERSION=${PY_VERSION:-3.12}
WORK=${DRY_RUN_WORK:-${STORAGE_TOOLING:?}/cloud_dry_run}
VENV=$WORK/venv
CLONE=$WORK/clone
PORT=${PORT:-8502}
FAKE_REPO=dry-run/lusc-atlas-db
FAKE_REV=0123456789abcdef0123456789abcdef01234567   # any 40-hex commit id
CEILING_MB=${CEILING_MB:-2700}                        # Community Cloud max memory (docs, Feb 2024 figures)

[[ -f $DB ]] || { echo "missing $DB: run 'pixi run -e dashboard export' first"; exit 1; }
SHA=$(sha256sum "$DB" | cut -d' ' -f1)

echo "== 1. clean clone of what git would push → $CLONE"
rm -rf "$CLONE"; mkdir -p "$CLONE"
git ls-files -z --cached --others --exclude-standard | while IFS= read -r -d '' f; do
  [[ -e $f ]] && { mkdir -p "$CLONE/$(dirname "$f")"; cp -p "$f" "$CLONE/$f"; }
done
[[ ! -e $CLONE/results/dashboard ]] || { echo "DB leaked into the clone"; exit 1; }
[[ -f $CLONE/app/requirements.txt && -f $CLONE/app/streamlit_app.py ]] || { echo "app missing from clone"; exit 1; }

echo "== 2. fresh venv (Python $PY_VERSION) from app/requirements.txt"
PY=$(pixi run -e dashboard python -c "import sys;print(sys.executable)")
"$PY" -c "import sys; assert sys.version.startswith('$PY_VERSION'), sys.version" \
  || { echo "need a Python $PY_VERSION interpreter"; exit 1; }
if [[ ! -x $VENV/bin/python ]] || ! cmp -s "$CLONE/app/requirements.txt" "$VENV/.requirements.txt"; then
  rm -rf "$VENV"
  "$PY" -m venv "$VENV"
  "$VENV/bin/pip" install -q --upgrade pip
  "$VENV/bin/pip" install -q -r "$CLONE/app/requirements.txt"
  cp "$CLONE/app/requirements.txt" "$VENV/.requirements.txt"
fi

echo "== 3. secrets.toml (Community Cloud 'Secrets' box)"
mkdir -p "$CLONE/.streamlit"
cat > "$CLONE/.streamlit/secrets.toml" <<EOF
LUSC_DB_HF_REPO = "$FAKE_REPO"
LUSC_DB_HF_REVISION = "$FAKE_REV"
LUSC_DB_SHA256 = "$SHA"
EOF

echo "== 4. offline HF cache seeded with the DB at revision ${FAKE_REV:0:8}"
export HF_HOME=$WORK/hf_home HF_HUB_OFFLINE=1
SNAP=$HF_HOME/hub/datasets--${FAKE_REPO//\//--}/snapshots/$FAKE_REV
rm -rf "$HF_HOME"; mkdir -p "$SNAP" "$HF_HOME/hub/datasets--${FAKE_REPO//\//--}/refs"
ln -s "$DB" "$SNAP/lusc_atlas.duckdb"
echo -n "$FAKE_REV" > "$HF_HOME/hub/datasets--${FAKE_REPO//\//--}/refs/main"

echo "== 5. render every page from the clone (cwd = repo root, as on Community Cloud)"
cd "$CLONE"
export HOME=$WORK/home; mkdir -p "$HOME"
env -u PYTHONPATH -u LUSC_DB -u LUSC_DB_HF_REPO -u LUSC_DB_HF_REVISION -u LUSC_DB_SHA256 \
  CEILING_MB="$CEILING_MB" FAKE_REPO="$FAKE_REPO" "$VENV/bin/python" -I - <<'PY'
import os, resource, sys, tomllib
from pathlib import Path
from streamlit.testing.v1 import AppTest

secrets = tomllib.loads(Path(".streamlit/secrets.toml").read_text())
pages = sorted(p.name for p in Path("app/views").glob("*.py"))
bad = []
for page in pages:
    at = AppTest.from_file("app/streamlit_app.py", default_timeout=240)
    for k, v in secrets.items():
        at.secrets[k] = v
    at.run()
    at.switch_page(f"views/{page}")
    at.run()
    if at.exception:
        bad.append((page, [e.value for e in at.exception]))
    print(f"  {'FAIL' if at.exception else 'ok':4s} {page}")

# the DB must come through hf_hub_download (HF cache snapshot) + secrets, not a local fallback
sys.path.insert(0, "app")
import lib.db as db  # noqa: E402

path = db._locate()
print(f"  DB resolved via secrets + hf_hub_download → {path}")
assert os.environ["FAKE_REPO"].replace("/", "--") in str(path), "DB not fetched via the HF cache"
assert db.sha256(path) == secrets["LUSC_DB_SHA256"], "sha256 pin mismatch"
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
ceiling = float(os.environ["CEILING_MB"])
print(f"  peak RSS {peak:.0f} MB (all {len(pages)} pages, one process) = {peak / ceiling:.0%} of {ceiling:.0f} MB")
if bad:
    print(bad)
    sys.exit(1)
PY

echo "== 5b. negative control: a wrong sha256 pin must stop the app (fresh process, no cached resource)"
env -u PYTHONPATH -u LUSC_DB -u LUSC_DB_HF_REPO -u LUSC_DB_HF_REVISION -u LUSC_DB_SHA256 \
  "$VENV/bin/python" -I - <<'PY'
import tomllib
from pathlib import Path
from streamlit.testing.v1 import AppTest

secrets = tomllib.loads(Path(".streamlit/secrets.toml").read_text())
secrets["LUSC_DB_SHA256"] = "0" * 64
at = AppTest.from_file("app/streamlit_app.py", default_timeout=240)
for k, v in secrets.items():
    at.secrets[k] = v
at.run()
msgs = " ".join(str(e.value) for e in at.exception)
assert "checksum mismatch" in msgs, f"wrong pin was not refused: {msgs[:200]!r}"
print("  ok   wrong sha256 refused:", msgs[:90], "…")
PY

echo "== 6. serve the clone on :$PORT"
env -u PYTHONPATH -u LUSC_DB "$VENV/bin/streamlit" run app/streamlit_app.py --server.port "$PORT" \
  --server.headless true --browser.gatherUsageStats false > "$WORK/serve.log" 2>&1 &
PID=$!
for _ in $(seq 1 60); do curl -fs "http://localhost:$PORT/_stcore/health" >/dev/null 2>&1 && break; sleep 1; done
curl -fs "http://localhost:$PORT/_stcore/health" && echo "  ← health OK (pid $PID)"
if [[ "${KEEP_SERVING:-0}" != 1 ]]; then kill $PID; echo "stopped dry-run server"; fi
echo "CLOUD DRY RUN PASSED"
