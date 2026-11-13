#!/usr/bin/env bash
# Publish the dashboard's DuckDB snapshot as a PUBLIC Hugging Face Dataset (free) — the only upload the
# free hosting path needs. The app itself is deployed by you on Streamlit Community Cloud from the public
# GitHub repo (docs/v3/C_dashboard_hosting.md §3.4). Run it yourself; nothing in the repo runs it.
#
#   export HF_TOKEN=hf_...                 # write token from https://huggingface.co/settings/tokens
#   scripts/publish_dataset.sh             # shows what it will do, asks, then uploads
#   scripts/publish_dataset.sh --yes       # no prompt
#   DATASET_NAME=lusc-atlas-db HF_USER=<you> scripts/publish_dataset.sh
#
# What it does: creates (if missing) the public dataset repo <HF_USER>/<DATASET_NAME>, uploads
# docs/v3/dataset_card.md as README.md and results/dashboard/lusc_atlas.duckdb, then prints the exact
# Streamlit secrets TOML (repo, the upload's commit sha as revision, the file's sha256) to paste into
# Community Cloud → Advanced settings → Secrets.
#
# Zero-cost guarantees (DEVELOPMENT_LOG §3.11, docs/v3/C_dashboard_hosting.md "Cost: $0"):
#   * it NEVER creates a Space (HF Docker/Gradio Spaces need a paid plan) and never sets hardware;
#   * the dataset is always public (public storage is free; private storage has a quota);
#   * it refuses files over MAX_MB (default 1024) so a mistake can't push a multi-GB object;
#   * the token comes from the environment, is passed only to huggingface_hub, never printed or written.
set -euo pipefail
cd "$(dirname "$0")/.."
source ../storage_env.sh 2>/dev/null || true

: "${HF_TOKEN:?set HF_TOKEN (a write token from https://huggingface.co/settings/tokens)}"
DATASET_NAME=${DATASET_NAME:-lusc-atlas-db}
MAX_MB=${MAX_MB:-1024}
DB=results/dashboard/lusc_atlas.duckdb
CARD=docs/v3/dataset_card.md
YES=0; [[ "${1:-}" == "--yes" ]] && YES=1
py() { pixi run -e dashboard python "$@"; }

[[ -f $DB ]] || { echo "missing $DB: run 'pixi run -e dashboard export' first"; exit 1; }
MB=$(( $(stat -c %s "$DB") / 1048576 ))
(( MB <= MAX_MB )) || { echo "refusing: $DB is ${MB} MB > MAX_MB=${MAX_MB}"; exit 1; }
SHA=$(sha256sum "$DB" | cut -d' ' -f1)

HF_USER=${HF_USER:-$(py - <<'PY'
from huggingface_hub import HfApi
print(HfApi().whoami()["name"])
PY
)}
DATASET="$HF_USER/$DATASET_NAME"

echo "About to publish (PUBLIC, free):"
echo "  dataset  https://huggingface.co/datasets/$DATASET"
echo "  files    README.md ← $CARD ;  lusc_atlas.duckdb ← $DB (${MB} MB, sha256 ${SHA:0:12}…)"
echo "  no Space, no hardware, no paid feature is created or changed"
if [[ $YES == 0 ]]; then read -r -p "Proceed? [y/N] " ok; [[ $ok == [yY]* ]] || { echo "aborted"; exit 1; }; fi

REV=$(DATASET="$DATASET" DB="$DB" CARD="$CARD" py - <<'PY'
import os
from huggingface_hub import HfApi

api = HfApi()
repo = os.environ["DATASET"]
api.create_repo(repo, repo_type="dataset", private=False, exist_ok=True)
if api.repo_info(repo, repo_type="dataset").private:
    raise SystemExit(f"{repo} exists and is PRIVATE; make it public or choose another DATASET_NAME")
api.upload_file(path_or_fileobj=os.environ["CARD"], path_in_repo="README.md", repo_id=repo,
                repo_type="dataset", commit_message="dataset card")
info = api.upload_file(path_or_fileobj=os.environ["DB"], path_in_repo="lusc_atlas.duckdb", repo_id=repo,
                       repo_type="dataset", commit_message="lusc_atlas.duckdb snapshot")
print(info.oid)
PY
)

cat <<EOF

Published: https://huggingface.co/datasets/$DATASET  (revision $REV)

Paste into Streamlit Community Cloud → your app → Settings → Secrets (none of these is a credential):

LUSC_DB_HF_REPO = "$DATASET"
LUSC_DB_HF_REVISION = "$REV"
LUSC_DB_SHA256 = "$SHA"
EOF
