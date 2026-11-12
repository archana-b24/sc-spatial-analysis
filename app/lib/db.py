"""Read-only access to the dashboard snapshot (built by `lusc-atlas export`).

Where the DuckDB file comes from, in order:
1. `LUSC_DB` environment variable (a local path);
2. `results/dashboard/lusc_atlas.duckdb` in this repo (local development);
3. a public Hugging Face Dataset named by `LUSC_DB_HF_REPO` (the deployed app on Streamlit Community
   Cloud), downloaded once and cached, at revision `LUSC_DB_HF_REVISION` (default: main; a commit sha
   is preferred so the download is immutable).

Each `LUSC_DB_*` setting is read from the environment first, then from `st.secrets` (Community Cloud's
*Advanced settings → Secrets* TOML; locally `.streamlit/secrets.toml`, git-ignored). None of them is a
credential: the dataset is public and no token is used.

If `LUSC_DB_SHA256` is set (scripts/publish_dataset.sh prints the value), the file's sha256 must match or
the app refuses to start: the site then serves exactly the snapshot that was tested, never a silently
replaced file.
This DuckDB file is the app's single serving source of truth; nothing else is read at runtime.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

import duckdb
import pandas as pd
import streamlit as st

REPO = Path(__file__).resolve().parents[2]
LOCAL_DB = REPO / "results" / "dashboard" / "lusc_atlas.duckdb"


def setting(name: str) -> str:
    """`name` from the environment, else from Streamlit secrets (top-level key), else ''."""
    val = os.environ.get(name, "")
    if not val:
        try:
            val = str(st.secrets.get(name, "") or "")
        except Exception:  # noqa: BLE001 — no secrets file at all is the normal local case
            val = ""
    return val.strip()


def _locate() -> Path:
    if setting("LUSC_DB"):
        return Path(setting("LUSC_DB"))
    if LOCAL_DB.exists():
        return LOCAL_DB
    repo_id = setting("LUSC_DB_HF_REPO")
    if repo_id:
        from huggingface_hub import hf_hub_download

        return Path(
            hf_hub_download(
                repo_id=repo_id,
                filename="lusc_atlas.duckdb",
                repo_type="dataset",
                revision=setting("LUSC_DB_HF_REVISION") or None,
            )
        )
    raise FileNotFoundError(
        "No dashboard database found. Run `pixi run -e dashboard export`, or set LUSC_DB / LUSC_DB_HF_REPO."
    )


def sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:  # the DB file itself  noqa: source-of-truth
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


@st.cache_resource(show_spinner="Locating the results database…")
def _db_path() -> Path:
    path = _locate()
    want = setting("LUSC_DB_SHA256").lower()
    if want:
        got = sha256(path)
        if got != want:
            raise RuntimeError(
                f"Database checksum mismatch for {path.name}: expected {want[:12]}…, got {got[:12]}…. "
                "Refusing to serve an unverified snapshot; republish the dataset or fix LUSC_DB_SHA256."
            )
    return path


@st.cache_resource(show_spinner="Opening results database…")
def connection() -> duckdb.DuckDBPyConnection:
    return duckdb.connect(str(_db_path()), read_only=True)


@st.cache_data(show_spinner=False, max_entries=256)
def q(sql: str, params: tuple | None = None) -> pd.DataFrame:
    """Run a read-only query; results are cached per (sql, params)."""
    cur = connection().cursor()
    return cur.execute(sql, list(params) if params else None).df()


def table(name: str) -> pd.DataFrame:
    return q(f'SELECT * FROM "{name}"')


@st.cache_data(show_spinner=False)
def meta() -> dict[str, str]:
    df = q("SELECT key, value FROM meta")
    return dict(zip(df["key"], df["value"], strict=True))


def has_table(name: str) -> bool:
    return not q("SELECT 1 FROM information_schema.tables WHERE table_name = ?", (name,)).empty
