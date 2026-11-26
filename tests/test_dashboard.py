"""Every dashboard page renders without an exception against the built snapshot.

Needs the `dashboard` pixi env and a built DB (`pixi run -e dashboard export`); skipped otherwise, so the
default-env CI job is unaffected.
"""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
pytest.importorskip("duckdb")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).resolve().parents[1] / "app"
DB = Path(__file__).resolve().parents[1] / "results" / "dashboard" / "lusc_atlas.duckdb"
PAGES = sorted(p.name for p in (APP / "views").glob("*.py"))

pytestmark = pytest.mark.skipif(not DB.exists(), reason="dashboard DB not built")


def test_schema():
    import duckdb

    con = duckdb.connect(str(DB), read_only=True)
    names = set(con.execute("SELECT table_name FROM information_schema.tables").df().table_name)
    for t in [
        "cells",
        "expr_pseudobulk",
        "cell_expr",
        "decisions",
        "known_issues",
        "meta",
        "step_runs",
        "visium_spots",
    ]:
        assert t in names, t
    # every decision cites the log; every known issue has a severity the UI understands
    assert con.execute("SELECT count(*) FROM decisions WHERE log_ref IS NULL").fetchone()[0] == 0
    assert set(con.execute("SELECT DISTINCT severity FROM known_issues").df().severity) <= {
        "high",
        "medium",
        "low",
    }


@pytest.mark.parametrize("page", PAGES)
def test_page_renders(page):
    at = AppTest.from_file(str(APP / "streamlit_app.py"), default_timeout=180)
    at.run()
    at.switch_page(f"views/{page}")
    at.run()
    assert not at.exception, [e.value for e in at.exception]


def test_schema_contract():
    """Every table/column the app and agent rely on (export.SCHEMA_CONTRACT) exists in the built DB, and
    the DB records the same contract it was built against."""
    import duckdb

    from lusc_atlas.export import SCHEMA_CONTRACT, SCHEMA_VERSION

    con = duckdb.connect(str(DB), read_only=True)
    cols = con.execute("SELECT table_name, column_name FROM information_schema.columns").df()
    have = set(zip(cols.table_name, cols.column_name, strict=True))
    missing = [(t, c) for t, cs in SCHEMA_CONTRACT.items() for c in cs if (t, c) not in have]
    assert not missing, (
        f"DB violates schema contract (rebuild with `pixi run -e dashboard export`): {missing}"
    )
    stored = con.execute("SELECT table_name, column_name FROM schema_contract").fetchall()
    assert set(stored) == {(t, c) for t, cs in SCHEMA_CONTRACT.items() for c in cs}
    meta = dict(con.execute("SELECT key, value FROM meta").fetchall())
    assert meta["schema_version"] == str(SCHEMA_VERSION)
    for k in ("built_at", "git_sha", "config_sha256", "content_sha256", "role", "run_id", "run_label"):
        assert meta.get(k), k
    # tables_info must describe every table actually present
    info = set(con.execute("SELECT table_name FROM tables_info").df().table_name)
    assert set(cols.table_name) <= info | {"tables_info"}


def test_app_reads_only_duckdb():
    """Single serving source of truth: no page or helper reads csv/h5ad/yaml/json files at runtime."""
    import re

    forbidden = re.compile(
        r"read_csv|read_h5ad|read_parquet|yaml\.|open\(|\.read_text\(|json\.load\(|h5py|anndata"
    )
    offenders = []
    for p in sorted((APP).rglob("*.py")):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if forbidden.search(line) and "noqa: source-of-truth" not in line:
                offenders.append(f"{p.relative_to(APP)}:{i}: {line.strip()}")
    assert not offenders, offenders


def test_guides_compute_interpretations():
    """Every stage has a reader guide whose 'What it means' sentences compute without falling back."""
    import sys

    sys.path.insert(0, str(APP))
    from lib.explain import FINDINGS, findings

    for stage in FINDINGS:
        out = findings(stage)
        assert out, stage
        assert not any("interpretation unavailable" in s for s in out), (stage, out)


def test_db_checksum_pin(monkeypatch):
    """With LUSC_DB_SHA256 set (deploy bundle), a mismatching snapshot is refused; a matching one served."""
    import hashlib

    h = hashlib.sha256()
    with DB.open("rb") as f:
        while b := f.read(1 << 22):
            h.update(b)
    import streamlit as st

    st.cache_resource.clear()  # earlier tests in this process cached the unpinned DB path/connection
    st.cache_data.clear()
    monkeypatch.setenv("LUSC_DB", str(DB))
    monkeypatch.setenv("LUSC_DB_SHA256", "0" * 64)
    at = AppTest.from_file(str(APP / "streamlit_app.py"), default_timeout=120)
    at.run()
    assert at.exception and "checksum mismatch" in at.exception[0].value
    monkeypatch.setenv("LUSC_DB_SHA256", h.hexdigest())
    at = AppTest.from_file(str(APP / "streamlit_app.py"), default_timeout=120)
    at.run()
    assert not at.exception
    st.cache_resource.clear()


def test_series_color_never_runs_out():
    """Run 2 produced 9 CellCharter niches, one more than SERIES; categorical colours must not IndexError."""
    import sys

    sys.path.insert(0, str(APP))
    from lib.ui import SERIES, series_color

    cols = [series_color(i) for i in range(12)]
    assert cols[: len(SERIES)] == SERIES
    assert len(set(cols)) == 12
    assert series_color(12) == series_color(0)


def test_every_page_names_its_part():
    """§3.13: every stage page belongs to exactly one part (built / Evidence / Context), the sidebar groups
    say so, and page titles carry no stage numbers (order no longer equals group)."""
    import re
    import sys

    import duckdb

    sys.path.insert(0, str(APP))
    from lib.ui import PART, STAGE_PART

    con = duckdb.connect(str(DB), read_only=True)
    stages = set(con.execute("SELECT stage FROM page_guide").df().stage)
    assert stages <= set(STAGE_PART), stages - set(STAGE_PART)
    assert set(STAGE_PART.values()) <= set(PART)
    assert {s for s, p in STAGE_PART.items() if p == "evidence"} == {"de", "composition"}
    entry = (APP / "streamlit_app.py").read_text()
    assert "Evidence · same-patient comparisons" in entry and "Context · descriptive" in entry
    assert not re.search(r'title="\d+ · ', entry)
    for stage in STAGE_PART:
        src = (APP / "views" / f"{stage}.py").read_text()
        assert f'part="{STAGE_PART[stage]}"' in src, stage
