"""Read-only tools over the dashboard DuckDB — the only way the agent can see data.

Guardrails (each one closes a failure mode seen in published agent evaluations, see
docs/v3/A_agentic_research.md §3):

1. **Read-only at three levels.** The file is opened ``read_only=True``; external file/network access is
   disabled and the configuration locked (so ``read_csv('/etc/...')`` or ``ATTACH`` cannot escape the
   snapshot); and ``run_sql`` accepts exactly one statement whose parsed type is SELECT.
2. **Bounded results.** Every result is capped (``MAX_ROWS``) and flagged ``truncated`` so the model
   knows it is looking at a slice; a watchdog interrupts queries that run longer than ``TIMEOUT_S``.
3. **Provenance on every call.** Each result carries the exact SQL executed and the tables it read, so the
   UI can show where every number came from and the eval can check it.
4. **Caveats travel with the data.** Semantic tools attach the open known issues that affect the rows
   they return (e.g. the scCODA reference-type problem), so an answer cannot quote a number without
   having been handed its caveat.
"""

from __future__ import annotations

import json
import math
import os
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import duckdb
import pandas as pd

MAX_ROWS = 200
DEFAULT_ROWS = 50
TIMEOUT_S = 20.0

REPO = Path(__file__).resolve().parents[3]
LOCAL_DB = REPO / "results" / "dashboard" / "lusc_atlas.duckdb"

CONTRASTS = {
    "paired_tumour_vs_adjacent": "tumour vs adjacent normal, same patients (primary; ~donor + condition)",
    "tumour_vs_healthy": "tumour vs healthy donors (secondary; confounded with study by design)",
    "malignant_vs_healthy_basal": "CNV-called malignant cells vs healthy basal cells (epithelial question)",
}
CONDITIONS = ["healthy", "adjacent_normal", "tumour", "tumour_region_unknown"]
STAGES = [
    "cohort",
    "qc",
    "integration",
    "malignant",
    "de",
    "composition",
    "trajectory",
    "spatial",
    "ccc",
    "provenance",
]


def default_db_path() -> Path:
    """Same resolution order as the app (app/lib/db.py): $LUSC_DB, then the repo's built snapshot."""
    if os.environ.get("LUSC_DB"):
        return Path(os.environ["LUSC_DB"])
    return LOCAL_DB


class ToolError(Exception):
    """A tool refused or failed; the message is returned to the model as an error tool_result."""


@dataclass
class ToolResult:
    tool: str
    args: dict
    rows: list[dict] = field(default_factory=list)
    columns: list[str] = field(default_factory=list)
    n_rows: int = 0
    truncated: bool = False
    sql: list[str] = field(default_factory=list)
    tables: list[str] = field(default_factory=list)
    caveats: list[dict] = field(default_factory=list)
    note: str = ""
    error: str | None = None
    duration_s: float = 0.0

    def to_model_json(self, max_chars: int = 12_000) -> str:
        """Compact JSON for the LLM. Long results are cut by rows (never mid-row)."""
        payload = {
            "tool": self.tool,
            "n_rows": self.n_rows,
            "truncated": self.truncated,
            "tables": self.tables,
            "columns": self.columns,
            "rows": self.rows,
        }
        if self.note:
            payload["note"] = self.note
        if self.caveats:
            payload["caveats"] = self.caveats
        if self.error:
            payload = {"tool": self.tool, "error": self.error}
        txt = json.dumps(payload, default=_json_default)
        rows = list(self.rows)
        while len(txt) > max_chars and len(rows) > 1:
            rows = rows[: max(1, len(rows) // 2)]
            payload["rows"] = rows
            payload["truncated"] = True
            payload["note"] = (
                payload.get("note", "") + f" Only the first {len(rows)} rows are shown."
            ).strip()
            txt = json.dumps(payload, default=_json_default)
        return txt

    def to_dict(self) -> dict:
        return {
            "tool": self.tool,
            "args": self.args,
            "n_rows": self.n_rows,
            "truncated": self.truncated,
            "sql": self.sql,
            "tables": self.tables,
            "columns": self.columns,
            "rows": self.rows,
            "caveats": self.caveats,
            "note": self.note,
            "error": self.error,
            "duration_s": round(self.duration_s, 3),
        }


def _json_default(o: Any):
    if hasattr(o, "isoformat"):
        return o.isoformat()
    if isinstance(o, bytes):
        return f"<{len(o)} bytes>"
    try:
        return o.item()  # numpy scalars
    except AttributeError:
        return str(o)


def _clean(v: Any) -> Any:
    if isinstance(v, float):
        if math.isnan(v):
            return None
        return float(f"{v:.6g}")  # 6 significant digits: enough to quote, small for the context window
    if isinstance(v, bytes):
        return f"<{len(v)} bytes>"
    if hasattr(v, "item"):
        return _clean(v.item())
    return v


def _records(df: pd.DataFrame) -> list[dict]:
    return [{k: _clean(v) for k, v in r.items()} for r in df.to_dict(orient="records")]


# ------------------------------------------------------------------------------------- tool specs
# JSON-schema specs shared by every provider (llm.py converts the envelope) and by the MCP server.
def _spec(name: str, description: str, properties: dict, required: list[str] | None = None) -> dict:
    return {
        "name": name,
        "description": description,
        "input_schema": {
            "type": "object",
            "properties": properties,
            "required": required or [],
            "additionalProperties": False,
        },
    }


_STAGE = {"type": "string", "enum": STAGES, "description": "pipeline stage"}
_CONTRAST = {"type": "string", "enum": list(CONTRASTS), "description": "DE contrast; see describe_contrasts"}
_N = {"type": "integer", "minimum": 1, "maximum": MAX_ROWS, "description": "number of rows (default 10)"}

TOOL_SPECS: list[dict] = [
    _spec(
        "list_tables",
        "List every table in the atlas database with its row count, the pipeline output it came from, and "
        "how it was reduced. Call this first if you do not know where an answer lives.",
        {},
    ),
    _spec(
        "describe_table",
        "Columns, types and 3 example rows of one table. Use before writing SQL against it.",
        {"table": {"type": "string"}},
        ["table"],
    ),
    _spec(
        "run_sql",
        "Run ONE read-only DuckDB SELECT (or WITH … SELECT) against the atlas. Results are capped at "
        f"{MAX_ROWS} rows (default {DEFAULT_ROWS}); aggregate in SQL instead of fetching raw rows. Quote "
        'identifiers with spaces or capitals, e.g. "Cell Type", "group". Prefer a curated tool when one '
        "fits.",
        {
            "sql": {"type": "string"},
            "max_rows": {"type": "integer", "minimum": 1, "maximum": MAX_ROWS},
        },
        ["sql"],
    ),
    _spec(
        "cohort_summary",
        "Headline counts after QC: cells, donors, samples, studies, paired tumour/adjacent patients, "
        "CNV-called malignant cells, Visium spots, and cells/donors per atlas × condition.",
        {},
    ),
    _spec(
        "describe_contrasts",
        "The differential-expression contrasts, what each compares, and how many cell types / significant "
        "genes each has.",
        {},
    ),
    _spec(
        "top_de_genes",
        "Top differentially expressed genes for one contrast and cell type (pseudobulk PyDESeq2), ranked by "
        "adjusted p then |log2FC|. direction: up = higher in the case group (tumour / malignant).",
        {
            "contrast": _CONTRAST,
            "cell_type": {"type": "string", "description": "DE group, e.g. 'Macrophage', 'CD8 T'"},
            "direction": {"type": "string", "enum": ["up", "down", "both"]},
            "n": _N,
            "padj_max": {"type": "number", "minimum": 0, "maximum": 1},
        },
        ["contrast", "cell_type"],
    ),
    _spec(
        "gene_de_lookup",
        "One gene's DE result in every contrast × cell type where it was tested.",
        {"gene": {"type": "string", "description": "HGNC symbol, e.g. SPP1"}},
        ["gene"],
    ),
    _spec(
        "gene_expression",
        "Pseudobulk mean log-normalised expression and fraction of expressing cells for one gene in every "
        "cell type (optionally one condition), from all post-QC cells.",
        {
            "gene": {"type": "string"},
            "condition": {
                "type": "string",
                "enum": ["healthy", "adjacent_normal", "tumour", "tumour_region_unknown"],
            },
        },
        ["gene"],
    ),
    _spec(
        "composition_change",
        "scCODA compositional analysis: change in each cell type's proportion in one condition vs adjacent "
        "normal (log2 fold change, inclusion probability; credible if ≥ 0.8). Includes the reference-type "
        "caveat.",
        {"condition": {"type": "string", "enum": ["tumour", "healthy", "tumour_region_unknown"]}},
    ),
    _spec(
        "ccc_top_interactions",
        "Ligand–receptor interactions that change most between paired tumour and adjacent tissue (LIANA+ "
        "differential, from the paired DE). Optionally filter by sender and/or receiver cell type or gene.",
        {
            "source": {"type": "string", "description": "sender cell type"},
            "target": {"type": "string", "description": "receiver cell type"},
            "gene": {"type": "string", "description": "ligand or receptor gene"},
            "direction": {"type": "string", "enum": ["up", "down"]},
            "n": _N,
        },
    ),
    _spec(
        "caveats_and_decisions",
        "Known issues (with severity/status) and documented design decisions, optionally for one stage. "
        "ALWAYS check this before interpreting composition, Unknown cells, spatial epithelium or tensor "
        "CCC results.",
        {"stage": _STAGE},
    ),
]
TOOL_NAMES = [s["name"] for s in TOOL_SPECS]


# ------------------------------------------------------------------------------------------- tools
class AtlasTools:
    """All tools, bound to one read-only connection. Thread-safe enough for Streamlit: every call uses
    its own cursor."""

    def __init__(self, db_path: str | Path | None = None):
        path = Path(db_path) if db_path else default_db_path()
        if not path.exists():
            raise FileNotFoundError(f"atlas DB not found at {path}; run `pixi run -e dashboard export`")
        self.db_path = path
        self.con = duckdb.connect(str(path), read_only=True)
        # no file system / network / extension access from SQL, and nobody can switch it back on
        # (DuckDB security guide: external access off, no extension auto-install/-load, then lock).
        # DuckDB shares one instance per file within a process, so a second AtlasTools (or the app's own
        # connection) may find it already locked: then verify the lockdown instead of failing.
        try:
            for stmt in (
                "SET enable_external_access = false",
                "SET autoinstall_known_extensions = false",
                "SET autoload_known_extensions = false",
                "SET lock_configuration = true",
            ):
                self.con.execute(stmt)
        except duckdb.InvalidInputException:
            ext = self.con.execute("SELECT current_setting('enable_external_access')").fetchone()[0]
            if ext not in (False, "false"):
                raise RuntimeError("DuckDB configuration is locked with external access enabled") from None
        self._tables = set(
            self.con.execute("SELECT table_name FROM information_schema.tables").df().table_name
        )
        self._issues = self._query("SELECT * FROM known_issues ORDER BY id")

    # ---- plumbing
    def _query(self, sql: str, params: list | None = None) -> pd.DataFrame:
        cur = self.con.cursor()
        timer = threading.Timer(TIMEOUT_S, cur.interrupt)
        timer.start()
        try:
            return cur.execute(sql, params or []).df()
        except duckdb.InterruptException as e:
            raise ToolError(
                f"query exceeded {TIMEOUT_S:.0f}s and was stopped; aggregate or filter more"
            ) from e
        finally:
            timer.cancel()
            cur.close()

    def _tables_in(self, sql: str) -> list[str]:
        try:
            names = self.con.get_table_names(sql)
        except duckdb.Error:
            names = set()
        return sorted(n for n in names if n in self._tables)

    def _result(self, tool: str, args: dict, df: pd.DataFrame, sql: str, n: int | None = None, **kw):
        n = n or MAX_ROWS
        truncated = len(df) > n
        df = df.head(n)
        return ToolResult(
            tool=tool,
            args=args,
            rows=_records(df),
            columns=[str(c) for c in df.columns],
            n_rows=len(df),
            truncated=truncated,
            sql=[" ".join(sql.split())],
            tables=self._tables_in(sql),
            **kw,
        )

    def _issues_for(self, *, stages: list[str] = (), ids: list[int] = ()) -> list[dict]:
        d = self._issues
        d = d[(d.stage.isin(list(stages))) | (d.id.isin(list(ids)))]
        return [
            {
                "issue": int(r.id),
                "severity": r.severity,
                "status": r.status,
                "title": r.title,
                "detail": r.detail,
            }
            for r in d.itertuples()
        ]

    # ---- dispatch
    def call(self, name: str, args: dict | None = None) -> ToolResult:
        args = dict(args or {})
        t0 = time.time()
        fn = getattr(self, f"t_{name}", None)
        if name not in TOOL_NAMES or fn is None:
            res = ToolResult(tool=name, args=args, error=f"unknown tool '{name}'; available: {TOOL_NAMES}")
        else:
            try:
                res = fn(**args)
            except ToolError as e:
                res = ToolResult(tool=name, args=args, error=str(e))
            except TypeError as e:  # bad/missing arguments from the model
                res = ToolResult(tool=name, args=args, error=f"bad arguments: {e}")
            except duckdb.Error as e:
                res = ToolResult(tool=name, args=args, error=f"DuckDB error: {e}")
        res.duration_s = time.time() - t0
        return res

    # ---- generic tools
    def t_list_tables(self) -> ToolResult:
        sql = "SELECT table_name, n_rows, source, note FROM tables_info ORDER BY table_name"
        return self._result("list_tables", {}, self._query(sql), sql)

    def t_describe_table(self, table: str) -> ToolResult:
        if table not in self._tables:
            raise ToolError(f"no table '{table}'. Call list_tables.")
        sql = (
            "SELECT column_name, data_type FROM information_schema.columns "
            "WHERE table_name = ? ORDER BY ordinal_position"
        )
        cols = self._query(sql, [table])
        sample_sql = f'SELECT * FROM "{table}" LIMIT 3'
        sample = self._query(sample_sql)
        info = self._query("SELECT n_rows, source, note FROM tables_info WHERE table_name = ?", [table])
        note = ""
        if not info.empty:
            r = info.iloc[0]
            note = f"{int(r.n_rows):,} rows; source {r.source}" + (f"; {r.note}" if r.note else "")
        res = self._result("describe_table", {"table": table}, cols, sql.replace("?", f"'{table}'"))
        res.sql.append(sample_sql)
        res.tables = [table]
        res.note = note + f". Example rows: {json.dumps(_records(sample), default=_json_default)[:3000]}"
        return res

    def t_run_sql(self, sql: str, max_rows: int = DEFAULT_ROWS) -> ToolResult:
        sql = sql.strip().rstrip(";").strip()
        try:
            stmts = duckdb.extract_statements(sql)
        except duckdb.Error as e:
            raise ToolError(f"could not parse SQL: {e}") from e
        if len(stmts) != 1:
            raise ToolError("exactly one statement is allowed")
        if stmts[0].type != duckdb.StatementType.SELECT:
            raise ToolError(f"only SELECT queries are allowed (got {stmts[0].type.name})")
        n = max(1, min(int(max_rows), MAX_ROWS))
        wrapped = f"SELECT * FROM ({sql}) AS q LIMIT {n + 1}"
        df = self._query(wrapped)
        res = self._result("run_sql", {"sql": sql, "max_rows": n}, df, sql, n)
        res.sql = [" ".join(wrapped.split())]
        if res.truncated:
            res.note = f"More than {n} rows; only the first {n} returned. Aggregate or filter in SQL."
        return res

    # ---- curated semantic tools
    def t_cohort_summary(self) -> ToolResult:
        sql = """
        WITH c AS (SELECT count(*) cells, count(DISTINCT donor_id) donors, count(DISTINCT sample) samples,
                          count(DISTINCT study) studies, count(*) FILTER (WHERE malignant) malignant_cells
                   FROM cells),
             p AS (SELECT count(*) paired_patients FROM (SELECT donor_id FROM cohort_donors GROUP BY donor_id
                   HAVING bool_or(condition = 'tumour') AND bool_or(condition = 'adjacent_normal'))),
             v AS (SELECT count(*) visium_spots FROM visium_spots)
        SELECT * FROM c, p, v"""
        res = self._result("cohort_summary", {}, self._query(sql), sql)
        by_sql = """SELECT atlas, condition, count(DISTINCT donor_id) donors, count(*) cells
                    FROM cells GROUP BY 1, 2 ORDER BY 1, 2"""
        res.rows.append({"breakdown_by_atlas_condition": _records(self._query(by_sql))})
        res.sql.append(" ".join(by_sql.split()))
        res.note = (
            "Counts are post-QC (sc_annotated.h5ad). The pre-QC cohort (cohort_donors) has 160 donors / "
            "249,402 cells; 26 samples with < 200 cells were dropped before QC, removing 9 donors. "
            "paired_patients is counted on the pre-QC cohort_donors table."
        )
        return res

    def t_describe_contrasts(self) -> ToolResult:
        sql = """SELECT contrast, any_value(design) design, count(DISTINCT "group") n_cell_types,
                        count(*) FILTER (WHERE padj < 0.05) n_sig_gene_tests, count(*) n_tests
                 FROM de_pseudobulk GROUP BY 1 ORDER BY 1"""
        res = self._result("describe_contrasts", {}, self._query(sql), sql)
        for r in res.rows:
            r["meaning"] = CONTRASTS.get(r["contrast"], "")
        res.caveats = self._issues_for(ids=[2])
        return res

    def _check_contrast(self, contrast: str) -> None:
        if contrast not in CONTRASTS:
            raise ToolError(f"unknown contrast '{contrast}'; one of {list(CONTRASTS)}")

    def t_top_de_genes(
        self,
        contrast: str,
        cell_type: str,
        direction: str = "both",
        n: int = 10,
        padj_max: float = 0.05,
    ) -> ToolResult:
        self._check_contrast(contrast)
        groups = self._query('SELECT DISTINCT "group" FROM de_pseudobulk WHERE contrast = ?', [contrast])
        groups = groups["group"].tolist()
        if cell_type not in groups:
            raise ToolError(f"'{cell_type}' was not tested in {contrast}. Tested groups: {sorted(groups)}")
        if direction not in ("up", "down", "both"):
            raise ToolError("direction must be up, down or both")
        sign = {"up": "AND log2FoldChange > 0", "down": "AND log2FoldChange < 0", "both": ""}[direction]
        n = max(1, min(int(n), MAX_ROWS))
        sql = f"""SELECT gene, log2FoldChange, padj, baseMean, stat, n_case, n_control
                  FROM de_pseudobulk
                  WHERE contrast = '{contrast}' AND "group" = '{cell_type.replace("'", "''")}'
                    AND padj < {float(padj_max)} {sign}
                  ORDER BY padj ASC, abs(log2FoldChange) DESC"""
        df = self._query(sql)
        tot = len(df)
        args = {
            "contrast": contrast,
            "cell_type": cell_type,
            "direction": direction,
            "n": n,
            "padj_max": padj_max,
        }
        res = self._result("top_de_genes", args, df, sql, n)
        res.note = (
            f"{tot} genes pass padj < {padj_max} ({direction}). log2FoldChange > 0 = higher in case group."
        )
        if cell_type == "Unknown" or contrast == "malignant_vs_healthy_basal":
            res.caveats = self._issues_for(ids=[2])
        return res

    def t_gene_de_lookup(self, gene: str) -> ToolResult:
        gene = gene.strip().upper()
        sql = f"""SELECT contrast, "group" AS cell_type, log2FoldChange, padj, baseMean, n_case, n_control
                  FROM de_pseudobulk WHERE gene = '{gene.replace("'", "''")}' ORDER BY padj NULLS LAST"""
        df = self._query(sql)
        res = self._result("gene_de_lookup", {"gene": gene}, df, sql)
        if df.empty:
            res.note = f"{gene} was not tested (not in the gene space or too lowly expressed in every group)."
        elif (df.cell_type == "Unknown").any():
            res.caveats = self._issues_for(ids=[2])
        return res

    def t_gene_expression(self, gene: str, condition: str | None = None) -> ToolResult:
        gene = gene.strip().upper()
        where = f"gene = '{gene.replace(chr(39), chr(39) * 2)}'"
        if condition and condition not in CONDITIONS:
            raise ToolError(f"unknown condition; one of {CONDITIONS}")
        if condition:
            where += f" AND condition = '{condition}'"
        sql = f"""SELECT cell_type, condition, mean_lognorm, frac_expressing, n_cells
                  FROM expr_pseudobulk WHERE {where} ORDER BY mean_lognorm DESC"""
        df = self._query(sql)
        res = self._result("gene_expression", {"gene": gene, "condition": condition}, df, sql)
        res.note = (
            "Descriptive pseudobulk over all post-QC cells (no statistics). Cell-type/condition groups "
            "where no "
            "cell expresses the gene are omitted."
        )
        if df.empty:
            res.note = f"{gene} not found or not expressed."
        return res

    def t_composition_change(self, condition: str = "tumour") -> ToolResult:
        if condition not in CONDITIONS or condition == "adjacent_normal":
            raise ToolError(
                "condition must be tumour, healthy or tumour_region_unknown (reference = adjacent_normal)"
            )
        sql = f"""SELECT "Cell Type" AS cell_type, "log2-fold change" AS log2_fold_change,
                         "Inclusion probability" AS inclusion_probability,
                         "Inclusion probability" >= 0.8 AS credible
                  FROM composition_sccoda
                  WHERE Covariate LIKE '%T.{condition}]%'
                  ORDER BY credible DESC, abs("log2-fold change") DESC"""
        df = self._query(sql)
        if df.empty:
            raise ToolError(f"no scCODA covariate for condition '{condition}'")
        res = self._result("composition_change", {"condition": condition}, df, sql)
        res.note = (
            f"vs adjacent_normal. Credible = inclusion probability ≥ 0.8 ({int(df.credible.sum())} of "
            f"{len(df)} types). Effects are relative to the automatically chosen reference type 'Unknown'."
        )
        res.caveats = self._issues_for(ids=[1, 2])
        return res

    def t_ccc_top_interactions(
        self,
        source: str | None = None,
        target: str | None = None,
        gene: str | None = None,
        direction: str = "up",
        n: int = 10,
    ) -> ToolResult:
        conds = []
        for col, v in [("source", source), ("target", target)]:
            if v:
                conds.append(f"{col} = '{v.replace(chr(39), chr(39) * 2)}'")
        if gene:
            g = gene.strip().upper().replace("'", "''")
            conds.append(f"(ligand_complex LIKE '%{g}%' OR receptor_complex LIKE '%{g}%')")
        where = ("WHERE " + " AND ".join(conds)) if conds else ""
        order = "DESC" if direction == "up" else "ASC"
        n = max(1, min(int(n), MAX_ROWS))
        sql = f"""SELECT source, target, ligand_complex, receptor_complex, interaction_stat, interaction_padj,
                         ligand_stat, receptor_stat
                  FROM ccc_differential {where} ORDER BY interaction_stat {order}"""
        df = self._query(sql)
        args = {"source": source, "target": target, "gene": gene, "direction": direction, "n": n}
        res = self._result("ccc_top_interactions", args, df, sql, n)
        res.note = (
            f"{len(df)} matching interactions. interaction_stat > 0 = ligand (sender) and receptor "
            "(receiver) "
            "both higher in tumour vs adjacent (paired DE). Scores are a ranking, not an independent test."
        )
        if df.empty and (source or target):
            types = self._query("SELECT DISTINCT source FROM ccc_differential ORDER BY 1").source.tolist()
            res.note += f" Cell types available: {types}"
        res.caveats = self._issues_for(ids=[2] if "Unknown" in {source, target} else [])
        return res

    def t_caveats_and_decisions(self, stage: str | None = None) -> ToolResult:
        if stage and stage not in STAGES:
            raise ToolError(f"unknown stage; one of {STAGES}")
        w = f"WHERE stage = '{stage}'" if stage else ""
        sql = f"SELECT id, stage, severity, status, title, detail, log_ref FROM known_issues {w} ORDER BY id"
        issues = self._query(sql)
        res = self._result("caveats_and_decisions", {"stage": stage}, issues, sql)
        dsql = f"SELECT stage, question, decision, why, log_ref FROM decisions {w}"
        res.rows = [{"known_issues": res.rows, "decisions": _records(self._query(dsql))}]
        res.sql.append(dsql)
        res.tables = ["decisions", "known_issues"]
        res.note = "log_ref points to DEVELOPMENT_LOG.md sections in the repository."
        return res
