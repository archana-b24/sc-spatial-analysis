"""The atlas tools over the Model Context Protocol (stdio), so any MCP client (Claude Desktop, Claude Code,
other agent frameworks) can query the same read-only snapshot with the same guardrails.

    pixi run -e dashboard mcp                      # stdio server (what MCP clients launch)
    LUSC_DB=/path/lusc_atlas.duckdb pixi run -e dashboard mcp

Claude Code:   claude mcp add lusc-atlas -- pixi run --manifest-path <repo>/pixi.toml -e dashboard mcp
Claude Desktop (claude_desktop_config.json):
    {"mcpServers": {"lusc-atlas": {"command": "pixi",
        "args": ["run", "--manifest-path", "<repo>/pixi.toml", "-e", "dashboard", "mcp"]}}}

Written for the `mcp` 2.x Python SDK (``MCPServer``; v1 called it ``FastMCP``). Every tool is annotated
read-only / idempotent / closed-world. The results are the JSON the in-app agent sees, plus the SQL and
table provenance.
"""

from __future__ import annotations

import functools
import inspect
import json

from .tools import TOOL_SPECS, AtlasTools

INSTRUCTIONS = (
    "Read-only tools over a lung squamous cell carcinoma single-cell + spatial atlas (DuckDB snapshot). "
    "Quote numbers only from tool results; every result lists the SQL and tables it came from. Check "
    "caveats_and_decisions before interpreting composition, `Unknown` cells, spatial epithelium or "
    "tensor CCC."
)


def build_server(tools: AtlasTools | None = None):
    from mcp.server.mcpserver import MCPServer
    from mcp.types import ToolAnnotations

    tools = tools or AtlasTools()
    server = MCPServer(name="lusc-atlas", instructions=INSTRUCTIONS)
    ro = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False)
    for spec in TOOL_SPECS:
        method = getattr(tools, f"t_{spec['name']}")

        def make(name=spec["name"], method=method):
            @functools.wraps(method)
            def fn(**kwargs) -> str:
                return json.dumps(tools.call(name, kwargs).to_dict(), default=str)

            # keep the typed parameters (schema inference) but declare the JSON-string return
            fn.__signature__ = inspect.signature(method).replace(return_annotation=str)
            fn.__annotations__ = {**getattr(method, "__annotations__", {}), "return": str}
            fn.__name__ = name
            return fn

        server.add_tool(make(), name=spec["name"], description=spec["description"], annotations=ro)
    return server


def main() -> None:
    build_server().run("stdio")


if __name__ == "__main__":
    main()
