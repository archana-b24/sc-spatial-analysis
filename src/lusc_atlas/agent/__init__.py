"""v3 "Ask the atlas" agent: an LLM that answers questions about the LUSC atlas ONLY through read-only
tools over the dashboard DuckDB, and returns every SQL statement and table it used as provenance.

Layout (design and rationale: docs/v3/A_agentic_research.md):

* ``tools``  - the tool layer: schema listing, guarded read-only SQL, and curated semantic tools. Pure
               DuckDB + pandas; no LLM dependency, so the same functions back the agent, the MCP server
               and the tests.
* ``llm``    - provider-agnostic chat sessions: Claude (Anthropic SDK, default), any OpenAI-compatible
               endpoint (NVIDIA NIM / build.nvidia.com for Nemotron), and a scripted fake for tests.
* ``loop``   - the tool-use loop, the system prompt, and the post-hoc number-grounding check.
* ``mcp_server`` - the same tools over the Model Context Protocol (stdio), for Claude Desktop / Claude Code.

Why it lives in ``src/lusc_atlas`` rather than ``app/``: four consumers (Streamlit page, MCP server, eval
runner, CLI) import it, and it depends only on duckdb/pandas (+ optional provider SDKs), never on the
analysis stack.
"""

from .loop import AgentResult, ask
from .tools import AtlasTools

__all__ = ["AgentResult", "AtlasTools", "ask"]
