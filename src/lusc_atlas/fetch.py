"""Data acquisition with provenance.

Every input file is either downloaded from a pinned, public source or extracted
from the immutable git tag ``v1.0.0``, and is always verified against a sha256
checksum declared in ``config/config.yaml``. A run never silently uses data it
cannot prove is identical to what the analysis was developed on.
"""

from __future__ import annotations

import hashlib
import logging
import subprocess
import urllib.request
from pathlib import Path

log = logging.getLogger(__name__)


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _try_download(url: str, dest: Path, timeout: int = 60) -> bool:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "lusc-atlas-fetch"})
        with urllib.request.urlopen(req, timeout=timeout) as r, open(dest, "wb") as out:
            while block := r.read(1 << 20):
                out.write(block)
        return True
    except Exception as exc:  # noqa: BLE001 — any network failure falls through to the next source
        log.warning("download failed (%s): %s", url, exc)
        dest.unlink(missing_ok=True)
        return False


def _try_git_tag(tag: str, repo_path: str, dest: Path, repo_root: Path) -> bool:
    try:
        blob = subprocess.run(
            ["git", "-C", str(repo_root), "show", f"{tag}:{repo_path}"],
            check=True,
            capture_output=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        log.warning("git fallback failed (%s:%s): %s", tag, repo_path, exc)
        return False
    dest.write_bytes(blob)
    return True


def fetch_visium(cfg: dict, out_dir: Path, repo_root: Path) -> Path:
    """Materialise the Visium Space Ranger outputs under ``out_dir``.

    Source order per file: existing verified copy → 10x URL → git tag fallback.
    Raises if no source yields a file matching the declared checksum.
    """
    vcfg = cfg["data"]["visium"]
    out_dir.mkdir(parents=True, exist_ok=True)
    for rel, meta in vcfg["files"].items():
        dest = out_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        expected = meta["sha256"]
        if dest.exists() and sha256(dest) == expected:
            log.info("ok (cached)  %s", rel)
            continue
        sources = []
        if meta.get("url"):
            sources.append(("10x", lambda d=dest, u=meta["url"]: _try_download(u, d)))
        if meta.get("git_path"):
            sources.append(
                (
                    "git:" + vcfg["git_tag"],
                    lambda d=dest, p=meta["git_path"]: _try_git_tag(vcfg["git_tag"], p, d, repo_root),
                )
            )
        for name, fn in sources:
            if fn() and sha256(dest) == expected:
                log.info("ok (%s)  %s", name, rel)
                break
            dest.unlink(missing_ok=True)
        else:
            raise RuntimeError(f"Could not obtain {rel} with sha256 {expected} from any source")
    return out_dir
