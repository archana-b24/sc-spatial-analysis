"""Config loading and global reproducibility settings."""

from __future__ import annotations

import logging
import os
import random
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path | None = None) -> dict:
    path = Path(path) if path else REPO_ROOT / "config" / "config.yaml"
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["_repo_root"] = str(REPO_ROOT)
    return cfg


def resolve(cfg: dict, key: str, *parts: str) -> Path:
    """Absolute path for ``cfg['paths'][key] / parts`` (relative paths are repo-relative)."""
    base = Path(cfg["paths"][key])
    if not base.is_absolute():
        base = Path(cfg["_repo_root"]) / base
    return base.joinpath(*parts)


def set_seed(seed: int) -> None:
    """Seed every RNG the pipeline touches so reruns are bit-for-bit comparable where possible."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    import numpy as np

    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
    except ImportError:
        pass
    try:
        import scvi

        scvi.settings.seed = seed
    except ImportError:
        pass


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)-7s %(name)s | %(message)s")
