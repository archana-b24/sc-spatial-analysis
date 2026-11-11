"""Safe loading of third-party pickles (HLCA reference model, Geneformer dictionaries).

Unpickling can execute arbitrary code, so files downloaded from the internet are never passed to
``pickle.load`` / ``torch.load(weights_only=False)``. Instead:

* every external file is **sha256-pinned** in config and verified before loading;
* plain pickles go through :class:`RestrictedUnpickler`, which can only construct an explicit whitelist of
  data types (established by a static ``pickletools`` scan of the files — see DEVELOPMENT_LOG §3.3 G);
* torch checkpoints are read with ``torch.load(weights_only=True)`` (torch's own restricted unpickler).

Anything outside the whitelist raises ``pickle.UnpicklingError`` instead of being imported.
"""

from __future__ import annotations

import io
import pickle
from pathlib import Path

from .fetch import sha256

# (module, name) pairs that only *construct data*; nothing here can run a shell, import code or eval.
SAFE_GLOBALS: set[tuple[str, str]] = {
    ("builtins", "str"), ("builtins", "list"), ("builtins", "dict"), ("builtins", "tuple"),
    ("builtins", "set"), ("builtins", "int"), ("builtins", "float"), ("builtins", "bool"),
    ("__builtin__", "list"), ("__builtin__", "unicode"), ("__builtin__", "dict"),  # py2-era aliases
    ("collections", "OrderedDict"), ("collections", "defaultdict"),
    ("numpy", "ndarray"), ("numpy", "dtype"),
    ("numpy.core.multiarray", "_reconstruct"), ("numpy.core.multiarray", "scalar"),
    ("numpy._core.multiarray", "_reconstruct"), ("numpy._core.multiarray", "scalar"),
    ("_codecs", "encode"),  # numpy uses this to re-encode raw byte buffers (latin-1); returns bytes only
}


class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str):
        if (module, name) not in SAFE_GLOBALS:
            raise pickle.UnpicklingError(f"refusing to load non-whitelisted global {module}.{name}")
        if module == "__builtin__":
            module = "builtins"
            name = {"unicode": "str"}.get(name, name)
        if module.startswith("numpy.core"):  # numpy ≥2 moved numpy.core → numpy._core
            module = module.replace("numpy.core", "numpy._core")
        return super().find_class(module, name)


def verify(path: Path, expected_sha256: str) -> None:
    got = sha256(path)
    if got != expected_sha256:
        raise RuntimeError(f"sha256 mismatch for {path}: expected {expected_sha256}, got {got}")


def load_pickle(path: Path, expected_sha256: str):
    verify(path, expected_sha256)
    with open(path, "rb") as fh:
        return RestrictedUnpickler(io.BytesIO(fh.read())).load()


def load_torch_weights(path: Path, expected_sha256: str):
    import torch

    verify(path, expected_sha256)
    return torch.load(path, map_location="cpu", weights_only=True)
