"""Fast unit tests on synthetic data (run in CI; no downloads, no model training)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lusc_atlas import steps


# ------------------------------------------------------------------ pipeline wiring
def test_every_input_is_produced_upstream():
    """A step may only consume files that an *earlier* step declares as output (the v1 failure mode R1)."""
    produced: set[str] = set()
    for s in steps.ORDER:
        missing = [i for i in steps.INPUTS[s] if i not in produced]
        assert not missing, f"step '{s}' reads {missing} which no earlier step writes"
        produced.update(steps.OUTPUTS[s])


def test_steps_registry_matches_order():
    assert list(steps.STEPS) == steps.ORDER
    assert set(steps.ENV) <= set(steps.ORDER)


# ------------------------------------------------------------------ QC
def test_mad_outlier_flags_extremes_only():
    pytest.importorskip("scanpy")  # pipeline env only
    from lusc_atlas.qc import _mad_outlier

    x = pd.Series(np.r_[np.random.default_rng(0).normal(10, 1, 500), [100.0, -80.0]])
    flags = _mad_outlier(x, 5)
    assert flags.iloc[-2:].all()
    assert flags.iloc[:-2].mean() < 0.01


def test_qc_is_per_sample():
    """A shallow sample must not be wiped out by a deep sample's distribution (v1 used global cutoffs)."""
    pytest.importorskip("scanpy")  # pipeline env only
    import anndata as ad

    from lusc_atlas.qc import flag_outliers

    rng = np.random.default_rng(0)
    n = 400
    obs = pd.DataFrame(
        {
            "sample": pd.Categorical(["deep"] * n + ["shallow"] * n),
            "log1p_total_counts": np.r_[rng.normal(9, 0.3, n), rng.normal(7, 0.3, n)],
            "log1p_n_genes_by_counts": np.r_[rng.normal(8, 0.3, n), rng.normal(6.5, 0.3, n)],
            "pct_counts_in_top_20_genes": rng.normal(20, 2, 2 * n),
            "pct_counts_mt": rng.normal(5, 1, 2 * n).clip(0),
            "n_genes_by_counts": np.full(2 * n, 1000),
            "total_counts": np.full(2 * n, 3000),
        },
        index=[f"c{i}" for i in range(2 * n)],
    )
    a = ad.AnnData(obs=obs)
    cfg = {
        "mad_log_counts": 5,
        "mad_log_genes": 5,
        "mad_pct_top20": 5,
        "mad_pct_mt": 3,
        "max_pct_mt": 20,
        "min_genes": 200,
    }
    flag_outliers(a, cfg)
    assert a.obs.loc[a.obs["sample"] == "shallow", "qc_outlier"].mean() < 0.02


# ------------------------------------------------------------------ DE design
def test_design_adds_study_only_when_not_collinear():
    from lusc_atlas.de import _design

    cfg = {"design": "~condition", "condition_key": "condition"}
    collinear = pd.DataFrame({"study": ["a", "a", "b", "b"], "condition": ["t", "t", "n", "n"]})
    crossed = pd.DataFrame({"study": ["a", "a", "b", "b"], "condition": ["t", "n", "t", "n"]})
    assert _design(collinear, cfg) == "~condition"
    assert _design(crossed, cfg) == "~study + condition"


# ------------------------------------------------------------------ integration guard
def test_knn_purity_detects_overcorrection():
    pytest.importorskip("scanpy")  # pipeline env only
    from lusc_atlas.integrate import knn_label_purity

    rng = np.random.default_rng(0)
    labels = pd.Series(["A"] * 200 + ["B"] * 200)
    separated = np.r_[rng.normal(0, 1, (200, 5)), rng.normal(8, 1, (200, 5))]
    collapsed = rng.normal(0, 1, (400, 5))  # an "over-corrected" embedding that merged A and B
    assert knn_label_purity(separated, labels) > 0.95
    assert knn_label_purity(collapsed, labels) < 0.65


def test_choose_rejects_methods_that_fail_the_guard():
    pytest.importorskip("scanpy")  # pipeline env only
    from lusc_atlas.integrate import choose

    scores = pd.DataFrame(
        {
            "embedding": ["X_pca", "X_scVI", "X_harmony"],
            "Total": [0.40, 0.60, 0.70],
            "knn_label_purity": [0.90, 0.88, 0.70],  # harmony scores higher but merges cell types
        }
    )
    assert choose(scores, tolerance=0.05) == "X_scVI"


# ------------------------------------------------------------------ CNV chr3q mapping
def test_chr3q_signal_uses_only_q_arm_windows():
    import anndata as ad

    from lusc_atlas.cnv import CHR3_CENTROMERE_END_GRCH38, chr3q_signal

    n_genes, window, step = 200, 20, 10
    starts = np.linspace(1e6, 190e6, n_genes)
    var = pd.DataFrame({"chromosome": "chr3", "start": starts}, index=[f"g{i}" for i in range(n_genes)])
    n_win = (n_genes - window) // step + 1
    centres = np.minimum(np.arange(n_win) * step + window // 2, n_genes - 1)
    is_q = starts[centres] > CHR3_CENTROMERE_END_GRCH38
    X = np.where(is_q, 1.0, -1.0)[None, :].repeat(3, axis=0)  # +1 on q windows, -1 on p windows
    a = ad.AnnData(np.zeros((3, n_genes)), var=var)
    a.obsm["X_cnv"] = X
    a.uns["cnv"] = {"chr_pos": {"chr3": 0}}
    np.testing.assert_allclose(chr3q_signal(a, window, step), 1.0)


# ------------------------------------------------------------------ provenance
def test_sha256_matches_hashlib(tmp_path):
    import hashlib

    from lusc_atlas.fetch import sha256

    p = tmp_path / "x.bin"
    payload = np.random.default_rng(0).bytes(3_000_000)  # > one 1 MiB read chunk
    p.write_bytes(payload)
    assert sha256(p) == hashlib.sha256(payload).hexdigest()


if __name__ == "__main__":
    pytest.main([__file__, "-q"])
