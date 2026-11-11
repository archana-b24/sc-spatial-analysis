# Changelog

All notable changes. Format: [Keep a Changelog](https://keepachangelog.com); versions: [SemVer](https://semver.org).
Detailed rationale for every item: [DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md).

## [2.0.0] — unreleased

### Added
- Importable package `lusc_atlas`, Snakemake workflow and `scripts/run_all.py` sharing one step registry.
- pixi environments (`default`, `c2l`, `workflow`) with lockfile; per-step provenance records.
- LuCA atlas (45 LUSC donors) alongside HLCA; cross-atlas deduplication, gene-space and label harmonisation.
- Integration benchmark (PCA, Harmony, scVI, scANVI, sysVI) with scib-metrics and an over-correction guard.
- Malignant-cell calling (infercnvpy), pseudobulk DE (PyDESeq2), decoupler 2 pathway/TF activity,
  scCODA/Milo abundance, Palantir + CellRank 2 trajectories, cell2location deconvolution,
  CellCharter niches, LIANA+ per-sample/differential/spatial communication.
- Unit tests, CI, pre-commit (nbstripout, ruff, large-file guard), CITATION.cff.

### Changed
- "SCLC" → **LUSC** throughout (SCLC denotes small-cell lung cancer).
- Visium input is the official 10x 2.0.1 release (the v1 matrix was not byte-identical to it).
- QC is per sample (MAD) and applied identically to every sample; doublets require Scrublet ∧ SOLO.

### Fixed
- v1 notebooks did not chain (doublet removal never reached downstream analysis).
- Control cells were filtered with the tumour QC threshold.
- Cell-level tests and phenotype-permutation GSEA treated cells as independent replicates.

### Removed
- Committed data, Space Ranger CSVs, `.DS_Store`, notebook outputs (all preserved in tag `v1.0.0`).

## [1.0.0] — 2024-05-29
Original BENG207 course project notebooks (tagged retrospectively on 2026-11-09).
