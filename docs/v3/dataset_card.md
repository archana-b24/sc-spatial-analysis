---
pretty_name: LUSC single-cell & spatial atlas (DuckDB results snapshot)
tags:
  - single-cell
  - spatial-transcriptomics
  - lung-cancer
  - duckdb
size_categories:
  - 1M<n<10M
---

# LUSC atlas results snapshot (`lusc_atlas.duckdb`)

A single read-only **DuckDB** file with the results of a reproducible analysis of **lung squamous cell carcinoma
(LUSC)**: HLCA + LuCA single-cell RNA-seq (231,504 cells after QC, 151 donors, 20 studies) and one 10x Visium
section. It is the serving snapshot behind the project's dashboard, rebuilt from
pipeline outputs by the `export` step and never edited by hand.

* Code, pipeline, decision log: https://github.com/archana-b24/sc-spatial-analysis
* Inputs: public CELLxGENE releases of the Human Lung Cell Atlas and the Lung Cancer Atlas (LuCA), and the 10x
  Genomics FFPE LUSC Visium dataset; see the repository for exact versions and checksums.

## What is inside
`SELECT * FROM tables_info` lists every table, its row count, the pipeline output it came from and how it was
reduced (e.g. per-cell expression only for a 50k-cell stratified subsample and a marker panel; pseudobulk for all
genes). Main tables: `cells`, `expr_pseudobulk`, `de_pseudobulk`, `composition_sccoda`, `cnv_*`, `lineage_cells`,
`visium_spots`, `c2l_abundance`, `spatial_*`, `ccc_*`, plus `decisions`, `known_issues`, `glossary`, `meta`
(build time, git commit, config hash) and `schema_contract`.

**Read `known_issues` before quoting numbers.** Three open issues (Milo contrast, the `Unknown` cell type holding
un-called tumour cells, a cell2location artefact) change how some results should be read.

## Use
```bash
duckdb -readonly lusc_atlas.duckdb "SELECT * FROM tables_info"
```
```python
import duckdb
con = duckdb.connect("lusc_atlas.duckdb", read_only=True)
con.sql("""SELECT gene, log2FoldChange, padj FROM de_pseudobulk
           WHERE contrast = 'paired_tumour_vs_adjacent' AND "group" = 'Macrophage'
           ORDER BY padj LIMIT 10""").show()
```

Licensing of the underlying data follows the source atlases (CELLxGENE / 10x Genomics terms).
