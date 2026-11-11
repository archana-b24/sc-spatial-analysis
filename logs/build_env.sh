#!/usr/bin/env bash
set -euo pipefail
ENV=lusc-v2
conda create -y -n $ENV -c conda-forge python=3.12 pip
PY=$(conda run -n $ENV which python)
$PY -m pip install --upgrade pip
$PY -m pip install torch --index-url https://download.pytorch.org/whl/cpu
$PY -m pip install "scanpy[leiden]" anndata squidpy scvi-tools cellxgene-census \
  infercnvpy decoupler pydeseq2 liana cellrank palantir pertpy gseapy \
  snakemake papermill jupyterlab ipykernel nbstripout pre-commit pytest ruff pyyaml pydantic
$PY -m pip freeze > logs/pip_freeze_initial.txt
echo BUILD_DONE
