"""Foundation-model candidate: zero-shot Geneformer cell embeddings (Theodoris 2023; V2 models 2024–25).

Integration-benchmark candidate ``X_geneformer`` (DEVELOPMENT_LOG §3.3 G). Question it answers: does a
transcriptome foundation model, used **zero-shot**, integrate a cross-atlas tumour cohort as well as
purpose-trained generative models? Published benchmarks (Kedzierska 2025, Genome Biol) report that
zero-shot FMs often trail scVI and even HVG-PCA; we test whether that holds here. The default checkpoint is
``Geneformer-V2-104M_CLcancer`` (continual pre-training on cancer cells) — the variant most likely to help
on tumour data.

Implementation follows Geneformer's own tokenizer/EmbExtractor defaults, without its package:
1. per cell: counts / total × 10 000, divided by the gene's corpus median (``gene_median_dictionary``);
2. rank genes by that value (non-zero only), map Ensembl ID → token, keep the top ``input_size − 2``;
3. ``<cls>`` + ranked tokens + ``<eos>``;
4. embedding = mean of the **second-to-last** hidden layer over gene tokens (special tokens excluded).

All files are sha256-pinned; weights are ``safetensors`` (no pickle); dictionaries go through the restricted
unpickler. CPU inference is the bottleneck, so the candidate is computed on the benchmark subsample only
and is *benchmark-only* (not eligible to become the downstream embedding).
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

from .safe_load import load_pickle, verify

log = logging.getLogger(__name__)


def _download(fcfg: dict, rel: str, dest_dir: Path) -> Path:
    from huggingface_hub import hf_hub_download

    p = Path(hf_hub_download(fcfg["repo_id"], rel, revision=fcfg["revision"], local_dir=str(dest_dir)))
    return p


def load_geneformer(fcfg: dict, repo: Path):
    import torch
    from transformers import BertConfig, BertModel

    d = repo / fcfg["local_dir"]
    files = {}
    for rel, sha in fcfg["files"].items():
        p = d / rel
        if not p.exists():
            p = _download(fcfg, rel, d)
        verify(p, sha)
        files[rel] = p
    tok = load_pickle(files[fcfg["token_dict"]], fcfg["files"][fcfg["token_dict"]])
    med = load_pickle(files[fcfg["median_dict"]], fcfg["files"][fcfg["median_dict"]])
    cfg = BertConfig.from_json_file(str(files[f"{fcfg['model']}/config.json"]))
    model = BertModel.from_pretrained(
        str(d / fcfg["model"]), config=cfg, use_safetensors=True, add_pooling_layer=False
    )
    model.eval()
    torch.set_grad_enabled(False)
    return model, tok, med, cfg


def tokenize(
    X: sp.csr_matrix, ensembl: np.ndarray, tok: dict, med: dict, input_size: int
) -> list[np.ndarray]:
    """Rank-value encoding (Geneformer). X = raw counts, rows = cells."""
    keep = np.array([g in tok and g in med for g in ensembl])
    gidx = np.where(keep)[0]
    tokens = np.array([tok[g] for g in ensembl[gidx]], dtype=np.int64)
    medians = np.array([med[g] for g in ensembl[gidx]], dtype=np.float64)
    Xk = sp.csr_matrix(X)[:, gidx]
    totals = np.asarray(X.sum(axis=1)).ravel()
    cls_id, eos_id = tok["<cls>"], tok["<eos>"]
    out = []
    for i in range(Xk.shape[0]):
        row = Xk.getrow(i)
        if row.nnz == 0:
            out.append(np.array([cls_id, eos_id]))
            continue
        vals = row.data / max(totals[i], 1.0) * 10_000 / medians[row.indices]
        order = np.argsort(-vals, kind="stable")[: input_size - 2]
        out.append(np.concatenate([[cls_id], tokens[row.indices[order]], [eos_id]]))
    return out


def embed(adata: ad.AnnData, fcfg: dict, repo: Path, n_threads: int = 16) -> tuple[np.ndarray, dict]:
    import torch

    torch.set_num_threads(n_threads)
    model, tok, med, cfg = load_geneformer(fcfg, repo)
    input_size = min(fcfg["input_size"], cfg.max_position_embeddings)
    seqs = tokenize(
        adata.layers["counts"],
        adata.var["ensembl_id"].astype(str).str.split(".").str[0].to_numpy(),
        tok,
        med,
        input_size,
    )
    pad = tok["<pad>"]
    order = np.argsort([len(s) for s in seqs])  # length-bucketed batches → minimal padding
    emb = np.zeros((len(seqs), cfg.hidden_size), dtype=np.float32)
    bs, t0 = fcfg["batch_size"], time.time()
    for b in range(0, len(order), bs):
        idx = order[b : b + bs]
        L = max(len(seqs[i]) for i in idx)
        ids = np.full((len(idx), L), pad, dtype=np.int64)
        for r, i in enumerate(idx):
            ids[r, : len(seqs[i])] = seqs[i]
        ids_t = torch.from_numpy(ids)
        attn = (ids_t != pad).long()
        hs = model(input_ids=ids_t, attention_mask=attn, output_hidden_states=True).hidden_states[-2]
        gene_mask = attn.clone()
        gene_mask[:, 0] = 0  # drop <cls>
        lengths = attn.sum(1)
        gene_mask[torch.arange(len(idx)), lengths - 1] = 0  # drop <eos>
        m = gene_mask.unsqueeze(-1).float()
        emb[idx] = ((hs * m).sum(1) / m.sum(1).clamp(min=1)).numpy()
        if (b // bs) % 50 == 0:
            done = b + len(idx)
            log.info("geneformer %d/%d cells (%.2f s/cell)", done, len(order), (time.time() - t0) / done)
    info = {
        "model": fcfg["model"],
        "revision": fcfg["revision"],
        "input_size": input_size,
        "n_cells": len(seqs),
        "median_tokens": int(np.median([len(s) for s in seqs])),
        "seconds": round(time.time() - t0, 1),
    }
    return emb, info
