"""Two-atlas cohort construction: HLCA + LuCA → one deduplicated, harmonised AnnData.

Implements DEVELOPMENT_LOG §3.3 A–E. Every decision is written to a table so it can be audited:

* ``cohort_donors.csv``          one row per donor: atlas, study, condition, cells, selected?, reason
* ``merge_study_crosswalk.csv``  one row per (atlas, study): overlap key, kept/dropped, count check
* ``merge_gene_space.csv``       genes per atlas, intersection used, genes lost

Steps
1. **Primary cells only** (CELLxGENE ``is_primary_data``) — removes cell-level duplicates that one atlas
   re-distributes from another (A).
2. **Study crosswalk** — normalised study keys (surname tokens + year); any study still present in both
   atlases after step 1 is kept from ``prefer_atlas_on_overlap`` only (A).
3. **Condition** from disease + sample origin: ``tumour`` / ``adjacent_normal`` / ``healthy`` /
   ``tumour_region_unknown`` (D). Origin missing in one atlas is imported by sample ID from the other.
4. **Raw-count check** per study: integer, non-negative (B). Failing studies are dropped.
5. **Gene space**: Ensembl IDs without version, intersection across atlases (C).
6. **Labels**: Cell Ontology → harmonised label + compartment via ``config/label_crosswalk.csv`` (E).
7. Stratified per-donor caps (separate for disease vs healthy donors) bound CPU cost.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger(__name__)


# ----------------------------------------------------------------------------- helpers
def study_key(name: str) -> str:
    """Normalise a study label for cross-atlas matching: sorted surname tokens + year.

    'Thienpont_2018' and 'Lambrechts_Thienpont_2018_6653' share {'thienpont', '2018'}; we match on any shared
    surname *and* the same year, which is how the same deposit is named differently by two atlases.
    """
    s = str(name).lower()
    years = re.findall(r"(?:19|20)\d{2}", s)
    words = [w for w in re.findall(r"[a-z]+", s) if len(w) > 2 and w not in {"unpubl", "unpublished", "et", "al"}]
    return "|".join(sorted(set(words))) + "#" + (years[0] if years else "")


def keys_overlap(k1: str, k2: str) -> bool:
    w1, y1 = k1.split("#")
    w2, y2 = k2.split("#")
    return bool(set(w1.split("|")) & set(w2.split("|"))) and y1 == y2 and y1 != ""


def age_years(obs: pd.DataFrame) -> pd.Series:
    if "age_or_mean_of_age_range" in obs:
        a = pd.to_numeric(obs["age_or_mean_of_age_range"], errors="coerce")
        if a.notna().any():
            return a
    stage = obs.get("development_stage", pd.Series(index=obs.index, dtype=str)).astype(str)
    return stage.str.extract(r"(\d+)-year-old")[0].astype(float)


def sample_core_id(sample: str) -> str:
    """Numeric core of a sample ID so the same physical sample matches across atlases.

    HLCA 'BT1247' ↔ LuCA 'Lambrechts_Thienpont_2018_6149v1_1247'; HLCA 'scrBT1425' ↔ LuCA
    '..._scrBT1425_hg19'. Uses the ≥3-digit number in the final '_' token (so 'BT1A' → '' : no match).
    """
    last = str(sample).replace("_hg19", "").split("_")[-1]  # final token only: never batch/study numbers
    nums = re.findall(r"\d{3,}", last)
    return nums[-1] if nums else ""


def is_count_matrix(X, n: int = 2000) -> bool:
    X = X[: min(n, X.shape[0])]
    v = X.data if sp.issparse(X) else np.asarray(X).ravel()
    return bool(v.size) and bool((v >= 0).all()) and bool(np.allclose(v, np.round(v)))


def stratified_cap(obs: pd.DataFrame, max_per_donor: int, label: str, seed: int) -> np.ndarray:
    """Per donor ≤ max cells, proportional across `label` (≥1 per present label). Returns positional idx."""
    rng = np.random.default_rng(seed)
    keep = []
    pos = pd.Series(np.arange(len(obs)), index=obs.index)
    for _, d in obs.groupby("donor_id", observed=True):
        if len(d) <= max_per_donor:
            keep.append(pos[d.index].to_numpy())
            continue
        frac = max_per_donor / len(d)
        for _, g in d.groupby(label, observed=True):
            n = min(len(g), max(1, int(round(len(g) * frac))))
            keep.append(pos[rng.choice(g.index, n, replace=False)].to_numpy())
    return np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)


# ----------------------------------------------------------------------------- per-atlas
def assign_condition(obs: pd.DataFrame, acfg: dict, ccfg: dict, origin_lookup: dict[str, str]) -> pd.Series:
    lusc = obs["disease_ontology_term_id"].astype(str) == ccfg["tumour_disease_ontology"]
    healthy = obs["disease"].astype(str) == ccfg["normal_disease"]
    origin_col = acfg.get("origin_col")
    if origin_col and origin_col in obs:
        origin = obs[origin_col].astype(str)
    else:  # import from the other atlas: lookup keyed by (study_key, sample core id)
        sk = obs[acfg["study_col"]].astype(str).map(study_key)
        sid = obs[acfg["sample_col"]].astype(str).map(sample_core_id)
        origin = pd.Series(
            [next((v for (k, i), v in origin_lookup.items() if i == b and keys_overlap(k, a)), "unknown")
             for a, b in zip(sk, sid, strict=True)], index=obs.index)
    omap = ccfg["origin_map"]
    cond = pd.Series("excluded", index=obs.index, dtype=object)
    cond[lusc & origin.isin(omap["tumour"])] = "tumour"
    cond[lusc & origin.isin(omap["adjacent_normal"])] = "adjacent_normal"
    cond[lusc & ~origin.isin(omap["tumour"] + omap["adjacent_normal"] + omap["exclude"])] = "tumour_region_unknown"
    cond[healthy & ~origin.isin(omap["exclude"])] = "healthy"
    return cond


def load_atlas_obs(name: str, acfg: dict, repo: Path) -> tuple[ad.AnnData, pd.DataFrame]:
    a = ad.read_h5ad(repo / acfg["file"], backed="r")
    o = a.obs.copy()
    o["atlas"] = name
    o["study"] = o[acfg["study_col"]].astype(str)
    o["sample"] = o[acfg["sample_col"]].astype(str)
    return a, o


def build_merged_cohort(cfg: dict) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    repo = Path(cfg["_repo_root"])
    ccfg = cfg["data"]["cohort"]
    atlases = cfg["data"]["atlases"]
    seed = cfg["project"]["seed"]
    xw = pd.read_csv(repo / "config" / "label_crosswalk.csv")
    cl2lab = dict(zip(xw.cell_type_ontology_term_id, xw.harmonised_label, strict=True))
    cl2comp = dict(zip(xw.cell_type_ontology_term_id, xw.compartment, strict=True))

    handles, obs_all = {}, {}
    for name, acfg in atlases.items():
        handles[name], obs_all[name] = load_atlas_obs(name, acfg, repo)

    # origin lookup across atlases (D): LuCA annotates tumour/adjacent for Lambrechts samples HLCA lacks.
    # Keyed by (study_key, numeric sample id) so IDs only match within the same deposited study.
    origin_lookup: dict[tuple[str, str], str] = {}
    for name, acfg in atlases.items():
        oc = acfg.get("origin_col")
        if oc and oc in obs_all[name]:
            o = obs_all[name]
            m = o.groupby(["study", "sample"], observed=True)[oc].agg(lambda x: x.astype(str).mode().iat[0])
            for (st, sa), v in m.items():
                if sample_core_id(sa):
                    origin_lookup[(study_key(st), sample_core_id(sa))] = v

    # ---- per-atlas filtering + condition; every filter's effect is recorded
    sel, excl = {}, []
    rel = {ccfg["tumour_disease_ontology"]}
    for name, acfg in atlases.items():
        o = obs_all[name]
        o = o[o["disease_ontology_term_id"].astype(str).isin(rel) | (o["disease"].astype(str) == ccfg["normal_disease"])]
        steps_ = [
            ("not primary data (duplicate of a cell in another CELLxGENE dataset)", ~o["is_primary_data"].astype(bool)),
            ("tissue not primary lung", ~o["tissue"].astype(str).isin(ccfg["tissues"])),
            ("non-UMI assay (read counts; gene-length bias)", o["assay"].astype(str).isin(ccfg["exclude_assays"])),
        ]
        for reason, m in steps_:
            for (st, cond), n in o[m].groupby(["study", "disease"], observed=True).size().items():
                if n:
                    excl.append({"atlas": name, "study": st, "disease": cond, "reason": reason, "n_cells": int(n)})
            o = o[~m]
        o = o.copy()
        o["condition"] = assign_condition(o, acfg, ccfg, origin_lookup).values
        o["age_years"] = age_years(o).values
        drop = (o["condition"] == "excluded") | ((o["condition"] == "healthy") & (o["age_years"] < ccfg["healthy_min_age"]))
        for (st, cond), n in o[drop].groupby(["study", "condition"], observed=True).size().items():
            if n:
                reason = "excluded origin (metastasis/effusion)" if cond == "excluded" else f"healthy donor age < {ccfg['healthy_min_age']}"
                excl.append({"atlas": name, "study": st, "disease": cond, "reason": reason, "n_cells": int(n)})
        sel[name] = o[~drop]
    exclusions = pd.DataFrame(excl)

    # ---- study crosswalk + overlap resolution (A)
    rows = []
    for name, o in sel.items():
        for st, g in o.groupby("study", observed=True):
            rows.append({"atlas": name, "study": st, "study_key": study_key(st), "n_cells": len(g),
                         "n_donors": g["donor_id"].nunique(),
                         **{f"n_{c}": int((g["condition"] == c).sum()) for c in
                            ["tumour", "adjacent_normal", "tumour_region_unknown", "healthy"]}})
    cw = pd.DataFrame(rows)
    cw["overlaps_with"] = ""
    cw["decision"] = "keep"
    prefer = ccfg["prefer_atlas_on_overlap"]
    for i, r in cw.iterrows():
        others = cw[(cw.atlas != r.atlas)]
        hits = [o.study for _, o in others.iterrows() if keys_overlap(r.study_key, o.study_key)]
        if hits:
            cw.at[i, "overlaps_with"] = ";".join(hits)
            if r.atlas != prefer:
                cw.at[i, "decision"] = f"drop (duplicate study; kept from {prefer})"

    # ---- raw-count check per kept study (B)
    cw["counts_ok"] = np.nan
    for i, r in cw[cw.decision == "keep"].iterrows():
        a, o = handles[r.atlas], sel[r.atlas]
        pos = np.sort(obs_all[r.atlas].index.get_indexer(o.index[o["study"] == r.study][:2000]))
        X = (a.raw.X if a.raw is not None else a.X)[pos]
        ok = is_count_matrix(X)
        cw.at[i, "counts_ok"] = ok
        if not ok:
            cw.at[i, "decision"] = "drop (matrix is not raw integer counts)"
    keep_pairs = set(zip(cw.loc[cw.decision == "keep", "atlas"], cw.loc[cw.decision == "keep", "study"], strict=True))

    # ---- donor table + caps
    donor_rows, parts = [], {}
    for name, o in sel.items():
        o = o[[(name, s) in keep_pairs for s in o["study"]]]
        lab_col = "cell_type_ontology_term_id"
        is_healthy = o["condition"] == "healthy"
        idx_h = stratified_cap(o[is_healthy], ccfg["healthy_max_cells_per_donor"], lab_col, seed)
        idx_d = stratified_cap(o[~is_healthy], ccfg["disease_max_cells_per_donor"], lab_col, seed)
        chosen = pd.Index(o[is_healthy].index[idx_h]).append(o[~is_healthy].index[idx_d])
        parts[name] = o.loc[chosen]
        for (donor, cond), g in o.groupby(["donor_id", "condition"], observed=True):
            donor_rows.append({"atlas": name, "donor_id": f"{name}:{donor}", "study": g["study"].iat[0],
                               "condition": cond, "n_cells_available": len(g),
                               "n_cells": int(g.index.isin(chosen).sum()), "age_years": g["age_years"].iat[0]})
    donors = pd.DataFrame(donor_rows)

    # ---- gene space (C)
    gene_sets = {}
    for name, a in handles.items():
        var = a.raw.var if a.raw is not None else a.var
        ids = pd.Index(var.index.astype(str).str.split(".").str[0])
        filt = var.get("feature_is_filtered", pd.Series(False, index=var.index)).astype(bool).to_numpy()
        gene_sets[name] = set(ids[~filt])
    shared = sorted(set.intersection(*gene_sets.values()))
    genes = pd.DataFrame([{"atlas": n, "n_genes": len(g), "n_shared": len(shared), "n_lost": len(g) - len(shared)}
                          for n, g in gene_sets.items()])

    # ---- materialise counts on the shared gene space
    mats, obs_list, symbols = [], [], None
    for name, a in handles.items():
        o = parts[name]
        pos = np.sort(obs_all[name].index.get_indexer(o.index))
        var = a.raw.var if a.raw is not None else a.var
        ids = pd.Index(var.index.astype(str).str.split(".").str[0])
        col = ids.get_indexer(shared)
        X = (a.raw.X if a.raw is not None else a.X)[pos]
        X = sp.csr_matrix(X)[:, col].astype(np.float32)
        mats.append(X)
        ob = obs_all[name].iloc[pos]
        keep_cols = ["donor_id", "sample", "study", "assay", "disease", "disease_ontology_term_id", "tissue",
                     "sex", "development_stage", "cell_type", "cell_type_ontology_term_id", "suspension_type"]
        ob = ob[[c for c in keep_cols if c in ob]].astype(str).copy()
        ob["atlas"] = name
        ob["condition"] = o.loc[ob.index, "condition"].values
        ob["age_years"] = o.loc[ob.index, "age_years"].values
        fine = atlases[name].get("fine_label_col")
        ob["fine_label"] = obs_all[name].iloc[pos][fine].astype(str).values if fine and fine in obs_all[name] else ob["cell_type"]
        for c in ("donor_id", "sample"):
            ob[c] = name + ":" + ob[c]
        ob.index = name + ":" + ob.index.astype(str)
        obs_list.append(ob)
        if symbols is None and "feature_name" in var:
            symbols = pd.Series(var["feature_name"].astype(str).values, index=ids).groupby(level=0).first().reindex(shared)
    for a in handles.values():
        a.file.close()

    adata = ad.AnnData(X=sp.vstack(mats).tocsr(), obs=pd.concat(obs_list), var=pd.DataFrame(index=shared))
    adata.var["ensembl_id"] = adata.var_names
    adata.var["gene_symbol"] = symbols.values if symbols is not None else adata.var_names
    adata.var_names = pd.Index(adata.var["gene_symbol"].astype(str))
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    adata.obs["harmonised_label"] = adata.obs["cell_type_ontology_term_id"].map(cl2lab).fillna("Unknown")
    adata.obs["compartment"] = adata.obs["cell_type_ontology_term_id"].map(cl2comp).fillna("Unknown")
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype("category")
    log.info("merged cohort: %d cells × %d genes; conditions %s", adata.n_obs, adata.n_vars,
             adata.obs["condition"].value_counts().to_dict())
    return adata, donors, cw, genes, exclusions
