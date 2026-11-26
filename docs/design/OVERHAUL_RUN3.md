# Overhaul + run 3: every page and every pipeline step, reviewed

*Started 2026-11-14 18:00 PST, run unattended in the `sc-deploy` worktree (branch `v2-dev`). Nothing here is
committed, pushed, deployed or uploaded. The engineering record is DEVELOPMENT_LOG §3.12; this page is the
summary and the plan. Study notes for the methods: [LEARNING_NOTES.md](LEARNING_NOTES.md).*

> **Status (2026-11-17): Run 4 DONE.** Run 3 complete (2026-11-16); then the control-group decision (§2.3,
> DEVELOPMENT_LOG §3.13) changed the healthy cohort, so run 4 reran everything from `select` (finished 2026-11-17
> 05:21). Run 4 changed only Context results, by design: the paired Evidence is identical to run 3. The dashboard
> is now split into **Evidence** (same-patient comparisons) and **Context** (descriptive). The app in tmux
> `lusc-v2:app` serves the run-4 database. Nothing is committed, pushed, published or deployed — redeploy steps
> in §7. Run-4 numbers: DEVELOPMENT_LOG §3.13 "Run 4 results".

## TL;DR — results (runs 3 and 4)

* **The design question is answered by the data.** Random halves of one study differ in ≈ 1.4 % of genes (the
  test is calibrated), two sets of healthy studies in 6.4–60.8 % (run 4, lung tissue only; run 3 had 14–73 % with
  a lavage study in one set): cross-study differences swamp biology, so only the within-patient comparison (15
  patients, 4 studies, immune cells only) is used as evidence. Run 4's label-swap null shows the paired test
  itself is calibrated (0.51 % of genes at p < 0.01, none of 10 types flagged).
* **Within patients**, 126 robust immune gene changes (*CXCL13*/*ENTPD1*/*HAVCR2* in CD8 T, *SPP1*/*MMP12* in
  macrophages) and a consistent immune shift (Treg, B, plasma up; alveolar macrophages, non-classical monocytes
  down).
* **Malignant calls failed their pre-registered false-positive check** (7.9 % > 5 %, one study); kept as
  provisional and flagged everywhere; dropping that study leaves the dependent results unchanged (§5.1).
* **The spatial section is now squamous carcinoma**; runs 1–2 analysed a neuroendocrine tumour.
* What is still open: §6. What to study: [LEARNING_NOTES.md](LEARNING_NOTES.md).

## TL;DR (plan stage, kept for the record)

Two independent audits (pipeline statistics; dashboard pages) found three problems that change conclusions,
and many that weaken them:

1. **The Visium section is not LUSC.** `CytAssist_11mm_FFPE_Human_Lung_Cancer` (6,195 spots) is a
   neuroendocrine carcinoma: CHGA is ~10 % of every spot's library (1,067 counts per 10k), SYP/NCAM1 in > 90 %
   of spots, while KRT5/TP63/SOX2 are essentially absent (≤ 0.04 per 10k). Every spatial result so far
   describes the wrong tumour type. 10x publishes a squamous section, `CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma`
   (3,858 spots, Space Ranger 2.0.0; the section v1 used): KRT5, TP63, SOX2 detected in 72–83 % of spots, CHGA in 3 %.
2. **Only the within-patient paired contrast is protected from study/atlas confounding.** Healthy tissue is
   100 % HLCA, tumour tissue is 95 % LuCA, and no study contains both; sorting protocols differ by study
   (CD45-enriched, CD31-enriched, epithelium-rich). `tumour_vs_healthy` DE cannot include study at all
   (aliased), scCODA pooled all 215 samples unpaired, and the CNV null came from the other atlas.
3. **"Malignant" calls are a cluster-level artefact resting on ~3 patients.** infercnvpy's `cnv_score` is the
   *cluster mean*, so the "95th percentile of healthy epithelium" threshold was a percentile over ~29 cluster
   values; 87 % of the 3,942 calls come from 3 donors.

Plus: the mito filter never ran (§3.10 #9), the healthy-age filter was bypassed for unknown ages (3 donors aged
22–38 kept), SOLO is miscalibrated (17 %), multi-region samples were counted as independent replicates in the
"paired" analyses, the scANVI benchmark is circular, and several pages state numbers the DB contradicts.

The answer to "is merging donors from here and there the right approach?" is **yes, for representation; no,
for inference** — see §2.1.

---

## 1. Caveat inventory

Severity: **C** = changes conclusions · **W** = weakens · **P** = presentation only. Source: A = statistics
audit, B = pages audit, U = user / earlier session.

### 1.1 Pipeline (by step)

| Step | # | Caveat | Sev | Src |
|---|---|---|---|---|
| select | S1 | Healthy min-age filter bypassed: `age_years()` falls back to `development_stage` only if the whole column is NaN; `NaN < 45` is False → 62/111 healthy donors have no age, 3 donors aged 22/24/38 kept | W | A |
| select/qc | S2 | Mito fraction never computed (LuCA ships no MT- genes → intersection drops them → 0 % everywhere) | W | U |
| qc | Q1 | Mito MAD rule two-sided (removes *low*-mito cells) and with MAD 0 flags every non-zero cell | W | A |
| qc | Q2 | SOLO flags 17 % (HLCA 20 %, some samples 95 %); `predict(soft=False)` = argmax of a balanced classifier, ignores expected rate. Scrublet threshold failed in HLCA:PX9 (972/1,001) | W | A, U |
| integrate | I1 | scANVI scored on the labels it was trained on → bio-conservation circular; it "won" partly by construction | W | A |
| integrate | I2 | Atlas-mixing metric compares healthy (HLCA) with tumour (LuCA) cells of the same type → "below random" is expected, not diagnostic; scIB batch = study is aliased with condition | W | A, B |
| cnv | C1 | `cnv_score` = cluster mean → threshold is a percentile over clusters; joint clustering of HLCA-healthy + LuCA-tumour epithelium | **C** | A |
| cnv | C2 | Null from healthy (HLCA) epithelium = other atlas/chemistry | W | A |
| cnv | C3 | 87 % of calls from 3 donors (Wu_P10 1,704, UKIM-V_P2 1,287, Wu_P4 451) → Malignant DE / trajectory are n ≈ 3 | **C** | A |
| cnv | C4 | 13 % "sensitivity" vs author labels (author labels are themselves inferCNV-informed, not truth) | W | A, B |
| de | D1 | Paired DE uses *samples* as units: multi-region tumour samples are pseudo-replicates within donor (Leader_Merad_706: 8 tumour + 8 adjacent) | W | A |
| de | D2 | `tumour_vs_healthy` is unidentifiable (study ≡ condition; assay ≈ atlas); 47–58 % of genes p < 0.01 in immune types | **C** | A |
| de | D3 | Groups tested with 3 donors; `Unknown` tested | W | A, B |
| de | D4 | Leave-one-study-out promised in §3.3 K, never implemented; no per-study replication | **C** | A, U |
| de | D5 | Malignant / CNV-low vs healthy basal: atlas + chemistry + state confounded; top hits MT pseudogenes | W | A |
| abundance | A1 | scCODA on all 215 samples, `~condition`, unpaired: study sorting protocols dominate; CNV-low +2.23 by construction; tumour-only types dilute everything else | **C** | A, B |
| abundance | A2 | Milo: samples, not donors (Leader 706 = 16/48 samples) | W | A |
| abundance | A3 | scCODA and Milo disagree (CD4 T) without reconciliation | W | B |
| trajectory | T1 | Root = HLCA healthy basal, endpoint = LuCA malignant from ~3 donors → pseudotime partly atlas distance; embedding label-guided | W | A, B |
| spatial | V1 | **Visium section is a neuroendocrine carcinoma, not LUSC** | **C** | A (verified) |
| deconvolve | V2 | Proportions from q05 (a lower bound), not posterior mean | W | A |
| deconvolve | V3 | Malignant vs CNV-low as separate signatures (differ by a CNV threshold, not expression) | W | A, B |
| ccc | X1 | Tensor: samples as units, outer-fill 0 encodes sorting protocol; Mann–Whitney ignores pairing → factors track protocol, not condition (§3.10 #4) | W | A |
| ccc | X2 | Malignant absent from tensor (§3.10 #8) | W | U |
| ccc | X3 | `ccc` waits on deconvolve only for the spatial part → needless 11 h serialisation | P (cost) | A |
| pathways | P1 | Pathways computed on confounded contrasts presented alongside paired ones | W | A |

### 1.2 Pages (dashboard)

| Page | Main caveats (P1 = misleading/wrong) | Src |
|---|---|---|
| All | P1: confound never stated plainly; no n donors / n studies on results; curated text contradicts DB (CNV reference, 5 vs 6 tensor factors, p ≥ 0.31 stale, "27 cell types"); decisions hidden in expanders; methods unexplained on pages 3–9; condition / compartment / niche colours collide; no run label | B |
| Start | CNV reference described wrongly; QC text silent on mito; hard-coded "(cDC)", "Seven approaches" | B |
| Overview | "Most healthy lung comes from HLCA" (it is 100 %); no studies column; 15 (pre-QC) vs 14 (post-QC) paired unexplained; healthy cap not mentioned; finding cards without n | B |
| Integration | circular label-transfer check; mixing chart uncaptioned and confounded; batch = sample vs benchmark batch = study unexplained; hard-coded numbers | B |
| Malignant | log axis hides 17/38 donors with 0 calls; 5 % false-positive rate by construction unstated; threshold line drawn against group means; violin on raw keys | B |
| DE | CNV-low contrast unreachable; n donors only in hover; confounded contrasts not flagged; no CIs (lfcSE unused); SATB1 tops macrophages (contamination?) | B |
| Composition | caption wrong ("log2FC 0 = reference": all non-credible types are +0.407); tumour-only types dilute; scCODA/Milo disagreement; "Healthy" radio confounded | B |
| Trajectory | atlas confound; promised CellRank outputs not shown (raw JSON); colours | B |
| Spatial | top SVGs are neuroendocrine (CHGB, CARTPT, TTR, GC) — symptom of V1; niches unnamed / over-split; co-localisation statistic undefined | B |
| CCC | factor 6 unselectable; per-sample counts without denominators; "top 1 %" misdescribes magnitude rank; only gains shown | B |
| Explorer / Glossary / Provenance | no guide/caveats; missing method terms; run not labelled | B |

---

## 2. Decisions (options considered, choice, why)

### 2.1 The central question: pooling HLCA + LuCA donors "from here and there"

> **Note (2026-11-16).** This section is the plan-stage text. Where it differs from what was run, it is superseded;
> the as-run design is in DEVELOPMENT_LOG §3.12 B. Main differences: the paired DE uses **15 patients**
> (Leader_Merad 5, UKIM-V 5, Goveia_Carmeliet 3, Thienpont 2; Thienpont tumour samples with no recorded region
> count as tumour tissue) and paired composition 14; a primary cell type needs ≥ 6 pairs and ≥ 2 studies with ≥ 2
> pairs; "robust" follows `de.paired_robustness`; composition uses a paired NB GLM plus CLR and arcsin-√ paired
> tests (not scCODA); the tensor has 11 patients. `adjacent_vs_healthy` is a **decomposition term, not a negative
> control**. The cohort numbers below are from run 2; run 3 after QC has 148 donors, 19 studies and 108 healthy
> donors, of which 17 are tumour-free surgery tissue from cancer patients and 10 are bronchoalveolar lavage (§2.3).
> Option 2 is worded too strongly: a study random effect can be fitted, but it would credit systematic study
> differences (atlas, sorting, chemistry) to the condition (DEVELOPMENT_LOG §3.13, decision 4c).

**Facts.** 151 donors post-QC from 20 studies. Healthy: 111 donors, 13 HLCA studies, no tumour tissue. Tumour:
38 LuCA donors (5 studies) + 2 HLCA (Thienpont). Paired tumour + adjacent: **14 donors from 4 studies** (Leader_Merad 5,
UKIM-V 5, Goveia_Carmeliet 3, Thienpont 1). Wu_Zhou (22 donors) is tumour-only. Assays are nearly aliased with atlas.
Sorting differs by study (Leader/UKIM-V ~90–97 % immune, Goveia ~50 % endothelial, Wu_Zhou 67 % epithelial).

**Options.**
1. *Pool everything, adjust for study as a covariate* — impossible for healthy-vs-tumour: study and condition are
   aliased (a covariate constant within condition is unidentifiable).
2. *Study as a random effect (mixed model)* — same aliasing; a random effect cannot separate a study effect from a
   condition effect when no study spans both conditions. Useful only where studies span conditions (the 4 paired studies).
3. *Restrict to one study* — e.g. Leader_Merad alone: 5 donors; too small, and results would not generalise.
4. **Within-donor inference, replicated across studies (chosen).** The pooled atlas is used for what pooling is good at
   — a shared embedding, a large reference for labels, CNV and deconvolution signatures — but every *inferential*
   claim rests on a contrast that is protected by design, and is checked for replication across the studies that
   contribute it:
   - **Primary**: paired tumour vs adjacent within the same patient (14 donors, 4 studies), pseudobulk summed per
     **donor × condition × cell type** (one unit per patient per tissue; multi-region samples no longer pseudo-replicates).
   - **Replication**: per-study paired fits (sign of log2FC per study) → "concordant in k of K studies", and
     **leave-one-study-out (LOSO)** refits (4) → a gene is called **robust** only if significant in the full fit, same
     sign in every LOSO fit, and not contradicted by any study with ≥ 2 donors. (Meta-analytic logic of Nguyen et al. 2023
     *Nat Commun* 14:1570 and Squair 2021: model donors, never pool cells; replication across batches > corrected data.)
   - **Eligibility**: a cell type is tested only with ≥ 5 paired donors from ≥ 2 studies; `Unknown` is never tested.
   - **Atlas-confounded contrasts are relabelled, not deleted**: `tumour_vs_healthy` →
     *exploratory (atlas-confounded)*, with a **negative-control contrast `adjacent_vs_healthy`** (adjacent-normal of
     LuCA patients vs healthy HLCA donors: no tumour, only atlas + study + chemistry + smoking/field effect). A gene in
     the exploratory contrast is shown as "supported" only if it is also significant with the same sign in the paired
     contrast. The malignant / CNV-low vs healthy-basal contrasts stay, labelled hypothesis-generating.
   - **Composition**: paired donors only, donor × condition units, proportions **within compartment** (immune subtypes
     as a share of immune cells; this neutralises CD45/CD31 enrichment differences between studies), paired model.
   - **Every page** shows n donors · n studies next to each result, and a study-grounded view (per-study replication).

**Why.** Pooling increases n only if the extra donors are exchangeable with respect to the contrast. They are not here:
every healthy donor differs from every tumour donor in atlas, study, protocol and often chemistry. Within-donor
contrasts cancel all of that exactly; replication across studies is the external-validity check that pooling was
supposed to provide. This is the standard of multi-cohort single-cell meta-analyses (e.g. Salcher 2022 used
patient-level tumour vs normal comparisons and study as a blocking factor; Squair 2021; Heumos 2023 §DE).

### 2.2 Other decisions

| # | Topic | Options | Choice | Why |
|---|---|---|---|---|
| D-V | Visium section | keep 6,195-spot NE section; drop spatial; **switch to 10x LUSC section (3,858 spots)** | switch + **marker gate** in `spatial_prep` (squamous markers must be detected; fails the step otherwise) | the current section is a different cancer; the squamous section is free (CC BY 4.0) and verified by markers. Reverses the 2026-11-09 decision, which assumed both files were the same tissue (they are not) |
| D-M | Mito QC | 3 MAD two-sided (configured); fixed 20 %; **upper-tail MAD within sample × compartment + ceiling** | upper tail, 5 MAD, MAD floor 1 pp, ceiling 25 %, MAD skipped where a sample has no mito signal | only high mito marks damage; both atlases were already author-filtered (LuCA ≤ 30 %), so a lenient rule; computing within compartment keeps metabolically active tumour cells (Yates 2025 *Genome Biol*); OSCA/sc-best-practices |
| D-A | Healthy age | exclude unknown ages; **row-wise fallback, exclude known < 45, keep unknown-age donors flagged** | as bold | dropping 59/111 healthy donors would gut the reference; healthy donors no longer enter any primary inference (§2.1), so unknown age only touches exploratory contrasts — stated on the page |
| D-D | Doublets | SOLO + Scrublet consensus; union; scDblFinder (R); **Scrublet only + per-sample failure guard** | Scrublet only; a sample where > 25 % of cells are called is a threshold failure → fall back to the top 5 % of scores; SOLO off | SOLO miscalibrated (argmax, balanced classes, rate falls with sample size); consensus was ≈ Scrublet anyway; scDblFinder (best in Xi & Li 2021) needs an R env — noted as future work |
| D-I | Integration benchmark | keep; drop label-guided methods; **mask the labels of the benchmark cells during scANVI training** | mask (held-out evaluation) + condition-matched mixing metric (tumour/adjacent cells only: Thienpont HLCA vs LuCA) | removes circularity without discarding scANVI; gives a fair label-transfer accuracy for the page. Geneformer kept (user-approved candidate) |
| D-C | Malignant calling | keep; lower quantile; per-sample runs; **per-cell score + two-axis rule, tumour/adjacent cells only** | inferCNV on tumour + adjacent samples only (LuCA + Thienpont; healthy HLCA excluded), reference = immune/endothelial/stroma cells *from the same samples*; per-cell score = mean squared CNV signal; per-cell correlation with its sample's consensus tumour profile; thresholds = 99th percentile of each axis in adjacent-normal epithelium (same atlas, same donors); malignant = epithelial ∧ tumour sample ∧ both axes above | Tirosh 2016 / Neftel 2019 two-axis rule; removes the atlas confound and the cluster-mean artefact; author labels become a validation, not a target |
| D-T | Tumour epithelium groups | keep Malignant + CNV-low separate everywhere; **merge for spatial / CCC / composition, keep separate for DE** | `Tumour epithelium` = Malignant ∪ CNV-low for cell2location, CCC and composition | the two differ by a threshold, not by expression; a deconvolution or communication model cannot separate them |
| D-P | Composition model | scCODA all samples; propeller (R); sccomp (R); **scCODA paired + CLR paired t-test** | both, on paired donors, donor × condition units, within compartment; scCODA with `condition + donor_id`; CLR + paired t-test (BH) as frequentist cross-check | R methods would need an R env; two independent models agreeing is stronger than one |
| D-Mi | Milo | per sample; **per donor × condition** | count neighbourhoods per donor × condition (28 units) | removes multi-region pseudo-replication |
| D-X | Tensor-cell2cell | 82 samples, outer fill; **14 paired donors × 2 conditions, paired test** | tensor on donor × condition LIANA (28 slices), cell types present in ≥ 90 % of slices, elbow rank search, paired Wilcoxon on loadings | condition contrast becomes testable; tumour epithelium enters via D-T. Elbow on a 28-slice tensor costs minutes, not 6.5 h |
| D-S | `ccc` step | one step; **split `ccc` (single-cell) + `ccc_spatial`** | split | the single-cell part no longer waits 11 h for deconvolution |
| D-c2l | cell2location outputs | q05; **posterior mean for proportions, q05 kept for presence** | as bold | q05 is a conservative lower bound; dividing lower bounds distorts proportions |
| D-Tr | Trajectory | keep; drop; **keep, labelled descriptive, per-donor summaries, scVI-latent sensitivity** | as bold | root/endpoint atlas confound cannot be removed (no LuCA basal cells); say so and show per-donor n |
| D-Pw | Pathways | all contrasts equal; **headline from paired only** | paired contrast headlines; exploratory contrasts flagged | inherits §2.1 |

---

### 2.3 Control groups (decided 2026-11-16, after run 3; revised after three reviews; run 4)

Full entry, with numbers and rejected options: DEVELOPMENT_LOG §3.13. In short:

* **No study that gives LUSC tumours also gives healthy donors.** So tumour vs healthy is always a comparison
  between studies, whatever the model.
* **The paired design already uses 15 of the 17 paired LUSC patients** in HLCA + LuCA. The other 2 are Smart-seq2
  (read counts) and excluded.
* **Main comparison unchanged: tumour vs the same patient's adjacent normal tissue, immune cell types only.**
  Pairing cancels what is the same for both samples of a patient (genotype, age, sex, smoking history, study, lab,
  chemistry). It does not cancel differences in how the two tissues were handled, tumour cells inside "adjacent"
  samples, or changes shared by both tissues. It measures tumour relative to adjacent tissue, not relative to
  healthy lung, and mostly in ever-smokers (14 of 15).
* **LUAD patients' adjacent tissue is not used**, neither as controls (unpaired, different disease, different sex,
  smoking, stage and site) nor for the LUSC- vs LUAD-adjacent "field-effect" check, which is dropped: neither of
  its outcomes could be interpreted or would change a result.
* **"Healthy" now means lung tissue from a person without cancer** (organ-donor, explant and autopsy lungs; run
  4: 87 donors selected, 81 after QC, 9 studies; the run-3 cohort had 116 selected, 108 after QC). Bronchoalveolar lavage (Janssen_2020, 10 donors; airway cells, not tissue) and
  tumour-free surgery tissue from cancer patients (Schiller_2020, Schiller_2021, Schultze_unpubl; 17 donors) are
  removed from the cohort by an allow-list in `select`, because neither has a job the remaining samples cannot do.
  Run 3's cross-study null had the lavage study in one set only, so part of its 14–73 % was a tissue difference;
  run 4 recomputed it on lung tissue only: 6.4–60.8 % across 12 cell types, so study effects still dominate. This changes the cohort, so run 4 is a complete rerun.
* **HLCA healthy donors are a reference** (shared embedding, label-transfer reference, CNV window mask, trajectory
  root, figure context). Every comparison against them is descriptive and carries the recomputed cross-study null.
* **Two checks of the main result, from existing data:** a sign-flip null for the paired test (labels swapped in a
  random half of patients; a cell type averaging > 3 % of genes at p < 0.01 is flagged), and a per-patient table of
  what each adjacent sample contains. Run 4: no type flagged (0.28–0.84 % per type); the adjacent samples are
  80–100 % immune in Leader_Merad and UKIM-V, 38–49 % endothelial in Goveia_Carmeliet, and UKIM-V's airway-rich
  adjacent epithelium has CNV call rates of 18–35 %. Both are on the DE page.
* **Limits stated next to every paired result:** 10–14 pairs per primary type; 2–4 studies with ≥ 2 pairs per type;
  Leader_Merad and UKIM-V give 10 of 15 patients and all Treg pairs; the pooled paired model is the test, per-study
  and leave-one-study-out results are checks; FDR of the "robust" subset is not formally controlled.
* **Hartung–Knapp random effects (§6) is closed, not done:** with 2–4 studies per cell type it cannot estimate
  between-study variance (with 2 studies the t value has 1 degree of freedom, 12.7).

## 3. Implementation + rerun plan

**Where it runs (hard rule 6).** Data and run stamps live in the main checkout
(`sc-spatial-analysis`, branch `v2.1-dashboard`). sc-deploy's `data/processed`, `results/tables`, `results/reports`
are already symlinks into it. Choice: **run the pipeline from the main checkout after syncing sc-deploy's changed
`src/`, `config/`, `scripts/`, `workflow/` into it** (a one-way rsync of tracked code files, recorded in the log);
results then appear in sc-deploy through the symlinks, and `export` is run in sc-deploy to write its own DB.
Rejected: symlinking `data/` wholesale into sc-deploy (would expose `data/raw` + `data/interim` under a second tree
and risk two run_all instances sharing stamps).

**Order.**
1. Archive run-2 outputs → `results/_archive_run2_2026-11-13/` (tables, reports, DB).
2. Code (sc-deploy): select (age, mito), qc (mito rule, Scrublet-only + guard), integrate (masked scANVI,
   matched mixing), cnv (two-axis per-cell), de (donor × condition, eligibility, LOSO, per-study, negative control,
   labels), abundance (paired, within compartment, CLR cross-check; Milo per donor × condition), trajectory
   (per-donor summary, scVI sensitivity), deconvolve (tumour-epithelium merge, posterior mean), spatial_prep
   (new section + marker gate), ccc split + paired tensor, report, export (new tables, run label).
3. Unit tests + subsample smoke test.
4. Senior reviewer on the diff → fixes.
5. Run: `fetch_visium` + `spatial_prep` → `select` → `qc` → `fm_embed` → `integrate` → `cnv` →
   in parallel {`de` → `ccc`} ‖ `abundance` ‖ `trajectory` ‖ {`deconvolve` → `spatial` → `ccc_spatial`} → `report` → `export`.
6. Pages rewritten while it runs (publication-standard figures, methods boxes, n donors · n studies, study-grounded views).
7. Final senior review of pages + numbers; figure/communication review; fixes; quality gates.

**ETA (estimates from run 1/2 timings):** select ~15 min, qc ~25 min (SOLO off), fm_embed ~2 h, integrate ~2 h,
cnv ~30 min, then the parallel block is set by cell2location (~2.3 h reference + ~5–6 h mapping on 3,858 spots) →
**~13–14 h wall time** from launch.

## 4. Senior review — verdicts and resolutions

### 4.1 Review #1 (on plan v1, 2026-11-14 evening)

**Verdict on the merge question:** "pool for representation, infer within donor" is the right framework; a study
random effect is *not* the answer (no study spans healthy and tumour; K = 4 studies cannot support a variance
component). The weak point was the replication machinery, and the plan under-stated scope: the paired design is
**15 patients from 4 studies** (Leader_Merad 5, UKIM-V 5, Goveia_Carmeliet 3, Thienpont 2 once region-unknown
tumour counts as tumour), and only **immune** cell types have enough paired cells — no epithelial or stromal type
can be tested within patient. Every tumour-epithelium claim stays atlas-confounded whatever the modelling.

| # | Required change | Resolution (code) |
|---|---|---|
| 1 | Region-unknown tumour is tumour tissue (15 pairs, not 14); say plainly that the paired design is immune-only | new obs `tissue` (`cnv.tissue_of`); all paired analyses (DE, composition, Milo, tensor) key on it; page text |
| 2 | Per-study replication from ONE joint model, not separate fits on 1–3 donors | `de.per_study_effects`: explicit design `~ donor + study:tissue` (shared dispersions) → per-study log2FC ± SE; `de.meta_summary`: inverse-variance fixed-effect estimate, Cochran's Q, I² (descriptive), k of K concordance; 1-pair studies are point estimates only |
| 3 | Redefine "robust" | `robust` = thresholded test (H0 \|log2FC\| ≤ log2 1.5, `padj_thr` < 0.05) ∧ no ≥ 2-pair study contradicting with a 95 % CI excluding 0 ∧ LOSO refits without each *major* study (≥ 25 % of pairs: Leader_Merad, UKIM-V) keep the sign with padj < 0.1; concordance is a graded column |
| 4 | Eligibility ≥ 6 pairs and ≥ 2 studies with ≥ 2 pairs each | `de.paired_eligibility`: tier `primary` / `descriptive` (≥ 3 pairs; fitted, shown, never headlined) |
| 5 | Subtype mix inside coarse labels (alveolar vs monocyte-derived macrophages; Treg vs conventional CD4) | `de.refine_groups` in config, applied in `make_de_groups` from author fine labels |
| 6 | Ambient/contamination flag | `de.ambient_flags`: ratio to the highest-expressing other type (≥ 10× → suspect) + fraction of other types where the gene also rises; suspects never headline |
| 7 | Drop scCODA with donor dummies; paired NB GLM within the immune compartment | `de.composition_paired`: NB GLM `count ~ donor + tissue`, offset log(immune total) (quasi-Poisson fallback) + arcsin-√ and CLR paired t-tests, BH each; ≥ 200 immune cells per unit; tumour epithelium / Unknown out of denominators; per-study direction count |
| 8 | Malignant calling: per-sample null, squamous-aware window mask, cross-fitted consensus | `cnv.run_infercnv` rewritten: per sample, same-sample reference with 20 % held out as null, window mask from normal epithelium (healthy basal + adjacent), cross-fitted consensus correlation, both axes > sample q99 of the null; "unresolved" < 50 epithelial cells; adjacent call rate + 3q/3p/5p arms + author labels as validation only |
| 9 | `adjacent_vs_healthy` is a decomposition term, not a negative control; add a real cross-study null | relabelled; `de.cross_study_null`: healthy HLCA studies of one chemistry split into two sets → share of genes p < 0.01 = cross-study floor |
| 10 | Mito code must match the plan; counts/genes MAD compartment-aware | all MAD rules per sample × compartment; MAD floor 1 pp; 5 MAD + 25 % ceiling in config; NaN count logged; `qc_per_compartment.csv` |
| 11 | Doublet fallback = expected rate, plus a marker-based heterotypic check | `qc.guard_scrublet` (0.8 % per 1,000 cells, cap 10 %); `qc.heterotypic_rate` (PTPRC + EPCAM/KRT8/18 or PECAM1) per sample |
| 12 | Tensor: no tumour-only types (they separate conditions by construction); drop slices, not types | `ccc.balanced_tensor_scope`: types present in ≥ 80 % of slices of each tissue; slices missing one are dropped with their pair; paired Wilcoxon + BH |
| 13 | c2l proportions = posterior mean of per-sample ratios; Basal + tumour sum; n = 1 stated | `deconvolve.posterior_proportions` (500 posterior samples) → `c2l_proportions.csv`; spatial page text |
| 14 | Run mechanics: stamps, wipe old Visium files, ccc split in Snakefile, memory/threads | `--force` for every step; per-section raw folder `data/raw/visium/<sample_id>/`; Snakefile is generated from `steps.py` (split is automatic); concurrent steps get capped threads |
| 15 | Global BH, shrunken LFC + CIs, power, sex-swap check | `padj_global`; `ci_low/ci_high`, `log2FoldChange_shrunk`; `de_power.csv` (Wald MDE at 80 % power); `de_sex_check.csv` (XIST vs RPS4Y1 per unit) |

### 4.2 Review #2 (on the code diff, before launch)

| Item | Finding | Resolution |
|---|---|---|
| B1 | Completion stamps + run-2 objects would be silently reused without `--force` | run launched with `--force`; run-2 tables archived first (`results/_archive_run2_2026-11-13/`) |
| B2 | Author-malignant cells in units too small to score were labelled "CNV-low" | `cnv_low` requires `cnv_status == "not called"`; unscored cells → `Tumour epithelial (CNV unresolved)` (in `Tumour epithelium` for c2l/CCC/composition, never tested in DE) |
| B3 | `robust` passed silently when the per-study model failed | rank + residual-df guard; `per_study_available`; `robust` requires it |
| B4 | CNV null/thresholds fragile; most adjacent samples unresolvable | CNV resolved per **patient × tissue** (regions pooled); held-out null ≥ 100 cells; `min_reference` 300; **pre-registered acceptance gate** in config (adjacent call rate ≤ 5 %, calls in ≥ 10 patients, called cells show 3q gain and 3p loss) — recorded in `cnv_summary.acceptance_passed`, shown on the pages |
| NB1 | LuCA `pct_counts_mito` is exactly 0 for 17/19 studies | studies with an all-zero mito share → NaN, `mt_source = "unavailable …"` (ceiling cannot apply; stated on the QC page) |
| NB2 | the new `tissue` column overwrote CELLxGENE's anatomical `tissue` | renamed `tissue_type` everywhere |
| NB3 | `padj_global` over all tiers, plain Wald p; `supported_by_paired` used `padj` | global BH over the primary tier, also for the thresholded test (`padj_thr_global`); support uses `padj_thr` |
| NB4 | LOSO/per-study for descriptive tiers | primary tier only |
| NB8 | c2l proportion fallback invisible | `proportion_method` column |
| NB9 | tensor could end with < 5 patients silently | dropped patients logged; raises below `tensor_min_donors` (5) |
| NB5/6/10/11 | thread oversubscription; by-sample LIANA on all samples; NB→quasi-Poisson fallback rate; Snakefile threads | accepted (slower, not wrong); fallback method is a column of `composition_paired`; per-sample LIANA kept for the dashboard browser |

Optional items taken: Milo merges region-unknown tumour (via `tissue`); per-donor trajectory summaries + scVI
sensitivity; label-masked scANVI (random benchmark cells — the reviewer notes whole-study masking would be
stricter; recorded as open). Not taken (recorded in §6): NMF meta-programmes (Gavish 2023), H&E nuclei counts for
`n_cells_per_location`, a log(n_cells) covariate sensitivity, Hartung–Knapp random-effects for top genes (later closed as uninformative at 2–4 studies, §2.3).

### 4.3 Final reviews (after run 3, on the rendered pages and the numbers)

Two independent reviewers read every rendered page against the database: the senior statistical reviewer (final
pass) and a figure / science-communication reviewer (Nature-figure standard). Findings and resolutions:

| Finding | Resolution |
|---|---|
| Failed CNV gate visible only on the Malignant page (both reviewers, blocker) | issue 16 also shown on DE and Trajectory (`ui.ALSO_SHOWN_ON`); headline sentences, the trajectory banner and the tumour-vs-basal badge say "provisional"; issue 16 text corrected (the paired DE holds no tumour epithelium) |
| Composition headline truncated to its header on Start/Overview | header folded into each result line |
| "Fixed in run 2" hard-coded for run-3 fixes | label derived from each issue's status (`ui.fixed_in`) |
| Generic palette reused tissue/atlas/direction hues; region-unknown tumour colour-only | new validated `SERIES` without semantic hues; `MALIGNANT_COLOR`/`BASAL_COLOR`; region-unknown tumour = lighter tint + own label; sex chart off the direction poles |
| Milo figure empty (0 significant) | strip plot of every neighbourhood's log2FC per type with the whole-type estimate overlaid; ρ computed on the page |
| "FDR ≤ 0" legend; legends that describe instead of state; numbering gaps; UMAPs without axis titles | fixed; Cohort, QC, Overview and DE figures numbered at render time; titles state the finding, computed from the data where they quote it |
| Robust scope mixed per-type and global FDR | one headline rule everywhere: robust ∧ `padj_thr_global` < 0.05 ∧ not a contamination suspect |
| IL2RA in macrophages headlined (Treg expansion) | export flags `ambient_expanding_source` (marker ratio ≥ 3 and the source type expands in the paired composition test) and folds it into `ambient_suspect` |
| Sex check counted unpaired donors; zero-zero units called "female" | page reports paired patients and recorded-vs-inferred agreement; `sex_check` returns "undetermined" when both markers are absent (test) |
| "No biology differs" in the cross-study null; is the unpaired test calibrated? | wording: "no tumour; study and cohort differences only"; new `de.within_study_null` (random halves of one study): ≈ 1.4 % at p < 0.01 → calibrated; shown in the DE figure |
| Mito filter silently absent for three LuCA studies | QC warning with the studies and cell count; issue 9 → partly resolved |
| CCC headline driven by neutrophils (descriptive, one study) | headline and default figure restricted to primary-tier sender/receiver pairs; others tagged "descriptive" |
| Malignant results may be driven by UKIM-V | sensitivity without UKIM-V (§5.1): unchanged |
| Composition 14 vs Milo 15 patients unexplained | dropped patient named on the page (< 200 immune cells in one tissue) |
| Power figure overstates sensitivity for the robust set | legend states the thresholded detection limit is ≈ 0.585 log2 higher |
| Smaller items (jargon, units, long donor IDs, "Unlabelled" vs "Unknown", study counts, glossary CI, `np.int64` reprs, retired issue 7) | fixed |

## 5. Before / after (run 2 → run 3)

Run-2 values from `results/_archive_run2_2026-11-13/`; run-3 values from the run-3 tables / DB.

| Quantity | Run 2 | Run 3 | Why it changed |
|---|---|---|---|
| Cells after `select` | 249,402 | 246,395 | 3 healthy donors aged < 45 (unknown-age bug) removed |
| Cells tested / kept by QC | 246,350 / 231,504 | 243,318 / 223,406 | mito filter active; compartment-aware MAD; Scrublet only |
| Cells flagged for high mito share | 0 (filter inactive) | 6,507 | §3.10 #9 fixed |
| Doublets removed | 672 (Scrublet ∧ SOLO) | 1,123 (Scrublet, no threshold failures) | SOLO dropped |
| Samples without mito data | 215 of 215 (all 0) | 33 of 212 (three LuCA studies) | computed before the merge |
| Malignant cells (patients with ≥ 10) | 3,942 (7 of 38) | 16,027 (28 patients); **pre-registered gate failed** (adjacent call rate 7.9 % > 5 %) | per-patient two-axis calling; see §5.1 |
| Paired DE: testable cell types | 15 (≥ 3 samples) | 17 fitted: 10 primary (all immune) + 7 descriptive | patient units, stricter eligibility |
| Paired DE: significant gene × type rows (padj < 0.05) | 8,474 | 5,363 (3,307 primary tier) | patient units; subtype split |
| … of which robust | — | 126 in 9 cell types (headline rule: robust ∧ global FDR ∧ not a contamination suspect; 174 robust before the last two filters, 26 of them contamination suspects) | new: thresholded + replicated |
| Composition | scCODA, 9 credible changes (all samples) | 5 of 14 immune types change consistently (14 patients): Treg, B, plasma up; alveolar macrophages, non-classical monocytes down | paired NB GLM within immune + 2 paired tests |
| Milo significant neighbourhoods | 1,671 of 6,225 | 0 of 6,208 (directions agree with composition, Spearman ρ = 0.96) | donor × tissue units, 15 pairs: neighbourhood-level power is low |
| Visium section | 6,195 spots, neuroendocrine | 3,858 spots (3,816 after QC), squamous | wrong section in runs 1–2 |
| cell2location signatures / niches | 25 / 9 | 26 / 8 (tumour epithelium 43 % of a spot on average) | new section, merged tumour epithelium |
| Tensor | 82 samples, rank 6, no factor differs (p ≥ 0.34, unpaired) | 11 patients × 2 tissues, 4 factors; none differs (smallest padj 0.17) | paired slices, tested within patient |
| scIB choice | scANVI (scored on its own training labels) | scANVI, total 0.668, kNN purity 0.92 (held-out labels) | label masking |
| Cross-study technical null | — | 14–73 % of genes at p < 0.01 between two sets of healthy studies | new: the floor under every cross-study comparison |

### 5.1 The CNV acceptance gate failed — what that means

The four criteria were written into `config.yaml` (`cnv.acceptance`) before run 3 started. Three pass: calls in
29 patients (≥ 10 required), mean chr3q gain in called cells, mean chr3p loss. The false-positive check fails:
**7.9 %** of adjacent-normal epithelium is called CNV-high against a limit of 5 %. Located post hoc, the excess
is one study: UKIM-V's three adjacent samples run 17.6–35.2 %, while the six adjacent samples of the other
studies pool to 3.2 % (27 of 845 cells). The threshold was **not** moved and UKIM-V was **not** dropped after the
fact. Calls in adjacent tissue can be real (tumour cells or field cancerisation in a sample taken near the
tumour) or artefact, and these data cannot separate the two.

Consequences, as shipped: the Malignant cells page shows the gate table and an error banner plus the post-hoc
per-study table (labelled as such); known issue 16 (open) is attached to every result that uses the calls
(malignant-vs-healthy-basal DE, trajectory endpoint, agent tools); issue 11 is "partly resolved". The **primary
paired DE does not use the calls** — tumour epithelium is pooled into one group whatever its CNV status — and
neither do composition, CCC or cell2location. **Sensitivity (final review):** refitting malignant-vs-healthy-basal DE without UKIM-V (25 instead of 28 malignant units) leaves the log2FC ranking unchanged (Spearman ≈ 1.00; 593 vs 597 significant genes, Jaccard 0.84; the top 200 genes keep their sign), and the malignant median pseudotime moves from 0.124 to 0.125 — the study that fails the gate does not drive the malignant-dependent results. Agreement with the LuCA authors' labels (tumour tissue): 15,723 of
29,151 author-malignant cells called (54 %), 257 of 2,101 author-non-malignant cells called (12 %).

### 5.2 Bugs found after the run, fixed, rerun

| Bug | Effect | Fix | Rerun |
|---|---|---|---|
| `composition_paired` built the design from a *categorical* `donor_id`: `get_dummies` emitted all-zero columns for donors outside the test | singular design → the NB GLM never converged and every type silently used the quasi-Poisson fallback; `k_studies` counted 19 categories instead of 4 | cast to strings; NB fitted by BFGS from Poisson start values (Newton returned NaN even on a valid design; BFGS matches Newton wherever Newton works); regression test | `abundance` |
| decoupler's pseudobulk drops obs columns that vary within a unit; two Leader_Merad patients pool 10x 3′ v2 and 5′ v1 libraries | `assay` vanished from every pseudobulk → the unpaired designs lost their chemistry covariate and the cross-study null returned nothing; export then failed on the empty table | restore design covariates as the "+"-joined set of values; regression test; export skips empty *optional* tables | `de` (paired contrast identical: max \|Δ log2FC\| = 0), `export` |
| an unquoted `:` in `dashboard_content.yaml` | export failed | quoted | `export` |
| `lineage_cells` has no `study` column | Trajectory page crashed on hover data | hover fields optional | — |
| agent tools still read scCODA and run-2 group names | MCP tools and eval set stale | `composition_change` → paired table; DE tools return tier/CI/robust and attach issues 16/17; eval set updated (+ q25 on the CNV gate) | — |

### 5.3 What run 3 says (all numbers from the DB)

* **Within-patient immune changes are robust and specific.** 126 gene × type changes pass the thresholded test
  (> 1.5-fold, global FDR 5 %), survive per-study and leave-one-study-out checks and are not contamination
  suspects, in 9 of 10 primary-tier immune types; the
  largest are classic tumour-microenvironment programmes — *CXCL13*, *HAVCR2*, *ENTPD1* in CD8 T cells
  (tumour-reactive, exhausted), *CXCL13* in conventional CD4 T cells, *SPP1* and *MMP12* in macrophages,
  *TGM2*/*SPP1* in classical monocytes. 26 robust rows are contamination suspects and never labelled:
  immunoglobulin genes in myeloid cells, and Treg-marker genes (IL2RA in macrophages; TIGIT, CTLA4, TNFRSF18 in
  conventional T cells) whose source type expands in tumour — for the T-cell genes real exhaustion biology is
  equally possible, so the flag means "cannot exclude carry-over", not "artefact".
* **Immune composition shifts within patients:** regulatory T, B and plasma cells up; alveolar macrophages and
  non-classical monocytes down — consistent across three models and in the same direction in 3–4 of 4 studies.
* **Cross-study comparisons are dominated by study effects:** two sets of *healthy* studies differ in 14–73 % of
  genes at p < 0.01, while random halves of one study differ in ≈ 1.4 % (≈ 1 % expected) — the unpaired test is
  calibrated, so the excess is study, not method (`de_within_study_null`, added after the final review). Tumour-vs-healthy and malignant-vs-healthy-basal results are therefore shown as
  hypothesis-generating only (issue 17).
* **Neighbourhood-level (Milo) and communication-tensor tests are underpowered** with 15 / 11 paired patients;
  both point the same way as the cell-type-level results but neither passes FDR.
* **The spatial section is now a squamous tumour** (KRT5/TP63/SOX2 in 72–83 % of spots, CHGA in 3 %); one section
  = one patient.

## 6. What is still open

| Item | Why it is open | Suggested next step |
|---|---|---|
| CNV false-positive check failed (issue 16) | adjacent call rate 7.9 % > 5 %, UKIM-V-driven | inspect UKIM-V adjacent cells (author labels, CNV profiles, 3q/3p arms); if they are tumour cells, document and keep; if artefact, pre-register a per-study reference rule for run 4 — do not tune on these data |
| Paired design is immune-only (issue 12) | adjacent tissue has few epithelial/stromal cells | add a dataset with matched tumour/normal epithelium (e.g. whole-tissue LUSC cohorts) before any epithelial claim |
| Whole-study label masking for scANVI | random benchmark cells are masked; whole-study masking is stricter (reviewer note) | rerun the benchmark with one study's labels masked at a time |
| NMF meta-programmes (Gavish 2023) | malignant heterogeneity is described only by CNV and pseudotime | per-patient NMF on tumour epithelium, cluster programmes across patients |
| scDblFinder as second doublet caller | SOLO dropped; Scrublet alone (with guard) | add scDblFinder (R) in a separate env and compare per sample |
| cell2location `N_cells_per_location` from H&E | set to a prior, not measured | nuclei segmentation of the H&E image → per-spot counts |
| Milo / tensor power | 15 / 11 patients | more paired patients; until then read both as descriptive |
| More Visium sections | n = 1 patient (issue 15) | add public LUSC Visium/Xenium sections |
| Trajectory stability (run 4) | scVI sensitivity ρ fell from 0.83 (run 3) to 0.59; the scANVI ordering itself reproduces run 3 (ρ 0.95), so the rerun scVI map orders cells differently (DEVELOPMENT_LOG §3.13 "Run 4 results") | repeat the sensitivity over several scVI seeds and roots to see which ρ is typical; until then read per-patient pseudotime as a rough picture; do not tune |

## 7. Redeploy steps (you run these; nothing was pushed, published or deployed)

Everything is uncommitted in the `sc-deploy` worktree (branch `v2-dev`). The database schema moved to **5**, so app
code and database must go live together.

1. **Review locally**: `ssh -L 8501:localhost:8501 <vm>` → http://localhost:8501 (served from tmux `lusc-v2:app`
   with the run-4 DB). The run label in every page footer must read "Run 4".
2. **Commit** in `/home/abalasubramanian/storage/sc_spatial/sc-deploy` (suggested split: pipeline + config + tests;
   app + content; docs). `results/`, `data/` and `logs/run_all/` are git-ignored.
3. **Rehearse the cloud install**: `pixi run -e dashboard cloud-dry-run` (clean clone, fresh Python 3.12 venv, all
   pages, wrong-sha negative control). It must print `CLOUD DRY RUN PASSED`.
4. **Publish the database**: `export HF_TOKEN=…; pixi run -e dashboard publish-dataset`. It prints the three
   secrets lines (repo, **new revision**, **new sha256**).
5. **Push**: `git push origin v2-dev`.
6. **Streamlit Community Cloud** → app → *Settings → Secrets*: replace `LUSC_DB_HF_REVISION` and
   `LUSC_DB_SHA256` with the printed values; *Reboot app* (the cached resource keeps the old DB until reboot).
7. Open https://lusc-atlas.streamlit.app/ and check the footer says Run 4, the sidebar shows the Evidence and
   Context groups, and the Spatial page shows the squamous identity check.
