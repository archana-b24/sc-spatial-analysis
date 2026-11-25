# Learning notes: single-cell + spatial analysis, taught through this project

*A study guide written alongside the run-3 overhaul ([OVERHAUL_RUN3.md](OVERHAUL_RUN3.md)). Each section follows
the pipeline: what the stage does, why it exists, the method and its alternatives, the statistical traps we fell
into or avoided, and what to read next. Numbers quoted here are illustrative; the dashboard holds the current
ones.*

---

## 0. The one idea that organises everything: what is the unit of replication?

A single-cell dataset of 200,000 cells from 40 patients has **40 independent observations of "a patient"**, not
200,000. Cells from one patient share genotype, smoking history, surgery, dissociation protocol and sequencing
batch. Treating cells as independent replicates makes p-values absurdly small — the classic *pseudoreplication*
error (Squair et al. 2021 *Nat Commun*; Zimmerman et al. 2021 *Nat Commun*).

Every inferential choice in this project follows from asking "what is the replicate here?":

| Question | Unit of replication | How the project implements it |
|---|---|---|
| Does gene X change in tumour macrophages? | patient (one tumour + one adjacent per patient) | pseudobulk = sum counts per **donor × tissue × cell type**, then a negative-binomial GLM |
| Does the share of Tregs change? | patient | counts per donor × tissue within the immune compartment, paired GLM |
| Is a cell malignant? | cell, but judged against a null from *the same sample* | per-sample inferCNV with held-out reference cells |
| Do niches exist in the tissue? | one section = one patient | descriptive only; no p-values across spots generalise to patients |

Run 2 still had a subtle version of the error: a patient with 8 tumour *regions* contributed 8 tumour "replicates".
Run 3 sums regions per patient (DEVELOPMENT_LOG §3.12).

---

## 1. Building the cohort: merging atlases without fooling yourself

**What.** Two curated atlases from CELLxGENE: HLCA (Sikkema 2023, healthy lung, ~2.3 M cells from many studies)
and LuCA (Salcher 2022, lung cancer, many studies). LUSC patients come almost entirely from LuCA; healthy donors
entirely from HLCA.

**Why merge at all.** HLCA alone has 3 LUSC patients from one study (n = 3, perfectly confounded with study).

**What can go wrong (and the checks):**
- *Duplicated cells/patients.* LuCA redistributes studies HLCA also holds. Use CELLxGENE `is_primary_data` and a
  study crosswalk; never match by donor ID strings (they differ between atlases).
- *Gene space.* Join on Ensembl IDs, not symbols. Taking the **intersection** of genes is safe for integration but
  silently deletes genes one atlas lacks — here LuCA ships no mitochondrial (`MT-`) genes, so the merged matrix had
  none and the mitochondrial QC filter computed 0 % for every cell for two runs. **Lesson:** compute per-cell QC
  metrics on each source's *full* gene set before intersecting, and write a test that fails if a QC metric is
  constant.
- *Metadata parsing.* A fallback that only triggers when a whole column is empty (`if col.notna().any(): return col`)
  hides row-level gaps; `NaN < 45` is `False`, so missing ages silently passed an age filter. Fall back row by row
  and make the policy for missing values explicit.
- *Confounding by design.* No study contains both healthy and tumour tissue. Any healthy-vs-tumour comparison is
  simultaneously HLCA-vs-LuCA, study-vs-study, chemistry-vs-chemistry and protocol-vs-protocol.

**The design answer (§2.1 of OVERHAUL_RUN3).** *Pool for representation, infer within donor.* The pooled atlas is
used for the shared embedding, label transfer, CNV references and deconvolution signatures. Inference uses only
contrasts protected by design — tumour vs adjacent normal **within the same patient** — and checks that results
replicate across the studies contributing patients.

Why not "just add study as a covariate"? A covariate that is constant within each condition is *aliased* with the
condition: the model cannot separate them (the design matrix is rank-deficient). Why not a random study effect?
It can be fitted (condition then acts like a study-level variable, as in a cluster-randomised trial), but it is only
right if healthy and tumour studies are otherwise alike. Here they differ systematically (atlas, sorting, chemistry),
so the model would credit study differences to the tumour. And with 2–4 studies per cell type a between-study
variance cannot be estimated.

*Read:* Luecken et al. 2022 *Nat Methods* (integration benchmark); Sikkema et al. 2023 *Nat Med* (HLCA); Salcher
et al. 2022 *Cancer Cell* (LuCA); Heumos et al. 2023 *Nat Rev Genet* (best practices).

---

## 2. Quality control

**What.** Remove barcodes that are not one intact living cell: damaged cells (low counts, high mitochondrial
share), empty droplets/debris (few genes, dominated by a few transcripts), doublets (two cells, one barcode).

**Per-sample, MAD-based thresholds.** Samples differ in depth by an order of magnitude across chemistries, so a
global cut-off removes whole shallow samples. A cell is an outlier if it is more than *k* median absolute
deviations (MAD) from *its own sample's* median (sc-best-practices; *k* = 5 is permissive).

**Refinements learned here:**
- Mitochondrial share is **one-sided**: only *high* mito indicates damage. A two-sided rule removes healthy
  low-mito cells. When a sample has no mito signal (MAD = 0 — the depositors removed MT reads), a MAD rule flags
  every non-zero cell; skip it and keep only a hard ceiling.
- Judge within **sample × compartment**: neutrophils have intrinsically low UMI counts and tumour/epithelial cells
  more mitochondrial RNA (Yates et al. 2025 *Genome Biol*); a sample-wide rule removes them preferentially.
- Floor the MAD (here 1 percentage point) so a tight distribution does not create many "outliers".

**Doublets.** Scrublet simulates doublets by adding random pairs of cells and scores each real cell by how many
simulated doublets surround it; a threshold is chosen per sample from the bimodal score histogram. SOLO (scVI-based)
trains a classifier on real vs simulated doublets. Here SOLO called 17 % of cells (argmax of a classifier trained on
balanced classes ignores the expected rate of ~1 %), so it was dropped. A Scrublet threshold can also fail (one
sample: 97 % called) — guard it with an expected-rate fallback. scDblFinder (R) is the benchmark winner (Xi & Li
2021 *Cell Syst*; Germain 2022).

*Read:* Heumos 2023; Osorio & Cai 2021 *Bioinformatics* (mito thresholds); McGinnis 2019 (DoubletFinder); Wolock 2019
(Scrublet); Bernstein 2020 (SOLO).

---

## 3. Integration: a shared space, not corrected data

**What.** Place cells from all samples in one low-dimensional space where the same cell type from different
studies sits together, so neighbours, clusters, UMAPs, label transfer and trajectories work across studies.

**Methods compared.** PCA (no correction, the floor); Harmony (iteratively adjusts PCA coordinates per batch);
scVI (variational autoencoder on raw counts with batch as a covariate); scANVI (scVI + cell-type labels,
semi-supervised); sysVI (scVI variant for *substantial* batch effects, VampPrior + cycle consistency); scArches
(map the query onto a published reference model — here HLCA's); Geneformer (a transformer foundation model used
zero-shot).

**How the choice is made.** scIB metrics (Luecken 2022): bio-conservation (do labels stay separate?) vs batch
removal (are batches mixed within a label?), plus an *over-correction guard* (kNN label purity must not drop below
unintegrated PCA).

**Traps.**
- *Circularity.* scANVI is trained with the labels used to score it. Run 3 hides the labels of the benchmark cells
  during training, so scANVI is scored on cells it never saw labelled.
- *Batch metrics on confounded data.* "Mix HLCA and LuCA cells of the same type" rewards erasing tumour-vs-healthy
  biology when atlas = condition. Compare like with like (tumour-patient tissue from both atlases).
- *Corrected values are not counts.* Never run DE on integrated values; DE uses raw counts with the design.

*Read:* Lopez 2018 *Nat Methods* (scVI); Xu 2021 *Mol Syst Biol* (scANVI); Korsunsky 2019 (Harmony); Hrovatin 2024
(sysVI); Lotfollahi 2022 (scArches); Theodoris 2023 *Nature* (Geneformer); Kedzierska 2025 (foundation-model
zero-shot limits).

---

## 4. Calling malignant cells from expression: inferCNV done carefully

**What.** Tumour cells carry copy-number changes (whole arms gained or lost). Averaging expression over ~100–250
neighbouring genes along a chromosome, relative to a diploid reference (immune/stromal cells), reveals those
changes as smooth bumps (Patel 2014; Tirosh 2016).

**Two-axis rule (Tirosh 2016; Neftel 2019).** A malignant cell has (a) a strong CNV signal *and* (b) a profile that
correlates with its tumour's consensus profile (clonal). Normal cells can have noise on (a) but not a clonal pattern.

**Traps found here.**
- `infercnvpy.tl.cnv_score` returns the **cluster mean**, not a per-cell score. A threshold on it is a threshold on
  ~30 cluster values. Compute per-cell scores yourself.
- The **null** must match the query in depth, chemistry and atlas: use held-out reference cells *of the same
  sample*, not healthy cells from another atlas.
- Expression programmes whose genes cluster on a chromosome (keratins on 12q13/17q21, the epidermal differentiation
  complex on 1q21) look like CNVs in squamous cells. **Mask** windows where *normal* epithelium deviates.
- The consensus profile must be **cross-fitted** (built from the other half of the cells), or each cell is compared
  with a profile it helped build.
- Author labels are not ground truth (they are themselves inferCNV-informed); use them, adjacent-tissue call rate
  (should be ~1 %) and LUSC hallmark arms (3q gain, 3p loss, 5p gain) as validation.

*Alternatives:* CopyKAT, SCEVAN (R; expression-only); Numbat (allele-aware, needs BAMs — best when available).

---

## 5. Differential expression: pseudobulk, pairing and replication

**Pseudobulk.** Sum raw counts per *donor × tissue × cell type*, then model with a negative binomial GLM (DESeq2 /
edgeR / limma-voom). This respects the unit of replication and has the best calibration in benchmarks (Squair 2021;
Crowell 2020 *muscat*; Junttila 2022).

**Paired design.** `~ donor + tissue`: each patient is its own control, so everything constant within a patient
(genotype, age, sex, smoking, study, atlas, chemistry) cancels. What differs *between the two tissues* does not:
harder digestion of dense tumour, blood in the adjacent lung, each sample's ambient RNA, and tumour cells inside an
"adjacent" sample (which shrinks real differences).

**Study-grounded replication (run 3).**
- One joint model `~ donor + study:tissue` gives each study its own tumour effect with shared dispersions →
  per-study log2FC ± SE → inverse-variance (fixed-effect) meta-estimate, Cochran's Q and I² for heterogeneity.
  (Do not fit separate DESeq2 models on 1–3 patients: dispersion cannot be estimated.)
- Leave-one-study-out refits: the result must survive dropping either of the two large studies.
- A thresholded test (H0: |log2FC| ≤ log2 1.5) instead of "different from exactly zero": with enough cells
  everything is "significant"; the question is whether the effect is big enough to matter.
- Eligibility: ≥ 6 pairs from ≥ 2 studies with ≥ 2 pairs each; otherwise "descriptive".

**Other traps.**
- *Composition inside a label.* "Macrophage" in adjacent lung is mostly alveolar macrophages; in tumour mostly
  monocyte-derived. DE on the coarse label measures the subtype shift. Split the label.
- *Ambient RNA.* Highly expressed transcripts of an abundant type (immunoglobulins from plasma cells, surfactant
  from AT2) leak into every droplet. A "tumour-up" gene in macrophages that is 50× higher in another cell type is a
  contamination suspect. CellBender/SoupX need raw droplets; without them, flag.
- *Multiple testing scope.* DESeq2 adjusts per cell type; report a global BH across all tests of the headline
  contrast as well.
- *Power.* With 6–15 pairs, only large effects are detectable; state the minimum detectable effect.

**Pathways.** decoupler (Badia-i-Mompel 2022) scores PROGENy (pathway footprints; Schubert 2018), CollecTRI
(TF regulons; Müller-Dott 2023) and MSigDB Hallmark gene sets on the DE statistic — inherit the DE design's
validity, so headline pathways come from the paired contrast.

---

## 6. Composition: proportions are compositional and sorting-dependent

**Why special.** Proportions sum to 1: if tumour cells expand, every other share falls even if nothing else
changed. And studies enrich different compartments (CD45+ sorting, CD31+ sorting), so cross-study proportions of
"all cells" compare protocols.

**Run 3 approach.** Paired patients; proportions *within* the immune compartment (immune subtype / all immune
cells); three views: NB GLM `count ~ donor + tissue` with offset log(immune total) (the edgeR-style model behind
Milo and propeller), arcsin-√ proportions + paired t-test (propeller; Phipson 2022 *Bioinformatics*), CLR +
paired t-test (compositional data analysis). A change is "consistent" only if all three agree.

**Why not scCODA here.** scCODA (Büttner 2021 *Nat Commun*) is a Bayesian compositional model with a reference cell
type and a spike-and-slab prior. With donor as 14 dummy covariates, the sparsity prior shrinks donor effects and
leaks donor variance into the condition effect.

**Milo** (Dann 2022 *Nat Biotechnol*) tests abundance changes in overlapping kNN *neighbourhoods*, so it sees states
within a type that move in opposite directions. Count per donor × tissue, design `~ donor + tissue`.

---

## 7. Trajectories: distance, not lineage

Palantir (Setty 2019) orders cells by diffusion distance from a root cell; CellRank 2 (Weiler 2024) combines
kernels (pseudotime, CytoTRACE developmental potential) into fate probabilities. In human tumours this is a
*transcriptional distance from normal basal cells*, not lineage tracing (malignant cells are patient-clonal; there
is no time course). Here the root (healthy basal) is HLCA and the endpoint (malignant) LuCA, so atlas distance is
part of the axis: report per patient, check sensitivity to the embedding (label-guided scANVI vs label-free scVI).

---

## 8. Spatial transcriptomics (Visium) and deconvolution

**Check tissue identity first.** The 10x section used in runs 1–2 was labelled "lung cancer"; markers showed a
neuroendocrine carcinoma (chromogranin A ≈ 10 % of every spot's RNA; no KRT5/TP63/SOX2). A one-line marker gate
would have caught it on day one. Run 3 uses 10x's squamous cell carcinoma section and the gate is in code.

**Deconvolution.** A Visium spot (55 µm) holds ~1–10 cells. cell2location (Kleshchevnikov 2022) learns reference
signatures per cell type from scRNA-seq (with batch covariates), then models each spot's counts as a mixture.
Report the posterior mean abundance; compute proportions per posterior sample and average them (a ratio of means is
not the mean of a ratio); use q05 for "confidently present". Signatures that differ only by a threshold (CNV-high vs
CNV-low tumour cells) cannot be separated — merge them.

**Niches.** CellCharter (Varrone 2024) clusters spots by their own and their neighbours' composition; squidpy
(Palla 2022) gives neighbourhood enrichment and Moran's I for spatially variable genes. One section = one patient:
everything is descriptive.

---

## 9. Cell–cell communication

LIANA+ (Dimitrov 2024 *Nat Cell Biol*) runs several ligand–receptor methods and aggregates their ranks per sample.
Differential CCC combines the paired DE statistics of ligand and receptor. Tensor-cell2cell (Armingol 2022
*Nat Commun*) decomposes a sample × LR pair × sender × receiver tensor into programmes.

**Traps.** Communication is *predicted from mRNA*, not measured binding. A tensor with zero-filled missing cell
types encodes sorting protocol; a tumour-only cell type separates conditions by construction. Run 3 builds the
tensor on paired patients only, with cell types present in both tissues, and tests factor loadings *within*
patient (paired Wilcoxon). CCC method agreement is modest (Dimitrov 2022 *Nat Commun*): treat single LR pairs as
hypotheses.

---

## 10. Checklist I would now apply to any atlas project

1. Write down the unit of replication for every question before computing anything.
2. Tabulate condition × study × atlas × chemistry; mark every aliased contrast as exploratory.
3. Compute QC metrics on each source's full gene set; assert no metric is constant.
4. Verify tissue identity of every external dataset with markers.
5. Read the method's source for what a returned "score" really is (cluster mean vs per cell).
6. Matched nulls (same sample, same chemistry), never a null from another batch.
7. Effect sizes with uncertainty, thresholded tests, global multiple-testing scope, power statement.
8. Replication across studies and leave-one-study-out before any headline.
9. Flag contamination (ambient) and composition-within-label explanations for DE hits.
10. Put n patients and n studies next to every number on every figure.
11. Pre-register acceptance criteria for any classifier you build, and keep the verdict when it fails.
12. Never trust a fallback silently: record which branch ran, and explain every empty output.

## 11. What run 3 taught, in the data themselves

Four lessons came out of the finished run rather than the plan. Each is worth remembering because each
was invisible until the outputs were read.

**A pre-registered check that fails is a result, not a nuisance.** The malignant-cell caller had four
acceptance criteria written into the config before the run. Three passed: calls in 29 patients, and called
cells show 3q gain and 3p loss on average, which is the LUSC signature. The false-positive check failed: 7.9 %
of adjacent-normal epithelium was called, against a 5 % limit. It is tempting to move the limit or drop the
study that causes it (UKIM-V; the others pool to 3.2 %). Either move turns a test into a description. The honest
response is to keep the verdict, show where the excess comes from labelled as post hoc, and make every
downstream result say "provisional". The design then pays off: the primary DE pools all tumour epithelium, so
it never depended on the calls.

**Measure the null instead of arguing about it.** Splitting *healthy* donors from one chemistry into two sets of
studies gives a comparison in which nothing biological differs. In this atlas 14–73 % of genes reach p < 0.01
in that comparison, depending on cell type. That one table is the quantitative case for the paired design:
any tumour-vs-healthy "hit" has to clear a floor that high before it means anything.

**Neighbourhood tests and cell-type tests answer different questions at different power.** With 15 paired
patients, the cell-type-level composition test finds five clear shifts. Milo, splitting the same cells into
about 6,000 neighbourhoods, finds none at FDR 0.1. Its per-type fold changes still track the cell-type test
(Spearman 0.96). Read Milo as *where within a cell type* a shift sits, not as an independent confirmation.

**Silent fallbacks hide bugs: check the method column, not just the exit code.** Two bugs passed every test
and finished without errors. A categorical `donor_id` made pandas emit all-zero dummy columns for patients
outside the test, so the NB model never converged and quietly used its fallback for every cell type. A
pseudobulk tool dropped the chemistry column because two patients pooled two chemistries, which emptied the
cross-study null. Both were found by asking "which method actually ran?" and "why is this table empty?".
Practical rules: cast keys to plain strings before building design matrices; log which branch a fallback took
and store it as a column; treat an empty optional output as a failure to explain.

## Key references

Heumos 2023 *Nat Rev Genet* · Squair 2021 *Nat Commun* · Zimmerman 2021 *Nat Commun* · Crowell 2020 *Nat Commun* ·
Nguyen 2023 *Nat Commun* 14:1570 · Luecken 2022 *Nat Methods* · Lopez 2018 *Nat Methods* · Xu 2021 *Mol Syst Biol* ·
Sikkema 2023 *Nat Med* · Salcher 2022 *Cancer Cell* · Tirosh 2016 *Science* · Neftel 2019 *Cell* · Patel 2014
*Science* · Gavish 2023 *Nature* · Xi & Li 2021 *Cell Syst* · Germain 2022 *F1000Res* · Yates 2025 *Genome Biol* ·
Büttner 2021 *Nat Commun* · Phipson 2022 *Bioinformatics* · Dann 2022 *Nat Biotechnol* · Setty 2019 *Nat Biotechnol* ·
Weiler 2024 *Nat Methods* · Kleshchevnikov 2022 *Nat Biotechnol* · Palla 2022 *Nat Methods* · Varrone 2024 *Nat
Genet* · Dimitrov 2022 *Nat Commun* · Dimitrov 2024 *Nat Cell Biol* · Armingol 2022 *Nat Commun* ·
Badia-i-Mompel 2022 *Bioinform Adv* · Schubert 2018 *Nat Commun* · Müller-Dott 2023 *NAR*.
