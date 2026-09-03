# Reproducibility package — RAENR-MM

**A Robust Adaptive Elastic-Net MM Estimator for Contaminated High-Dimensional
Genomic Regression: Alpha-Fetoprotein Prediction in Hepatocellular Carcinoma**

Submitted to the *Journal of Biostatistics and Epidemiology*.

This archive contains everything needed to reproduce every number, table and
figure in the manuscript and Supplementary Material from the raw GEO data.
No intermediate result is taken on trust: the pipeline runs end to end from
the series matrix file.

---

## 0. Version history

**v3.0.0, 2 September 2026 — the current version.** Corrects a penalty-scaling
defect in the proposed estimator and regenerates every reported result. The
coordinate-descent routine was receiving the cross-validated penalty parameter
without the factor n that the loss normalisation requires, so RAENR-MM was
penalised about n = 356 times too weakly. `analysis/run_shard.R` is corrected at
one line; the pre-correction file is retained under `superseded/` so the change
can be inspected by diff. Ridge, LASSO, Elastic Net and MM-RWAL are fitted
through glmnet, were unaffected, and reproduce their previous values exactly.
Adds `verification/`, carrying the code and raw output behind Supplementary
Tables S5 to S8.

**v2.0.0, 28 August 2026 — prepared but never published.** Added the first
`verification/` folder. Superseded by v3.0.0 before release.

## 1. What changed in this revision, and why this archive exists

The editorial letter raised two points that bear directly on reproducibility:

- **S4** asked whether preprocessing used information from the held-out data.
  On audit it did, at four separate points. The pipeline was rewritten and the
  entire analysis rerun. Section 4 below identifies each control by line number
  so the claim can be checked rather than believed.
- **S3** asked for a single authoritative results set. There is now exactly one.
  Every figure in the manuscript is generated programmatically from one stored
  run; nothing is transcribed by hand.

The superseded pipeline is included in `superseded/` rather than deleted, so
that the difference between the two is inspectable. See Section 8.

---

## 2. Contents

```
analysis/
  run_shard.R            Main leakage-free analysis. Runs standalone.
  cd_sweep.cpp           Compiled coordinate-descent kernel (written out by
                         run_shard.R at startup; included here for reading).
  run_all_shards.bat     Optional Windows launcher for a 4-way split.
  combine_shards.R       Merges shard outputs when the run is split.
  inspect_shard1.R       Produces the summary tables from a single run file.
  extract_stability.R    Selection frequencies and Jaccard indices (Tables S5–S8b).
  export_diagnostics.R   Collinearity diagnostics (Table A1).

results/
  shard1_summ.csv        Per-framework/dimensionality/estimator summaries.
  shard1_paired.csv      Paired RAENR-MM vs MM-RWAL comparisons.
  shard1_stability.csv   Selection stability measures.
  shard1_diagnostics.csv Collinearity diagnostics, including the superseded
                         legacy_kappa_top300 column (see section 8).

superseded/
  GSE14520_Real_Analysis_1.Rmd   The original pipeline. Not used for any
                                 reported result. Retained for comparison.
```

---

## 2a. Scope

This archive covers the real-data analysis, which is the whole of the study.
An earlier draft of the manuscript referred to a controlled Monte Carlo
simulation; that material has been withdrawn from the paper in this revision
and no claim in the manuscript now rests on it, so no simulation code is
included here. Everything in the paper is reproduced by the scripts below.

---

## 3. Data

Not redistributed here. The analysis reads one file:

| Item | Value |
|---|---|
| Accession | **GSE14520**, platform GPL3921 |
| File | `GSE14520-GPL3921_series_matrix.txt.gz` |
| Source | NCBI Gene Expression Omnibus |
| Samples used | n = 445 |
| Response | Alpha-fetoprotein (AFP) expression |
| Predictors | Gene expression, evaluated at p ∈ {1,000; 5,000; 10,000; 15,000; 20,000} |

Download the series matrix from GEO and set the path at the top of
`analysis/run_shard.R`:

```r
GEO_FILE <- ".../GSE14520-GPL3921_series_matrix.txt.gz"
```

The file is read as distributed by the depositors. There are **no missing
expression values** after their preprocessing, and no missing-data imputation
is performed anywhere in this pipeline. (The "median-imputed matrix" referred
to in the manuscript is the deliberate construction of the cellwise evaluation
set, which is a different operation — see Section 6.)

---

## 4. Leakage controls

Every data-dependent quantity is estimated on the training fold and then
applied to the held-out fold. The four controls are marked in the source with
the tags `L1`–`L4`:

| Tag | Line | Control |
|---|---|---|
| **L1** | 476 | Centring and scaling constants from training rows only (`colMeans`, `colSds` of `Xtr_raw`). |
| **L2** | 479 | Variance ranking for the top-*p* filter computed on training rows only; the resulting index set `sel` is then applied to the held-out rows. |
| **L3** | 487 | Cellwise medians and MADs (`med_tr`, `mad_tr`) computed on the standardised **training** matrix; used to flag and replace cells in the held-out predictors. |
| **L4** | 500 | Rowwise outlier rule (PCA + MCD, `alpha = 0.75`, 30 components) **fitted** on training rows and **scored** on held-out rows. |

Two further reproducibility gaps found in the same audit and closed here:

- **L5** — the split ratio and fold construction are fixed and seeded
  (`train_prop <- 0.80`, line 62; `set.seed(100000L + rep_id)`, line 521).
- **L6** — all tuning is by cross-validation inside the training fold
  (`cv.glmnet(..., nfolds = NFOLD)`), never spanning the split.

Because the per-replication seed depends only on `rep_id`, **all five
estimators see identical training and held-out sets within a replication.**
The comparisons are therefore paired by construction, which is why the
manuscript reports paired *t*-tests on the per-replication differences rather
than two-sample tests.

---

## 5. Parameters

All set at the top of `run_shard.R` (lines 60–70):

| Parameter | Value | Meaning |
|---|---|---|
| `set.seed` | 42 | Global seed |
| per-replication seed | `100000 + rep_id` | Guarantees identical splits across estimators |
| `P_GRID` | 1000, 5000, 10000, 15000, 20000 | Dimensionalities |
| `n_reps` | 1000 | Replications |
| `train_prop` | 0.80 | Training proportion |
| `gamma_raen` | 2 | Adaptive weight exponent γ |
| `alpha_raen` | 0.2 | RAENR-MM elastic-net mixing |
| `alpha_en` | 0.5 | Elastic Net mixing (comparator) |
| `c_val` | 4.685 | Tukey bisquare tuning constant |
| `delta` | 1e-6 | Adaptive weight stabiliser |
| `NLAM` / `NFOLD` | 20 / 3 | λ grid size, CV folds |
| `MAX_IRLS` / `MAX_CD` | 5 / 50 | Outer IRLS and inner coordinate-descent caps |

γ is fixed at **2** throughout. The previous version of the manuscript was
internally inconsistent on this point — the text specified γ = 2 while the
analysis code used γ = 1 — and the value is now the same in both.

The robust scale is **held fixed** across outer IRLS iterations, which brings
the implementation inside the hypotheses of the monotone-descent theorem
(editorial item S1).

---

## 6. The two contamination frameworks

Both are evaluated from the same fitted models in the same replication; only
the held-out set differs (`run_shard.R`, lines 541–560).

- **THCM (rowwise).** Held-out rows flagged by the L4 rule are excluded from
  evaluation. Mean held-out size ≈ 47.
- **ICM (cellwise).** No rows excluded; held-out cells flagged at |z| > 3
  against the training median/MAD are replaced by the training column median.
  Held-out size = 89.

**The two frameworks therefore do not evaluate on sets of the same size.**
This is a property of the design, not a bug, but it was not disclosed in the
previous version and it means the ICM/THCM error *ratio* is not a pure
contamination effect. It is now reported in Section 5.8 of the manuscript.
Comparisons **within** a framework — which is every comparison the conclusions
rest on — are unaffected.

Note that `Selected` (the number of genes chosen) is identical across the two
frameworks by construction: selection happens on the training fold, which is
shared, and the frameworks differ only in how the held-out set is built.

---

## 7. Running it

### Requirements

R ≥ 4.3.1 with a working C++ toolchain (Rtools on Windows), and:

```r
install.packages(c("glmnet", "robustbase", "matrixStats", "Rcpp",
                   "dplyr", "tidyr", "tibble"))
```

`run_shard.R` writes `cd_sweep.cpp` to its working directory and compiles it
with `Rcpp::sourceCpp` at startup. It then **verifies** that the compiled
kernel is bound and agrees with the pure-R reference implementation
(lines 137–143 and 464–472), and aborts if not. This check exists because a
silent fallback to the R kernel slows the run by roughly a factor of forty
without changing any result.

### Single process — the configuration used for the reported results

```
"C:\Program Files\R\R-4.3.1\bin\Rscript.exe" run_shard.R 1 1
```

Runtime **≈ 18 hours** (1,073 minutes measured) for 1,000 replications across
all five dimensionalities. Writes `C:/raenr_out/shard_01.rds`.

### Split across four processes

```
run_all_shards.bat            # or four manual calls: run_shard.R k 4, k = 1..4
```

Replications are interleaved (`seq(SHARD, n_reps, by = NSHARD)`), so each
process carries an equal load. ≈ 4.5 hours each. Then:

```r
source("combine_shards.R")    # refuses to run unless all four shards are present
```

### Producing the reported tables

```r
source("inspect_shard1.R")      # summaries, ranks, paired comparisons
source("extract_stability.R")   # selection frequencies and Jaccard indices
source("export_diagnostics.R")  # collinearity diagnostics
```

---

## 8. Mapping from manuscript to code

| Manuscript item | Produced by |
|---|---|
| Table 2 / A1 — diagnostics | `export_diagnostics.R` → `collinearity_diag()`, `run_shard.R` lines 562–577 |
| Table 3 / A2 — RMSE | `inspect_shard1.R` → `shard1_summ.csv` |
| Table 4 / A3 — R² | `inspect_shard1.R` → `shard1_summ.csv` |
| Table 5 / A4 — MAE | `inspect_shard1.R` → `shard1_summ.csv` |
| Table 6 / A5 — detail by framework | `inspect_shard1.R` → `shard1_summ.csv` |
| Table 7 / A6 — paired improvement | `inspect_shard1.R` → `shard1_paired.csv` |
| Table 8 / A7 — framework comparison | `inspect_shard1.R` → `shard1_summ.csv` |
| Table 9 / A8 — comprehensive ranking | `inspect_shard1.R` → `shard1_summ.csv` |
| Tables S5–S8b — selection stability | `extract_stability.R` → `shard1_stability.csv` |

### On the condition number (editorial item S6)

The letter correctly observed that a condition number cannot be identical at
p = 1,000 and p = 20,000. The cause is in `collinearity_diag()`:

```r
p_s   <- min(300L, ncol(Ztr))
R_sub <- cor(Ztr[, seq_len(p_s), drop = FALSE])
```

The diagnostic used the first 300 columns only. Because the matrix is sorted
by decreasing marginal variance, those are the same 300 genes at every
dimensionality — so the figure was constant because it was the same 300 × 300
matrix each time. The quantity reported alongside it as a "maximum condition
index" was the square root of that number, not a Belsley–Kuh–Welsch condition
index.

`collinearity_diag()` now returns the corrected quantities —
`rank`, `kappa_rank_restricted`, `max_condition_index`, `prop_sv_below_1pct` —
**and** retains `legacy_kappa_top300`, so that the superseded figure and its
replacement appear side by side in the same output rather than the change
being made silently.

---

## 9. The superseded pipeline

`superseded/GSE14520_Real_Analysis_1.Rmd` produced the results in the previous
submission. **No number in the current manuscript comes from it.** It is
included so that the four leakage points can be located in the original source
and compared against the controls listed in Section 4.

The estimator implementations in `run_shard.R` are ported **verbatim** from
this file — the blocks are marked as such in the source (`safe_coef`,
`cd_helpers`, `glmnet_estimators`, `mm_rwal`, `raenr_mm`). Only the
preprocessing and evaluation pipeline was rewritten. This was deliberate: it
keeps the change confined to the thing the editorial comment was about, so
that differences between the old and new results are attributable to the
leakage correction and not to an incidental reimplementation of the estimators.

---

## 10. Effect of the correction

Two corrections have been applied to this pipeline, in successive versions.

**The leakage correction (v2).** Centring, variance filtering, standardisation,
cellwise imputation and the rowwise flagging rule are now fitted on the training
fold only and applied to held-out rows, closing four points of train/test
leakage present in v1.0.0.

**The penalty correction (v3).** `cv.glmnet` minimises
(1/2n)||y - Xb||^2 + lambda[.], so a lambda selected by cross-validation belongs
with a penalty of n*lambda once the loss is written with a 1/2 factor. The
implementation passed lambda. The effect on the proposed estimator was large:

| Quantity, p = 15,000 | Before | After |
|---|---|---|
| Mean genes selected | 835.2 | 44.4 |
| Mean pairwise Jaccard | 0.125 | 0.4635 |
| Genes ever selected, p = 20,000 | 19,479 of 20,000 | 653 |
| Outer loop meets its tolerance | 0% of fits | 39% |
| Inner passes meet their tolerance | 0% | 75.9% |
| Mean test RMSE, THCM | 0.1683 | 0.1678 |

Prediction is essentially unchanged; sparsity and selection stability are
transformed. The four comparators are unaffected in both corrections.

To see the penalty change:

```
diff superseded/run_shard_PRE_PATCH_2026-09-01.R analysis/run_shard.R
```

They differ in one functional line.

## 11. Citation

If you use this code, please cite the paper. The code itself is citable as
release v3.0.0 of this repository.

Correspondence regarding the code: see the corresponding author's details on
the title page.
