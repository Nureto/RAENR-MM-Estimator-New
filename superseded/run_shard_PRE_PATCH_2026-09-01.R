## =====================================================================
##  GSE14520 LEAKAGE-FREE ANALYSIS -- ONE SHARD, STANDALONE
##
##  No future. No parallel workers. No serialisation of functions, pointers
##  or data between processes. Every failure so far came from that machinery.
##
##  Instead: launch N of these as SEPARATE R processes, each handling a slice
##  of the replications. Each process loads its own data, compiles its own
##  kernel, and writes its own .rds. Nothing is shared but the input file.
##
##  USAGE (from a Windows command prompt, one per shard):
##      "C:\Program Files\R\R-4.3.1\bin\Rscript.exe" run_shard.R 1 4
##      "C:\Program Files\R\R-4.3.1\bin\Rscript.exe" run_shard.R 2 4
##      "C:\Program Files\R\R-4.3.1\bin\Rscript.exe" run_shard.R 3 4
##      "C:\Program Files\R\R-4.3.1\bin\Rscript.exe" run_shard.R 4 4
##  or just double-click run_all_shards.bat
##
##  Then, in RStudio:  source("combine_shards.R")
## =====================================================================

args   <- commandArgs(trailingOnly = TRUE)
SHARD  <- if (length(args) >= 1) as.integer(args[1]) else 1L
NSHARD <- if (length(args) >= 2) as.integer(args[2]) else 1L
stopifnot(SHARD >= 1L, SHARD <= NSHARD)

## ------------------------------------------------------------------ SET THIS
## Path to the GEO series matrix. Download GSE14520 (platform GPL3921) from
## https://www.ncbi.nlm.nih.gov/geo/ and either place the .txt.gz beside this
## script or set GSE14520_PATH in the environment.
GEO_FILE  <- Sys.getenv("GSE14520_PATH",
                        "GSE14520-GPL3921_series_matrix.txt.gz")
if (!file.exists(GEO_FILE))
  stop("GEO series matrix not found at: ", GEO_FILE,
       "\n  Download GSE14520-GPL3921_series_matrix.txt.gz from GEO and place\n  it beside this script, or set the GSE14520_PATH environment variable.")
## Output and kernel-build directories. Defaults are relative so the script
## runs anywhere; override with environment variables if preferred.
## NOTE: on Windows, keep the build directory path SHORT. Compiling under a
## deeply nested path can exceed MAX_PATH and fail with "CreateProcess".
OUT_DIR   <- Sys.getenv("RAENR_OUT",   file.path(getwd(), "raenr_out"))
BUILD_DIR <- Sys.getenv("RAENR_BUILD", file.path(tempdir(), "raenr_kernel"))
## ---------------------------------------------------------------------------

dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(BUILD_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(glmnet); library(robustbase)
  library(matrixStats)
})

## single-threaded BLAS: N shards already use N cores
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
}

set.seed(42L); options(stringsAsFactors = FALSE)
RNGkind("Mersenne-Twister", "Inversion")

P_GRID     <- c(1000L, 5000L, 10000L, 15000L, 20000L)
n_reps     <- 1000L
train_prop <- 0.80
gamma_raen <- 2; alpha_en <- 0.5; alpha_raen <- 0.2
c_val <- 4.685; delta <- 1e-6
NLAM <- 20L; NFOLD <- 3L; MAX_IRLS <- 5L; MAX_CD <- 50L
TOL_IRLS <- 1e-3; TOL_CD <- 1e-4
METHODS <- c("Ridge","LASSO","ENET","MM-RWAL","RAENR-MM")
SPARSE  <- c("LASSO","ENET","MM-RWAL","RAENR-MM")

MY_REPS <- seq(SHARD, n_reps, by = NSHARD)      # interleaved: even load
cat(sprintf("shard %d/%d | %d replications | pid %d\n",
            SHARD, NSHARD, length(MY_REPS), Sys.getpid()))

## ===================== compiled coordinate-descent kernel ==================
CPP <- c(
'#include <Rcpp.h>','using namespace Rcpp;','',
'// [[Rcpp::export]]',
'NumericVector cd_sweep_cpp(NumericMatrix Xw, NumericVector yw,',
'                           NumericVector beta_in, NumericVector wp,',
'                           double lam, double alpha) {',
'  const int n = Xw.nrow(), p = Xw.ncol();',
'  NumericVector beta = clone(beta_in);',
'  std::vector<double> d(p), r(n);',
'  for (int j = 0; j < p; ++j) {',
'    const double* col = &Xw(0, j); double s = 0.0;',
'    for (int i = 0; i < n; ++i) s += col[i] * col[i];',
'    d[j] = s;',
'  }',
'  for (int i = 0; i < n; ++i) r[i] = yw[i];',
'  for (int j = 0; j < p; ++j) {',
'    const double bj = beta[j];',
'    if (bj != 0.0) { const double* col = &Xw(0, j);',
'      for (int i = 0; i < n; ++i) r[i] -= col[i] * bj; }',
'  }',
'  const double l2 = lam * (1.0 - alpha);',
'  for (int j = 0; j < p; ++j) {',
'    const double* col = &Xw(0, j); const double bj = beta[j];',
'    if (bj != 0.0) for (int i = 0; i < n; ++i) r[i] += col[i] * bj;',
'    double z = 0.0;',
'    for (int i = 0; i < n; ++i) z += col[i] * r[i];',
'    const double den = d[j] + l2; double bnew;',
'    if (std::fabs(den) < 1e-12) { bnew = 0.0; }',
'    else { const double zz = z / den;',
'           const double g  = lam * alpha * wp[j] / den;',
'           const double a  = std::fabs(zz) - g;',
'           bnew = (a > 0.0) ? (zz > 0 ? a : (zz < 0 ? -a : 0.0)) : 0.0; }',
'    beta[j] = bnew;',
'    if (bnew != 0.0) for (int i = 0; i < n; ++i) r[i] -= col[i] * bnew;',
'  }',
'  return beta;',
'}')

## R reference kernel, kept for the equivalence check
soft_threshold <- function(z, gamma) sign(z) * pmax(abs(z) - gamma, 0)
cd_sweep_R <- function(Xw, yw, beta, wp, lam, alpha) {
  p <- length(beta)
  d <- colSums(Xw^2)
  r <- as.vector(yw - Xw %*% beta)
  for (j in seq_len(p)) {
    r   <- r + Xw[, j] * beta[j]
    z_j <- sum(Xw[, j] * r)
    den <- d[j] + lam * (1 - alpha)
    beta[j] <- if (abs(den) < 1e-12) 0L else
      soft_threshold(z_j / den, lam * alpha * wp[j] / den)
    r <- r - Xw[, j] * beta[j]
  }
  beta
}

wd  <- file.path(BUILD_DIR, paste0("shard", SHARD, "_", Sys.getpid()))
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
src <- file.path(wd, "cd_sweep.cpp"); writeLines(CPP, src)
old <- setwd(wd)                       # short path: R CMD SHLIB needs it
Rcpp::sourceCpp(src)                   # compiles into THIS process only
setwd(old)

## prove it works here, and that it reproduces the R kernel
set.seed(99)
.Xw <- matrix(rnorm(50*20), 50, 20); .yw <- rnorm(50)
.b1 <- cd_sweep_cpp(.Xw, .yw, numeric(20), rep(1,20), 0.1, 0.2)
.b2 <- cd_sweep_R  (.Xw, .yw, numeric(20), rep(1,20), 0.1, 0.2)
stopifnot(all(is.finite(.b1)), max(abs(.b1 - .b2)) < 1e-10)
cat("kernel compiled and verified against the R version\n")

## NOTE: cd_sweep is bound to the compiled kernel AFTER the estimator block
## below, because that block contains the original pure-R cd_sweep and would
## otherwise overwrite it. Binding here would be silently undone.

## ===================== estimators (verbatim) ===============
## ===== ported verbatim from GSE14520_Real_Analysis_1.Rmd, chunk `safe_coef` =====
safe_coef_vec <- function(beta, p) {
  if (is.null(beta) || length(beta)!=p || any(!is.finite(beta)))
    return(rep(0, p))
  as.vector(beta)
}

## ===== ported verbatim from GSE14520_Real_Analysis_1.Rmd, chunk `cd_helpers` =====
## Soft-thresholding — equation (2.37)
## S(z, gamma) = sign(z) · max(|z| − gamma, 0)
soft_threshold <- function(z, gamma) sign(z) * pmax(abs(z) - gamma, 0)

## One full coordinate descent sweep — equations (2.34)–(2.38)
##
## At entry: beta is current coefficient vector
##           Xw  = diag(sqrt(w_i^{(r)})) X    (observation-weighted X)
##           yw  = sqrt(w_i^{(r)}) * y         (observation-weighted y)
##           wp  = adaptive penalty weights w_j^{(a)}
##           lam = fixed lambda
##           alpha = elastic-net mixing parameter
##
## For each predictor j = 1, …, p:
##   (2.34) r_{i,-j} = yw_i - sum_{k != j} Xw_ik * beta_k
##   (2.35) z_j      = sum_i Xw_ij * r_{i,-j}
##          d_j      = sum_i Xw_ij^2
##   (2.36) denom    = d_j + lambda * (1 - alpha)
##          beta_j  <- S( z_j / denom,  lambda * alpha * w_j^{(a)} / denom )
##   (2.38) residuals updated in-place after beta_j changes
cd_sweep <- function(Xw, yw, beta, wp, lam, alpha) {
  p <- length(beta)
  d <- colSums(Xw^2)                      # d_j = sum_i Xw_ij^2   -- eq (2.35)
  r <- as.vector(yw - Xw %*% beta)        # initialise full residual

  for (j in seq_len(p)) {

    r   <- r + Xw[, j] * beta[j]          # eq (2.34): restore partial residual
    z_j <- sum(Xw[, j] * r)               # eq (2.35): weighted partial corr

    den     <- d[j] + lam * (1 - alpha)   # eq (2.36): denominator
    beta[j] <- if (abs(den) < 1e-12) 0L else
      soft_threshold(z_j / den,           # eq (2.36): coordinate update
                     lam * alpha * wp[j] / den)

    r <- r - Xw[, j] * beta[j]            # eq (2.38): update residuals
  }
  return(beta)
}

## ===== ported verbatim from GSE14520_Real_Analysis_1.Rmd, chunk `glmnet_estimators` =====
fit_ridge <- function(X, y) {
  p <- ncol(X); if (sd(y)<1e-8) return(rep(0,p))
  f <- tryCatch(cv.glmnet(X,y,alpha=0,nfolds=NFOLD,nlambda=NLAM),
                error=function(e)NULL)
  if (is.null(f)) return(rep(0,p))
  safe_coef_vec(as.vector(coef(f,s="lambda.min"))[-1],p)
}
fit_lasso <- function(X, y) {
  p <- ncol(X); if (sd(y)<1e-8) return(rep(0,p))
  f <- tryCatch(cv.glmnet(X,y,alpha=1,nfolds=NFOLD,nlambda=NLAM),
                error=function(e)NULL)
  if (is.null(f)) return(rep(0,p))
  safe_coef_vec(as.vector(coef(f,s="lambda.min"))[-1],p)
}
fit_enet <- function(X, y, al=0.5) {
  p <- ncol(X); if (sd(y)<1e-8) return(rep(0,p))
  f <- tryCatch(cv.glmnet(X,y,alpha=al,nfolds=NFOLD,nlambda=NLAM),
                error=function(e)NULL)
  if (is.null(f)) return(rep(0,p))
  safe_coef_vec(as.vector(coef(f,s="lambda.min"))[-1],p)
}

## ===== ported verbatim from GSE14520_Real_Analysis_1.Rmd, chunk `mm_rwal` =====
## MM-RWAL (Machkour et al., 2020; Toka et al., 2021)
##
## STEP 1:  MM-LASSO pilot  →  beta_hat_mm_lasso  (equations 2.19-2.20)
##   p < n : lmrob (true MM) residuals → Tukey weights → weighted LASSO
##   p >= n: Huber-IRLS proxy → Tukey weights → weighted LASSO
##           (lmrob infeasible when p >= n)
##
## STEP 2:  Observation weights from MM-LASSO residuals
##   w(r_i) = (1 - (r_i / (c * s_n))^2)^2  for |r_i/(c*s_n)| <= 1, else 0
##
## STEP 3:  z_j normalisation — equation (2.20)
##   z_j = mean(w) * p / sum(mean(w)) ,  sum(z_j) = p
##
## STEP 4:  Adaptive weights — equation (2.20)
##   w_j = 1 / (z_j * |beta_j_mm_lasso| + delta)
##
## STEP 5:  One cv.glmnet to select fixed lambda
##
## STEP 6:  IRLS loop with fixed lambda and adaptive weights
fit_mm_rwal <- function(X, y, c_val=4.685, delta=1e-6,
                         max_iter=3L, tol=1e-3) {
  n <- nrow(X); p <- ncol(X)

  ## ── STEP 1: MM-LASSO pilot ──────────────────────────────────────────────
  beta_mm_lasso <- tryCatch({

    if (p < n) {
      ## True MM: lmrob residuals -> Tukey weights -> weighted LASSO
      mm_fit <- suppressWarnings(
        lmrob(y ~ X - 1, setting="KS2014", k.max=200, maxit.scale=200))
      r_mm  <- residuals(mm_fit); s_mm <- mad(r_mm)
      u_mm  <- r_mm / max(s_mm, 1e-10)
      rw_mm <- ifelse(abs(u_mm)<=c_val, (1-(u_mm/c_val)^2)^2, 0)
      Xw_mm <- (X * sqrt(rw_mm)); yw_mm <- y*sqrt(rw_mm)
      f0    <- cv.glmnet(Xw_mm, yw_mm, alpha=1, nfolds=NFOLD, nlambda=NLAM)
      safe_coef_vec(as.vector(coef(f0, s="lambda.1se"))[-1], p)

    } else {
      ## p >= n: Huber-IRLS proxy (3 passes, c=1.345)
      f0    <- cv.glmnet(X, y, alpha=1, nfolds=NFOLD, nlambda=NLAM)
      b_h   <- safe_coef_vec(as.vector(coef(f0,s="lambda.1se"))[-1], p)
      lam_h <- f0$lambda.1se
      for (hs in seq_len(3)) {
        r_h  <- as.vector(y - X %*% b_h)
        s_h  <- mad(r_h); if (!is.finite(s_h)||s_h<1e-10) break
        u_h  <- r_h/s_h; rw_h <- pmin(1, 1.345/abs(u_h))
        rw_h[!is.finite(rw_h)] <- 1
        Xw_h <- (X * sqrt(rw_h)); yw_h <- y*sqrt(rw_h)
        fg_h <- tryCatch(
          glmnet(Xw_h,yw_h,alpha=1,lambda=lam_h,standardize=FALSE),
          error=function(e)NULL)
        if (is.null(fg_h)) break
        bn   <- as.vector(coef(fg_h))[-1]
        if (length(bn)!=p||any(!is.finite(bn))) break
        if (sqrt(sum((bn-b_h)^2))<1e-3){b_h<-bn;break}; b_h <- bn
      }
      ## Tukey weights from Huber-IRLS residuals
      r_rob  <- as.vector(y - X %*% b_h)
      s_rob  <- mad(r_rob); if (!is.finite(s_rob)||s_rob<1e-10) s_rob <- mad(y)
      u_rob  <- r_rob/s_rob
      rw_rob <- ifelse(abs(u_rob)<=c_val,(1-(u_rob/c_val)^2)^2,0)
      Xw_rob <- (X * sqrt(rw_rob)); yw_rob <- y*sqrt(rw_rob)
      f1     <- tryCatch(
        cv.glmnet(Xw_rob,yw_rob,alpha=1,nfolds=NFOLD,nlambda=NLAM),
        error=function(e)f0)
      safe_coef_vec(as.vector(coef(f1,s="lambda.1se"))[-1], p)
    }

  }, error=function(e) rep(0,p))

  ## ── STEP 2: observation weights from MM-LASSO residuals ─────────────────
  r_p <- as.vector(y - X %*% beta_mm_lasso)
  s_n <- mad(r_p); if (!is.finite(s_n)||s_n<1e-10) s_n <- 1
  u_p <- r_p/s_n
  ow  <- ifelse(abs(u_p)<=c_val, (1-(u_p/c_val)^2)^2, 0)
  mw  <- mean(ow); if (mw<1e-12) mw <- 1

  ## ── STEP 3: z_j normalisation — equation (2.20),  sum(z_j) = p ─────────
  z_j <- rep(mw, p); z_j <- z_j * p / sum(z_j)

  ## ── STEP 4: adaptive weights — equation (2.20) ──────────────────────────
  ## w_j = 1 / (z_j * |beta_j_mm_lasso| + delta)
  wp  <- 1 / (z_j * abs(beta_mm_lasso) + delta)
  wp  <- pmin(wp, quantile(wp, 0.95)); wp <- wp/mean(wp)

  ## ── STEP 5: one cv.glmnet with adaptive weights to fix lambda ───────────
  pf  <- tryCatch(
    cv.glmnet(X,y,alpha=1,penalty.factor=wp,nfolds=NFOLD,nlambda=NLAM),
    error=function(e)NULL)
  if (is.null(pf)) return(beta_mm_lasso)
  lam  <- pf$lambda.1se
  beta <- safe_coef_vec(as.vector(coef(pf,s="lambda.1se"))[-1],p)

  ## ── STEP 6: IRLS loop — minimises equation (2.19) ───────────────────────
  for (i in seq_len(max_iter)) {
    r  <- as.vector(y - X %*% beta); sn <- mad(r)
    if (!is.finite(sn)||sn<1e-10) break
    u  <- r/sn; rw <- ifelse(abs(u)<=c_val,(1-(u/c_val)^2)^2,0)
    if (sum(rw)<1e-6) break
    Xw <- (X * sqrt(rw)); yw <- y*sqrt(rw)
    if (!is.finite(sd(yw))||sd(yw)<1e-8) break
    fg <- tryCatch(
      glmnet(Xw,yw,alpha=1,lambda=lam,penalty.factor=wp,standardize=FALSE),
      error=function(e)NULL)
    if (is.null(fg)) break
    bn <- as.vector(coef(fg))[-1]
    if (length(bn)!=p||any(!is.finite(bn))) break
    if (sum((bn-beta)^2)<tol^2){beta<-bn;break}; beta <- bn
  }
  return(beta)
}

## ===== ported verbatim from GSE14520_Real_Analysis_1.Rmd, chunk `raenr_mm` =====
## RAENR-MM (Alabi & Oyeyemi, 2025)
## Main novelty: adaptive weights from MM initial pilot, not from EN pilot.
##
## STEP 1: EN pilot       -> fixed lambda,  seed beta_0
## STEP 2: MM pilot       -> beta_hat_mm  for adaptive weights  (eq 2.22)
##   p < n  : lmrob  (true MM estimator)
##   p >= n : Huber-IRLS proxy  (lmrob infeasible)
## STEP 3: Adaptive weights  w_j^{(a)} = 1/(|beta_j_mm|^gamma + delta)  (eq 2.22)
##   NOTE: wp is computed ONCE from the MM pilot and stays FIXED for ALL
##         IRLS iterations. Only the observation weights w_i^{(r)} change.
## STEP 4: Outer IRLS + Inner CD loop
##   IRLS (outer, eq 2.23): recomputes Tukey observation weights each pass
##   CD   (inner, eq 2.34-2.38): cycles through all p predictors per pass
##   Convergence (eq 2.39): ||beta_new - beta_old|| < tol for both loops
fit_raen_mm <- function(X, y, gamma    = 2,
                               c_val   = 4.685,
                               delta   = 1e-6,
                               al      = 0.2,
                               max_irls = MAX_IRLS,
                               max_cd   = MAX_CD,
                               tol_irls = TOL_IRLS,
                               tol_cd   = TOL_CD) {
  n <- nrow(X); p <- ncol(X)
  if (!is.finite(sd(y))||sd(y)<1e-8) return(rep(0,p))

  ## ── STEP 1: EN pilot — selects lambda, seeds IRLS starting point ────────
  pilot_en <- tryCatch(
    cv.glmnet(X, y, alpha=al, nfolds=NFOLD, nlambda=NLAM),
    error=function(e)NULL)
  if (is.null(pilot_en)) return(rep(0,p))
  lambda_fixed <- pilot_en$lambda.1se
  beta         <- safe_coef_vec(as.vector(coef(pilot_en,s="lambda.1se"))[-1],p)

  ## ── STEP 2: MM pilot — beta_hat_mm for adaptive weights  eq (2.22) ──────
  ## This is the MAIN NOVELTY: adaptive weights must come from the MM
  ## initial estimator, not from the EN pilot. Using EN here would cause
  ## outlier-inflated coefficients to receive spuriously small penalties,
  ## defeating the purpose of the adaptive weighting.
  beta_mm <- tryCatch({

    if (p < n) {
      ## True MM estimator (50% breakdown point)
      mm_fit <- suppressWarnings(
        lmrob(y ~ X - 1, setting="KS2014", k.max=200, maxit.scale=200))
      out <- as.vector(coef(mm_fit))
      if (length(out)!=p||any(!is.finite(out))) rep(0,p) else out

    } else {
      ## Huber-IRLS proxy for p >= n (3 passes, Huber c=1.345)
      b_h   <- beta; lam_h <- lambda_fixed
      for (hs in seq_len(3)) {
        r_h  <- as.vector(y - X %*% b_h)
        s_h  <- mad(r_h); if (!is.finite(s_h)||s_h<1e-10) break
        u_h  <- r_h/s_h; rw_h <- pmin(1, 1.345/abs(u_h))
        rw_h[!is.finite(rw_h)] <- 1
        Xw_h <- (X * sqrt(rw_h)); yw_h <- y*sqrt(rw_h)
        fg_h <- tryCatch(
          glmnet(Xw_h,yw_h,alpha=al,lambda=lam_h,standardize=FALSE),
          error=function(e)NULL)
        if (is.null(fg_h)) break
        bn_h <- as.vector(coef(fg_h))[-1]
        if (length(bn_h)!=p||any(!is.finite(bn_h))) break
        if (sqrt(sum((bn_h-b_h)^2))<1e-3){b_h<-bn_h;break}; b_h <- bn_h
      }
      b_h
    }

  }, error=function(e) beta)   # fall back to EN pilot only if MM crashes

  ## ── STEP 3: Adaptive penalty weights — equation (2.22) ──────────────────
  ## w_j^{(a)} = 1 / ( |beta_j_mm|^gamma + delta )
  ## Large |beta_j_mm| -> small w_j -> weak penalty -> predictor survives
  ## Small |beta_j_mm| -> large w_j -> strong penalty -> predictor shrinks to 0
  ## wp is computed ONCE here and is NOT updated inside the IRLS loop.
  wp <- 1 / (abs(beta_mm)^gamma + delta)   # equation (2.22)
  wp <- pmin(wp, quantile(wp, 0.95))       # winsorise extreme weights
  wp <- wp / mean(wp)                      # rescale: mean(wp) = 1

  ## ── STEP 4: Outer IRLS + Inner CD — equations (2.23), (2.34)–(2.39) ────
  ##
  ## Outer IRLS loop:
  ##   - Compute Tukey biweight observation weights w_i^{(r)}  (eq 2.23)
  ##   - Form weighted matrices Xw = diag(sqrt(rw)) X, yw = sqrt(rw) y
  ##   - Run inner CD loop to convergence
  ##   - Check outer convergence on ||beta_new - beta_old||  (eq 2.39)
  ##
  ## Inner CD loop:
  ##   - Each call to cd_sweep cycles through ALL p predictors once
  ##   - Implements equations (2.34)-(2.38) exactly
  ##   - Check inner convergence on ||beta_new - beta_old||  (eq 2.39)
  ##
  ## The inner problem solved each IRLS step is:
  ##   min_beta  ||yw - Xw beta||^2
  ##             + lambda [ alpha  sum_j w_j^{(a)} |beta_j|
  ##                      + (1-alpha) sum_j beta_j^2 ]
  ## with FIXED Xw, yw (current IRLS weights), FIXED lambda, FIXED wp.
  for (irls in seq_len(max_irls)) {

    ## eq (2.23): Tukey biweight observation weights
    r         <- as.vector(y - X %*% beta)
    sigma_hat <- mad(r)
    if (!is.finite(sigma_hat)||sigma_hat<1e-10) break

    u  <- r / sigma_hat
    rw <- ifelse(abs(u)<=c_val, (1-(u/c_val)^2)^2, 0)   # eq (2.23)
    if (sum(rw)<1e-6) break

    Xw <- (X * sqrt(rw))         # Xw_ij = sqrt(w_i^{(r)}) * x_ij
    yw <- y * sqrt(rw)                       # yw_i  = sqrt(w_i^{(r)}) * y_i
    if (!is.finite(sd(yw))||sd(yw)<1e-8) break

    beta_pre_cd <- beta                      # save for outer convergence check

    ## Inner CD loop — equations (2.34)-(2.38)
    for (cd_iter in seq_len(max_cd)) {
      beta_old <- beta
      beta     <- cd_sweep(Xw, yw, beta, wp, lambda_fixed, al)

      ## Inner convergence — equation (2.39)
      if (sqrt(sum((beta-beta_old)^2)) < tol_cd) break
    }

    ## Outer convergence — equation (2.39)
    if (sqrt(sum((beta-beta_pre_cd)^2)) < tol_irls) break
  }

  return(beta)
}
## ---- bind the compiled kernel NOW, after the verbatim estimator block -----
## The block above defines the original pure-R cd_sweep(). Overriding it here
## is what makes the run fast. Verify the binding actually took effect.
stopifnot(exists("cd_sweep_cpp"))
cd_sweep <- cd_sweep_cpp
.chk <- cd_sweep(matrix(c(1,0,0,1,1,0), 3, 2), c(1,2,3), c(0,0), c(1,1), 0.1, 0.2)
stopifnot(is.numeric(.chk), length(.chk) == 2L, all(is.finite(.chk)))
if (!identical(body(cd_sweep), body(cd_sweep_cpp)))
  stop("cd_sweep is NOT the compiled kernel - the run would take ~40x longer")
cat("cd_sweep bound to the compiled kernel (verified)\n")

## ===================== preprocessing =======================
make_preprocessor <- function(Xtr_raw, p_keep) {
  ctr   <- colMeans(Xtr_raw)                       # L1
  sdv   <- matrixStats::colSds(Xtr_raw)
  keep0 <- which(sdv > 1e-8)
  v     <- sdv[keep0]^2                            # L2
  ord   <- order(v, decreasing = TRUE)
  sel   <- keep0[ord[seq_len(min(p_keep, length(ord)))]]
  ctr_s <- ctr[sel]; sdv_s <- sdv[sel]

  std <- function(Xraw)
    sweep(sweep(Xraw[, sel, drop=FALSE], 2, ctr_s, "-"), 2, sdv_s, "/")

  Ztr    <- std(Xtr_raw)                           # L3
  med_tr <- matrixStats::colMedians(Ztr)
  mad_tr <- pmax(matrixStats::colMads(Ztr), 1e-8)

  list(sel = sel, apply_std = std,
       impute = function(Z) {
         Zc <- Z
         bad <- abs(sweep(sweep(Z, 2, med_tr, "-"), 2, mad_tr, "/")) > 3
         for (j in seq_len(ncol(Z))) { b <- bad[, j]; if (any(b)) Zc[b, j] <- med_tr[j] }
         Zc
       })
}

## L4: rowwise rule FITTED on training, SCORED on held-out rows
train_rowwise_flagger <- function(Ztr, n_pc = 30L) {
  n_pc <- min(n_pc, nrow(Ztr) - 2L, ncol(Ztr))
  pca  <- prcomp(Ztr, center=FALSE, scale.=FALSE, rank.=n_pc)
  mcd  <- tryCatch(robustbase::covMcd(pca$x, alpha=0.75), error=function(e) NULL)
  if (is.null(mcd)) return(NULL)
  Si <- tryCatch(solve(mcd$cov), error=function(e) NULL)
  if (is.null(Si)) return(NULL)
  list(thr = sqrt(qchisq(0.975, df=n_pc)),
       score = function(Znew) {
         dd <- sweep(Znew %*% pca$rotation, 2, mcd$center, "-")
         sqrt(pmax(rowSums((dd %*% Si) * dd), 0))
       })
}

perf_metrics <- function(yt, yp) {
  ss_r <- sum((yt-yp)^2); ss_t <- sum((yt-mean(yt))^2)
  c(RMSE=sqrt(mean((yt-yp)^2)), MAE=mean(abs(yt-yp)),
    R2=if (ss_t<1e-12) NA_real_ else 1-ss_r/ss_t)
}
run_rep <- function(rep_id, pp) {
  set.seed(100000L + rep_id)
  n_te <- n_all - floor(train_prop * n_all)
  te   <- sort(sample(seq_len(n_all), n_te))
  tr   <- setdiff(seq_len(n_all), te)

  pre <- make_preprocessor(X_raw[tr, , drop=FALSE], pp)
  Ztr <- pre$apply_std(X_raw[tr, , drop=FALSE])
  Zte <- pre$apply_std(X_raw[te, , drop=FALSE])
  ytr <- y_full[tr]; yte <- y_full[te]; mu <- mean(ytr)
  yc  <- ytr - mu

  fits <- list(
    Ridge      = fit_ridge(Ztr, yc),
    LASSO      = fit_lasso(Ztr, yc),
    ENET       = fit_enet (Ztr, yc, al = alpha_en),
    `MM-RWAL`  = fit_mm_rwal(Ztr, yc, c_val = c_val, delta = delta),
    `RAENR-MM` = fit_raen_mm(Ztr, yc, gamma = gamma_raen, al = alpha_raen,
                             c_val = c_val, delta = delta))

  flag    <- train_rowwise_flagger(Ztr)                       # L4
  keep_te <- if (is.null(flag)) seq_along(te) else {
    k <- which(flag$score(Zte) <= flag$thr)
    if (length(k) < 5) seq_along(te) else k
  }
  Zte_imp <- pre$impute(Zte)

  rows <- list(); sets <- list()
  for (m in METHODS) {
    b <- as.vector(fits[[m]])
    rows[[length(rows)+1]] <- data.frame(
      Rep=rep_id, p=pp, Framework="THCM", Method=m,
      t(perf_metrics(yte[keep_te], mu + as.vector(Zte[keep_te, , drop=FALSE] %*% b))),
      Selected=sum(b!=0), n_test=length(keep_te), row.names=NULL)
    rows[[length(rows)+1]] <- data.frame(
      Rep=rep_id, p=pp, Framework="ICM", Method=m,
      t(perf_metrics(yte, mu + as.vector(Zte_imp %*% b))),
      Selected=sum(b!=0), n_test=length(te), row.names=NULL)
    if (m %in% SPARSE) sets[[m]] <- pre$sel[which(b != 0)]
  }
  list(metrics = do.call(rbind, rows), sets = sets)
}
collinearity_diag <- function(pp, rep_id = 1L, tol_ratio = 1e-10) {
  set.seed(100000L + rep_id)
  te  <- sort(sample(seq_len(n_all), n_all - floor(train_prop * n_all)))
  tr  <- setdiff(seq_len(n_all), te)
  pre <- make_preprocessor(X_raw[tr, , drop=FALSE], pp)
  Ztr <- pre$apply_std(X_raw[tr, , drop=FALSE])
  sv  <- svd(Ztr, nu=0, nv=0)$d
  r   <- sum(sv > sv[1] * tol_ratio)
  p_s   <- min(300L, ncol(Ztr))
  R_sub <- cor(Ztr[, seq_len(p_s), drop=FALSE])
  ev    <- pmax(eigen(R_sub, symmetric=TRUE, only.values=TRUE)$values, 1e-12)
  data.frame(p=pp, n_train=nrow(Ztr), p_cols=ncol(Ztr), rank=r,
             kappa_rank_restricted=round((sv[1]/sv[r])^2, 2),
             max_condition_index=round(sv[1]/sv[r], 2),
             prop_sv_below_1pct=round(mean(sv < 0.01*sv[1]), 4),
             legacy_kappa_top300=round(max(ev)/min(ev), 2), row.names=NULL)
}
## ===================== data =================================================
gse_raw <- getGEO(filename = GEO_FILE, getGPL = FALSE)
eset <- if (is(gse_raw,"ExpressionSet")) gse_raw else {
  np <- sapply(gse_raw, function(x) if (is(x,"ExpressionSet")) nrow(exprs(x)) else 0L)
  gse_raw[[which.max(np)]]
}
expr_mat <- exprs(eset); feature_info <- fData(eset)

probe_ids <- function(fi) if ("ID" %in% colnames(fi)) fi$ID else rownames(fi)
safe_chr  <- function(x) if (is.list(x))
  vapply(x, function(v) paste(unlist(v), collapse=" "), character(1)) else as.character(x)
search_col <- function(fi, col, pat) {
  v <- tryCatch(safe_chr(fi[[col]]), error=function(e) character(0))
  if (!length(v)) return(character(0))
  probe_ids(fi)[grepl(pat, v, ignore.case=TRUE)]
}
cand <- c("Gene Symbol","Gene.Symbol","gene_assignment","symbol")
sym_col <- cand[cand %in% colnames(feature_info)][1]
afp_probes <- if (!is.na(sym_col)) search_col(feature_info, sym_col, "\\bAFP\\b") else character(0)
if (!length(afp_probes)) {
  known <- c("208076_s_at","208075_s_at","208074_s_at","211368_s_at")
  afp_probes <- known[known %in% rownames(expr_mat)]
}
valid <- afp_probes[afp_probes %in% rownames(expr_mat)]
if (length(valid) > 1)
  valid <- names(which.max(rowMeans(expr_mat[valid,,drop=FALSE])))
response_row <- which(rownames(expr_mat) == valid)

y_full <- as.numeric(expr_mat[response_row, ])
X_raw  <- t(expr_mat[-response_row, ])
n_all  <- nrow(X_raw)
rm(gse_raw, eset, expr_mat); gc(FALSE)
cat("data ready: n =", n_all, "genes =", ncol(X_raw), "\n")

## ===================== run ==================================================
metrics <- list(); sets <- list(); diags <- list()
t_start <- Sys.time()

for (pp in P_GRID) {
  t0 <- Sys.time()
  if (SHARD == 1L) diags[[as.character(pp)]] <- collinearity_diag(pp)
  out <- vector("list", length(MY_REPS))
  st  <- vector("list", length(MY_REPS))
  for (k in seq_along(MY_REPS)) {
    r <- run_rep(MY_REPS[k], pp)
    out[[k]] <- r$metrics; st[[k]] <- r$sets
    if (k %% 25L == 0L)
      cat(sprintf("  p=%d  %d/%d  (%.1f min elapsed)\n", pp, k, length(MY_REPS),
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
  }
  metrics[[as.character(pp)]] <- do.call(rbind, out)
  sets[[as.character(pp)]]    <- st
  saveRDS(list(shard=SHARD, nshard=NSHARD, reps=MY_REPS,
               metrics=metrics, sets=sets, diags=diags),
          file.path(OUT_DIR, sprintf("shard_%02d.rds", SHARD)))
  cat(sprintf("p = %d done in %.1f min\n", pp,
              as.numeric(difftime(Sys.time(), t0, units="mins"))))
}

cat(sprintf("\nSHARD %d COMPLETE in %.1f min -> %s\n", SHARD,
            as.numeric(difftime(Sys.time(), t_start, units="mins")),
            file.path(OUT_DIR, sprintf("shard_%02d.rds", SHARD))))
