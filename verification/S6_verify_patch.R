## =====================================================================
##  S6. Verify the one-line correction to fit_raen_mm, end to end.
##
##  S5 established, on a reconstruction of the fitting path, that
##  cd_sweep is handed a lambda selected by cv.glmnet without the factor
##  n that the loss scaling requires. This script checks the correction
##  against the PRODUCTION function in run_shard.R, not a reconstruction,
##  and on the full fitting path including its own cv.glmnet pilot.
##
##  THE CORRECTION
##  --------------
##  In run_shard.R, inside fit_raen_mm, the inner loop currently reads
##
##      beta <- cd_sweep(Xw, yw, beta, wp, lambda_fixed, al)
##
##  and should read
##
##      beta <- cd_sweep(Xw, yw, beta, wp, lambda_fixed * n, al)
##
##  n is already in scope (line 2 of the function). Nothing else changes.
##
##  WHY n AND NOT sum(rw)
##  ---------------------
##  lambda_fixed comes from cv.glmnet(X, y, ...), which minimises
##  1/(2n)||y - Xb||^2 + lambda[...] over the n = 356 training rows.
##  Rewriting that with a 1/2 loss gives a penalty of n*lambda. The
##  weighting that produces Xw does not change the row count, so the
##  factor is n, not sum(rw).
##
##  WHAT THIS SCRIPT CHECKS
##  -----------------------
##  Both arms are given identical CV folds by seeding immediately before
##  each call, so the pilot, lambda and adaptive weights are shared and
##  the penalty scaling is the only difference.
##
##  Expect, from S5 at p = 15,000: about 950 genes selected in the
##  current arm and about 47 in the corrected one, with test RMSE
##  statistically indistinguishable. If this script reproduces that on
##  the production path, the diagnosis is settled and the patch is safe
##  to apply.
##
##  USAGE
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(GSE14520_PATH = "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
##    source("C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R")
##    s1_load_pipeline()
##    source("C:/Users/alabi/Desktop/JBE_Revision_3/S6_verify_patch.R")
##
##    V <- s6_verify(p = 15000L, reps = 1:20)
##    s6_summary(V)
##    s6_write(V)
##
##  Then, if it checks out, the staged rerun plan is in
##  PATCH_AND_IMPACT.md. Do not start the full rerun before this passes.
## =====================================================================

stopifnot(exists("fit_raen_mm"), exists("s1_split_data"),
          exists("cd_sweep"), exists("safe_coef_vec"))


## ---------------------------------------------------------------------
##  Production fit_raen_mm, copied verbatim from run_shard.R, with the
##  single change marked <<< PATCH >>>. Kept as a copy so the original
##  file is untouched until you decide to apply it.
## ---------------------------------------------------------------------
fit_raen_mm_patched <- function(X, y, gamma = 2, c_val = 4.685,
                                delta = 1e-6, al = 0.2,
                                max_irls = MAX_IRLS, max_cd = MAX_CD,
                                tol_irls = TOL_IRLS, tol_cd = TOL_CD) {
  n <- nrow(X); p <- ncol(X)
  if (!is.finite(sd(y)) || sd(y) < 1e-8) return(rep(0, p))

  pilot_en <- tryCatch(cv.glmnet(X, y, alpha = al, nfolds = NFOLD,
                                 nlambda = NLAM),
                       error = function(e) NULL)
  if (is.null(pilot_en)) return(rep(0, p))
  lambda_fixed <- pilot_en$lambda.1se
  beta <- safe_coef_vec(as.vector(coef(pilot_en, s = "lambda.1se"))[-1], p)

  beta_mm <- tryCatch({
    if (p < n) {
      mm_fit <- suppressWarnings(
        lmrob(y ~ X - 1, setting = "KS2014", k.max = 200,
              maxit.scale = 200))
      out <- as.vector(coef(mm_fit))
      if (length(out) != p || any(!is.finite(out))) rep(0, p) else out
    } else {
      b_h <- beta; lam_h <- lambda_fixed
      for (hs in seq_len(3)) {
        r_h <- as.vector(y - X %*% b_h)
        s_h <- mad(r_h); if (!is.finite(s_h) || s_h < 1e-10) break
        u_h <- r_h / s_h; rw_h <- pmin(1, 1.345 / abs(u_h))
        rw_h[!is.finite(rw_h)] <- 1
        Xw_h <- (X * sqrt(rw_h)); yw_h <- y * sqrt(rw_h)
        fg_h <- tryCatch(glmnet(Xw_h, yw_h, alpha = al, lambda = lam_h,
                                standardize = FALSE),
                         error = function(e) NULL)
        if (is.null(fg_h)) break
        bn_h <- as.vector(coef(fg_h))[-1]
        if (length(bn_h) != p || any(!is.finite(bn_h))) break
        if (sqrt(sum((bn_h - b_h)^2)) < 1e-3) { b_h <- bn_h; break }
        b_h <- bn_h
      }
      b_h
    }
  }, error = function(e) beta)

  wp <- 1 / (abs(beta_mm)^gamma + delta)
  wp <- pmin(wp, quantile(wp, 0.95))
  wp <- wp / mean(wp)

  for (irls in seq_len(max_irls)) {
    r <- as.vector(y - X %*% beta)
    sigma_hat <- mad(r)
    if (!is.finite(sigma_hat) || sigma_hat < 1e-10) break
    u  <- r / sigma_hat
    rw <- ifelse(abs(u) <= c_val, (1 - (u / c_val)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- (X * sqrt(rw)); yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break

    beta_pre_cd <- beta
    for (cd_iter in seq_len(max_cd)) {
      beta_old <- beta
      ## <<< PATCH >>> lambda_fixed -> lambda_fixed * n
      beta <- cd_sweep(Xw, yw, beta, wp, lambda_fixed * n, al)
      if (sqrt(sum((beta - beta_old)^2)) < tol_cd) break
    }
    if (sqrt(sum((beta - beta_pre_cd)^2)) < tol_irls) break
  }
  beta
}


## ---------------------------------------------------------------------
##  Paired comparison on the production path.
## ---------------------------------------------------------------------
s6_verify <- function(p = 15000L, reps = 1:20) {
  res <- data.frame()
  cat("S6: production-path patch verification, p =", p,
      "| reps", length(reps), "\n")

  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next

    ## identical CV folds in both arms
    set.seed(500000L + reps[k])
    b_cur <- try(fit_raen_mm(d$Ztr, d$yc, gamma = gamma_raen,
                             al = alpha_raen, c_val = c_val,
                             delta = delta), silent = TRUE)
    set.seed(500000L + reps[k])
    b_new <- try(fit_raen_mm_patched(d$Ztr, d$yc, gamma = gamma_raen,
                                     al = alpha_raen, c_val = c_val,
                                     delta = delta), silent = TRUE)
    if (inherits(b_cur, "try-error") || inherits(b_new, "try-error")) next

    flag <- train_rowwise_flagger(d$Ztr)
    keep <- if (is.null(flag)) {
      seq_along(d$te)
    } else {
      kk <- which(flag$score(d$Zte) <= flag$thr)
      if (length(kk) < 5) seq_along(d$te) else kk
    }
    rmse <- function(b) {
      pr <- d$mu + as.vector(d$Zte[keep, , drop = FALSE] %*% b)
      sqrt(mean((d$yte[keep] - pr)^2))
    }

    sa <- which(b_cur != 0); sb <- which(b_new != 0)
    res <- rbind(res, data.frame(
      rep = reps[k], p = p,
      sel_current = length(sa), sel_patched = length(sb),
      rmse_current = rmse(b_cur), rmse_patched = rmse(b_new),
      jaccard = length(intersect(sa, sb)) / max(1, length(union(sa, sb))),
      patched_subset_frac = if (length(sb)) {
        length(intersect(sa, sb)) / length(sb)
      } else NA_real_))

    if (k %% 5L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  res
}


s6_summary <- function(df) {
  dd <- df$rmse_patched - df$rmse_current
  tt <- t.test(dd)
  cat("\n--- S6 result: production path, current vs patched ---\n")
  cat(sprintf("  replications                    : %d\n", nrow(df)))
  cat(sprintf("  genes selected, current         : %.1f\n",
              mean(df$sel_current)))
  cat(sprintf("  genes selected, patched         : %.1f\n",
              mean(df$sel_patched)))
  cat(sprintf("  mean RMSE, current              : %.4f\n",
              mean(df$rmse_current)))
  cat(sprintf("  mean RMSE, patched              : %.4f\n",
              mean(df$rmse_patched)))
  cat(sprintf("  paired difference (patched-cur) : %+.5f\n", mean(dd)))
  cat(sprintf("  95%% CI                          : %+.5f to %+.5f\n",
              tt$conf.int[1], tt$conf.int[2]))
  cat(sprintf("  p-value                         : %.3f\n", tt$p.value))
  cat(sprintf("  support Jaccard                 : %.3f\n",
              mean(df$jaccard)))
  cat(sprintf("  patched panel inside current    : %.1f%%\n",
              100 * mean(df$patched_subset_frac, na.rm = TRUE)))

  cat("\n  Expected from S5 at p = 15,000: about 950 vs about 47 genes,\n")
  cat("  RMSE indistinguishable. If that is what you see, the diagnosis\n")
  cat("  is settled on the production path and the patch is safe.\n")
  cat("  The last line is a sanity check: the corrected panel should be\n")
  cat("  largely a SUBSET of the over-dense one, not a different set.\n")
  invisible(df)
}


s6_write <- function(V, dir = ".") {
  write.csv(V, file.path(dir, "S6_patch_verification_raw.csv"),
            row.names = FALSE)
  cat("wrote S6_patch_verification_raw.csv\n")
}


cat("Loaded S6 patch verification.\n")
cat("  V <- s6_verify(p = 15000L, reps = 1:20)\n")
cat("  s6_summary(V)\n")
cat("  s6_write(V)\n")
