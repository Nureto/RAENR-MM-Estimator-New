## =====================================================================
##  S8. The S2 diagnostics, redone on the CORRECTED estimator.
##
##  This replaces S2_convergence.R and S4_fixed_scale_convergence.R. Both
##  of those measured the defective fit, and both scored the objective
##  with a penalty n times heavier than the one the algorithm used. The
##  numbers they produced (0% convergence, 0 of 50 runs descending, a
##  median net change of +13.32) are artefacts and must not be reported.
##
##  Two things are fixed here:
##
##  1. THE PENALTY. cd_sweep is given fe$lambda * n, matching the patched
##     production fit_raen_mm.
##  2. THE OBJECTIVE. Q_fixed multiplies its penalty by n internally, so
##     it is called with qlam = fe$lambda, giving an objective whose
##     penalty is n * fe$lambda -- the same one the inner loop targets.
##     Objective and optimiser now refer to one function, which is what
##     makes the descent test meaningful.
##
##  WHAT IT PRODUCES, mapped to the Editor's four bullets in item S2:
##    - proportion of fits meeting the outer IRLS tolerance
##    - proportion reaching the outer-iteration cap
##    - the same for the inner coordinate-descent loop
##    - the distribution of iterations used
##  plus a descent record for the fixed-scale variant, which is a direct
##  empirical test of Theorem 2b, and a cap-sensitivity analysis to run
##  only if the cap is still commonly reached.
##
##  RUN THIS AFTER THE SHARDS FINISH. It competes for the same cores.
##
##  USAGE
##    ## fresh session, after run_all_shards_CORRECTED.bat has completed
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(GSE14520_PATH = "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
##    source("C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R")
##    s1_load_pipeline()
##    source("C:/Users/alabi/Desktop/JBE_Revision_3/S8_corrected_diagnostics.R")
##
##    R8 <- s8_run_all(p = 15000L, reps = 1:100)
##
##  Roughly the cost of the original S2 audit. Send me the printout.
## =====================================================================

stopifnot(exists("s1_front_end"), exists("s1_split_data"),
          exists("cd_sweep"), exists("Q_fixed"))


## ---------------------------------------------------------------------
##  One instrumented fit at the corrected penalty.
##
##  adaptive = TRUE  reproduces the production estimator (scale
##                   recomputed each outer iteration).
##  adaptive = FALSE holds the scale at its pilot value, which is the
##                   hypothesis of Theorem 2b. Only in this case is
##                   Q_fixed the same function at every iterate, so only
##                   here is the descent record a test of the theorem.
## ---------------------------------------------------------------------
s8_fit_record <- function(X, y, fe, adaptive = TRUE,
                          max_irls = MAX_IRLS, max_cd = MAX_CD) {
  beta <- fe$beta0
  wp   <- fe$wp
  n    <- nrow(X)
  lam  <- fe$lambda * n          # corrected penalty, as patched
  qlam <- fe$lambda              # Q_fixed multiplies by n internally
  sig0 <- fe$sigma_pilot
  al   <- alpha_raen
  cc   <- c_val

  trace <- Q_fixed(beta, X, y, sig0, wp, qlam, al, cc)
  inner_used <- integer(0)
  inner_ok   <- logical(0)
  outer_ok   <- FALSE
  used       <- 0L
  last_delta <- NA_real_
  scales     <- numeric(0)

  for (irls in seq_len(max_irls)) {
    r <- as.vector(y - X %*% beta)
    s <- if (adaptive) mad(r) else sig0
    if (!is.finite(s) || s < 1e-10) break
    scales <- c(scales, s)
    u  <- r / s
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- X * sqrt(rw)
    yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break

    b_pre   <- beta
    cd_used <- max_cd
    cd_ok   <- FALSE
    for (k in seq_len(max_cd)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) {
        cd_used <- k
        cd_ok   <- TRUE
        break
      }
    }
    inner_used <- c(inner_used, cd_used)
    inner_ok   <- c(inner_ok, cd_ok)

    trace <- c(trace, Q_fixed(beta, X, y, sig0, wp, qlam, al, cc))
    last_delta <- sqrt(sum((beta - b_pre)^2))
    used <- irls
    if (last_delta < TOL_IRLS) {
      outer_ok <- TRUE
      break
    }
  }

  steps <- diff(trace)
  list(beta = beta,
       outer_iters = used,
       outer_converged = outer_ok,
       outer_hit_cap = (!outer_ok && used >= max_irls),
       final_outer_delta = last_delta,
       inner_mean = if (length(inner_used)) mean(inner_used) else NA_real_,
       inner_max = if (length(inner_used)) max(inner_used) else NA_integer_,
       inner_conv_frac = if (length(inner_ok)) mean(inner_ok) else NA_real_,
       inner_cap_frac = if (length(inner_used)) mean(inner_used >= max_cd)
                        else NA_real_,
       n_steps = length(steps),
       n_desc = sum(steps <= 0),
       worst_increase = if (length(steps)) max(steps) else NA_real_,
       obj_change = trace[length(trace)] - trace[1],
       scale_ratio = if (length(scales) > 1)
                       max(scales) / min(scales) else 1,
       n_selected = sum(beta != 0))
}


s8_audit <- function(p = 15000L, reps = 1:100, adaptive = TRUE,
                     max_irls = MAX_IRLS) {
  out <- list()
  lab <- if (adaptive) "adaptive-scale (production)" else "fixed-scale"
  cat("S8:", lab, "audit, p =", p, "| reps", length(reps), "\n")
  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next
    r <- try(s8_fit_record(d$Ztr, d$yc, fe, adaptive = adaptive,
                           max_irls = max_irls), silent = TRUE)
    if (inherits(r, "try-error")) next
    out[[length(out) + 1L]] <- data.frame(
      rep = reps[k], p = p, cap = max_irls, adaptive = adaptive,
      outer_iters = r$outer_iters, outer_converged = r$outer_converged,
      outer_hit_cap = r$outer_hit_cap,
      final_outer_delta = r$final_outer_delta,
      inner_mean = r$inner_mean, inner_max = r$inner_max,
      inner_conv_frac = r$inner_conv_frac, inner_cap_frac = r$inner_cap_frac,
      n_steps = r$n_steps, n_desc = r$n_desc,
      worst_increase = r$worst_increase, obj_change = r$obj_change,
      scale_ratio = r$scale_ratio, n_selected = r$n_selected)
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  do.call(rbind, out)
}


s8_summary <- function(df) {
  lab <- if (df$adaptive[1]) "ADAPTIVE-SCALE (the estimator of record)" else
    "FIXED-SCALE (the object Theorem 2b covers)"
  cat("\n---", lab, "---\n")
  cat(sprintf("  fits audited                          : %d\n", nrow(df)))
  cat(sprintf("  outer loop met tolerance (%.0e)      : %d (%.1f%%)\n",
              TOL_IRLS, sum(df$outer_converged),
              100 * mean(df$outer_converged)))
  cat(sprintf("  outer loop stopped at the cap (%d)    : %d (%.1f%%)\n",
              df$cap[1], sum(df$outer_hit_cap),
              100 * mean(df$outer_hit_cap)))
  cat(sprintf("  outer iterations: median %.1f, range %d to %d\n",
              median(df$outer_iters), min(df$outer_iters),
              max(df$outer_iters)))
  cat(sprintf("  final outer change: median %.4g, max %.4g\n",
              median(df$final_outer_delta, na.rm = TRUE),
              max(df$final_outer_delta, na.rm = TRUE)))
  cat(sprintf("  inner passes meeting tolerance (%.0e): %.1f%%\n",
              TOL_CD, 100 * mean(df$inner_conv_frac, na.rm = TRUE)))
  cat(sprintf("  inner passes reaching the %d-sweep cap: %.1f%%\n",
              MAX_CD, 100 * mean(df$inner_cap_frac, na.rm = TRUE)))
  cat(sprintf("  inner sweeps per outer pass: mean %.1f, max %d\n",
              mean(df$inner_mean, na.rm = TRUE),
              max(df$inner_max, na.rm = TRUE)))
  cat(sprintf("  genes selected: mean %.1f\n", mean(df$n_selected)))
  cat(sprintf("  within-fit scale ratio: median %.2f\n",
              median(df$scale_ratio)))

  cat("\n  OBJECTIVE (matched penalty; negative change = descent)\n")
  cat(sprintf("    outer steps evaluated               : %d\n",
              sum(df$n_steps)))
  cat(sprintf("    steps at which the objective fell   : %d (%.1f%%)\n",
              sum(df$n_desc),
              100 * sum(df$n_desc) / max(1, sum(df$n_steps))))
  cat(sprintf("    runs ending below their start       : %d of %d\n",
              sum(df$obj_change < 0), nrow(df)))
  cat(sprintf("    median net change                   : %+.4g\n",
              median(df$obj_change)))
  cat(sprintf("    largest single increase             : %+.4g\n",
              max(df$worst_increase, na.rm = TRUE)))
  if (!df$adaptive[1]) {
    cat("    (for the fixed-scale variant this is a direct test of\n")
    cat("     Theorem 2b: descent should be 100%.)\n")
  } else {
    cat("    (for the adaptive variant the objective is evaluated at the\n")
    cat("     pilot scale for comparability across iterations; it is not\n")
    cat("     a function the iteration is guaranteed to decrease.)\n")
  }
  invisible(df)
}


s8_cap_sensitivity <- function(p = 15000L, reps = 1:50, caps = c(5L, 20L)) {
  res <- data.frame()
  cat("S8: cap sensitivity at the corrected penalty, p =", p,
      "| caps", paste(caps, collapse = " vs "), "\n")
  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next
    fa <- try(s8_fit_record(d$Ztr, d$yc, fe, TRUE, max_irls = caps[1]),
              silent = TRUE)
    fb <- try(s8_fit_record(d$Ztr, d$yc, fe, TRUE, max_irls = caps[2]),
              silent = TRUE)
    if (inherits(fa, "try-error") || inherits(fb, "try-error")) next

    flag <- train_rowwise_flagger(d$Ztr)
    keep <- if (is.null(flag)) seq_along(d$te) else {
      kk <- which(flag$score(d$Zte) <= flag$thr)
      if (length(kk) < 5) seq_along(d$te) else kk
    }
    rmse <- function(b) {
      pr <- d$mu + as.vector(d$Zte[keep, , drop = FALSE] %*% b)
      sqrt(mean((d$yte[keep] - pr)^2))
    }
    sa <- which(fa$beta != 0); sb <- which(fb$beta != 0)
    res <- rbind(res, data.frame(
      rep = reps[k], p = p,
      rmse_lo = rmse(fa$beta), rmse_hi = rmse(fb$beta),
      sel_lo = length(sa), sel_hi = length(sb),
      jaccard = length(intersect(sa, sb)) / max(1, length(union(sa, sb))),
      beta_rel_change = sqrt(sum((fb$beta - fa$beta)^2)) /
        max(1e-12, sqrt(sum(fa$beta^2)))))
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  dd <- res$rmse_hi - res$rmse_lo
  tt <- t.test(dd)
  cat("\n--- cap sensitivity, corrected penalty ---\n")
  cat(sprintf("  replications                  : %d\n", nrow(res)))
  cat(sprintf("  mean RMSE at cap %-3d          : %.4f\n",
              caps[1], mean(res$rmse_lo)))
  cat(sprintf("  mean RMSE at cap %-3d          : %.4f\n",
              caps[2], mean(res$rmse_hi)))
  cat(sprintf("  paired difference             : %+.5f\n", mean(dd)))
  cat(sprintf("  95%% CI                        : %+.5f to %+.5f\n",
              tt$conf.int[1], tt$conf.int[2]))
  cat(sprintf("  p-value                       : %.3f\n", tt$p.value))
  cat(sprintf("  selected: %.1f vs %.1f | Jaccard %.3f\n",
              mean(res$sel_lo), mean(res$sel_hi), mean(res$jaccard)))
  cat(sprintf("  median relative change in beta: %.4f\n",
              median(res$beta_rel_change)))
  invisible(res)
}


s8_write <- function(R8, dir = ".") {
  if (!is.null(R8$adaptive))
    write.csv(R8$adaptive, file.path(dir, "S8_adaptive_audit.csv"),
              row.names = FALSE)
  if (!is.null(R8$fixed))
    write.csv(R8$fixed, file.path(dir, "S8_fixed_audit.csv"),
              row.names = FALSE)
  if (!is.null(R8$cap))
    write.csv(R8$cap, file.path(dir, "S8_cap_sensitivity.csv"),
              row.names = FALSE)
  cat("wrote the S8 CSVs\n")
}


s8_run_all <- function(p = 15000L, reps = 1:100, cap_reps = 1:50,
                       dir = ".") {
  A <- s8_audit(p = p, reps = reps, adaptive = TRUE)
  s8_summary(A)
  F <- s8_audit(p = p, reps = reps, adaptive = FALSE)
  s8_summary(F)

  cap <- NULL
  hit <- mean(A$outer_hit_cap)
  cat(sprintf("\n  Outer cap reached in %.1f%% of fits.\n", 100 * hit))
  if (hit > 0.20) {
    cat("  Above 20%, so running the cap-sensitivity analysis the Editor\n")
    cat("  asks for in item S2.\n")
    cap <- s8_cap_sensitivity(p = p, reps = cap_reps)
  } else {
    cat("  Below 20%, so no cap-sensitivity analysis is required; we say\n")
    cat("  so and give this proportion.\n")
  }

  R8 <- list(adaptive = A, fixed = F, cap = cap)
  s8_write(R8, dir = dir)
  invisible(R8)
}


cat("Loaded S8 corrected diagnostics.\n")
cat("Run AFTER the shards finish:\n\n")
cat("  R8 <- s8_run_all(p = 15000L, reps = 1:100)\n")
