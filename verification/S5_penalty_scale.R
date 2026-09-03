## =====================================================================
##  S5. Is the penalty in the inner solve scaled by n or not?
##
##  WHY THIS EXISTS
##  ---------------
##  The S4 fixed-scale audit returned almost exactly the S5-table numbers
##  for the adaptive algorithm: 79.8% vs 80.0% of steps descending, 0 of
##  100 vs 0 of 50 runs ending below their start, and a largest single
##  increase of +56.04 in BOTH. Identical to four figures is not a
##  coincidence.
##
##  It is not one. At the first outer iteration the two variants are the
##  same algorithm: the adaptive rule sets s = mad(y - X beta0) and the
##  fixed rule uses sigma_pilot, which is the same quantity. Table S5
##  also records that every observed increase occurs at iteration 1. So
##  both audits are measuring one event, the first inner solve, twice.
##
##  That points at the inner solve, and there the arithmetic does not
##  line up. The coordinate update in cd_sweep is
##
##      b_j = S( z_j / den , lam * alpha * wp_j / den ),
##      den = d_j + lam * (1 - alpha)
##
##  which is the stationarity condition of
##
##      F(b) = 1/2 ||yw - Xw b||^2
##             + lam * alpha * SUM_j wp_j |b_j|
##             + lam * (1 - alpha) / 2 * ||b||^2 .
##
##  The manuscript's objective, equation (19), carries a factor n on both
##  penalty terms, and Q_fixed() in S1_scale_verification.R was written
##  from equation (19). The MM surrogate of that objective's loss is
##  exactly 1/2||yw - Xw b||^2, so the loss scalings agree and only the
##  penalty differs -- by a factor of n.
##
##  This matters twice over:
##
##  1. lam comes from cv.glmnet, which minimises
##     1/(2n)||y - Xb||^2 + lam[...]. Written with a 1/2 loss that is
##     n*lam[...]. So the lambda selected by cross-validation belongs
##     with a penalty of n*lam. As implemented the penalty is lam, i.e.
##     about n = 356 times weaker than the tuning implies.
##
##  2. The descent diagnostics in Table S5 and in the S4 run above
##     evaluate Q_fixed, which uses n*lam. If the algorithm optimises the
##     lam version, those tables compare an iterate against an objective
##     it never targeted, and an apparent "increase" is expected rather
##     than anomalous.
##
##  TWO SEPARATE CLAIMS. THIS SCRIPT TESTS THEM SEPARATELY.
##  -------------------------------------------------------
##  CLAIM 1 -- the descent diagnostic is mis-specified. Q_fixed multiplies
##  its penalty by n internally, while cd_sweep uses lam. So Table S5 and
##  the S4 run compared each iterate against an objective the algorithm
##  never optimised. Majorize-minimize theory says that when the scale is
##  held fixed and the inner loop decreases the surrogate OF THE SAME
##  objective, that objective cannot rise. So evaluating the MATCHED
##  objective should show descent. This claim is arithmetic, and the run
##  either confirms it or refutes it outright.
##
##  CLAIM 2 -- the implementation may under-penalise. lam comes from
##  cv.glmnet, whose lambda belongs with a penalty of n*lam once the loss
##  is written as 1/2||.||^2. If so the fit is penalised about n = 356
##  times too weakly, which would predict implausibly dense "sparse"
##  solutions -- and the paper reports 678-855 genes against 3.4-8.8 for
##  MM-RWAL. This claim is about intent, not arithmetic, and the run
##  informs it rather than settling it.
##
##  Each fit is therefore scored against BOTH conventions:
##    obj_matched     penalty in Q equals the penalty cd_sweep used
##    obj_asreported  penalty in Q is n*lam, the Table S5 convention
##
##  READING THE OUTPUT
##  ------------------
##  (a) mult = 1 descends under obj_matched but not under obj_asreported
##      -> Claim 1 confirmed. Table S5 and the S4 result above are
##      artefacts of a mis-specified diagnostic, and both must be
##      withdrawn and recomputed. The S2 convergence and cap-sensitivity
##      findings are untouched: neither uses Q_fixed.
##  (b) mult = 1 fails to descend under obj_matched too -> Claim 1 is
##      wrong, Table S5 stands, and the inner solve is at fault.
##  Compare sel_code against sel_eq19 for Claim 2.
##
##  USAGE
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(GSE14520_PATH = "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
##    source("C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R")
##    s1_load_pipeline()
##    source("C:/Users/alabi/Desktop/JBE_Revision_3/S5_penalty_scale.R")
##
##    s5_thresholds(p = 15000L, reps = 1:3)     # arithmetic, seconds
##    PS <- s5_compare(p = 15000L, reps = 1:20) # the decisive run
##    s5_summary(PS)
##    s5_write(PS)
##
##  COST. Cheap. Two fits per replication at the normal 5 x 50 budget.
## =====================================================================

stopifnot(exists("s1_front_end"), exists("s1_split_data"),
          exists("cd_sweep"), exists("Q_fixed"))


## ---------------------------------------------------------------------
##  Report the raw magnitudes. No fitting, no interpretation.
## ---------------------------------------------------------------------
s5_thresholds <- function(p = 15000L, reps = 1:3) {
  cat("S5: penalty magnitudes at p =", p, "\n\n")
  for (k in reps) {
    d <- try(s1_split_data(k, p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next

    X <- d$Ztr
    n <- nrow(X)
    lam <- fe$lambda
    al <- alpha_raen

    r <- as.vector(d$yc - X %*% fe$beta0)
    s <- mad(r)
    u <- r / s
    rw <- ifelse(abs(u) <= c_val, (1 - (u / c_val)^2)^2, 0)
    Xw <- X * sqrt(rw)
    dj <- colSums(Xw^2)

    den <- median(dj) + lam * (1 - al)
    thr_code <- lam * al * 1 / den
    thr_eq19 <- n * lam * al * 1 / den

    cat(sprintf("  rep %-3d  n = %d  lambda = %.6g\n", k, n, lam))
    cat(sprintf("           median d_j = %.4g   sum(rw) = %.1f\n",
                median(dj), sum(rw)))
    cat(sprintf("           soft threshold as implemented : %.3e\n",
                thr_code))
    cat(sprintf("           soft threshold under Eq. (19) : %.3e\n",
                thr_eq19))
    cat(sprintf("           ratio                         : %.1f\n\n",
                thr_eq19 / max(thr_code, 1e-300)))
  }
  cat("  A threshold of order 1e-4 on standardized columns removes almost\n")
  cat("  nothing, which is the signature of an under-penalised fit.\n")
}


## ---------------------------------------------------------------------
##  One fixed-scale fit at a given penalty multiplier, recording the
##  Q_fixed trajectory. mult = 1 is the implementation; mult = n is
##  equation (19).
## ---------------------------------------------------------------------
s5_fit_record <- function(X, y, fe, mult, max_irls = MAX_IRLS,
                          max_cd = MAX_CD) {
  beta <- fe$beta0
  wp   <- fe$wp
  lam  <- fe$lambda * mult
  sig  <- fe$sigma_pilot
  al   <- alpha_raen
  cc   <- c_val

  ## Two scorings of every iterate.
  ##   qlam_m : Q's internal factor n cancels, so Q's penalty equals the
  ##            penalty cd_sweep actually used. Objective and optimiser
  ##            agree; MM predicts monotone descent.
  ##   qlam_o : the Table S5 convention, Q penalty = n * fe$lambda,
  ##            regardless of what the inner loop used.
  qlam_m <- fe$lambda * mult / nrow(X)
  qlam_o <- fe$lambda
  trace_m <- Q_fixed(beta, X, y, sig, wp, qlam_m, al, cc)
  trace_o <- Q_fixed(beta, X, y, sig, wp, qlam_o, al, cc)
  inner_ok <- logical(0)
  inner_used <- integer(0)
  outer_ok <- FALSE
  used <- 0L
  last_delta <- NA_real_

  for (irls in seq_len(max_irls)) {
    r  <- as.vector(y - X %*% beta)
    u  <- r / sig
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- X * sqrt(rw)
    yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break

    b_pre <- beta
    cd_used <- max_cd
    cd_ok <- FALSE
    for (k in seq_len(max_cd)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) {
        cd_used <- k
        cd_ok <- TRUE
        break
      }
    }
    inner_used <- c(inner_used, cd_used)
    inner_ok <- c(inner_ok, cd_ok)

    trace_m <- c(trace_m, Q_fixed(beta, X, y, sig, wp, qlam_m, al, cc))
    trace_o <- c(trace_o, Q_fixed(beta, X, y, sig, wp, qlam_o, al, cc))
    last_delta <- sqrt(sum((beta - b_pre)^2))
    used <- irls
    if (last_delta < TOL_IRLS) {
      outer_ok <- TRUE
      break
    }
  }

  sm <- diff(trace_m)
  so <- diff(trace_o)
  list(beta = beta,
       n_steps = length(sm),
       n_desc = sum(sm <= 0),
       n_desc_asrep = sum(so <= 0),
       worst_increase = if (length(sm)) max(sm) else NA_real_,
       worst_increase_asrep = if (length(so)) max(so) else NA_real_,
       obj_change = trace_m[length(trace_m)] - trace_m[1],
       obj_change_asrep = trace_o[length(trace_o)] - trace_o[1],
       outer_iters = used,
       outer_converged = outer_ok,
       final_outer_delta = last_delta,
       inner_conv_frac = mean(inner_ok),
       inner_mean = mean(inner_used),
       n_selected = sum(beta != 0))
}


## ---------------------------------------------------------------------
##  Paired comparison: penalty as implemented vs penalty as defined.
## ---------------------------------------------------------------------
s5_compare <- function(p = 15000L, reps = 1:20) {
  res <- data.frame()
  cat("S5: penalty-scale comparison, p =", p,
      "| reps", length(reps), "\n")

  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next
    nn <- nrow(d$Ztr)

    fa <- try(s5_fit_record(d$Ztr, d$yc, fe, mult = 1), silent = TRUE)
    fb <- try(s5_fit_record(d$Ztr, d$yc, fe, mult = nn), silent = TRUE)
    if (inherits(fa, "try-error") || inherits(fb, "try-error")) next

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

    res <- rbind(res, data.frame(
      rep = reps[k], p = p, n = nn,
      desc_frac_code = fa$n_desc / max(1, fa$n_steps),
      desc_frac_eq19 = fb$n_desc / max(1, fb$n_steps),
      desc_frac_asrep = fa$n_desc_asrep / max(1, fa$n_steps),
      obj_change_code = fa$obj_change,
      obj_change_eq19 = fb$obj_change,
      obj_change_asrep = fa$obj_change_asrep,
      worst_inc_code = fa$worst_increase,
      worst_inc_eq19 = fb$worst_increase,
      worst_inc_asrep = fa$worst_increase_asrep,
      outer_conv_code = fa$outer_converged,
      outer_conv_eq19 = fb$outer_converged,
      outer_iters_code = fa$outer_iters,
      outer_iters_eq19 = fb$outer_iters,
      inner_conv_code = fa$inner_conv_frac,
      inner_conv_eq19 = fb$inner_conv_frac,
      sel_code = fa$n_selected, sel_eq19 = fb$n_selected,
      rmse_code = rmse(fa$beta), rmse_eq19 = rmse(fb$beta)))

    if (k %% 5L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  res
}


s5_summary <- function(df) {
  cat("\n--- S5 result: penalty as implemented vs penalty as defined ---\n")
  cat(sprintf("  replications : %d   (n = %d)\n", nrow(df), df$n[1]))

  cat("\n  CLAIM 1. THE SAME FITS, SCORED TWO WAYS.\n")
  cat("  Both rows below are the penalty AS IMPLEMENTED (lam). Only the\n")
  cat("  objective used to score them differs.\n")
  cat(sprintf("    steps descending, matched objective   : %.1f%%\n",
              100 * mean(df$desc_frac_code)))
  cat(sprintf("    steps descending, Table S5 convention : %.1f%%\n",
              100 * mean(df$desc_frac_asrep)))
  cat(sprintf("    runs ending below start, matched      : %d of %d\n",
              sum(df$obj_change_code < 0), nrow(df)))
  cat(sprintf("    runs ending below start, Table S5     : %d of %d\n",
              sum(df$obj_change_asrep < 0), nrow(df)))
  cat(sprintf("    largest increase, matched             : %+.4g\n",
              max(df$worst_inc_code, na.rm = TRUE)))
  cat(sprintf("    largest increase, Table S5            : %+.4g\n",
              max(df$worst_inc_asrep, na.rm = TRUE)))
  cat("    (the Table S5 row should reproduce ~80%% and 0 runs below\n")
  cat("     start. If it does, the diagnostic is confirmed as the\n")
  cat("     source of that number.)\n")

  cat("\n  CLAIM 2. PENALTY AS DEFINED BY EQ. (19), n*lam.\n")
  cat(sprintf("    steps descending, matched objective   : %.1f%%\n",
              100 * mean(df$desc_frac_eq19)))
  cat(sprintf("    runs ending below start                : %d of %d\n",
              sum(df$obj_change_eq19 < 0), nrow(df)))
  cat(sprintf("    largest increase                       : %+.4g\n",
              max(df$worst_inc_eq19, na.rm = TRUE)))

  cat("\n  CONVERGENCE\n")
  cat(sprintf("    outer met tolerance, lam   : %.1f%%\n",
              100 * mean(df$outer_conv_code)))
  cat(sprintf("    outer met tolerance, n*lam : %.1f%%\n",
              100 * mean(df$outer_conv_eq19)))
  cat(sprintf("    inner met tolerance, lam   : %.1f%%\n",
              100 * mean(df$inner_conv_code)))
  cat(sprintf("    inner met tolerance, n*lam : %.1f%%\n",
              100 * mean(df$inner_conv_eq19)))

  cat("\n  SPARSITY AND PREDICTION\n")
  cat(sprintf("    genes selected, lam        : %.1f\n", mean(df$sel_code)))
  cat(sprintf("    genes selected, n*lam      : %.1f\n", mean(df$sel_eq19)))
  cat(sprintf("    mean RMSE, lam             : %.4f\n", mean(df$rmse_code)))
  cat(sprintf("    mean RMSE, n*lam           : %.4f\n", mean(df$rmse_eq19)))
  dd <- df$rmse_eq19 - df$rmse_code
  tt <- t.test(dd)
  cat(sprintf("    paired difference          : %+.5f\n", mean(dd)))
  cat(sprintf("    95%% CI                     : %+.5f to %+.5f\n",
              tt$conf.int[1], tt$conf.int[2]))
  cat(sprintf("    p-value                    : %.3f\n", tt$p.value))

  cat("\n  If the n*lam arm descends at every step and the lam arm does\n")
  cat("  not, the diagnosis is confirmed and Table S5 is measuring the\n")
  cat("  wrong objective. Send me the output before changing anything.\n")
  invisible(df)
}


s5_write <- function(PS, dir = ".") {
  write.csv(PS, file.path(dir, "S5_penalty_scale_raw.csv"), row.names = FALSE)
  cat("wrote S5_penalty_scale_raw.csv\n")
}


cat("Loaded S5 penalty-scale check.\n")
cat("  s5_thresholds(p = 15000L, reps = 1:3)\n")
cat("  PS <- s5_compare(p = 15000L, reps = 1:20)\n")
cat("  s5_summary(PS)\n")
cat("  s5_write(PS)\n")
