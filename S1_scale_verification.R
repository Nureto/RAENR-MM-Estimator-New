## =====================================================================
##  S1: scale handling -- descent record and fixed-scale sensitivity
##
##  Editor's finding, confirmed: analysis/run_shard.R line 435 recomputes
##  the scale inside the outer IRLS loop, so the implemented estimator is
##  ADAPTIVE-scale and the fixed-scale monotone-descent theorem does not
##  apply to it. This script does not argue otherwise. It produces:
##
##    PART A  an empirical descent record for the algorithm as implemented
##    PART B  a paired fixed-scale vs adaptive-scale sensitivity analysis
##
##  Neither is a rerun of the 1,000-replication study.
##
##  ---------------------------------------------------------------
##  WHY THIS SCRIPT LOADS run_shard.R THE WAY IT DOES
##  ---------------------------------------------------------------
##  Sourcing run_shard.R directly would launch the whole 1,000-rep
##  analysis (the loop at the bottom of that file). s1_load_pipeline()
##  therefore evaluates run_shard.R only as far as the "## === run ==="
##  marker: data are loaded, the Rcpp kernel is compiled and bound, and
##  every function and constant is defined -- but nothing runs.
##
##  Everything below reuses YOUR machinery: the same split rule
##  (set.seed(100000L + rep_id)), the same make_preprocessor(), the same
##  train_rowwise_flagger(), the same perf_metrics(), and the same
##  constants (NFOLD = 3, NLAM = 20, MAX_IRLS = 5, MAX_CD = 50,
##  TOL_CD = 1e-4, alpha_raen = 0.2, gamma_raen = 2). The only thing that
##  changes between the two arms is where mad() is called.
##
##  ---------------------------------------------------------------
##  USAGE
##  ---------------------------------------------------------------
##    setwd("<the analysis/ directory>")
##    source("S1_scale_verification.R")
##    s1_load_pipeline()                       # a few minutes, one-off
##
##    A <- s1_descent_summary(p = 15000L, reps = 1:50)
##    B <- s1_sensitivity(p = 15000L, reps = 1:100)
##    B2 <- s1_sensitivity(p = 20000L, reps = 1:100)
##    s1_write_supplement(A, B, B2)
## =====================================================================


## ---------------------------------------------------------------------
##  Load definitions and data from run_shard.R WITHOUT running it
## ---------------------------------------------------------------------
s1_load_pipeline <- function(path = "run_shard.R") {
  stopifnot(file.exists(path))
  txt <- readLines(path, warn = FALSE)
  cut <- grep("^## =+ run =+", txt)
  if (!length(cut)) cut <- grep("^## =====================\\s*run", txt)
  if (!length(cut)) stop("could not find the '## === run ===' marker in ", path)
  head_txt <- txt[seq_len(cut[1] - 1L)]
  tf <- tempfile(fileext = ".R"); writeLines(head_txt, tf)
  cat("Evaluating run_shard.R definitions only (", cut[1] - 1L,
      " of ", length(txt), " lines)...\n", sep = "")
  sys.source(tf, envir = globalenv())
  unlink(tf)
  need <- c("fit_raen_mm", "cd_sweep", "safe_coef_vec", "make_preprocessor",
            "train_rowwise_flagger", "perf_metrics", "X_raw", "y_full",
            "n_all", "alpha_raen", "gamma_raen", "c_val", "delta",
            "NLAM", "NFOLD", "MAX_IRLS", "MAX_CD", "TOL_CD", "TOL_IRLS",
            "train_prop")
  miss <- need[!vapply(need, exists, logical(1), envir = globalenv())]
  if (length(miss))
    stop("missing after load: ", paste(miss, collapse = ", "))
  cat("Pipeline ready. n =", n_all, " genes =", ncol(X_raw), "\n")
  invisible(TRUE)
}


## ---------------------------------------------------------------------
##  The objective, on the manuscript's parametrisation
##    Q(b) = s^2 * sum_i rho_c(r_i/s) + n*lam*al*sum_j wp_j|b_j|
##                                     + n*lam*(1-al)/2 * sum_j b_j^2
##  s is supplied and HELD FIXED, so values are comparable across
##  iterations. That is the whole point of the diagnostic.
## ---------------------------------------------------------------------
rho_bisquare <- function(u, cc = 4.685)
  ifelse(abs(u) <= cc, (cc^2 / 6) * (1 - (1 - (u / cc)^2)^3), cc^2 / 6)

Q_fixed <- function(beta, X, y, s, wp, lam, al, cc = 4.685) {
  n <- nrow(X); r <- as.vector(y - X %*% beta)
  s^2 * sum(rho_bisquare(r / s, cc)) +
    n * lam * al * sum(wp * abs(beta)) +
    n * lam * (1 - al) / 2 * sum(beta^2)
}


## ---------------------------------------------------------------------
##  Shared front end: STEPS 1-3 of fit_raen_mm, verbatim
##  Returns lambda, the seed beta, the adaptive weights, and the pilot
##  scale (the sigma a fixed-scale implementation would use).
## ---------------------------------------------------------------------
s1_front_end <- function(X, y, gamma = gamma_raen, al = alpha_raen,
                         cc = c_val, dl = delta) {
  p <- ncol(X)
  pilot_en <- cv.glmnet(X, y, alpha = al, nfolds = NFOLD, nlambda = NLAM)
  lam  <- pilot_en$lambda.1se
  beta <- safe_coef_vec(as.vector(coef(pilot_en, s = "lambda.1se"))[-1], p)

  b_h <- beta
  for (hs in seq_len(3L)) {
    r_h <- as.vector(y - X %*% b_h)
    s_h <- mad(r_h); if (!is.finite(s_h) || s_h < 1e-10) break
    u_h <- r_h / s_h; rw_h <- pmin(1, 1.345 / abs(u_h)); rw_h[!is.finite(rw_h)] <- 1
    fg <- tryCatch(glmnet(X * sqrt(rw_h), y * sqrt(rw_h), alpha = al,
                          lambda = lam, standardize = FALSE),
                   error = function(e) NULL)
    if (is.null(fg)) break
    bn <- as.vector(coef(fg))[-1]
    if (length(bn) != p || any(!is.finite(bn))) break
    if (sqrt(sum((bn - b_h)^2)) < 1e-3) { b_h <- bn; break }
    b_h <- bn
  }
  wp <- 1 / (abs(b_h)^gamma + dl)
  wp <- pmin(wp, quantile(wp, 0.95)); wp <- wp / mean(wp)

  s0 <- mad(as.vector(y - X %*% beta))
  if (!is.finite(s0) || s0 < 1e-10) s0 <- mad(y)
  if (!is.finite(s0) || s0 < 1e-10) s0 <- 1

  list(lambda = lam, beta0 = beta, wp = wp, sigma_pilot = s0)
}


## =====================================================================
##  PART A -- descent record, algorithm exactly as implemented
## =====================================================================
s1_descent_record <- function(X, y, fe = NULL) {
  if (is.null(fe)) fe <- s1_front_end(X, y)
  beta <- fe$beta0; wp <- fe$wp; lam <- fe$lambda; s0 <- fe$sigma_pilot
  al <- alpha_raen; cc <- c_val

  rec <- data.frame(iter = 0L, sigma_used = NA_real_,
                    Q = Q_fixed(beta, X, y, s0, wp, lam, al, cc),
                    n_sel = sum(beta != 0))

  for (irls in seq_len(MAX_IRLS)) {
    r <- as.vector(y - X %*% beta)
    sigma_hat <- mad(r)                      # ADAPTIVE -- as in run_shard.R
    if (!is.finite(sigma_hat) || sigma_hat < 1e-10) break
    u  <- r / sigma_hat
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- X * sqrt(rw); yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break
    for (cd_iter in seq_len(MAX_CD)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) break
    }
    rec <- rbind(rec, data.frame(
      iter = irls, sigma_used = sigma_hat,
      Q = Q_fixed(beta, X, y, s0, wp, lam, al, cc),
      n_sel = sum(beta != 0)))
  }
  rec$dQ <- c(NA, diff(rec$Q))
  rec$decreased <- rec$dQ <= 1e-8
  rec
}


s1_descent_summary <- function(p = 15000L, reps = 1:50) {
  out <- list()
  cat("PART A: descent record, p =", p, ", reps =", length(reps), "\n")
  for (k in seq_along(reps)) {
    d <- s1_split_data(reps[k], p)
    r <- try(s1_descent_record(d$Ztr, d$yc), silent = TRUE)
    if (inherits(r, "try-error")) next
    r$rep <- reps[k]; out[[length(out) + 1L]] <- r
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  df <- do.call(rbind, out)
  st <- df[!is.na(df$dQ), ]
  cat("\n--- PART A result ---\n")
  cat(sprintf("  replications          : %d\n", length(unique(df$rep))))
  cat(sprintf("  outer IRLS steps      : %d\n", nrow(st)))
  cat(sprintf("  steps decreasing Q    : %d (%.2f%%)\n",
              sum(st$decreased), 100 * mean(st$decreased)))
  cat(sprintf("  largest increase in Q : %.6g\n",
              if (any(st$dQ > 0)) max(st$dQ[st$dQ > 0]) else 0))
  sg <- df$sigma_used[!is.na(df$sigma_used)]
  cat(sprintf("  sigma range across iterations: %.4g to %.4g\n",
              min(sg), max(sg)))
  cat("\n  Q is evaluated at the PILOT scale at every iterate, so the\n")
  cat("  values are comparable. This describes observed behaviour; it is\n")
  cat("  not a guarantee, and the manuscript must not present it as one.\n")
  invisible(df)
}


## =====================================================================
##  PART B -- fixed-scale variant
##  Identical to fit_raen_mm except mad() is called once, before the loop
## =====================================================================
fit_raen_mm_fixed_scale <- function(X, y, gamma = gamma_raen,
                                    al = alpha_raen, cc = c_val, dl = delta,
                                    fe = NULL) {
  if (is.null(fe)) fe <- s1_front_end(X, y, gamma, al, cc, dl)
  beta <- fe$beta0; wp <- fe$wp; lam <- fe$lambda
  sigma_fixed <- fe$sigma_pilot            # <-- computed ONCE

  for (irls in seq_len(MAX_IRLS)) {
    r  <- as.vector(y - X %*% beta)
    u  <- r / sigma_fixed                  # <-- not recomputed
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- X * sqrt(rw); yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break
    b_pre <- beta
    for (cd_iter in seq_len(MAX_CD)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) break
    }
    if (sqrt(sum((beta - b_pre)^2)) < TOL_IRLS) break
  }
  beta
}


## ---------------------------------------------------------------------
##  Reproduce ONE replication's split and preprocessing exactly as run_rep
## ---------------------------------------------------------------------
s1_split_data <- function(rep_id, pp) {
  set.seed(100000L + rep_id)
  n_te <- n_all - floor(train_prop * n_all)
  te   <- sort(sample(seq_len(n_all), n_te))
  tr   <- setdiff(seq_len(n_all), te)
  pre  <- make_preprocessor(X_raw[tr, , drop = FALSE], pp)
  Ztr  <- pre$apply_std(X_raw[tr, , drop = FALSE])
  Zte  <- pre$apply_std(X_raw[te, , drop = FALSE])
  ytr  <- y_full[tr]; yte <- y_full[te]; mu <- mean(ytr)
  list(tr = tr, te = te, pre = pre, Ztr = Ztr, Zte = Zte,
       yc = ytr - mu, yte = yte, mu = mu)
}


## ---------------------------------------------------------------------
##  Paired sensitivity. Both arms share the split, the preprocessing, the
##  pilot, lambda and the adaptive weights -- so the ONLY difference is
##  the scale rule. Scored on both frameworks, exactly as run_rep does.
## ---------------------------------------------------------------------
s1_sensitivity <- function(p = 15000L, reps = 1:100) {
  res <- data.frame()
  cat("PART B: fixed vs adaptive scale, p =", p, ", reps =", length(reps), "\n")
  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next

    ba <- try(fit_raen_mm(d$Ztr, d$yc, gamma = gamma_raen, al = alpha_raen,
                          c_val = c_val, delta = delta), silent = TRUE)
    bb <- try(fit_raen_mm_fixed_scale(d$Ztr, d$yc, fe = fe), silent = TRUE)
    if (inherits(ba, "try-error") || inherits(bb, "try-error")) next
    ba <- as.vector(ba); bb <- as.vector(bb)

    flag    <- train_rowwise_flagger(d$Ztr)
    keep_te <- if (is.null(flag)) seq_along(d$te) else {
      kk <- which(flag$score(d$Zte) <= flag$thr)
      if (length(kk) < 5) seq_along(d$te) else kk
    }
    Zte_imp <- d$pre$impute(d$Zte)

    mA_t <- perf_metrics(d$yte[keep_te],
              d$mu + as.vector(d$Zte[keep_te, , drop = FALSE] %*% ba))
    mB_t <- perf_metrics(d$yte[keep_te],
              d$mu + as.vector(d$Zte[keep_te, , drop = FALSE] %*% bb))
    mA_i <- perf_metrics(d$yte, d$mu + as.vector(Zte_imp %*% ba))
    mB_i <- perf_metrics(d$yte, d$mu + as.vector(Zte_imp %*% bb))

    sa <- which(ba != 0); sb <- which(bb != 0)
    res <- rbind(res, data.frame(
      rep = reps[k], p = p,
      THCM_RMSE_adaptive = mA_t["RMSE"], THCM_RMSE_fixed = mB_t["RMSE"],
      THCM_R2_adaptive   = mA_t["R2"],   THCM_R2_fixed   = mB_t["R2"],
      ICM_RMSE_adaptive  = mA_i["RMSE"], ICM_RMSE_fixed  = mB_i["RMSE"],
      ICM_R2_adaptive    = mA_i["R2"],   ICM_R2_fixed    = mB_i["R2"],
      sel_adaptive = length(sa), sel_fixed = length(sb),
      jaccard = length(intersect(sa, sb)) / max(1, length(union(sa, sb))),
      row.names = NULL))
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }

  rep_pair <- function(a, b, lab) {
    d <- b - a; tt <- t.test(d)
    cat(sprintf(
      "  %-11s adaptive %.4f | fixed %.4f | diff %+.5f (95%% CI %+.5f, %+.5f)  p = %.3g\n",
      lab, mean(a), mean(b), mean(d),
      tt$conf.int[1], tt$conf.int[2], tt$p.value))
    invisible(c(mean(a), mean(b), mean(d), tt$conf.int, tt$p.value))
  }
  cat("\n--- PART B result, p =", p, ", n =", nrow(res), "replications ---\n")
  rep_pair(res$THCM_RMSE_adaptive, res$THCM_RMSE_fixed, "THCM RMSE")
  rep_pair(res$THCM_R2_adaptive,   res$THCM_R2_fixed,   "THCM R2  ")
  rep_pair(res$ICM_RMSE_adaptive,  res$ICM_RMSE_fixed,  "ICM RMSE ")
  rep_pair(res$ICM_R2_adaptive,    res$ICM_R2_fixed,    "ICM R2   ")
  cat(sprintf("  selected: adaptive %.1f | fixed %.1f | support Jaccard %.3f\n",
              mean(res$sel_adaptive), mean(res$sel_fixed), mean(res$jaccard)))
  cat("\n  Report this whatever it shows. A negligible difference is the\n")
  cat("  useful finding; a material one must be stated in the manuscript.\n")
  invisible(res)
}


## =====================================================================
##  PART A2 -- effective sample size along the iterations
##
##  Part A showed the scale collapses by a median factor of 85 WITHIN
##  each fit, in one step. The bisquare gives weight zero to any
##  observation with |r| > c*sigma, so a collapsing sigma means a
##  collapsing rejection threshold. This records how much data is
##  actually carrying the fit at each iteration:
##
##    n_pos  observations with strictly positive weight
##    ess    Kish effective sample size, (sum w)^2 / sum(w^2), which
##           also penalises a few observations carrying most of the
##           weight even when n_pos looks healthy
##
##  n_pos and ess are the numbers that decide whether the reported
##  results rest on a fit supported by the data or by a handful of
##  points. Report them whatever they show.
## =====================================================================
s1_ess_record <- function(X, y, fe = NULL) {
  if (is.null(fe)) fe <- s1_front_end(X, y)
  beta <- fe$beta0; wp <- fe$wp; lam <- fe$lambda
  al <- alpha_raen; cc <- c_val
  n <- nrow(X)

  rec <- data.frame()
  for (irls in seq_len(MAX_IRLS)) {
    r <- as.vector(y - X %*% beta)
    sigma_hat <- mad(r)                      # ADAPTIVE, as implemented
    if (!is.finite(sigma_hat) || sigma_hat < 1e-10) break
    u  <- r / sigma_hat
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break

    sw <- sum(rw); sw2 <- sum(rw^2)
    rec <- rbind(rec, data.frame(
      iter    = irls,
      n_train = n,
      sigma   = sigma_hat,
      cutoff  = cc * sigma_hat,
      n_pos   = sum(rw > 0),
      pct_pos = 100 * sum(rw > 0) / n,
      ess     = if (sw2 > 0) sw^2 / sw2 else 0,
      max_w_share = if (sw > 0) max(rw) / sw else NA_real_,
      n_sel   = sum(beta != 0)))

    Xw <- X * sqrt(rw); yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break
    for (cd_iter in seq_len(MAX_CD)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) break
    }
  }
  rec
}


s1_ess_check <- function(p = 15000L, reps = 1:30) {
  out <- list()
  cat("PART A2: effective sample size, p =", p, ", reps =", length(reps), "\n")
  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    r <- try(s1_ess_record(d$Ztr, d$yc), silent = TRUE)
    if (inherits(r, "try-error") || !nrow(r)) next
    r$rep <- reps[k]; out[[length(out) + 1L]] <- r
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  df <- do.call(rbind, out)

  cat("\n--- PART A2 result, averaged over replications ---\n")
  ag <- aggregate(cbind(sigma, cutoff, n_pos, pct_pos, ess, n_sel) ~ iter,
                  df, mean)
  print(ag, row.names = FALSE, digits = 4)

  last <- df[df$iter == max(df$iter), ]
  cat("\n  at the FINAL iteration:\n")
  cat(sprintf("    n_train                     : %d\n", last$n_train[1]))
  cat(sprintf("    observations with weight > 0: %.1f  (%.1f%%)\n",
              mean(last$n_pos), mean(last$pct_pos)))
  cat(sprintf("    effective sample size       : %.1f\n", mean(last$ess)))
  cat(sprintf("    range of ess across reps    : %.1f to %.1f\n",
              min(last$ess), max(last$ess)))
  cat(sprintf("    largest single weight share : %.3f\n",
              mean(last$max_w_share)))
  cat("\n  Interpretation, decided in advance:\n")
  cat("    ess > 100 of 356  -- collapse is ugly but fits remain supported\n")
  cat("    ess 30-100        -- marginal; must be reported prominently\n")
  cat("    ess < 30          -- the reported fits are not defensible as they\n")
  cat("                         stand and the fixed-scale rerun is required\n")
  invisible(df)
}


## ---------------------------------------------------------------------
##  Write both supplementary tables to CSV, ready to paste
## ---------------------------------------------------------------------
s1_write_supplement <- function(A, ..., dir = ".") {
  B <- list(...)
  st <- A[!is.na(A$dQ), ]
  sa <- data.frame(
    quantity = c("replications", "outer IRLS steps", "steps decreasing Q",
                 "percent decreasing", "largest increase in Q"),
    value = c(length(unique(A$rep)), nrow(st), sum(st$decreased),
              round(100 * mean(st$decreased), 2),
              signif(if (any(st$dQ > 0)) max(st$dQ[st$dQ > 0]) else 0, 6)))
  write.csv(sa, file.path(dir, "S1_partA_descent.csv"), row.names = FALSE)
  if (length(B)) {
    bb <- do.call(rbind, B)
    write.csv(bb, file.path(dir, "S1_partB_sensitivity_raw.csv"),
              row.names = FALSE)
  }
  cat("wrote S1_partA_descent.csv",
      if (length(B)) "and S1_partB_sensitivity_raw.csv" else "", "\n")
}


cat("Loaded S1 verification (driver version).\n")
cat("  1. setwd() to the analysis/ directory containing run_shard.R\n")
cat("  2. s1_load_pipeline()                    # loads data, compiles kernel\n")
cat("  3. A  <- s1_descent_summary(p = 15000L, reps = 1:50)\n")
cat("  4. B1 <- s1_sensitivity(p = 15000L, reps = 1:100)\n")
cat("     B2 <- s1_sensitivity(p = 20000L, reps = 1:100)\n")
cat("  5. s1_write_supplement(A, B1, B2)\n")
