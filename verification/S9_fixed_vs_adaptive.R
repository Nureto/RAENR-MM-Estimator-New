## =====================================================================
##  S9. Replacement for Supplementary Table S6.
##
##  WHY IT CANNOT BE PATCHED
##  Table S6 compares the fixed-scale and adaptive-scale variants on
##  prediction at p = 15,000 and p = 20,000. BOTH arms were run under the
##  defective penalty, so every number in it describes fits that are no
##  longer the estimator of record. Unlike Tables A2 to A8, there is no
##  comparator column that survives: the whole table must be recomputed.
##
##  This mirrors the original design exactly -- same two dimensionalities,
##  same R = 100, same paired splits, same scoring under both frameworks --
##  with the corrected penalty in both arms.
##
##  The Editor also asked (item S2) that the conclusion drawn from S6 be
##  restricted to the settings actually assessed. That restriction still
##  applies and is written into the note.
##
##  USAGE, fresh session
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(GSE14520_PATH = "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
##    source("C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R")
##    s1_load_pipeline()
##    source("C:/Users/alabi/Desktop/JBE_Revision_3/S9_fixed_vs_adaptive.R")
##
##    S6 <- s9_run(reps = 1:100)
##
##  Two fits per replication at two dimensionalities: comparable to the
##  original S6 run.
## =====================================================================

stopifnot(exists("s1_front_end"), exists("s1_split_data"),
          exists("cd_sweep"), exists("train_rowwise_flagger"))


s9_fit <- function(X, y, fe, adaptive, max_irls = MAX_IRLS,
                   max_cd = MAX_CD) {
  beta <- fe$beta0
  wp   <- fe$wp
  lam  <- fe$lambda * nrow(X)        # corrected penalty, both arms
  sig0 <- fe$sigma_pilot
  al   <- alpha_raen
  cc   <- c_val

  for (irls in seq_len(max_irls)) {
    r <- as.vector(y - X %*% beta)
    s <- if (adaptive) mad(r) else sig0
    if (!is.finite(s) || s < 1e-10) break
    u  <- r / s
    rw <- ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
    if (sum(rw) < 1e-6) break
    Xw <- X * sqrt(rw); yw <- y * sqrt(rw)
    if (!is.finite(sd(yw)) || sd(yw) < 1e-8) break
    b_pre <- beta
    for (k in seq_len(max_cd)) {
      b_old <- beta
      beta  <- cd_sweep(Xw, yw, beta, wp, lam, al)
      if (sqrt(sum((beta - b_old)^2)) < TOL_CD) break
    }
    if (sqrt(sum((beta - b_pre)^2)) < TOL_IRLS) break
  }
  beta
}


s9_one_p <- function(p, reps) {
  res <- data.frame()
  cat("S9: fixed vs adaptive at the corrected penalty, p =", p,
      "| reps", length(reps), "\n")
  for (k in seq_along(reps)) {
    d <- try(s1_split_data(reps[k], p), silent = TRUE)
    if (inherits(d, "try-error")) next
    fe <- try(s1_front_end(d$Ztr, d$yc), silent = TRUE)
    if (inherits(fe, "try-error")) next

    ba <- try(s9_fit(d$Ztr, d$yc, fe, TRUE), silent = TRUE)
    bf <- try(s9_fit(d$Ztr, d$yc, fe, FALSE), silent = TRUE)
    if (inherits(ba, "try-error") || inherits(bf, "try-error")) next

    ## THCM: drop held-out rows flagged by the training-fitted rule.
    flag <- train_rowwise_flagger(d$Ztr)
    keep <- if (is.null(flag)) seq_along(d$te) else {
      kk <- which(flag$score(d$Zte) <= flag$thr)
      if (length(kk) < 5) seq_along(d$te) else kk
    }
    ## ICM: no rows dropped; cells imputed by the stated rule.
    Zte_imp <- d$pre$impute(d$Zte)

    score <- function(b, Z, idx, yt) {
      pr <- d$mu + as.vector(Z[idx, , drop = FALSE] %*% b)
      ss_r <- sum((yt[idx] - pr)^2)
      ss_t <- sum((yt[idx] - mean(yt[idx]))^2)
      c(rmse = sqrt(mean((yt[idx] - pr)^2)),
        r2 = if (ss_t < 1e-12) NA_real_ else 1 - ss_r / ss_t)
    }
    all_te <- seq_along(d$te)
    ta <- score(ba, d$Zte, keep, d$yte); tf <- score(bf, d$Zte, keep, d$yte)
    ia <- score(ba, Zte_imp, all_te, d$yte)
    iff <- score(bf, Zte_imp, all_te, d$yte)

    res <- rbind(res, data.frame(
      rep = reps[k], p = p,
      thcm_rmse_adaptive = ta["rmse"], thcm_rmse_fixed = tf["rmse"],
      thcm_r2_adaptive = ta["r2"],     thcm_r2_fixed = tf["r2"],
      icm_rmse_adaptive = ia["rmse"],  icm_rmse_fixed = iff["rmse"],
      icm_r2_adaptive = ia["r2"],      icm_r2_fixed = iff["r2"],
      sel_adaptive = sum(ba != 0), sel_fixed = sum(bf != 0),
      jaccard = length(intersect(which(ba != 0), which(bf != 0))) /
        max(1, length(union(which(ba != 0), which(bf != 0)))),
      row.names = NULL))
    if (k %% 10L == 0L) cat("  ", k, "/", length(reps), "\n")
  }
  res
}


s9_report <- function(res) {
  cat("\n--- Table S6 replacement, p =", res$p[1], "---\n")
  cat("  Differences are FIXED minus ADAPTIVE, matching the original.\n")
  for (fw in c("thcm", "icm")) {
    for (met in c("rmse", "r2")) {
      a <- res[[paste0(fw, "_", met, "_adaptive")]]
      f <- res[[paste0(fw, "_", met, "_fixed")]]
      ok <- is.finite(a) & is.finite(f)
      d <- f[ok] - a[ok]
      tt <- t.test(d)
      pv <- if (tt$p.value < .001) "< .001" else sprintf("%.3f", tt$p.value)
      ## NB: R does not concatenate adjacent string literals, so the
      ## format must be built with paste0 rather than split across lines.
      fmt <- paste0("  %-4s %-4s adaptive %.4f | fixed %.4f | ",
                    "diff %+.5f (95%% CI %+.5f, %+.5f) p %s\n")
      cat(sprintf(fmt,
                  toupper(fw), toupper(met), mean(a[ok]), mean(f[ok]),
                  mean(d), tt$conf.int[1], tt$conf.int[2], pv))
    }
  }
  fmt2 <- paste0("  genes selected: adaptive %.1f | fixed %.1f | ",
                 "support Jaccard between variants %.3f\n")
  cat(sprintf(fmt2,
              mean(res$sel_adaptive), mean(res$sel_fixed),
              mean(res$jaccard)))
}


s9_run <- function(reps = 1:100, ps = c(15000L, 20000L), dir = ".") {
  out <- list()
  for (p in ps) {
    r <- s9_one_p(p, reps)
    s9_report(r)
    out[[as.character(p)]] <- r
  }
  all <- do.call(rbind, out)
  write.csv(all, file.path(dir, "S9_fixed_vs_adaptive.csv"),
            row.names = FALSE)
  cat("\nwrote S9_fixed_vs_adaptive.csv\n")
  invisible(out)
}


cat("Loaded S9 fixed-vs-adaptive comparison (Table S6 replacement).\n")
cat("  S6 <- s9_run(reps = 1:100)\n")
