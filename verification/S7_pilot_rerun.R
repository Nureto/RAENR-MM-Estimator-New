## =====================================================================
##  S7. Pilot rerun on the corrected pipeline.
##
##  run_shard.R has been patched: fit_raen_mm now passes lambda_fixed * n
##  to cd_sweep. Verified end to end by S6_verify_patch.R (831.0 -> 45.3
##  genes at p = 15,000, RMSE difference +0.00009, p = .968, corrected
##  panel 99.9% a subset of the previous one).
##
##  This runs the EXISTING run_rep() over a reduced grid so you can see
##  the direction and size of every change before committing compute to
##  5 dimensionalities x 1,000 replications x 5 estimators.
##
##  IMPORTANT. Start a FRESH R session before running this. The patched
##  run_shard.R must be re-evaluated; a session that already has the old
##  fit_raen_mm in memory will silently reproduce the old results.
##
##  USAGE
##    ## fresh session
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(GSE14520_PATH = "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
##    source("C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R")
##    s1_load_pipeline()
##    source("C:/Users/alabi/Desktop/JBE_Revision_3/S7_pilot_rerun.R")
##
##    s7_check_patch()                                  # do this first
##    PILOT <- s7_pilot(p_grid = c(1000L, 15000L), reps = 1:100)
##    s7_summary(PILOT)
##    s7_write(PILOT)
##
##  Send me the two CSVs and I will tell you which claims in the
##  manuscript survive before you start the full grid.
## =====================================================================

stopifnot(exists("run_rep"), exists("METHODS"), exists("SPARSE"))


## ---------------------------------------------------------------------
##  Refuse to run against an unpatched pipeline. This is the single most
##  likely way to waste a long run.
## ---------------------------------------------------------------------
s7_check_patch <- function(verbose = TRUE) {
  ## deparse() wraps at width.cutoff = 60, so the call of interest can be
  ## split across elements and any regex containing spaces may fail on the
  ## join. Flatten all whitespace first, then match a fixed string.
  src  <- paste(deparse(fit_raen_mm, width.cutoff = 500L), collapse = " ")
  flat <- gsub("[[:space:]]+", "", src)
  ok   <- grepl("lambda_fixed*n,", flat, fixed = TRUE)
  old  <- grepl("cd_sweep(Xw,yw,beta,wp,lambda_fixed,al)", flat, fixed = TRUE)

  if (verbose) {
    ## Show the evidence rather than asking you to trust the matcher.
    lines <- deparse(fit_raen_mm, width.cutoff = 500L)
    hit <- grep("cd_sweep", lines, value = TRUE)
    hit <- hit[!grepl("^\\s*#", hit)]
    cat("cd_sweep call(s) in the loaded fit_raen_mm:\n")
    if (length(hit)) {
      for (h in hit) cat("   ", trimws(h), "\n")
    } else {
      cat("    (none found -- is fit_raen_mm the function you think it is?)\n")
    }
  }

  if (ok) {
    cat("PATCH PRESENT: fit_raen_mm passes lambda_fixed * n. Safe to run.\n")
  } else if (old) {
    cat("*** PATCH ABSENT ***\n")
    cat("This session holds the pre-patch fit_raen_mm. The file on disk is\n")
    cat("patched, but a running session keeps whatever it loaded. Restart R\n")
    cat("and repeat the four setup lines, then check again.\n")
  } else {
    cat("*** CANNOT TELL ***\n")
    cat("Neither the old nor the new call was recognised. Read the line\n")
    cat("printed above and tell me what it says.\n")
  }
  invisible(ok)
}


## ---------------------------------------------------------------------
##  Pilot grid. Uses run_rep() unchanged, so the metrics, the flagging
##  rule, the imputation and the scoring are all exactly as in the
##  reported study.
## ---------------------------------------------------------------------
s7_pilot <- function(p_grid = c(1000L, 15000L), reps = 1:100) {
  if (!s7_check_patch()) {
    stop("refusing to run against an unpatched pipeline")
  }
  met <- list()
  sel <- list()
  tim <- list()
  t0 <- Sys.time()

  for (pp in p_grid) {
    cat("\n=== p =", pp, "===\n")
    for (k in seq_along(reps)) {
      tr <- Sys.time()
      out <- try(run_rep(reps[k], pp), silent = TRUE)
      if (inherits(out, "try-error")) next
      tim[[length(tim) + 1L]] <- data.frame(
        p = pp, rep = reps[k],
        secs = as.numeric(difftime(Sys.time(), tr, units = "secs")))
      met[[length(met) + 1L]] <- out$metrics
      for (m in names(out$sets)) {
        sel[[length(sel) + 1L]] <- list(p = pp, rep = reps[k],
                                        method = m, set = out$sets[[m]])
      }
      if (k %% 10L == 0L) {
        el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
        cat("  ", k, "/", length(reps), " (", el, " min)\n", sep = "")
      }
    }
  }
  list(metrics = do.call(rbind, met), sets = sel,
       timing = do.call(rbind, tim),
       minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}


## ---------------------------------------------------------------------
##  Mean pairwise Jaccard within a method, matching Table S1.
## ---------------------------------------------------------------------
s7_jaccard <- function(sets, pp, method, max_pairs = 2000L) {
  ss <- lapply(Filter(function(z) z$p == pp && z$method == method, sets),
               function(z) z$set)
  if (length(ss) < 2L) return(NA_real_)
  n <- length(ss)
  tot <- 0
  cnt <- 0L
  pairs <- min(max_pairs, n * (n - 1L) / 2L)
  set.seed(1L)
  for (t in seq_len(pairs)) {
    i <- sample.int(n, 1L)
    j <- sample.int(n, 1L)
    if (i == j) next
    a <- ss[[i]]; b <- ss[[j]]
    u <- length(union(a, b))
    if (u == 0L) next
    tot <- tot + length(intersect(a, b)) / u
    cnt <- cnt + 1L
  }
  if (cnt == 0L) NA_real_ else tot / cnt
}


s7_summary <- function(P) {
  M <- P$metrics
  cat("\n--- S7 pilot on the corrected pipeline ---\n")
  cat(sprintf("  rows: %d | elapsed: %.1f min\n", nrow(M), P$minutes))

  for (pp in sort(unique(M$p))) {
    for (fw in c("THCM", "ICM")) {
      cat(sprintf("\n  p = %d, %s\n", pp, fw))
      cat("    method      RMSE      R2       MAE     Selected\n")
      for (m in METHODS) {
        s <- M[M$p == pp & M$Framework == fw & M$Method == m, ]
        if (!nrow(s)) next
        cat(sprintf("    %-11s %.4f   %.4f   %.4f   %8.1f\n", m,
                    mean(s$RMSE, na.rm = TRUE),
                    mean(s$R2, na.rm = TRUE),
                    mean(s$MAE, na.rm = TRUE),
                    mean(s$Selected, na.rm = TRUE)))
      }
    }
  }

  cat("\n  SELECTION STABILITY (mean pairwise Jaccard, cf. Table S1)\n")
  cat("    p        method        Jaccard\n")
  for (pp in sort(unique(M$p))) {
    for (m in SPARSE) {
      j <- s7_jaccard(P$sets, pp, m)
      cat(sprintf("    %-8d %-13s %.4f\n", pp, m, j))
    }
  }

  cat("\n  Compare against the published values at p = 15,000:\n")
  cat("    RAENR-MM selected 835.2 genes, Jaccard 0.125.\n")
  cat("    A corrected panel of tens with a higher Jaccard is the\n")
  cat("    expected outcome. A LOWER Jaccard would be a surprise and\n")
  cat("    worth stopping for.\n")
  invisible(M)
}


s7_write <- function(P, dir = ".") {
  write.csv(P$metrics, file.path(dir, "S7_pilot_metrics.csv"),
            row.names = FALSE)
  rows <- do.call(rbind, lapply(P$sets, function(z) {
    data.frame(p = z$p, rep = z$rep, method = z$method,
               n_sel = length(z$set),
               genes = paste(z$set, collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  write.csv(rows, file.path(dir, "S7_pilot_sets.csv"), row.names = FALSE)
  if (!is.null(P$timing)) {
    write.csv(P$timing, file.path(dir, "S7_pilot_timing.csv"),
              row.names = FALSE)
  }
  cat("wrote S7_pilot_metrics.csv, S7_pilot_sets.csv",
      "and S7_pilot_timing.csv\n")
}


## ---------------------------------------------------------------------
##  One call that does everything: check the patch, run the pilot,
##  print the summary, write both CSVs. Sourcing this file only DEFINES
##  these functions -- it does not start anything -- so this is the
##  entry point.
## ---------------------------------------------------------------------
## ---------------------------------------------------------------------
##  Project the wall-clock cost of the FULL grid from the pilot timings.
##
##  Cost per replication is close to linear in p, because the inner loop
##  is O(np) per sweep. Two pilot dimensionalities therefore give a line
##  through which the three unmeasured ones can be interpolated. This is
##  an estimate, not a promise: it assumes the same machine, the same
##  number of shards, and no contention.
## ---------------------------------------------------------------------
s7_project <- function(P, full_reps = 1000L,
                       full_grid = c(1000L, 5000L, 10000L, 15000L, 20000L),
                       shards = 4L) {
  tt <- P$timing
  if (is.null(tt) || length(unique(tt$p)) < 2L) {
    cat("Need at least two pilot dimensionalities to project.\n")
    return(invisible(NULL))
  }
  agg <- aggregate(secs ~ p, data = tt, FUN = mean)
  fit <- lm(secs ~ p, data = agg)
  pred <- predict(fit, newdata = data.frame(p = full_grid))
  pred[pred < 0] <- min(agg$secs)

  cat("\n--- projected cost of the full rerun ---\n")
  cat("  measured, seconds per replication (all five estimators):\n")
  for (i in seq_len(nrow(agg))) {
    cat(sprintf("    p = %-7d %.1f s\n", agg$p[i], agg$secs[i]))
  }
  cat("  projected:\n")
  for (i in seq_along(full_grid)) {
    tag <- if (full_grid[i] %in% agg$p) " (measured)" else ""
    cat(sprintf("    p = %-7d %.1f s%s\n", full_grid[i], pred[i], tag))
  }
  total_s <- sum(pred) * full_reps
  cat(sprintf("\n  total single-threaded : %.1f hours (%.1f days)\n",
              total_s / 3600, total_s / 86400))
  cat(sprintf("  across %d shards      : %.1f hours (%.1f days)\n",
              shards, total_s / 3600 / shards, total_s / 86400 / shards))
  cat("\n  Assumes the same machine and no contention between shards.\n")
  cat("  Add the rebuild on top: tables, figures, the corrected S2 and\n")
  cat("  descent diagnostics, the rewrite, and your own review.\n")
  invisible(data.frame(p = full_grid, secs_per_rep = pred))
}


s7_run_all <- function(p_grid = c(1000L, 15000L), reps = 1:100,
                       dir = ".") {
  P <- s7_pilot(p_grid = p_grid, reps = reps)
  s7_summary(P)
  s7_project(P)
  s7_write(P, dir = dir)
  cat("\nDone in", round(P$minutes, 1), "minutes.\n")
  cat("Send S7_pilot_metrics.csv and S7_pilot_sets.csv.\n")
  invisible(P)
}


cat("Loaded S7 pilot rerun. Sourcing only defines the functions.\n")
cat("To run it, copy this one line:\n\n")
cat("  PILOT <- s7_run_all()\n\n")
cat("Or, step by step:\n")
cat("  s7_check_patch()\n")
cat("  PILOT <- s7_pilot(p_grid = c(1000L, 15000L), reps = 1:100)\n")
cat("  s7_summary(PILOT)\n")
cat("  s7_write(PILOT)\n")
