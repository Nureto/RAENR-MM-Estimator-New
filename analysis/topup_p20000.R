## =====================================================================
##  Top-up: run ONLY p = 20,000 for a shard whose file stopped at 15,000.
##
##  WHAT HAPPENED
##  Shards 1, 2 and 3 terminated during the p = 20,000 stage. run_shard.R
##  saves the shard file after every dimensionality, so their results for
##  p = 1,000 to 15,000 are complete and valid on disk -- 2,500 rows each,
##  no NULL entries. Only the last dimensionality is missing. Shard 4 ran
##  to completion and has all five.
##
##  combine_shards.R took its P_GRID from the first file it listed, which
##  was the incomplete shard 1, and therefore reported shard 4 -- the one
##  correct file -- as the odd one out.
##
##  This rebuilds only the missing piece and merges it into the existing
##  file, so the hours already spent on p = 1,000 to 15,000 are not
##  repeated.
##
##  TWO WAYS TO RUN IT
##
##  (a) In RStudio, one shard at a time:
##        setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##        TOPUP_SHARD <- 1
##        source("topup_p20000.R")
##      then repeat with 2 and 3. Running in RStudio is fine and lets you
##      see any error rather than losing it with the window.
##
##  (b) From a command prompt, or via run_topup_p20000.bat:
##        Rscript topup_p20000.R 1
##
##  RUN THE SHARDS ONE AT A TIME, not concurrently. Three of the four
##  died at exactly the dimensionality with the largest working set,
##  while the survivor sped up once the others had gone. Do not recreate
##  that contention.
## =====================================================================

ANALYSIS <- "C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis"
S1FILE   <- "C:/Users/alabi/Desktop/JBE_Revision_2/S1_scale_verification.R"
PP       <- 20000L

## Deliberately no NSHARD constant here. run_shard.R defines SHARD and
## NSHARD, and s1_load_pipeline() evaluates it, so any such name in this
## file would be silently overwritten with 1 the moment the pipeline
## loads. The shard count is read from each shard file instead.


topup_p20000 <- function(SHARD) {
  SHARD <- as.integer(SHARD)
  stopifnot(!is.na(SHARD), SHARD >= 1L, SHARD <= 99L)

  setwd(ANALYSIS)
  if (Sys.getenv("GSE14520_PATH") == "")
    Sys.setenv(GSE14520_PATH =
      "C:/Users/alabi/Downloads/GSE14520-GPL3921_series_matrix.txt.gz")
  OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(ANALYSIS, "raenr_out"))

  f <- file.path(OUT_DIR, sprintf("shard_%02d.rds", SHARD))
  if (!file.exists(f)) stop("no existing shard file at ", f)
  S <- readRDS(f)

  cat(sprintf("shard %d | existing dimensionalities: %s\n", SHARD,
              paste(names(S$metrics), collapse = ", ")))
  if (as.character(PP) %in% names(S$metrics)) {
    cat("p = 20000 is already present. Nothing to do for this shard.\n")
    return(invisible(FALSE))
  }

  ## Load the pipeline definitions only. This picks up the PATCHED
  ## fit_raen_mm, which is what the existing four dimensionalities were
  ## also run under. Skip if it is already loaded in this session.
  if (!exists("run_rep", envir = globalenv())) {
    source(S1FILE)
    s1_load_pipeline()
  } else {
    cat("pipeline already loaded in this session\n")
  }

  ## Refuse to proceed against unpatched code -- mixing a corrected
  ## p <= 15000 with a defective p = 20000 would be worse than either.
  src <- paste(deparse(get("fit_raen_mm", envir = globalenv()),
                       width.cutoff = 500L), collapse = " ")
  if (!grepl("lambda_fixed*n,", gsub("[[:space:]]+", "", src), fixed = TRUE))
    stop("run_shard.R is not patched: fit_raen_mm passes lambda_fixed ",
         "unscaled. Refusing to produce results inconsistent with the ",
         "rest of this shard.")
  cat("patch check: PASSED\n")

  ## Take the shard count from the FILE, not from a global. Sourcing
  ## run_shard.R for its definitions also re-runs its argument parsing,
  ## which with no command-line arguments sets SHARD <- 1, NSHARD <- 1
  ## in the global environment ("shard 1/1 | 1000 replications"). Reading
  ## NSHARD from there would build seq(1, 1000, by = 1) and never match.
  NS <- as.integer(S$nshard)
  stopifnot(!is.na(NS), NS >= 1L)

  MY_REPS <- seq(SHARD, 1000L, by = NS)
  if (!identical(as.integer(MY_REPS), as.integer(S$reps)))
    stop(sprintf(paste("replication list does not match the stored shard.",
                       "Computed %d reps (%d to %d, by %d); stored %d reps",
                       "(%d to %d). Aborting."),
                 length(MY_REPS), min(MY_REPS), max(MY_REPS), NS,
                 length(S$reps), min(S$reps), max(S$reps)))
  cat(sprintf("replications: %d (%d to %d, every %d)\n",
              length(MY_REPS), min(MY_REPS), max(MY_REPS), NS))

  t0  <- Sys.time()
  out <- vector("list", length(MY_REPS))
  st  <- vector("list", length(MY_REPS))
  for (k in seq_along(MY_REPS)) {
    r <- run_rep(MY_REPS[k], PP)
    out[[k]] <- r$metrics
    st[[k]]  <- r$sets
    if (k %% 25L == 0L)
      cat(sprintf("  p=%d  %d/%d  (%.1f min elapsed)\n", PP, k,
                  length(MY_REPS),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  S$metrics[[as.character(PP)]] <- do.call(rbind, out)
  S$sets[[as.character(PP)]]    <- st

  ## Shard 1 also carries the collinearity diagnostics.
  if (SHARD == 1L && !is.null(S$diags) &&
      !(as.character(PP) %in% names(S$diags))) {
    cat("computing the collinearity diagnostic for p = 20000\n")
    S$diags[[as.character(PP)]] <- collinearity_diag(PP)
  }

  ## Keep dimensionalities ascending, so combine_shards.R's identical()
  ## check passes however it happens to list the files.
  S$metrics <- S$metrics[order(as.integer(names(S$metrics)))]
  S$sets    <- S$sets[order(as.integer(names(S$sets)))]
  if (!is.null(S$diags) && length(S$diags))
    S$diags <- S$diags[order(as.integer(names(S$diags)))]

  stopifnot(nrow(S$metrics[[as.character(PP)]]) == length(MY_REPS) * 2L * 5L)

  saveRDS(S, f)
  cat(sprintf("\nshard %d topped up in %.1f min. Dimensionalities now: %s\n",
              SHARD, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              paste(names(S$metrics), collapse = ", ")))
  invisible(TRUE)
}


## --- resolve which shard to run --------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.shard <- suppressWarnings(as.integer(.args[1]))
if (is.na(.shard) && exists("TOPUP_SHARD")) {
  .shard <- suppressWarnings(as.integer(TOPUP_SHARD))
}

if (is.na(.shard)) {
  cat("\nNo shard specified. Choose one of:\n\n")
  cat("  In RStudio:\n")
  cat("    TOPUP_SHARD <- 1\n")
  cat("    source(\"topup_p20000.R\")\n\n")
  cat("  From a command prompt:\n")
  cat("    Rscript topup_p20000.R 1\n\n")
  cat("Shards 1, 2 and 3 each need it. Shard 4 is already complete.\n")
  cat("Run them ONE AT A TIME.\n")
} else {
  topup_p20000(.shard)
}
