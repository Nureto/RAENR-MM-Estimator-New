## =====================================================================
##  Why did combine_shards.R reject shard 4?
##
##  The check that failed is
##      identical(as.integer(names(s$metrics)), P_GRID)
##  where P_GRID is taken from the FIRST shard file. It fails if a shard's
##  metrics list has different names, a different ORDER, a different
##  length, or a NULL entry -- not only if a dimensionality is missing.
##
##  Shard 4's console log shows all five dimensionalities completing, so a
##  genuinely missing p is unlikely. This prints enough to tell which of
##  the other possibilities it is.
##
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    source("diagnose_shards.R")
## =====================================================================

OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
files <- sort(list.files(OUT_DIR, "^shard_\\d+\\.rds$", full.names = TRUE))
cat("shard files:", length(files), "in", OUT_DIR, "\n\n")

sh <- lapply(files, readRDS)

for (s in sh) {
  cat(sprintf("--- shard %d of %d ---\n", s$shard, s$nshard))
  cat("  reps            :", length(s$reps),
      sprintf("(%d to %d)\n", min(s$reps), max(s$reps)))
  cat("  metrics length  :", length(s$metrics), "\n")
  cat("  metrics names   :", paste(names(s$metrics), collapse = ", "), "\n")
  cat("  as.integer()    :", paste(as.integer(names(s$metrics)),
                                   collapse = ", "), "\n")
  nulls <- vapply(s$metrics, is.null, logical(1))
  cat("  NULL entries    :",
      if (any(nulls)) paste(names(s$metrics)[nulls], collapse = ", ")
      else "none", "\n")
  rows <- vapply(s$metrics, function(x) if (is.null(x)) NA_integer_
                 else nrow(x), integer(1))
  cat("  rows per p      :", paste(rows, collapse = ", "),
      " (expect", length(s$reps) * 2 * 5, "each)\n")
  cat("  sets length     :", length(s$sets),
      "| names:", paste(names(s$sets), collapse = ", "), "\n\n")
}

## --- the comparison combine_shards.R actually makes -------------------
P1 <- as.integer(names(sh[[1]]$metrics))
cat("P_GRID taken from the first file:", paste(P1, collapse = ", "), "\n\n")
for (s in sh) {
  Pi <- as.integer(names(s$metrics))
  ok <- identical(Pi, P1)
  cat(sprintf("  shard %d identical to P_GRID: %s\n", s$shard, ok))
  if (!ok) {
    cat("      same set of values, different order? ",
        setequal(Pi, P1) && !identical(Pi, P1), "\n")
    cat("      values in shard not in P_GRID     : ",
        paste(setdiff(Pi, P1), collapse = ", "), "\n")
    cat("      values in P_GRID not in shard     : ",
        paste(setdiff(P1, Pi), collapse = ", "), "\n")
    cat("      any NA from as.integer()          : ", anyNA(Pi), "\n")
  }
}

cat("\nIf the only difference is ORDER, the run is fine and the check is\n")
cat("too strict; combine_shards_SAFE.R reorders and proceeds.\n")
cat("If a dimensionality is genuinely absent, that p must be rerun for\n")
cat("that shard before combining. Send me this output either way.\n")
