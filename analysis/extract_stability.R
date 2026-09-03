## =====================================================================
##  SELECTION STABILITY from the authoritative run
##
##  Answers editorial comment S3(e): the manuscript states that selection
##  frequencies and mean pairwise Jaccard indices are reported in the
##  Supplementary Material, while the Supplementary Material says they
##  "will be added when that file is available". This computes them.
##
##  With the 80-100 gene sparsity claim withdrawn, reproducibility of the
##  selected panel is now the paper's principal practical claim, so these
##  numbers matter more than they did before.
##
##  source("extract_stability.R")   -- about 2 minutes
## =====================================================================

## Output directory. Defaults to ./raenr_out; override with the RAENR_OUT
## environment variable. Must match the directory used by run_shard.R.
OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
f <- file.path(OUT_DIR, "shard_01.rds")
stopifnot(file.exists(f))
s <- readRDS(f)

SPARSE <- c("LASSO", "ENET", "MM-RWAL", "RAENR-MM")
P      <- as.integer(names(s$sets))
NPAIR  <- 2000L          # random replication pairs per method x p

jaccard <- function(A, B) {
  u <- length(union(A, B))
  if (!u) return(NA_real_)
  length(intersect(A, B)) / u
}

set.seed(2026)
out <- list()

for (pk in names(s$sets)) {
  reps <- s$sets[[pk]]                 # list over replications
  R <- length(reps)
  for (m in SPARSE) {
    S <- lapply(reps, function(x) x[[m]])
    S <- S[!vapply(S, is.null, logical(1))]
    R2 <- length(S)

    ## mean pairwise Jaccard over random replication pairs
    idx <- replicate(NPAIR, sample(R2, 2))
    jac <- apply(idx, 2, function(k) jaccard(S[[k[1]]], S[[k[2]]]))

    ## per-gene selection frequency
    tab  <- table(unlist(S))
    freq <- as.numeric(tab) / R2

    sizes <- lengths(S)
    out[[length(out) + 1]] <- data.frame(
      p = as.integer(pk), Method = m, R = R2,
      Mean_Jaccard = mean(jac, na.rm = TRUE),
      SE_Jaccard   = sd(jac, na.rm = TRUE) / sqrt(NPAIR),
      Median_Jaccard = median(jac, na.rm = TRUE),
      N_ever_selected = length(freq),
      N_freq_ge_50 = sum(freq >= 0.50),
      N_freq_ge_80 = sum(freq >= 0.80),
      N_freq_ge_95 = sum(freq >= 0.95),
      Mean_set_size = mean(sizes),
      SD_set_size   = sd(sizes),
      CV_set_size   = sd(sizes) / mean(sizes),
      row.names = NULL)
  }
  cat("p =", pk, "done\n")
}

stability <- do.call(rbind, out)
stability <- stability[order(stability$p,
                             match(stability$Method, SPARSE)), ]

cat("\n=========== MEAN PAIRWISE JACCARD ===========\n")
print(reshape(stability[, c("p","Method","Mean_Jaccard")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE, digits = 4)

cat("\n=========== GENES SELECTED IN >= 80% OF REPLICATIONS ===========\n")
print(reshape(stability[, c("p","Method","N_freq_ge_80")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE)

cat("\n=========== RELATIVE DISPERSION OF PANEL SIZE (SD/mean) ===========\n")
cat("lower = more consistent panel size\n")
print(reshape(stability[, c("p","Method","CV_set_size")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE, digits = 3)

cat("\n=========== FULL TABLE ===========\n")
print(stability, row.names = FALSE, digits = 4)

write.csv(stability, file.path(OUT_DIR, "shard1_stability.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT_DIR, "shard1_stability.csv"), "\n")

