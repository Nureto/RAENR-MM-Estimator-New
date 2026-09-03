## =====================================================================
##  Selection stability from the CORRECTED authoritative run, all four
##  shards, all 1,000 replications.
##
##  WHY THIS EXISTS
##  extract_stability.R reads shard_01.rds only, so it summarises 250
##  replications, not 1,000. combine_shards.R does use all four shards
##  but omits the 95% selection-frequency threshold and the coefficient
##  of variation of panel size, both of which Supplementary Tables S3 and
##  S4 report. This computes everything Tables S1 to S4b need, from the
##  full run.
##
##    setwd("C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis")
##    Sys.setenv(RAENR_OUT = "C:/Users/alabi/Desktop/RAENR-MM_GitHub_release_v2/analysis/raenr_out")
##    source("extract_stability_v2.R")
##
##  Two to three minutes. No fitting, only summarising stored output.
## =====================================================================

OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
files <- sort(list.files(OUT_DIR, "^shard_\\d+\\.rds$", full.names = TRUE))
stopifnot(length(files) > 0)
cat("reading", length(files), "shard files\n")
sh <- lapply(files, readRDS)

SPARSE <- c("LASSO", "ENET", "MM-RWAL", "RAENR-MM")
NPAIR  <- 2000L
P      <- sort(as.integer(names(sh[[1]]$sets)))

jaccard <- function(A, B) {
  u <- length(union(A, B))
  if (!u) return(NA_real_)
  length(intersect(A, B)) / u
}

set.seed(2026)
out <- list()

for (pk in as.character(P)) {
  ## pool the replication-level selection sets across ALL shards
  reps <- unlist(lapply(sh, function(s) s$sets[[pk]]), recursive = FALSE)
  R <- length(reps)

  for (m in SPARSE) {
    S <- lapply(reps, function(x) x[[m]])
    S <- S[!vapply(S, is.null, logical(1))]
    R2 <- length(S)

    idx <- replicate(NPAIR, sample(R2, 2))
    jac <- apply(idx, 2, function(k) jaccard(S[[k[1]]], S[[k[2]]]))

    tab   <- table(unlist(S))
    freq  <- as.numeric(tab) / R2
    sizes <- lengths(S)
    k_bar <- mean(sizes)
    pp    <- as.integer(pk)

    ## Expected Jaccard for two independent uniform subsets of size k
    ## drawn from p candidates is k / (2p - k). Table S1b reports the
    ## observed index as a multiple of this.
    exp_j <- k_bar / (2 * pp - k_bar)

    out[[length(out) + 1L]] <- data.frame(
      p = pp, Method = m, R = R2,
      Mean_Jaccard   = mean(jac, na.rm = TRUE),
      SE_Jaccard     = sd(jac, na.rm = TRUE) / sqrt(NPAIR),
      Median_Jaccard = median(jac, na.rm = TRUE),
      Expected_Jaccard_random = exp_j,
      Jaccard_multiple_of_random = mean(jac, na.rm = TRUE) / exp_j,
      N_ever_selected = length(freq),
      N_freq_ge_50 = sum(freq >= 0.50),
      N_freq_ge_80 = sum(freq >= 0.80),
      N_freq_ge_95 = sum(freq >= 0.95),
      Mean_set_size = k_bar,
      SD_set_size   = sd(sizes),
      CV_set_size   = sd(sizes) / k_bar,
      row.names = NULL)
  }
  cat("p =", pk, "done (", R, "replications )\n")
}

stability <- do.call(rbind, out)
stability <- stability[order(stability$p,
                             match(stability$Method, SPARSE)), ]

cat("\n=========== MEAN PAIRWISE JACCARD (Table S1) ===========\n")
print(reshape(stability[, c("p", "Method", "Mean_Jaccard")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE, digits = 4)

cat("\n===== AS A MULTIPLE OF RANDOM SELECTION (Table S1b) =====\n")
print(reshape(stability[, c("p", "Method", "Jaccard_multiple_of_random")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE, digits = 4)

cat("\n=========== GENES AT >= 80% (Table S2) ===========\n")
print(reshape(stability[, c("p", "Method", "N_freq_ge_80")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE)

cat("\n=========== GENES AT >= 95% (Table S3) ===========\n")
print(reshape(stability[, c("p", "Method", "N_freq_ge_95")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE)

cat("\n===== PANEL SIZE DISPERSION, SD/mean (Table S4) =====\n")
print(reshape(stability[, c("p", "Method", "CV_set_size")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE, digits = 3)

cat("\n===== DISTINCT GENES EVER SELECTED (Table S4b) =====\n")
print(reshape(stability[, c("p", "Method", "N_ever_selected")],
              idvar = "Method", timevar = "p", direction = "wide"),
      row.names = FALSE)

f <- file.path(OUT_DIR, "GSE14520_v2_stability_full.csv")
write.csv(stability, f, row.names = FALSE)
cat("\nwritten:", f, "\n")
cat("Send me that file and I will rebuild Supplementary Tables S1 to S4b.\n")
