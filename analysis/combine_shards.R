## =====================================================================
##  Combine the shard outputs into the single authoritative run.
##  Run in RStudio once all four windows say SHARD COMPLETE.
## =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## Output directory. Defaults to ./raenr_out; override with the RAENR_OUT
## environment variable. Must match the directory used by run_shard.R.
OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
files   <- list.files(OUT_DIR, "^shard_\\d+\\.rds$", full.names = TRUE)
cat("found", length(files), "shard files\n")
stopifnot(length(files) > 0)

sh <- lapply(files, readRDS)

## ---- integrity checks BEFORE combining ------------------------------------
nsh <- unique(vapply(sh, function(x) x$nshard, integer(1)))
ids <- sort(vapply(sh, function(x) x$shard, integer(1)))
stopifnot(length(nsh) == 1L)
if (!identical(ids, seq_len(nsh)))
  stop("missing shards: expected 1..", nsh, ", have ", paste(ids, collapse = ","))

reps_all <- sort(unlist(lapply(sh, `[[`, "reps")))
if (anyDuplicated(reps_all)) stop("shards overlap - replications counted twice")
cat("replications covered:", length(reps_all),
    "| contiguous 1..n:", identical(reps_all, seq_along(reps_all)), "\n")

P_GRID <- as.integer(names(sh[[1]]$metrics))
for (s in sh)
  if (!identical(as.integer(names(s$metrics)), P_GRID))
    stop("shard ", s$shard, " did not finish every dimensionality")

## ---- combine ---------------------------------------------------------------
results <- do.call(rbind, lapply(sh, function(s) do.call(rbind, s$metrics)))
results <- results[order(results$p, results$Framework, results$Method, results$Rep), ]
rownames(results) <- NULL

diagnostics <- do.call(rbind, sh[[which(ids == 1L)]]$diags)

cat("\nrows:", nrow(results), "| expected:",
    length(reps_all) * length(P_GRID) * 2 * 5, "\n")
stopifnot(nrow(results) == length(reps_all) * length(P_GRID) * 2 * 5)
stopifnot(!anyNA(results$RMSE))

## ---- summaries --------------------------------------------------------------
summ <- results |>
  group_by(Framework, p, Method) |>
  summarise(R = n(),
            Mean_RMSE = mean(RMSE), SD_RMSE = sd(RMSE),
            SE_RMSE = sd(RMSE)/sqrt(n()), Median_RMSE = median(RMSE),
            Mean_R2 = mean(R2, na.rm = TRUE), SE_R2 = sd(R2, na.rm=TRUE)/sqrt(n()),
            Mean_MAE = mean(MAE), SE_MAE = sd(MAE)/sqrt(n()),
            Mean_Sel = mean(Selected), SD_Sel = sd(Selected),
            Mean_ntest = mean(n_test), .groups = "drop") |>
  arrange(Framework, p, Mean_RMSE) |>
  group_by(Framework, p) |> mutate(Rank = row_number()) |> ungroup()

paired <- results |>
  select(Rep, p, Framework, Method, RMSE) |>
  pivot_wider(names_from = Method, values_from = RMSE) |>
  group_by(Framework, p) |>
  summarise(R = n(),
            d_RWAL   = mean(`MM-RWAL` - `RAENR-MM`),
            se_RWAL  = sd(`MM-RWAL` - `RAENR-MM`)/sqrt(n()),
            pct_RWAL = 100*mean((`MM-RWAL` - `RAENR-MM`)/`MM-RWAL`),
            p_RWAL   = 2*pt(-abs(d_RWAL/se_RWAL), df = n()-1),
            d_ENET   = mean(ENET - `RAENR-MM`),
            se_ENET  = sd(ENET - `RAENR-MM`)/sqrt(n()),
            p_ENET   = 2*pt(-abs(d_ENET/se_ENET), df = n()-1), .groups = "drop")

## ---- selection stability ----------------------------------------------------
SPARSE <- c("LASSO","ENET","MM-RWAL","RAENR-MM")
jaccard <- function(A,B){u <- length(union(A,B)); if(!u) NA_real_ else length(intersect(A,B))/u}

stability <- do.call(rbind, lapply(as.character(P_GRID), function(pk) {
  S_all <- unlist(lapply(sh, function(s) s$sets[[pk]]), recursive = FALSE)
  do.call(rbind, lapply(SPARSE, function(m) {
    S <- lapply(S_all, function(x) x[[m]]); R <- length(S)
    np <- min(500L, choose(R, 2)); idx <- replicate(np, sample(R, 2))
    jac <- apply(idx, 2, function(k) jaccard(S[[k[1]]], S[[k[2]]]))
    freq <- as.numeric(table(unlist(S)))/R
    data.frame(p = as.integer(pk), Method = m, R = R,
               Mean_Jaccard = mean(jac, na.rm=TRUE),
               SE_Jaccard = sd(jac, na.rm=TRUE)/sqrt(np),
               N_ever_selected = length(freq),
               N_freq_ge_50 = sum(freq >= .50), N_freq_ge_80 = sum(freq >= .80),
               Mean_set_size = mean(lengths(S)), SD_set_size = sd(lengths(S)),
               row.names = NULL)
  }))
}))

## ---- save -------------------------------------------------------------------
saveRDS(list(results=results, summ=summ, paired=paired,
             diagnostics=diagnostics, stability=stability,
             n_shards=nsh, sessionInfo=sessionInfo()),
        file.path(OUT_DIR, "GSE14520_AUTHORITATIVE_RUN_v2.rds"))
for (nm in c("results","summ","paired","diagnostics","stability"))
  write.csv(get(nm), file.path(OUT_DIR, paste0("GSE14520_v2_", nm, ".csv")),
            row.names = FALSE)

cat("\n=========== RAENR-MM SUMMARY ===========\n")
print(as.data.frame(summ |> filter(Method == "RAENR-MM") |>
        select(Framework, p, R, Mean_RMSE, SE_RMSE, Rank)))
cat("\n=========== vs MM-RWAL ===========\n")
print(as.data.frame(paired |> select(Framework, p, pct_RWAL, p_RWAL)))
cat("\nwritten to", OUT_DIR, "\n")
cat("Send GSE14520_AUTHORITATIVE_RUN_v2.rds (or the five CSVs).\n")
