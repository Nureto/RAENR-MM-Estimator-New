## =====================================================================
##  Analyse shard 1 on its own (250 replications, all five dimensionalities).
##
##  combine_shards.R deliberately refuses to run without all four shards.
##  This script is the single-shard equivalent: it labels the result R = 250
##  everywhere so it can never be mistaken for the 1,000-replication run.
##
##  250 replications is ample to settle the sign and significance of a ~5%
##  paired difference. It is a legitimate analysis, just with wider Monte
##  Carlo error (2x that of R = 1000).
##
##  source("inspect_shard1.R")
## =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## Output directory. Defaults to ./raenr_out; override with the RAENR_OUT
## environment variable. Must match the directory used by run_shard.R.
OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
f <- file.path(OUT_DIR, "shard_01.rds")
stopifnot(file.exists(f))
s <- readRDS(f)

cat("=========== COMPLETENESS ===========\n")
cat("shard", s$shard, "of", s$nshard, "| replications:", length(s$reps), "\n")
cat("dimensionalities finished:", paste(names(s$metrics), collapse = ", "), "\n")

results <- do.call(rbind, s$metrics)
rownames(results) <- NULL
P <- sort(unique(results$p))
R <- length(unique(results$Rep))

cat("rows:", nrow(results), "| expected:", R * length(P) * 2 * 5, "\n")
cat("NAs in RMSE:", sum(is.na(results$RMSE)), "\n")
ok <- nrow(results) == R * length(P) * 2 * 5 && !anyNA(results$RMSE)
cat("complete and clean:", ok, "\n")
if (!ok) cat("  -> shard is PARTIAL; treat what follows as provisional\n")

## ---------------------------------------------------------------- summaries
summ <- results |>
  group_by(Framework, p, Method) |>
  summarise(R = n(),
            Mean_RMSE = mean(RMSE), SD_RMSE = sd(RMSE),
            SE_RMSE = sd(RMSE)/sqrt(n()), Median_RMSE = median(RMSE),
            Mean_R2 = mean(R2, na.rm = TRUE),
            Mean_MAE = mean(MAE),
            Mean_Sel = mean(Selected), SD_Sel = sd(Selected),
            Mean_ntest = mean(n_test), .groups = "drop") |>
  arrange(Framework, p, Mean_RMSE) |>
  group_by(Framework, p) |> mutate(Rank = row_number()) |> ungroup()

ord <- c("Ridge","LASSO","ENET","MM-RWAL","RAENR-MM")
for (fw in c("THCM","ICM")) {
  cat("\n=========== ", fw, ": mean RMSE ===========\n", sep = "")
  print(round(as.data.frame(
    summ |> filter(Framework == fw) |>
      select(Method, p, Mean_RMSE) |>
      pivot_wider(names_from = p, values_from = Mean_RMSE) |>
      slice(match(ord, Method)) |> tibble::column_to_rownames("Method")), 4))
  cat("\nrank (1 = lowest RMSE):\n")
  print(as.data.frame(
    summ |> filter(Framework == fw) |>
      select(Method, p, Rank) |>
      pivot_wider(names_from = p, values_from = Rank) |>
      slice(match(ord, Method)) |> tibble::column_to_rownames("Method")))
}

## ------------------------------------------------- THE HEADLINE COMPARISON
paired <- results |>
  select(Rep, p, Framework, Method, RMSE) |>
  pivot_wider(names_from = Method, values_from = RMSE) |>
  group_by(Framework, p) |>
  summarise(R = n(),
            pct_vs_RWAL = 100*mean((`MM-RWAL` - `RAENR-MM`)/`MM-RWAL`),
            d = mean(`MM-RWAL` - `RAENR-MM`),
            se = sd(`MM-RWAL` - `RAENR-MM`)/sqrt(n()),
            p_value = 2*pt(-abs(d/se), df = n()-1),
            pct_vs_ENET = 100*mean((ENET - `RAENR-MM`)/ENET),
            .groups = "drop")

cat("\n=========== RAENR-MM vs MM-RWAL (paired) ===========\n")
cat("positive pct = RAENR-MM BETTER. Manuscript claims +4.78% to +13.81% (THCM).\n\n")
print(as.data.frame(paired |>
        mutate(across(c(pct_vs_RWAL, pct_vs_ENET), ~round(., 2)),
               d = round(d, 5), se = round(se, 5),
               p_value = signif(p_value, 3))))

cat("\n=========== selection counts (THCM) ===========\n")
print(round(as.data.frame(
  summ |> filter(Framework == "THCM") |>
    select(Method, p, Mean_Sel) |>
    pivot_wider(names_from = p, values_from = Mean_Sel) |>
    slice(match(ord, Method)) |> tibble::column_to_rownames("Method")), 1))

cat("\n=========== collinearity (S6) ===========\n")
print(do.call(rbind, s$diags), row.names = FALSE)

saveRDS(list(results = results, summ = summ, paired = paired,
             diagnostics = do.call(rbind, s$diags), R = R, shards_used = 1L),
        file.path(OUT_DIR, "GSE14520_SHARD1_ONLY.rds"))
write.csv(summ,   file.path(OUT_DIR, "shard1_summ.csv"),   row.names = FALSE)
write.csv(paired, file.path(OUT_DIR, "shard1_paired.csv"), row.names = FALSE)
cat("\nwritten: GSE14520_SHARD1_ONLY.rds, shard1_summ.csv, shard1_paired.csv\n")
cat("NOTE: R =", R, "replications, not 1,000. Label it as such everywhere.\n")
