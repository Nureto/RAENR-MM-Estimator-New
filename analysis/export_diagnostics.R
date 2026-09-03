## =====================================================================
##  Export the collinearity / contamination diagnostics from the
##  authoritative run so that Table A1 can be rebuilt.
##
##  This is the ONE table in the manuscript that has not been regenerated.
##  It still shows the artefact values (condition number 49060.170,
##  "max condition index" 221.500, identical at every p) that comment S6
##  objected to, because those numbers live in shard_01.rds and have not
##  been exported. Nothing was written in their place: an unverified
##  number in a table is worse than a stale one that is known to be stale.
##
##  source("export_diagnostics.R")   -- a few seconds
## =====================================================================

## Output directory. Defaults to ./raenr_out; override with the RAENR_OUT
## environment variable. Must match the directory used by run_shard.R.
OUT_DIR <- Sys.getenv("RAENR_OUT", file.path(getwd(), "raenr_out"))
f <- file.path(OUT_DIR, "shard_01.rds")
stopifnot(file.exists(f))
s <- readRDS(f)

diags <- do.call(rbind, s$diags)
print(diags, row.names = FALSE)

cat("\ncolumns available:\n"); print(names(diags))

write.csv(diags, file.path(OUT_DIR, "shard1_diagnostics.csv"), row.names = FALSE)
cat("\nwritten:", file.path(OUT_DIR, "shard1_diagnostics.csv"), "\n")

