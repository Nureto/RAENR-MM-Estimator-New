# verification/

Materials supporting Supplementary Tables S5 and S6 of the manuscript, which
characterise the behaviour of the implemented (adaptive-scale) algorithm and
compare it against the fixed-scale variant.

| File | Contents |
|---|---|
| `S1_scale_verification.R` | All code. Loads the definitions in `analysis/run_shard.R` without running the 1,000-replication study, then reproduces both analyses. |
| `S1_partA_descent.csv` | Every quantity reported in Supplementary Table S5, with the function that produced it. |
| `S1_partB_sensitivity_raw.csv` | Replication-level output behind Supplementary Table S6: 200 rows, 100 replications at each of p = 15,000 and p = 20,000. |

## Reproducing

```r
setwd("<the analysis/ directory>")
Sys.setenv(GSE14520_PATH = "<path to GSE14520-GPL3921_series_matrix.txt.gz>")
source("../verification/S1_scale_verification.R")
s1_load_pipeline()

A  <- s1_descent_summary(p = 15000L, reps = 1:50)
E  <- s1_ess_check(p = 15000L, reps = 1:30)
B1 <- s1_sensitivity(p = 15000L, reps = 1:100)
B2 <- s1_sensitivity(p = 20000L, reps = 1:100)
s1_write_supplement(A, B1, B2)
```

Both arms of Part B share the training/test split, the preprocessing, the
elastic-net pilot, the penalty parameter and the adaptive penalty weights, so
the only difference between them is where the robust scale is computed.

The gene-expression data are not redistributed. They are available from the
NCBI Gene Expression Omnibus under accession GSE14520, platform GPL3921; the
pipeline expects `GSE14520-GPL3921_series_matrix.txt.gz`, located via the
`GSE14520_PATH` environment variable.
