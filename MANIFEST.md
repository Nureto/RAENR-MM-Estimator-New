# Manifest

Version 3.0.0, 2026-09-02. SHA-256 checksums for every file.

**What changed since v2.0.0.** The coordinate-descent routine was receiving the cross-validated penalty parameter without the factor n that the loss normalisation requires, so the proposed estimator was penalised about n = 356 times too weakly. `analysis/run_shard.R` is corrected at one line. Every result in `results/` is regenerated from a fresh run of 1,000 replications. The comparators are unaffected and reproduce their previous values exactly. The pre-correction pipeline is retained under `superseded/` so the change can be inspected by diff.

| File | Bytes | SHA-256 (first 16) | Description |
|---|---|---|---|
| `README.md` | 14,186 | `50b114d81d8589fc` | This document. |
| `analysis/cd_sweep.cpp` | 1,519 | `14ee0a5e6d1e64b6` | Compiled coordinate-descent kernel. |
| `analysis/combine_shards.R` | 5,271 | `9ff58a42065a8c28` | Merges shard outputs; refuses an incomplete set. |
| `analysis/diagnose_shards.R` | 3,016 | `9e5d9573bae5ea09` | Checks shard files for completeness. |
| `analysis/export_diagnostics.R` | 1,327 | `7a76e16b03879c02` | Collinearity diagnostics (Table A1). |
| `analysis/extract_stability.R` | 3,572 | `3f12c799639ac7e1` | Selection stability, shard 1 only (superseded by extract_stability_v2.R). |
| `analysis/extract_stability_v2.R` | 4,652 | `84eea2caa5a718df` | Selection stability over all four shards, all thresholds (Tables S1-S4b). |
| `analysis/inspect_shard1.R` | 4,763 | `ccaf31a3c0546db1` | Summaries and paired comparisons from one run file. |
| `analysis/run_all_shards.bat` | 1,568 | `b4f5299fcfddffa3` | Original 4-way launcher. |
| `analysis/run_all_shards_CORRECTED.bat` | 3,034 | `3e3cec42c77ae799` | Launcher that sets GSE14520_PATH and refuses to start against unpatched code. |
| `analysis/run_shard.R` | 28,680 | `a79a5edb06d027ce` | Main analysis pipeline, CORRECTED. Line 459 passes lambda_fixed * n to cd_sweep. |
| `analysis/run_topup_p20000.bat` | 1,983 | `94f4a08472ceb4b1` | Sequential launcher for the above. |
| `analysis/topup_p20000.R` | 6,833 | `33a1c790820c4acc` | Rebuilds a single dimensionality for a shard. |
| `results/GSE14520_v2_diagnostics.csv` | 380 | `76cf7976a6e9b783` | Collinearity diagnostics. Table A1. |
| `results/GSE14520_v2_paired.csv` | 1,548 | `edcd8ad9d74328ad` | Paired RAENR-MM comparisons. Tables A6, A6b. |
| `results/GSE14520_v2_results.csv` | 4,296,959 | `08955f3a8958bd6b` | Replication-level output, 50,000 rows. |
| `results/GSE14520_v2_stability.csv` | 1,957 | `1d4f1de11d8a1bf0` | Selection stability, 500 pairs. |
| `results/GSE14520_v2_stability_full.csv` | 3,521 | `beaf06a4a36302a9` | Selection stability, 2,000 pairs, all thresholds. Tables S1-S4b. |
| `results/GSE14520_v2_summ.csv` | 10,418 | `0ecd820d93195bef` | Per framework/dimensionality/estimator summaries, R = 1,000. Tables A2-A5, A8. |
| `superseded/GSE14520_Real_Analysis_1.Rmd` | 40,367 | `22a91504d75cd972` | Original pipeline with the four leakage points. NOT used for any reported result. |
| `superseded/S2_cap_sensitivity_raw.csv` | 4,522 | `0198b5c9867685d2` | Cap sensitivity of the DEFECTIVE fit. Retained for the same reason. |
| `superseded/S2_convergence_raw.csv` | 5,948 | `7c86b55031a2a9d2` | Convergence audit of the DEFECTIVE fit. Retained to document what was found, not for citation. |
| `superseded/run_shard_PRE_PATCH_2026-09-01.R` | 28,120 | `50ca728885a916b7` | The pipeline as it stood before the penalty correction. Diff against analysis/run_shard.R to see the defect and the fix; they differ in one functional line. |
| `verification/S1_scale_verification.R` | 19,017 | `55902d9be72e8de2` | Shared harness: pilot, splits, fixed-scale variant. |
| `verification/S5_penalty_scale.R` | 14,765 | `fd75070bd1fa7800` | Establishes the penalty-scaling defect: same fits scored two ways. |
| `verification/S5_penalty_scale_raw.csv` | 4,296 | `653347e3e56028bd` |  |
| `verification/S6_patch_verification_raw.csv` | 1,208 | `927465598b5a3e9a` |  |
| `verification/S6_verify_patch.R` | 9,131 | `db8f7e7823e5cf6f` | Verifies the one-line correction against the production routine. |
| `verification/S7_pilot_metrics.csv` | 169,712 | `52552adaa20515f4` |  |
| `verification/S7_pilot_rerun.R` | 10,839 | `cf5ff653f23dc71f` | Pilot rerun at two dimensionalities. |
| `verification/S8_adaptive_audit.csv` | 13,020 | `ff5f49ee31edd7c9` |  |
| `verification/S8_cap_sensitivity.csv` | 3,534 | `ede464b6249dcda7` |  |
| `verification/S8_corrected_diagnostics.R` | 12,970 | `e824df61a9e0eb14` | Convergence, descent and cap sensitivity. Tables S5, S7, S8. |
| `verification/S8_fixed_audit.csv` | 11,757 | `ad67c1e09e5f667e` |  |
| `verification/S9_fixed_vs_adaptive.R` | 6,186 | `c40d0485ebac499a` | Fixed against adaptive scale. Table S6. |
| `verification/S9_fixed_vs_adaptive.csv` | 34,511 | `b552f8647ee5dbcf` |  |

## Full checksums

```
50b114d81d8589fc388e771c2117603ce57373ec9ad7e4776b873356ec9a36d1  README.md
14ee0a5e6d1e64b62c86cda86b41e5d74e5d2bdfaf83bf33dfcd302649c41318  analysis/cd_sweep.cpp
9ff58a42065a8c284c6f0f15dd0c8773ea88ee0c3fe70bfe7d0b8ce029baaba8  analysis/combine_shards.R
9e5d9573bae5ea09e6e79817055737d7a51ef36ca5b1c1f00f7e2b80091d9e16  analysis/diagnose_shards.R
7a76e16b03879c024a52256922a0903cd1a33d061f7c793d2158086272321465  analysis/export_diagnostics.R
3f12c799639ac7e15b4e2cc924ee1bdf062f61d385510d32833461d3eaf0bd4e  analysis/extract_stability.R
84eea2caa5a718df23b51f2d714bda3da8ed7f26ab96740d56e785a2a5a3c931  analysis/extract_stability_v2.R
ccaf31a3c0546db12a122619c999d6fd2a68c4cea27fab70692f59475642cc31  analysis/inspect_shard1.R
b4f5299fcfddffa38489b78b4fa4ba28d31b6bd561d4008bf137a8d703233e6c  analysis/run_all_shards.bat
3e3cec42c77ae799e443922e72f99539fb5e70331ef164c4f9bb0540779a62a0  analysis/run_all_shards_CORRECTED.bat
a79a5edb06d027ce8e3447819e371ae36b25eb86ef3f414456c322d1d5734532  analysis/run_shard.R
94f4a08472ceb4b16cdccdb88b8d557b09f16617697d5e5d7d52d515cab5bd8e  analysis/run_topup_p20000.bat
33a1c790820c4accea97d1b48bf9dbaf47a9c22a2b04b40a9f463e21f0d42393  analysis/topup_p20000.R
76cf7976a6e9b783e141ae31e145c11f76277f774b97eb9f6427f7238013fb9a  results/GSE14520_v2_diagnostics.csv
edcd8ad9d74328ad76f347e62848430e8a791e70c8f32510cb2130f0f3d0a666  results/GSE14520_v2_paired.csv
08955f3a8958bd6b0d93eb1732500214e2d8aa30b9b72e1a3da07d94d216e6da  results/GSE14520_v2_results.csv
1d4f1de11d8a1bf06cbd633c20857868aaa4f89c2256190e7227de98dd670140  results/GSE14520_v2_stability.csv
beaf06a4a36302a90a3e35f0a69135766b5dfcc3a0f4c13d84cae53057e1d464  results/GSE14520_v2_stability_full.csv
0ecd820d93195bef54e234008efdd6cc072a5e9047ecadc7b72801b0ed7548d1  results/GSE14520_v2_summ.csv
22a91504d75cd9720632d0a1e447534f33aea274e8c99b5cf7d0e0f68da109ab  superseded/GSE14520_Real_Analysis_1.Rmd
0198b5c9867685d23150279be2ede3bb7dbbf992acb65541bbd90f6e58d24397  superseded/S2_cap_sensitivity_raw.csv
7c86b55031a2a9d245dea72088e4f3fde894b45d35c7c65a35120aa403daddeb  superseded/S2_convergence_raw.csv
50ca728885a916b768ce444621e866c374b1b817c7e9104cb37da3e1c58090db  superseded/run_shard_PRE_PATCH_2026-09-01.R
55902d9be72e8de24cc94e8980b49d7733ec4a82231aa38146c107bfb71e2664  verification/S1_scale_verification.R
fd75070bd1fa78008281685e0d7790fb84f5e22d12c60c26d9880458ed14ceb0  verification/S5_penalty_scale.R
653347e3e56028bda455d5e6853d4d79f99766bae52459e63e42843458a86c00  verification/S5_penalty_scale_raw.csv
927465598b5a3e9a059435b7bb44968e6ea497969ac71e5f062f3676a9f5e5a2  verification/S6_patch_verification_raw.csv
db8f7e7823e5cf6f3ab3d3ccb1750f3274ac10737aef4b0f53cb365f21c3bd7b  verification/S6_verify_patch.R
52552adaa20515f4ba60c9981746b94481c5ff4a764041e60599bea7d45366f7  verification/S7_pilot_metrics.csv
cf5ff653f23dc71f81c250fe129e8b5bad016cafeec63806b5142107a38b5e3c  verification/S7_pilot_rerun.R
ff5f49ee31edd7c983953055ccdc8cf839d2032f1bf2dc88dd299248c51e5491  verification/S8_adaptive_audit.csv
ede464b6249dcda7a7f0131e4f27e85d0cde6ea2783c2cdfddb51fd85a3129af  verification/S8_cap_sensitivity.csv
e824df61a9e0eb14982295c874fa64cb93c57e2f15ae0cf6eb482d150f0029a9  verification/S8_corrected_diagnostics.R
ad67c1e09e5f667e090e60b874b96ec9080390abb761fbad5188a274a742997d  verification/S8_fixed_audit.csv
c40d0485ebac499a0cf39fad3063c5aa51b364b89c8c5d062129dfabaafc21b9  verification/S9_fixed_vs_adaptive.R
b552f8647ee5dbcf20c4927417ea201729d0c0db1d60275a0c9879c921f50b97  verification/S9_fixed_vs_adaptive.csv
```