# Filter comparison

Admixture and family simulations share `inference.R`: DDMM with a separate
BIC-selected nodewise penalty for each SNP.

## Install

```r
install.packages(c("pacman", "bnpsd", "popkin", "simfam", "kinship2", "gaston",
                   "glmnet", "ggmix", "plmmr", "dplyr", "RSpectra",
                   "ggplot2", "tidyr", "patchwork", "ragg"))
```

## Run

From `filter_comparison/`, run one state across all five values of `p`:

```sh
(cd admixture && Rscript compare_new_extend.R 1)
(cd family && Rscript simfam.R)
(cd family && Rscript compare_family_extend.R 1)
Rscript plot_figure2.R
```

Generate the family panels once. Repeat the association commands with state
numbers `1` through `50`; their seed is `2026 + state`.

**Output:** `<experiment>/results_bic/result-state<state>/<p>-results.csv`;
PDF/PNG/TIFF figures in `figure/`. Each CSV reports FDR, power, MSPE, and L1 error.

## Parameters

| Setting | Admixture | Family |
|---|---|---|
| Individuals | 1,000 | 100 per generation × 3 generations = 300 |
| Tested variants `p` | 1,000; 2,000; 4,000; 6,000; 8,000 | Same grid |
| Causal variants | 30 | 30 |
| GRM markers | `10 * p` external BN-PSD markers | Second 10,000-marker panel, before MAF filtering |
| Random-effect / iid-noise variance | 0.1 / 0.9 | 0.5 / 0.5, before phenotype standardization |
| Training / test split | 80% / 20% | 80% / 20% |
| PCA directions removed | 5, from tested SNPs (O) or external markers (E) | 10 |
| TRIM `rho` | 0.3 | 0.3 |
| DDMM FDR threshold | 0.1 | 0.1 |

Edit `admixture/compare_new_extend.R` or `family/compare_family_extend.R`
for the association settings; edit `family/simfam.R` for the pedigree and
panel sizes. The family generator uses seed `2026` and reuses `X_1` for the
second panel, exactly as in the local source.

Admixture compares ggmix, plmmr, PCA-Lasso-O/E, DDMM-PCA-O/E, and DDMM-Trim.
Family compares LMM-BH, ggmix, plmmr, LMM-Lasso, DDMM-LMM, DDMM-PCA, and DDMM-Trim.
`LMM-Lasso` is Lasso on the LMM-filtered data. DDMM uses `multiple_testing()` for thresholding and then post-selection estimation.
Set its FDR level in `inference.R` (`run_ddmm_variant`, `alpha = 0.1`).
