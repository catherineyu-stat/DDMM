# Filter comparison

Scripts copied unchanged from `code_official/synthetic data/admixture` and
`family`, with their original DDMM cores. These replace the earlier adapted runners.

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

**Output:** `<experiment>/results/result-state<state>/<p>-results.csv`;
PDF/PNG/TIFF figures in `figure/`. Each CSV reports FDR, power, MSPE, and L1 error.

## Settings in the local scripts

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
`LMM-Lasso` is Lasso on the LMM-filtered data. DDMM uses the original
`multiple_testing()` threshold rule and post-selection estimation.

## Local fixed-core variant

The supplied `*_fixed.R` drivers retain the same simulation settings and use
`inference_short_fixed.R`, which selects the nodewise penalty separately for
each SNP by BIC. Outputs go to `results_fixed/`.

```sh
(cd admixture && RUN_ORIGINAL=1 Rscript compare_new_extend_fixed.R 1)
(cd family && RUN_ORIGINAL=1 Rscript compare_family_extend_fixed.R 1)
python3 plot_fix_comparison.py
```

`RUN_ORIGINAL=1` also evaluates the original core on the same data. The paired
plot requires Python packages `numpy`, `pandas`, and `matplotlib` and writes to
`figure_fixed/`. Change settings in the corresponding `*_fixed.R` driver when
using this variant.
