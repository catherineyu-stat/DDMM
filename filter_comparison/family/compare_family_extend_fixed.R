## ============================================================================
## Section 5.2 — Family-data PCA / DDMM comparison driver
##
## Sister script to `compare_new_extend.R` in the Admixture folder. The two
## files share the same orchestration, the same per-method block layout,
## the same global-env handoff to `run_ddmm_variant()`, and the same
## per-(state, p) CSV output convention. The only differences are
##
##   (a) the data source -- BED genotype files via `gaston`, instead of
##       `gen_ddmm_data()` simulated genotypes;
##   (b) two extra methods that only make sense with real / synthetic-real
##       genotype data: BH (gaston's marginal LMM score test) and LMM
##       (cv.glmnet on the LMM-whitened design);
##   (c) MSPE and L1 are now reported for BH as well, via a post-selection
##       OLS on the discovered SNPs inside the same LMM whitening geometry.
##
## Methods compared (one CSV row each, per (state, p)):
##
##   * BH        : marginal LMM-score test + Benjamini-Hochberg
##   * ggmix     : penalised LMM with HDBIC
##   * plmmr     : preconditioned Lasso with 5-fold CV
##   * DDMM-Trim : DDMM with the TRIM spectral filter
##   * LMM       : Lasso applied to the LMM-whitened design (no debiasing)
##   * DDMM-LMM  : DDMM with the LMM whitening filter
##   * DDMM-PCA  : DDMM with the PCA-projection filter (q = 10)
##
## All four DDMM variants and PCA-Lasso share the same single API
##
##     Q     <- build_filter(method, ...)
##     stat  <- test_stat_decorrelation_short(Y, X, transform_matrix = Q)
##     out   <- run_ddmm_variant("DDMM-X", stat,
##                               list(transform_matrix = Q, centralize = TRUE))
##
## which lives in `inference_short_refactored.R`.
##
## Run as:    Rscript compare_family_extend.R <state>
##
## ---------------------------------------------------------------------------
## FIXED VERSION (compare_family_extend_fixed.R)
##
## Identical to compare_family_extend.R except for three lines of plumbing:
##   * it sources inference_short_fixed.R, which applies FIX 1 (per-variant
##     BIC nodewise lambda) to test_stat_decorrelation_short();
##   * it writes to results_fixed/ so the original results are not touched;
##   * an optional second argument runs a single p and checkpoints the RNG
##     state, so a state can be split across several short jobs and still
##     draw exactly the same random numbers as a single uninterrupted run.
##
##   Rscript compare_family_extend_fixed.R <state>            # all five p
##   Rscript compare_family_extend_fixed.R <state> <p_input>  # one p, resumable
##
## The seed, the data generation and every non-DDMM method are unchanged.
## ---------------------------------------------------------------------------
## ============================================================================

if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  ggmix, plmmr, dplyr, RSpectra, gaston, glmnet
)

Args  <- commandArgs(trailingOnly = TRUE)
state <- as.integer(Args[1])
if (is.na(state)) state <- 0L
p_only <- if (length(Args) >= 2) as.integer(Args[2]) else NA_integer_

## ---- helpers + DDMM core with FIX 1 --------------------------------------
source("helper.R")                                   # calculate_l1_error / fdr_power
source("inference_short_fixed.R")
options(ddmm.cores = as.integer(Sys.getenv("DDMM_CORES", "1")))
## Optional, for short-lived sandbox sessions only (see sandbox_checkpoint_hooks.R).
## Unset on a cluster, where it has no effect.
if (nzchar(Sys.getenv("DDMM_HOOKS"))) source(Sys.getenv("DDMM_HOOKS"))
OUT_ROOT <- file.path(".", "results_fixed")
## RUN_ORIGINAL=1 also runs the original DDMM core on the same data and writes
## it to <p>-results_original_core.csv, for a paired before/after comparison.
## It does not change any number in <p>-results.csv.
RUN_ORIGINAL <- Sys.getenv("RUN_ORIGINAL", "0") == "1"
## NB: sourcing the refactored core *after* helper.R intentionally overrides
## helper.R's legacy `run_sun_variant` with the new `run_ddmm_variant` alias.

## run_ddmm_variant() expects calculate_mspe() in scope; define a local one
## so the family driver doesn't depend on utils_extend.R.
calculate_mspe <- function(beta_est, X_new, y_true) {
  mean((as.numeric(y_true) -
        as.numeric(X_new %*% as.numeric(beta_est)))^2)
}

set.seed(2026 + state)


# ============================================================================
# 1. LOAD DATA  (this is the *only* block that meaningfully differs from
#    compare_new_extend.R; everything below the data block is shared logic)
# ============================================================================
file_prefix <- "family_data_synthetic/family_data_syn"

message("Loading full genotype data...")
x_full     <- read.bed.matrix(file_prefix)
x_full     <- select.snps(x_full, maf > 0.05)
current_p  <- ncol(x_full)

message("Loading genotype data (GRM)...")
file_prefix_grm <- "family_data_synthetic/family_data_syn_chr2"
x_grm     <- read.bed.matrix(file_prefix_grm)
x_grm     <- select.snps(x_grm, maf > 0.05)

message("Computing true GRM...")
K_true   <- GRM(x_grm)
h2       <- 0.5
eig_true <- eigen(K_true)
eig_true$values[eig_true$values < 0] <- 0
K_sqrt   <- eig_true$vectors %*% diag(sqrt(eig_true$values))


# ============================================================================
# 2. SIMULATION PARAMETERS
# ============================================================================
n_causal     <- 30
fixed_scale  <- 10
train_frac   <- 0.8

p_seq <- c(1000, 2000, 4000, 6000, 8000)


# ============================================================================
# 3. PER-(p_input) RUNNER  (mirrors the Admixture compare_new_extend.R layout)
# ============================================================================
run_one_p <- function(p_target, state) {

  cat(sprintf("\n=========== p = %d, state = %d ===========\n", p_target, state))

  ## ---- A. Sample p SNPs and standardise (full sample) -----------------
  sample_p <- sample(seq_len(current_p), p_target)
  x_all    <- x_full[, sample_p]
  p_vec    <- x_all@p
  X_all    <- as.matrix(x_all)

  for (j in seq_len(ncol(x_all))) {
    mu_j    <- 2 * p_vec[j]
    sigma_j <- sqrt(2 * p_vec[j] * (1 - p_vec[j]))
    if (sigma_j == 0) sigma_j <- 1
    X_all[, j] <- (X_all[, j] - mu_j) / sigma_j
    X_all[is.na(X_all[, j]), j] <- 0
  }
  zero_idx <- which(apply(X_all, 2, var) < 1e-8)
  if (length(zero_idx) > 0) {
    X_all <- X_all[, -zero_idx]
    x_all <- x_all[, -zero_idx]
  }
  N_total <- nrow(x_all)
  M       <- ncol(x_all)

  ## ---- B. Generate phenotype on the full sample -----------------------
  u <- as.numeric(K_sqrt %*% rnorm(N_total)) * sqrt(h2)
  e <- rnorm(N_total) * sqrt(1 - h2)

  common_snps <- which(x_all@snps$maf > 0.3)
  if (length(common_snps) < n_causal) common_snps <- seq_len(M)
  common_snps    <- sample(common_snps, n_causal)
  causal_indices <- common_snps

  betas             <- rnorm(n_causal, mean = 0,
                             sd = fixed_scale * sqrt(log(M) / N_total))
  true_beta         <- numeric(M)
  true_beta[common_snps] <- betas

  genetic_signal <- as.numeric(X_all %*% true_beta)
  y_full         <- as.vector(scale(u + e + genetic_signal))

  ## ---- C. 80/20 train-test split --------------------------------------
  n_train   <- floor(train_frac * N_total)
  train_idx <- sort(sample(seq_len(N_total), n_train))
  test_idx  <- setdiff(seq_len(N_total), train_idx)

  X_train <- X_all[train_idx, , drop = FALSE]
  X_test  <- X_all[test_idx,  , drop = FALSE]
  Y_train <- y_full[train_idx]
  Y_test  <- y_full[test_idx]

  x_train <- x_all[train_idx, ]
  x_train@ped$pheno <- Y_train

  K_train <- K_true[train_idx, train_idx]
  p       <- M

  results_list <- list()
  orig_list    <- list()   # original-core DDMM on the same data (paired)

  ## ---- expose the variables that run_ddmm_variant() reads from --------
  ## (the function uses lexical scoping: causal_indices, true_beta,
  ## X_train, Y_train, X_test, Y_test, calculate_mspe, calculate_l1_error)
  assign("causal_indices", causal_indices, envir = .GlobalEnv)
  assign("true_beta",      true_beta,      envir = .GlobalEnv)
  assign("X_train",        X_train,        envir = .GlobalEnv)
  assign("Y_train",        Y_train,        envir = .GlobalEnv)
  assign("X_test",         X_test,         envir = .GlobalEnv)
  assign("Y_test",         Y_test,         envir = .GlobalEnv)
  assign("p",              p,              envir = .GlobalEnv)

  ## ---- Pre-build shared spectral filters (re-used downstream) --------
  Q_trim <- build_filter("trim",          X = X_train,
                         rho = 0.3,        centralize = TRUE)
  Q_pca  <- build_filter("decorrelation", X = X_train,
                         q   = 10,         centralize = TRUE)
  Q_lmm  <- build_filter("lmm",           kinship = K_train, h2 = h2)


  ## ===================================================================
  ## BH  (marginal LMM score test + post-selection OLS for MSPE / L1)
  ## ===================================================================
  cat("--- BH ---\n")
  null_model  <- lmm.aireml(Y = Y_train, K = K_train, verbose = FALSE)
  res         <- association.test(x_train, Y = Y_train, K = K_train,
                                  method = "lmm", response = "quantitative",
                                  test = "score")
  res$padj    <- p.adjust(res$p, method = "BH")
  sig_indices <- which(res$padj <= 0.1)
  tp_count    <- sum(causal_indices %in% sig_indices)
  fdp_bh      <- if (length(sig_indices) == 0) 0 else
                 (length(sig_indices) - tp_count) / length(sig_indices)
  power_bh    <- tp_count / n_causal

  ## post-selection OLS in the same LMM whitening geometry so MSPE and L1
  ## are directly comparable to DDMM-LMM and to LMM Lasso.
  beta_bh_full <- numeric(M)
  if (length(sig_indices) > 0) {
    post_bh <- post_estimation(
      Y                = Y_train,
      X                = X_train,
      index            = sig_indices,
      transform_matrix = Q_lmm,
      centralize       = TRUE
    )
    beta_bh_full[sig_indices] <- post_bh$beta
  }
  results_list$BH <- list(
    Method   = "BH",
    MSPE     = calculate_mspe(beta_bh_full, X_test, Y_test),
    FDR      = fdp_bh,
    Power    = power_bh,
    L1_Error = calculate_l1_error(beta_bh_full, true_beta)
  )


  ## ===================================================================
  ## ggmix  (penalised LMM, no spectral filter)
  ## ===================================================================
  cat("--- ggmix ---\n")
  fit_ggmix       <- ggmix(x = X_train, y = Y_train, kinship = K_train)
  beta_ggmix      <- as.numeric(coef(gic(fit_ggmix)))[-1]
  results_list$ggmix <- list(
    Method   = "ggmix",
    MSPE     = calculate_mspe(beta_ggmix, X_test, Y_test),
    FDR      = calculate_fdr_power(beta_ggmix, true_beta)$FDR,
    Power    = calculate_fdr_power(beta_ggmix, true_beta)$Power,
    L1_Error = calculate_l1_error(beta_ggmix, true_beta)
  )


  ## ===================================================================
  ## plmmr
  ## ===================================================================
  cat("--- plmmr ---\n")
  admix_design <- create_design(X = X_train, y = Y_train)
  cv_plmmr     <- cv_plmm(design = admix_design)
  fit_plmmr    <- plmm(design = admix_design)
  beta_plmmr   <- as.numeric(coef(fit_plmmr, lambda = cv_plmmr$lambda_min))[-1]
  results_list$plmmr <- list(
    Method   = "plmmr",
    MSPE     = calculate_mspe(beta_plmmr, X_test, Y_test),
    FDR      = calculate_fdr_power(beta_plmmr, true_beta)$FDR,
    Power    = calculate_fdr_power(beta_plmmr, true_beta)$Power,
    L1_Error = calculate_l1_error(beta_plmmr, true_beta)
  )


  ## ===================================================================
  ## DDMM-Trim   (uniform API; TRIM filter)
  ## ===================================================================
  cat("--- DDMM-Trim ---\n")
  if (RUN_ORIGINAL) orig_list$DDMM_Trim <- ddmm_original_variant("DDMM-Trim", Y_train, X_train, Q_trim)
  ddmm_trim_stat <- test_stat_decorrelation_short(
    Y                = Y_train, X = X_train,
    transform_matrix = Q_trim,  centralize = TRUE
  )
  ddmm_trim_var <- run_ddmm_variant(
    "DDMM-Trim", ddmm_trim_stat,
    list(transform_matrix = Q_trim, centralize = TRUE)
  )
  results_list$DDMM_Trim <- ddmm_trim_var$post


  ## ===================================================================
  ## LMM Lasso  (cv.glmnet on the LMM-whitened design, no debiasing)
  ## ===================================================================
  cat("--- LMM Lasso ---\n")
  F_X_lmm   <- as.matrix(Q_lmm %*% X_train)
  F_Y_lmm   <- as.numeric(Q_lmm %*% Y_train)
  fit_lmm   <- cv.glmnet(F_X_lmm, F_Y_lmm,
                         intercept = FALSE, standardize = FALSE)
  coefs_lmm <- as.numeric(coef(fit_lmm, s = "lambda.min"))[-1]
  sig_lmm   <- which(coefs_lmm != 0)
  tp_lmm    <- sum(causal_indices %in% sig_lmm)
  results_list$LMM <- list(
    Method   = "LMM",
    MSPE     = calculate_mspe(coefs_lmm, X_test, Y_test),
    FDR      = if (length(sig_lmm) == 0) 0 else
               (length(sig_lmm) - tp_lmm) / length(sig_lmm),
    Power    = tp_lmm / n_causal,
    L1_Error = calculate_l1_error(coefs_lmm, true_beta)
  )


  ## ===================================================================
  ## DDMM-LMM   (uniform API; LMM whitening filter -- same Q_lmm as above)
  ## ===================================================================
  cat("--- DDMM-LMM ---\n")
  if (RUN_ORIGINAL) orig_list$DDMM_LMM <- ddmm_original_variant("DDMM-LMM", Y_train, X_train, Q_lmm)
  ddmm_lmm_stat <- test_stat_decorrelation_short(
    Y                = Y_train, X = X_train,
    transform_matrix = Q_lmm,   centralize = TRUE
  )
  ddmm_lmm_var <- run_ddmm_variant(
    "DDMM-LMM", ddmm_lmm_stat,
    list(transform_matrix = Q_lmm, centralize = TRUE)
  )
  results_list$DDMM_LMM <- ddmm_lmm_var$post


  ## ===================================================================
  ## DDMM-PCA   (uniform API; PCA projection filter, q = 10)
  ## ===================================================================
  cat("--- DDMM-PCA ---\n")
  if (RUN_ORIGINAL) orig_list$DDMM_PCA <- ddmm_original_variant("DDMM-PCA", Y_train, X_train, Q_pca)
  ddmm_pca_stat <- test_stat_decorrelation_short(
    Y                = Y_train, X = X_train,
    transform_matrix = Q_pca,   centralize = TRUE
  )
  ddmm_pca_var <- run_ddmm_variant(
    "DDMM-PCA", ddmm_pca_stat,
    list(transform_matrix = Q_pca, centralize = TRUE)
  )
  results_list$DDMM_PCA <- ddmm_pca_var$post


  ## ===================================================================
  ## summarise + persist (same CSV layout as the Admixture driver)
  ## ===================================================================
  results_df <- bind_rows(results_list)
  cat(sprintf("\n--- Final Performance Comparison (p = %d) ---\n", p_target))
  print(results_df, digits = 4)

  out_dir <- file.path(OUT_ROOT,
                       paste0("result-state", state))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df,
            file = file.path(out_dir, paste0(p_target, "-results.csv")),
            row.names = FALSE)
  if (RUN_ORIGINAL)
    write.csv(bind_rows(orig_list),
              file = file.path(out_dir, paste0(p_target, "-results_original_core.csv")),
              row.names = FALSE)
  invisible(results_df)
}


# ============================================================================
# 4. MAIN LOOP
# ============================================================================
rng_file <- function(p) file.path(OUT_ROOT, paste0("result-state", state),
                                  paste0(".rng_after_", p, ".rds"))
save_rng <- function(p) {
  dir.create(dirname(rng_file(p)), recursive = TRUE, showWarnings = FALSE)
  saveRDS(get(".Random.seed", envir = globalenv()), rng_file(p))
}

if (is.na(p_only)) {
  for (p_input in p_seq) {
    run_one_p(p_input, state)
    save_rng(p_input)
  }
} else {
  k <- match(p_only, p_seq)
  if (is.na(k)) stop("p_input must be one of ", paste(p_seq, collapse = ", "))
  if (k > 1) {
    ## resume the random stream exactly where the previous p left it
    prev <- rng_file(p_seq[k - 1])
    if (!file.exists(prev)) stop("Run p = ", p_seq[k - 1], " for this state first.")
    assign(".Random.seed", readRDS(prev), envir = globalenv())
  }
  run_one_p(p_only, state)
  save_rng(p_only)
}
