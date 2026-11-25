## ============================================================================
## Section 5.1 — Admixture PCA-O / PCA-E comparison driver
##
## Compares, on simulated admixed data, six pipelines that share the same
## raw genotype + phenotype:
##
##   * ggmix       : penalised LMM (kinship-aware, no spectral pre-filter)
##   * plmmr       : another penalised LMM
##   * PCA-E       : PCA-Lasso, where Q is the projection off the top-k
##                   left singular vectors of an *external* marker set X_GRM
##   * DDMM-E      : debiased DDMM with the same Q (= PCA-E filter)
##   * PCA-O       : PCA-Lasso, where Q comes from the test SNPs themselves
##                   (no external markers)
##   * DDMM-O      : debiased DDMM with the PCA-O filter
##   * DDMM-Trim   : DDMM with the TRIM filter (high-rank baseline)
##
## All six DDMM variants now share a single API:
##
##     Q     <- build_filter("decorrelation", X = ..., q = 5)
##     stat  <- test_stat_decorrelation_short(Y, X, transform_matrix = Q)
##     out   <- run_ddmm_variant("DDMM-X", stat,
##                               list(transform_matrix = Q, centralize = TRUE))
##
## which lives in `inference_short_refactored.R`. The trace-Lasso machinery
## from the old `solve_trace_lasso_flexible` is gone; PCA-Lasso is now
## simply Lasso applied to (Q Y, Q X), implemented in `pca_lasso()` inside
## `utils_extend.R`. This makes "PCA-Lasso vs DDMM" a clean head-to-head
## using literally the same filter.
##
## Run as:    Rscript compare_new_extend.R <state>
##
## ---------------------------------------------------------------------------
## FIXED VERSION (compare_new_extend_fixed.R)
##
## Identical to compare_new_extend.R except for three lines of plumbing:
##   * it sources inference_short_fixed.R, which applies FIX 1 (per-variant
##     BIC nodewise lambda) to test_stat_decorrelation_short();
##   * it writes to results_fixed/ so the original results/ are not touched;
##   * an optional second argument runs a single p and checkpoints the RNG
##     state, so a state can be split across several short jobs and still
##     draw exactly the same random numbers as a single uninterrupted run.
##
##   Rscript compare_new_extend_fixed.R <state>            # all five p, as before
##   Rscript compare_new_extend_fixed.R <state> <p_input>  # one p, resumable
##
## The seed, the data generation and every non-DDMM method are unchanged.
## ---------------------------------------------------------------------------
## ============================================================================

Args  <- commandArgs(trailingOnly = TRUE)
state <- as.integer(Args[1])
if (is.na(state)) state <- 0L
p_only <- if (length(Args) >= 2) as.integer(Args[2]) else NA_integer_

source("utils_extend.R")
## DDMM core with FIX 1 applied on top of the unchanged refactored core
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

set.seed(2026 + state)


## ---------------------------------------------------------------------------
## helper : run the full per-(p_input) sweep and append into results_list
## ---------------------------------------------------------------------------
run_one_p <- function(p_input, state) {

  cat(sprintf("\n=========== p = %d, state = %d ===========\n", p_input, state))
  n = 1000
  
  beta_strength = 10 * sqrt(log(p_input) / n*0.8)
  ## ---- generate admixed data -------------------------------------------
  admixed <- gen_ddmm_data(
    n          = 1000,
    p_design   = p_input,
    p_kinship  = 10 * p_input,
    n_causal   = 30,            # absolute, not a fraction
    k = 5, s = 0.5, Fst = 0.1,
    eta = 0.1, sigma2 = 1,
    beta_strength = beta_strength,
    geography = "1d",
    train_frac = 0.8
  )
  
  

  X_train   <- admixed$xtrain
  Y_train   <- admixed$ytrain
  kin_train <- admixed$kin_train
  X_test    <- admixed$xtest
  Y_test    <- admixed$ytest
  X_GRM     <- admixed$Xkinship          # external marker set, used by PCA-E
  true_beta <- admixed$beta

  ## drop zero-variance SNPs from training & test
  zero_idx <- which(apply(X_train, 2, var) == 0)
  if (length(zero_idx) > 0) {
    X_train   <- X_train[, -zero_idx]
    X_test    <- X_test[, -zero_idx]
    true_beta <- true_beta[-zero_idx]
    if (!is.null(X_GRM) && ncol(X_GRM) >= max(zero_idx)) {
      ## only drop if the external set indexes overlap with X (defensive)
      keep <- setdiff(seq_len(ncol(X_GRM)), zero_idx)
      X_GRM <- X_GRM[, keep, drop = FALSE]
    }
  }

  p              <- ncol(X_train)
  causal_indices <- which(true_beta != 0)
  results_list   <- list()
  orig_list      <- list()   # original-core DDMM on the same data (paired)

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


  ## ===================================================================
  ## ggmix  (penalised LMM, no spectral filter)
  ## ===================================================================
  cat("--- ggmix ---\n")
  fit_ggmix       <- ggmix(x = X_train, y = Y_train, kinship = kin_train)
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
  ## PCA-E  (filter from external markers X_GRM)  +  DDMM-E
  ## ===================================================================
  cat("--- PCA-E / DDMM-E ---\n")
  Q_E <- build_filter("decorrelation", X = X_GRM, q = 5, centralize = TRUE)

  ## (a) PCA-Lasso with Q_E as the filter
  pca_e_fit <- pca_lasso(Y_train, X_train, Q = Q_E, centralize = TRUE)
  beta_pca_e <- pca_e_fit$beta
  results_list$PCA_E <- list(
    Method   = "PCA-E",
    MSPE     = calculate_mspe(beta_pca_e, X_test, Y_test),
    FDR      = calculate_fdr_power(beta_pca_e, true_beta)$FDR,
    Power    = calculate_fdr_power(beta_pca_e, true_beta)$Power,
    L1_Error = calculate_l1_error(beta_pca_e, true_beta)
  )

  ## (b) DDMM with the same Q_E
  if (RUN_ORIGINAL) orig_list$DDMM_E <- ddmm_original_variant("DDMM-E", Y_train, X_train, Q_E)
  ddmm_e_stat <- test_stat_decorrelation_short(
    Y = Y_train, X = X_train,
    transform_matrix = Q_E, centralize = TRUE
  )
  ddmm_e_var <- run_ddmm_variant(
    "DDMM-E", ddmm_e_stat,
    list(transform_matrix = Q_E, centralize = TRUE)
  )
  results_list$DDMM_E <- ddmm_e_var$post


  ## ===================================================================
  ## PCA-O  (filter from the test SNPs themselves)  +  DDMM-O
  ## ===================================================================
  cat("--- PCA-O / DDMM-O ---\n")
  Q_O <- build_filter("decorrelation", X = X_train, q = 5, centralize = TRUE)

  pca_o_fit <- pca_lasso(Y_train, X_train, Q = Q_O, centralize = TRUE)
  beta_pca_o <- pca_o_fit$beta
  results_list$PCA_O <- list(
    Method   = "PCA-O",
    MSPE     = calculate_mspe(beta_pca_o, X_test, Y_test),
    FDR      = calculate_fdr_power(beta_pca_o, true_beta)$FDR,
    Power    = calculate_fdr_power(beta_pca_o, true_beta)$Power,
    L1_Error = calculate_l1_error(beta_pca_o, true_beta)
  )

  if (RUN_ORIGINAL) orig_list$DDMM_O <- ddmm_original_variant("DDMM-O", Y_train, X_train, Q_O)
  ddmm_o_stat <- test_stat_decorrelation_short(
    Y = Y_train, X = X_train,
    transform_matrix = Q_O, centralize = TRUE
  )
  ddmm_o_var <- run_ddmm_variant(
    "DDMM-O", ddmm_o_stat,
    list(transform_matrix = Q_O, centralize = TRUE)
  )
  results_list$DDMM_O <- ddmm_o_var$post


  ## ===================================================================
  ## DDMM-Trim  (TRIM filter; high-rank baseline)
  ## ===================================================================
  cat("--- DDMM-Trim ---\n")
  Q_trim <- build_filter("trim", X = X_train, rho = 0.3, centralize = TRUE)
  if (RUN_ORIGINAL) orig_list$DDMM_Trim <- ddmm_original_variant("DDMM-Trim", Y_train, X_train, Q_trim)
  ddmm_trim_stat <- test_stat_decorrelation_short(
    Y = Y_train, X = X_train,
    transform_matrix = Q_trim, centralize = TRUE
  )
  ddmm_trim_var <- run_ddmm_variant(
    "DDMM-Trim", ddmm_trim_stat,
    list(transform_matrix = Q_trim, centralize = TRUE)
  )
  results_list$DDMM_Trim <- ddmm_trim_var$post


  ## ===================================================================
  ## summarise + persist
  ## ===================================================================
  results_df <- bind_rows(results_list)
  cat("\n--- Final Performance Comparison (p =", p_input, ") ---\n")
  print(results_df, digits = 4)

  out_dir <- file.path(OUT_ROOT,
                       paste0("result-state", state))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df,
            file = file.path(out_dir, paste0(p_input, "-results.csv")),
            row.names = FALSE)
  if (RUN_ORIGINAL)
    write.csv(bind_rows(orig_list),
              file = file.path(out_dir, paste0(p_input, "-results_original_core.csv")),
              row.names = FALSE)
  invisible(results_df)
}


## ---------------------------------------------------------------------------
## main loop
## ---------------------------------------------------------------------------
p_seq <- c(1000, 2000, 4000, 6000, 8000)
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
