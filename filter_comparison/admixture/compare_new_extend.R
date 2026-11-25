# Admixture simulation with per-SNP BIC nodewise inference.
# Run: Rscript compare_new_extend.R <state>

Args  <- commandArgs(trailingOnly = TRUE)
state <- as.integer(Args[1])
if (is.na(state)) state <- 0L

source("utils_extend.R")
source("../inference.R")

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

  out_dir <- file.path(".", "results_bic",
                       paste0("result-state", state))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df,
            file = file.path(out_dir, paste0(p_input, "-results.csv")),
            row.names = FALSE)
  invisible(results_df)
}


## ---------------------------------------------------------------------------
## main loop
## ---------------------------------------------------------------------------
p_seq <- c(1000, 2000, 4000, 6000, 8000)
for (p_input in p_seq) {
  run_one_p(p_input, state)
}
