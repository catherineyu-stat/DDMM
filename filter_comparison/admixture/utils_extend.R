# Admixture generation, PCA-Lasso, and performance metrics.
# Shared DDMM inference is in ../inference.R.

tryCatch({
  source("ggmix_utils.R")
}, warning = function(w) {
  message("Warning: Could not source helper files. Using placeholder ",
          "fallbacks. Please check file paths.")
})

if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  glmnet,
  ggmix,
  plmmr,
  dplyr,
  RSpectra
)


## ---------------------------------------------------------------------------
## Metric helpers
## ---------------------------------------------------------------------------

calculate_mspe <- function(beta_est, X_new, y_true) {
  y_pred <- X_new %*% as.numeric(beta_est)
  mean((as.numeric(y_true) - as.numeric(y_pred))^2)
}

calculate_l1_error <- function(beta_est, beta_true) {
  beta_est  <- as.numeric(beta_est)
  beta_true <- as.numeric(beta_true)
  if (length(beta_est) != length(beta_true))
    stop("Estimated and true beta vectors must have the same length.")
  sum(abs(beta_est - beta_true))
}

calculate_fdr_power <- function(beta_est, beta_true) {
  causal   <- which(beta_true != 0)
  selected <- which(beta_est  != 0)
  if (length(selected) == 0) return(list(FDR = 0, Power = 0))
  fp <- length(setdiff(selected, causal))
  tp <- length(intersect(selected, causal))
  list(FDR = fp / length(selected), Power = tp / length(causal))
}


## ---------------------------------------------------------------------------
## pca_lasso() : PCA-Lasso with a user-supplied spectral filter Q
##
## This replaces the trace-Lasso route used previously. The pipeline is now
## identical to what DDMM does in its initial step (filter Y, X by Q; run
## cv.glmnet on the filtered pair) -- so the PCA-O / PCA-E numbers and the
## downstream DDMM-O / DDMM-E numbers come from the *same* Q, which is what
## the comparison table in Section 5.1 is meant to show.
##
## Arguments
##   Y, X        : response and design (untransformed; centring is applied
##                 inside the function when `centralize = TRUE`).
##   Q           : n x n spectral filter (build via build_filter()). PCA-O
##                 uses Q built from X; PCA-E uses Q built from X_GRM.
##   centralize  : same flag used in build_filter() and DDMM, so all three
##                 stages agree on the centred design.
##   nlambda     : grid length for the cv.glmnet path.
##   nfolds      : K-fold CV.
## ---------------------------------------------------------------------------
pca_lasso <- function(Y, X, Q, centralize = TRUE,
                      nlambda = 100, nfolds = 10) {
  n <- nrow(X); p <- ncol(X)
  if (!all(dim(Q) == c(n, n)))
    stop("Q must have dimension nrow(X) by nrow(X).")

  Y <- as.numeric(scale(Y, centralize, FALSE))
  X <- scale(X, centralize, FALSE)

  F_Y <- as.numeric(Q %*% Y)
  F_X <- as.matrix(Q %*% X)

  ## same lambda grid scaling as the DDMM initial Lasso for consistency
  sequence <- (nlambda:1) / nlambda * 3
  lambda <- sqrt(sum(F_Y^2) / n) * sqrt(log(max(p, 2)) / n) * sequence

  cv_fit <- glmnet::cv.glmnet(
    x = F_X, y = F_Y, lambda = lambda,
    alpha = 1, intercept = FALSE, standardize = FALSE,
    family = "gaussian", nfolds = min(nfolds, n)
  )
  beta <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
  list(beta = beta, lambda = cv_fit$lambda.min, fit = cv_fit)
}
