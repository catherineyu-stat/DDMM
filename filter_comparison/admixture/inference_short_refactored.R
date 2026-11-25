## ============================================================================
## DDMM core: spectral-deconfounding debiased inference for structured GWAS
##
## Design choices that distinguish this implementation from the Sun et al.
## (2024) reference code:
##
##   1. The spectral filter Q is supplied externally as `transform_matrix`.
##      The earlier code dispatched on `initial = "trim" | "decorrelation"
##      | "LMM"` and built Q internally; that conflated filter construction
##      with inference and made it awkward to test arbitrary user-defined
##      filters. Q is now a first-class input, constructed once via
##      `build_filter()` and re-used by both the inference step and the
##      post-estimation step.
##
##   2. When the initial Lasso recovers fewer than `min_selected` non-zero
##      coordinates, the procedure falls back to a marginal-correlation
##      screening (`marginal_screen_signal()`). This is the same screening
##      principle that BASIL (snpnet, GhostKnockoff) uses to bootstrap its
##      working set, here adapted to keep DDMM well-defined when the active
##      set is too small for the nodewise step to be informative.
##
##   3. The nodewise λ on each coordinate is selected by BIC over a fixed
##      grid, removing the need for a separate cv.glmnet call inside the
##      inner loop.
##
## Naming: the post-estimation runner is now `run_ddmm_variant()` to reflect
## that the pipeline is no longer a wrapper around Sun's code path.
## `run_sun_variant` is kept as an alias for backward compatibility.
## ============================================================================


## ---------------------------------------------------------------------------
## build_filter() : construct one of the three spectral filters Q in O(n^2 r)
##
## Returns an n x n matrix Q to be applied to both Y and X. The same Q is
## passed to `test_stat_decorrelation_short()` and to `post_estimation()`.
##
## Arguments
##   method   : "decorrelation" (PCA), "trim", or "lmm"
##   X        : design matrix (required for decorrelation / trim)
##   q        : rank for the PCA projection (decorrelation only).
##              If NULL, chosen by the eigen-ratio rule on D[1:kmax].
##   rho      : trim cutoff fraction of singular values (trim only)
##   kmax     : maximum candidate rank for the eigen-ratio rule
##   kinship  : n x n GRM (required for lmm)
##   h2       : presumed heritability for the lmm filter, default 0.5
## ---------------------------------------------------------------------------
build_filter <- function(method = c("decorrelation", "trim", "lmm"),
                         X = NULL, q = NULL, rho = 0.3, kmax = 20,
                         kinship = NULL, h2 = 0.5, centralize = TRUE)
{
  method <- match.arg(tolower(method), c("decorrelation", "trim", "lmm"))
  
  if (method == "lmm") {
    if (is.null(kinship)) stop("LMM filter requires `kinship`.")
    n <- nrow(kinship)
    eig <- eigen(kinship, symmetric = TRUE)
    eig$values[eig$values < 0] <- 0
    sigma_g2 <- h2; sigma_e2 <- 1 - h2
    weights <- 1 / sqrt(sigma_g2 * eig$values + sigma_e2)
    return(eig$vectors %*% diag(weights) %*% t(eig$vectors))
  }
  
  if (is.null(X)) stop(sprintf("`%s` filter requires X.", method))
  if (isTRUE(centralize)) X <- scale(X, TRUE, FALSE)
  
  n <- nrow(X); p <- ncol(X)
  num_sv <- max(c(floor(rho * min(n, p)), kmax + 1,
                  if (is.null(q)) 0 else q))
  num_sv <- min(num_sv, min(n, p) - 1L)
  UDV <- RSpectra::svds(A = X, k = num_sv, nu = num_sv, nv = 0)
  U <- UDV$u; D <- UDV$d
  
  if (method == "trim") {
    r <- floor(rho * min(n, p))
    D_tilde <- sqrt(1 - D[r] / D[seq_len(r)])
    F_trim <- U[, seq_len(r), drop = FALSE] %*% diag(D_tilde)
    return(diag(n) - F_trim %*% t(F_trim))
  }
  
  ## decorrelation (PCA projection off the top-q left singular vectors)
  if (is.null(q)) {
    D1 <- D[seq_len(kmax)]; D2 <- D[seq_len(kmax) + 1L]
    q <- which.max(D1 / D2)
  }
  diag(n) - U[, seq_len(q), drop = FALSE] %*% t(U[, seq_len(q), drop = FALSE])
}

## ---------------------------------------------------------------------------
## nodewise_lasso_short() : nodewise residuals on the selected support
## ---------------------------------------------------------------------------
# nodewise_lasso_short = function(signal, X, nlambda = 100)
# {
#   p = ncol(X); n = nrow(X)
#   signal = sort(unique(as.integer(signal)))
#   signal = signal[signal >= 1 & signal <= p]
#   if (length(signal) == 0)
#     stop("nodewise_lasso_short requires at least one selected variable.")
# 
#   if (p == 1) {
#     gamma = Matrix::Matrix(0, nrow = 0, ncol = length(signal), sparse = TRUE)
#     denominator = colSums(X[, signal, drop = FALSE]^2)
#     return(list(gamma = gamma, denominator = denominator))
#   }
# 
#   gamma = Matrix::Matrix(0, nrow = p - 1, ncol = length(signal), sparse = TRUE)
#   denominator = rep(NA_real_, length(signal))
#   sequence = (nlambda:1) / nlambda * 3 * sqrt(log(max(p, 2)) / n)
# 
#   for (idx in seq_along(signal)) {
#     j = signal[idx]
#     X_minus_j = X[, -j, drop = FALSE]
#     X_j = X[, j]
#     lambda_j = sqrt(sum(X_j^2) / n) * sequence
# 
#     fit_j = glmnet::glmnet(
#       x = X_minus_j, y = X_j, family = "gaussian",
#       lambda = lambda_j, intercept = FALSE
#     )
# 
#     pred = as.matrix(predict(fit_j, newx = X_minus_j, s = lambda_j))
#     resid_mat = matrix(X_j, nrow = n, ncol = ncol(pred)) - pred
#     rss = colSums(resid_mat^2)
#     df  = Matrix::colSums(fit_j$beta != 0)
#     bic = log(pmax(rss / n, .Machine$double.eps)) + log(n) * as.numeric(df) / n
#     index = which.min(bic)
# 
#     gamma[, idx] = fit_j$beta[, index]
#     denominator[idx] = as.numeric(crossprod(resid_mat[, index], X_j))
#     if (!is.finite(denominator[idx]) ||
#         abs(denominator[idx]) <= .Machine$double.eps) {
#       stop("Non-finite or near-zero nodewise normalization for coordinate ", j, ".")
#     }
#   }
#   list(gamma = gamma, denominator = denominator)
# }


## ---------------------------------------------------------------------------
## multiple_testing() : BH-style threshold on a vector of test statistics
## ---------------------------------------------------------------------------
multiple_testing = function(statistics, signal = NULL, alpha = 0.1,
                            num_discrete = 100)
{
  p = length(statistics)
  t_p = sqrt(2 * log(p) - 2 * log(log(p)))
  t_grid = t_p * (1:num_discrete) / num_discrete

  fdp = numeric(num_discrete)
  for (i in seq_len(num_discrete)) {
    nd = sum(abs(statistics) >= t_grid[i])
    fdp[i] = 2 * p * (1 - pnorm(t_grid[i])) / max(nd, 1)
  }

  threshold = if (min(fdp) <= alpha) t_grid[min(which(fdp <= alpha))] else
              sqrt(2 * log(p))

  index_discovery = which(abs(statistics) >= threshold)
  num_discovery   = length(index_discovery)

  fdp_val = NULL; power_val = NULL
  index_Tdiscovery = NULL; num_Tdiscovery = NULL
  if (!is.null(signal)) {
    index_Tdiscovery = intersect(index_discovery, signal)
    num_Tdiscovery   = length(index_Tdiscovery)
    power_val = num_Tdiscovery / length(signal)
    fdp_val   = (num_discovery - num_Tdiscovery) / max(num_discovery, 1)
  }

  list(fdp = fdp_val, power = power_val, alpha = alpha,
       threshold = threshold, num_discovery = num_discovery,
       index_discovery = index_discovery,
       num_Tdiscovery = num_Tdiscovery, index_Tdiscovery = index_Tdiscovery)
}


## ---------------------------------------------------------------------------
## initial_lasso_fit() : cv.glmnet wrapper used inside the short pipeline
## ---------------------------------------------------------------------------
initial_lasso_fit = function(X, Y, nlambda = 100, nfolds = 10)
{
  n = nrow(X); p = ncol(X)
  Y = as.numeric(Y)
  sequence = (nlambda:1) / nlambda * 3
  lambda = sqrt(sum(Y^2) / n) * sqrt(log(max(p, 2)) / n) * sequence
  fit = glmnet::cv.glmnet(
    x = X, y = Y, lambda = lambda, intercept = FALSE,
    family = "gaussian", nfolds = min(nfolds, n)
  )
  beta = as.numeric(coef(fit, s = fit$lambda.min)[-1])
  list(beta = beta, lambda = fit$lambda.min, fit = fit)
}


## ---------------------------------------------------------------------------
## marginal_screen_signal() : BASIL-style fallback for tiny initial supports
## ---------------------------------------------------------------------------
marginal_screen_signal = function(Y, X, selected = integer(0),
                                  min_selected = 3, marginal_top = NULL)
{
  p = ncol(X)
  selected = sort(unique(as.integer(selected)))
  selected = selected[selected >= 1 & selected <= p]
  min_selected = min(max(1, min_selected), p)

  if (is.null(marginal_top)) {
    marginal_top = min(p, max(min_selected, 20))
  } else {
    marginal_top = min(p, max(min_selected, marginal_top))
  }

  if (length(selected) >= min_selected) {
    return(list(signal = selected, marginal_stats = NULL,
                used_fallback = FALSE))
  }

  ## Same screening principle as BASIL: rank by the marginal score
  ## |X_j' Y| / (||X_j|| ||Y||) and keep the top batch. This guarantees
  ## the nodewise step has at least min_selected candidates to work with.
  xty = abs(as.numeric(crossprod(X, as.numeric(Y))))
  x_norm = sqrt(colSums(X^2))
  y_norm = sqrt(sum(Y^2))
  marginal_stats = xty / (x_norm * y_norm + 1e-15)
  marginal_stats[!is.finite(marginal_stats)] = 0

  top_idx = order(marginal_stats, decreasing = TRUE)[seq_len(marginal_top)]
  signal  = sort(unique(c(selected, top_idx)))
  list(signal = signal, marginal_stats = marginal_stats,
       used_fallback = TRUE)
}


## ---------------------------------------------------------------------------
## estimate_filtered_sigma() : refit-after-selection style sigma estimator
## ---------------------------------------------------------------------------
estimate_filtered_sigma = function(Y, X, transform_matrix, beta_init,
                                   active_set = integer(0))
{
  n = nrow(X)
  rss = sum((Y - X %*% beta_init)^2)
  active_set = sort(unique(as.integer(active_set)))
  active_set = active_set[active_set >= 1 & active_set <= ncol(X)]

  if (length(active_set) == 0) {
    denom = sum(transform_matrix * transform_matrix)
  } else {
    svd_fit = svd(X[, active_set, drop = FALSE], nv = 0)
    rank_val = sum(svd_fit$d > 1e-10)
    if (rank_val == 0) {
      denom = sum(transform_matrix * transform_matrix)
    } else {
      projection = tcrossprod(svd_fit$u[, seq_len(rank_val), drop = FALSE])
      projection = (diag(1, n) - projection) %*% transform_matrix
      denom = sum(projection * projection)
    }
  }
  sqrt(rss / max(denom, .Machine$double.eps))
}


## ---------------------------------------------------------------------------
## test_stat_decorrelation_short()
## ---------------------------------------------------------------------------
##  Inputs
##    Y                : response, length n
##    X                : design, n x p
##    transform_matrix : n x n spectral filter Q (build via `build_filter()`)
##    centralize       : TRUE/FALSE
##    nlambda          : grid length for both initial and nodewise Lasso
##    min_selected     : minimum |hat S_0| needed; if smaller, marginal
##                       screening is invoked
##    marginal_top     : how many variables the marginal fallback adds
##    nfolds           : nfolds for cv.glmnet inside `initial_lasso_fit`
##
##  Returns a list with: beta_init, beta_debias, statistics_short, sigma_hat,
##  screen_idx, initial_idx, marginal_fallback, marginal_stats, time_short
##  (`statistics_short` is the canonical name; downstream `run_ddmm_variant`
##  uses it directly.)
## ---------------------------------------------------------------------------
test_stat_decorrelation_short = function(Y, X, transform_matrix,
                                         centralize = TRUE,
                                         nlambda = 100,
                                         min_selected = 3,
                                         marginal_top = NULL,
                                         nfolds = 10)
{
  n = nrow(X); p = ncol(X)
  if (missing(transform_matrix) || is.null(transform_matrix))
    stop("Please provide `transform_matrix` (an n x n spectral filter Q).")
  if (!all(dim(transform_matrix) == c(n, n)))
    stop("`transform_matrix` must have dimension nrow(X) by nrow(X).")

  Y = as.numeric(scale(Y, centralize, FALSE))
  X = scale(X, centralize, FALSE)
  transform_matrix = as.matrix(transform_matrix)

  F_Y = as.numeric(transform_matrix %*% Y)
  F_X = transform_matrix %*% X

  init_fit = initial_lasso_fit(F_X, F_Y, nlambda = nlambda, nfolds = nfolds)
  beta_init = init_fit$beta
  initial_signal = which(beta_init != 0)

  screen = marginal_screen_signal(
    Y = F_Y, X = F_X, selected = initial_signal,
    min_selected = min_selected, marginal_top = marginal_top
  )
  temp_signal = screen$signal

  sigma_hat = estimate_filtered_sigma(
    Y = F_Y, X = F_X, transform_matrix = transform_matrix,
    beta_init = beta_init, active_set = temp_signal
  )

  start_time = proc.time()[3]
  nodewise_fit = nodewise_lasso_short(temp_signal, F_X, nlambda)
  gamma_short = nodewise_fit$gamma
  denominator_short = nodewise_fit$denominator

  statistics_short  = rep(0, p)
  beta_debias_short = rep(0, p)
  residual = as.numeric(F_Y - F_X %*% beta_init)

  for (idx in seq_along(temp_signal)) {
    j = temp_signal[idx]
    z_j = if (p == 1) F_X[, j] else
          as.vector(F_X[, j] - F_X[, -j, drop = FALSE] %*% gamma_short[, idx])
    beta_debias_short[j] = beta_init[j] +
      sum(z_j * residual) / denominator_short[idx]
    z_j_norm = sqrt(sum(z_j^2))
    statistics_short[j] = beta_debias_short[j] * denominator_short[idx] /
      (sigma_hat * z_j_norm)
  }

  time_taken = proc.time()[3] - start_time

  list(beta_debias = beta_debias_short,
       beta_init = beta_init,
       statistics_short = statistics_short,
       sigma_hat = sigma_hat,
       time_short = as.numeric(time_taken),
       screen_idx = temp_signal,
       initial_idx = initial_signal,
       marginal_fallback = screen$used_fallback,
       marginal_stats = screen$marginal_stats,
       lambda = init_fit$lambda)
}


## ---------------------------------------------------------------------------
## post_estimation()
## ---------------------------------------------------------------------------
##  After multiple-testing has selected `index`, refit a small unpenalized
##  least-squares model on the filtered (Q Y, Q X[,index]) pair. The same
##  spectral filter Q used for inference is passed in here, which keeps the
##  noise model consistent across the two stages.
##
##  Inputs
##    Y, X             : original (untransformed) response and design
##    index            : indices of selected variables from multiple_testing
##    transform_matrix : the same Q passed to test_stat_decorrelation_short
##    centralize       : whether to centre Y and X (must match the inference
##                       call to keep the post estimator on the same scale)
##
##  Returns list(beta = ..., index = ...) where length(beta) = length(index).
## ---------------------------------------------------------------------------
post_estimation = function(Y, X, index, transform_matrix,
                           centralize = TRUE)
{
  if (missing(transform_matrix) || is.null(transform_matrix))
    stop("post_estimation requires the same `transform_matrix` Q ",
         "used in test_stat_decorrelation_short().")

  n = nrow(X)
  index = sort(unique(as.integer(index)))
  index = index[index >= 1 & index <= ncol(X)]
  if (length(index) == 0)
    return(list(beta = numeric(0), index = integer(0)))

  Y = as.numeric(scale(Y, centralize, FALSE))
  X = scale(X, centralize, FALSE)

  F_Y = as.numeric(transform_matrix %*% Y)
  F_X_post = (transform_matrix %*% X)[, index, drop = FALSE]

  ## small unpenalized OLS on the filtered residual, no intercept
  fit = lm.fit(x = F_X_post, y = F_Y)
  list(beta = as.numeric(fit$coefficients), index = index)
}


## ---------------------------------------------------------------------------
## run_ddmm_variant()
## ---------------------------------------------------------------------------
##  Glue function that takes the output of test_stat_decorrelation_short(),
##  applies the BH-style threshold, and runs post_estimation on the
##  discovered set. Used inside simulation drivers.
##
##  `post_est_args` should be `list(transform_matrix = Q, centralize = ...)`,
##  where Q is the same filter used to produce `ddmm_result`.
##
##  The function expects `causal_indices`, `true_beta`, `X_test`, `Y_test`,
##  `Y_train`, `X_train`, and the helpers `calculate_mspe`, `calculate_l1_error`
##  to exist in the calling environment, matching the previous code's
##  conventions.
## ---------------------------------------------------------------------------
run_ddmm_variant <- function(method_name, ddmm_result, post_est_args)
{
  cat(sprintf("  - Evaluating variant: %s\n", method_name))

  inference_res <- multiple_testing(
    statistics = ddmm_result$statistics_short,
    signal     = causal_indices,
    alpha      = 0.1
  )

  beta_debiased   <- ddmm_result$beta_debias
  p_full          <- length(beta_debiased)
  beta_sparse_est <- numeric(p_full)
  if (inference_res$num_discovery > 0) {
    beta_sparse_est[inference_res$index_discovery] <-
      beta_debiased[inference_res$index_discovery]
  }

  mspe_debiased <- calculate_mspe(beta_debiased, X_test, Y_test)
  l1_debiased   <- calculate_l1_error(beta_sparse_est, true_beta)

  variant_results <- list(
    Method   = sprintf("DDMM (%s)", method_name),
    MSPE     = mspe_debiased,
    FDR      = inference_res$fdp,
    Power    = inference_res$power,
    L1_Error = l1_debiased
  )

  if (inference_res$num_discovery > 0) {
    post_fit <- do.call(
      post_estimation,
      c(list(Y = Y_train, X = X_train,
             index = inference_res$index_discovery),
        post_est_args)
    )
    beta_post_est <- numeric(p_full)
    beta_post_est[inference_res$index_discovery] <- post_fit$beta

    mspe_post <- calculate_mspe(beta_post_est, X_test, Y_test)
    l1_post   <- calculate_l1_error(beta_post_est, true_beta)

    post_results <- list(
      Method   = sprintf("DDMM (%s, Post-Est)", method_name),
      MSPE     = mspe_post,
      FDR      = inference_res$fdp,
      Power    = inference_res$power,
      L1_Error = l1_post
    )
  } else {
    post_results <- list(
      Method   = sprintf("DDMM (%s, Post-Est)", method_name),
      MSPE     = NA,
      FDR      = 0,
      Power    = 0,
      L1_Error = sum(abs(true_beta))
    )
  }

  list(debiased = variant_results, post = post_results)
}




nodewise_lasso_short = function(signal, X, nlambda = 100)
{
  # X: matrix, size = n * p
  # nlambda: scalar, number of candidate parameters
  
  p = ncol(X)
  n = nrow(X)
  signal = as.integer(signal)
  if (length(signal) == 0 || anyNA(signal) || anyDuplicated(signal) ||
      any(signal < 1L | signal > p)) {
    stop("`signal` must contain unique, valid column indices.")
  }
  
  # nodewise lasso procedure
  Gamma = NULL
  denominator = NULL
  sequence = (nlambda:1) / nlambda * 3 * sqrt(log(p) / n)
  for (j in signal){
    lambda_j = sqrt(sum(X[, j]^2) / n) * sequence
    fit_j = glmnet::glmnet(x = X[, -j], y = X[, j], family = "gaussian", lambda = lambda_j, intercept = FALSE)
    Gamma = rbind(Gamma, fit_j$beta) # matrix, size = (p - 1) * nlambda
    temp = as.vector(crossprod(X[, j] - predict(fit_j, newx = X[, -j]), X[, j]))
    denominator = rbind(denominator, temp)  # vector, length = nlambda
  }
  rm(lambda_j, fit_j, temp)
  
  # tuning parameter selection
  # lambda: from large to small
  index = 1
  IC = Inf
  p_short = length(signal)
  for (i in 1:nlambda)
  {
    gamma = matrix(Gamma[, i], nrow = (p - 1), ncol = p_short)

    # Embed each nodewise coefficient vector in the original p coordinates,
    # then restrict it to the screened set. Triangular indexing is valid only
    # when signal is exactly 1:p in canonical order.
    gamma_full = matrix(0, nrow = p, ncol = p_short)
    for (h in seq_len(p_short)) gamma_full[-signal[h], h] = gamma[, h]
    Omega = -gamma_full[signal, , drop = FALSE]
    diag(Omega) = 1
    rm(gamma, gamma_full)
    
    Omega = t(t(Omega) / (denominator[, i] / n))
    
    # symmetrization
    Omega = Omega * (abs(Omega) <= abs(t(Omega))) + t(Omega) * (abs(Omega) > abs(t(Omega)))
    Omega = Matrix::Matrix(Omega, sparse = TRUE)
    
    # IC
    eig_min = RSpectra::eigs_sym(as(Omega, "matrix"), k = 1, which = "SA", opts = list(retvec = FALSE))$values
    if (length(eig_min) == 0) eig_min = -Inf # algorithm may not be converged
    
    if (eig_min <= 0)
    {
      next
    }else
    {
      IC_new = - determinant(as.matrix(Omega))$modulus[1] + log(n) / n * (sum(Omega != 0) - p_short) / 2 +
        sum(Omega * crossprod(X[,signal]) / n)
      if (IC_new <= IC)
      {
        IC = IC_new
        index = i
      }
    }
  }
  gamma = Matrix::Matrix(Gamma[, index], nrow = (p - 1), ncol = p_short, sparse = TRUE)
  denominator = denominator[, index]
  rm(Gamma, Omega)
  
  return(list(gamma = gamma, denominator = denominator))
}
## Backward-compatible alias for older simulation drivers that still call
## the function by its previous name.
run_sun_variant <- run_ddmm_variant
