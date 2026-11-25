# Shared DDMM inference for admixture and family simulations.
# Apply Q to X and y; fit the pilot Lasso; choose each SNP's nodewise
# penalty by BIC; form debiased statistics; threshold and refit.

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

nodewise_direct <- function(S, FX, r, nlambda = 30L, eps = 0.01,
                            cores = getOption("ddmm.cores", 1L)) {
  n <- nrow(FX); p <- ncol(FX)
  one <- function(j) {
    xj  <- FX[, j]
    Xmj <- FX[, -j, drop = FALSE]
    Gmj <- as.numeric(crossprod(Xmj, xj))
    lam_max <- max(abs(Gmj)) / p
    if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1e-3
    lam <- lam_max * exp(seq(0, log(eps), length.out = nlambda)) * sqrt(p / n)
    fit <- glmnet::glmnet(Xmj, xj, lambda = lam, intercept = FALSE)
    Z   <- xj - as.matrix(Xmj %*% fit$beta)          # n x nlambda residuals
    rss <- colSums(Z^2)
    b   <- which.min(n * log(pmax(rss, .Machine$double.eps) / n) + log(n) * fit$df)
    z   <- Z[, b]
    c(zr = sum(z * r), zx = sum(z * xj), zz = sum(z^2))
  }
  ## glmnet is deterministic, so running the variants in parallel changes
  ## nothing but the wall time. The caller's RNG stream is saved and restored
  ## so that forking cannot shift the random numbers used downstream.
  if (cores > 1L && length(S) > 1L) {
    had_seed <- exists(".Random.seed", envir = globalenv())
    if (had_seed) seed <- get(".Random.seed", envir = globalenv())
    out <- parallel::mclapply(S, one, mc.cores = cores, mc.set.seed = FALSE)
    if (had_seed) assign(".Random.seed", seed, envir = globalenv())
  } else {
    out <- lapply(S, one)
  }
  out <- do.call(rbind, out)
  list(zr = out[, "zr"], zx = out[, "zx"], zz = out[, "zz"])
}

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

  ## Residual variance: RSS / ||(I - P_S) Q||_F^2
  sigma_hat = estimate_filtered_sigma(
    Y = F_Y, X = F_X, transform_matrix = transform_matrix,
    beta_init = beta_init, active_set = temp_signal
  )

  start_time = proc.time()[3]
  residual = as.numeric(F_Y - F_X %*% beta_init)

  ## Nodewise lambda chosen per variant by BIC
  nw = nodewise_direct(temp_signal, F_X, residual)

  statistics_short  = rep(0, p)
  beta_debias_short = rep(0, p)
  beta_debias_short[temp_signal] = beta_init[temp_signal] + nw$zr / nw$zx
  ## exact studentisation, as in the original: beta_bar * z'x / (sigma ||z||)
  statistics_short[temp_signal] = beta_debias_short[temp_signal] * nw$zx /
    (sigma_hat * sqrt(nw$zz))

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

