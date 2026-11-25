## ============================================================================
## inference_short_fixed.R -- DDMM core with the cluster fixes applied
##
## This file sources the unchanged core `inference_short_refactored.R` and then
## overrides a single function, `test_stat_decorrelation_short()`. Everything
## else (build_filter, multiple_testing, post_estimation, run_ddmm_variant,
## marginal_screen_signal, estimate_filtered_sigma) is used as it is.
##
## Status of the three fixes in the single-variant pipeline
##
##   FIX 1  APPLIED. The nodewise lambda is chosen per variant by BIC.
##          `nodewise_lasso_short()` chose one lambda jointly for the whole
##          support, which is too large for tightly linked pairs: the nodewise
##          residual of a tag keeps a component along its causal neighbour and
##          the tag inherits the signal. `nodewise_direct()` below runs the same
##          regression as `nodewise_z()` in ddmm_group/fast_nodewise.R (30-point
##          path, BIC per variant) but works on the filtered design itself, so it
##          does not need the p x p eigendecomposition that fast_nodewise.R uses.
##          The lambda path is rescaled so that both return the same fit.
##
##   FIX 2  NOT NEEDED. It concerns the group tests. In the single-variant test
##          the pilot is already the CV lasso that sigma is estimated from.
##
##   FIX 3  ALREADY IN THE CORE. `estimate_filtered_sigma()` divides the RSS by
##          ||(I - P_S) Q||_F^2, the degrees-of-freedom corrected form.
##
## Usage: source this file instead of inference_short_refactored.R.
## ============================================================================

source("inference_short_refactored.R")

## keep a handle on the original statistic so a driver can run both versions
## on the same data and attribute any difference to the nodewise step alone
test_stat_decorrelation_short_original <- test_stat_decorrelation_short


## ---------------------------------------------------------------------------
## ddmm_original_variant(): the ORIGINAL DDMM on the same filtered data, used
## only for the paired before/after comparison. The caller's RNG state is
## saved before and restored after, so (a) the original and the fixed DDMM see
## the same cross-validation folds and hence the same initial lasso, and (b)
## the random numbers drawn by everything downstream are exactly those of a
## run without this comparison.
## ---------------------------------------------------------------------------
ddmm_original_variant <- function(name, Y, X, Q) {
  seed <- get(".Random.seed", envir = globalenv())
  st <- test_stat_decorrelation_short_original(
    Y = Y, X = X, transform_matrix = Q, centralize = TRUE)
  out <- run_ddmm_variant(name, st, list(transform_matrix = Q, centralize = TRUE))
  assign(".Random.seed", seed, envir = globalenv())
  out$post
}


## ---------------------------------------------------------------------------
## nodewise_direct(): per-variant nodewise lasso, lambda by BIC  (FIX 1)
##
##   S   : indices of the variants to debias
##   FX  : filtered design Q X, n x p
##   r   : filtered residual Q Y - Q X beta_init
##
## Returns, for every j in S, the three inner products the debiased statistic
## needs: zr = z_j' r, zx = z_j' x_j, zz = z_j' z_j.
##
## Equivalence with fast_nodewise.R::nodewise_z(). That function fits glmnet on
## R, a p x p square root of G = FX'FX, with lambda_max = max|G_{-j,j}| / p. With
## glmnet's standardisation the two problems have the same solution path when
## the lambda sequence is multiplied by sqrt(p / n). BIC is computed from the
## same RSS and the same df, so the selected fit is identical.
## ---------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------
## test_stat_decorrelation_short(): same inputs, outputs and statistic as the
## original; only the nodewise step is replaced.
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

  ## FIX 3 is already here: RSS / ||(I - P_S) Q||_F^2
  sigma_hat = estimate_filtered_sigma(
    Y = F_Y, X = F_X, transform_matrix = transform_matrix,
    beta_init = beta_init, active_set = temp_signal
  )

  start_time = proc.time()[3]
  residual = as.numeric(F_Y - F_X %*% beta_init)

  ## FIX 1: nodewise lambda chosen per variant by BIC
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
