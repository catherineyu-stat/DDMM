calculate_l1_error <- function(beta_est, beta_true) {
  if (length(beta_est) != length(beta_true)) {
    stop("Estimated and true beta vectors must have the same length.")
  }
  sum(abs(beta_est - beta_true))
}

calculate_fdr_power <- function(beta_est, beta_true) {
  causal_indices_true <- which(beta_true != 0)
  selected_indices_est <- which(beta_est != 0)
  
  if (length(selected_indices_est) == 0) {
    fdr <- 0
    power <- 0
  } else {
    fp <- length(setdiff(selected_indices_est, causal_indices_true))
    tp <- length(intersect(selected_indices_est, causal_indices_true))
    fdr <- fp / length(selected_indices_est)
    power <- tp / length(causal_indices_true)
  }
  return(list(FDR = fdr, Power = power))
}
