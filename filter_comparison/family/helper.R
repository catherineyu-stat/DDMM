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
run_sun_variant <- function(method_name, sun_result, post_est_args) {
  
  cat(paste("  - Evaluating variant:", method_name, "\n"))
  
  inference_res <- multiple_testing(
    statistics = sun_result$statistics,
    signal = causal_indices,
    alpha = 0.1
  )
  
  if (inference_res$num_discovery > 0) {
    post_fit <- do.call(post_estimation, c(list(
      Y = y, X = X_std, index = inference_res$index_discovery
    ), post_est_args))
    
    beta_post_est <- rep(0, M)
    beta_post_est[inference_res$index_discovery] <- post_fit$beta
    
    l1_post <- calculate_l1_error(beta_post_est, beta_true)
    
    post_results <- list(
      Method = paste0("Sun (", method_name, ", Post-Est)"),
      FDR = inference_res$fdp,
      Power = inference_res$power,
      L1_Error = l1_post
    )
  } else {
    
    post_results <- list(
      Method = paste0("Sun (", method_name, ", Post-Est)"),
      FDR = 0, Power = 0, L1_Error = sum(abs(beta_true))
    )
  }
  
  return(list(post = post_results))
}
