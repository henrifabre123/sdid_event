#' Bootstrap for SDID Event Study
#'
#' @param data Original data frame
#' @param core_results Results from core estimation
#' @param effects Number of effects
#' @param placebo Number of placebo periods
#' @param method Estimation method
#' @param vce Bootstrap type
#' @param brep Number of replications
#' @param combine Combinations
#' @param boot_ci Use empirical CI
#' @param verbose Show progress
#'
#' @return List with bootstrap results
#' @keywords internal
bootstrap_sdid_event <- function(data, core_results, effects, placebo, 
                                 method, vce, brep, combine, boot_ci, 
                                 verbose = TRUE) {
  
  # Initialize bootstrap results matrix
  n_estimates <- nrow(core_results$estimates)
  boot_matrix <- matrix(NA, nrow = brep, ncol = n_estimates)
  
  # Progress bar setup
  if (verbose) {
    cat("|0%", paste0(rep("-", 41), collapse = ""), "100%|\n|", sep = "")
    progress_step <- ceiling(brep / 50)
  }
  
  # Run bootstrap replications
  successful_reps <- 0
  failed_reps <- 0
  rep <- 1
  
  while (successful_reps < brep) {
    tryCatch({
      # Run core with sampling
      boot_results <- sdid_event_core(
        data = data,
        effects = effects, 
        placebo = placebo,
        method = method,
        sampling = vce,
        combine = combine
      )
      
      # Store results
      boot_matrix[successful_reps + 1, ] <- boot_results$estimates[, "Estimate"]
      successful_reps <- successful_reps + 1
      
      # Update progress
      if (verbose && successful_reps %% progress_step == 0) {
        cat(".")
      }
      
    }, error = function(e) {
      failed_reps <<- failed_reps + 1
    })
    
    rep <- rep + 1
    
    # Safety check to avoid infinite loop
    if (rep > brep * 3) {
      warning(paste("Too many failed bootstrap replications. Stopping with", 
                    successful_reps, "successful replications."))
      break
    }
  }
  
  if (verbose) cat("|\n")
  
  # Trim bootstrap matrix if needed
  if (successful_reps < brep) {
    boot_matrix <- boot_matrix[1:successful_reps, ]
    warning(paste("Only", successful_reps, "bootstrap replications succeeded."))
  }
  
  # Calculate standard errors and confidence intervals
  estimates_with_se <- calculate_se(
    estimates = core_results$estimates,
    boot_matrix = boot_matrix,
    boot_ci = boot_ci
  )
  
  # Handle combinations if specified
  H_combined <- NULL
  if (!is.null(combine)) {
    # Calculate bootstrap results for combinations
    n_combine <- length(combine)
    boot_matrix_cb <- matrix(NA, nrow = successful_reps, ncol = n_combine)
    
    for (b in 1:successful_reps) {
      for (k in 1:n_combine) {
        effects_to_combine <- as.numeric(strsplit(combine[k], " ")[[1]])
        
        weighted_sum <- 0
        total_weight <- 0
        
        for (s in effects_to_combine) {
          if (s > 0 && s < n_estimates) {
            weight <- core_results$estimates[s + 1, "Switchers"]
            if (!is.na(weight) && !is.na(boot_matrix[b, s + 1])) {
              weighted_sum <- weighted_sum + boot_matrix[b, s + 1] * weight
              total_weight <- total_weight + weight
            }
          }
        }
        
        if (total_weight > 0) {
          boot_matrix_cb[b, k] <- weighted_sum / total_weight
        }
      }
    }
    
    # Calculate SE for combinations
    H_combined <- calculate_se_combinations(
      H_cb = core_results$H_combined,
      boot_matrix_cb = boot_matrix_cb,
      boot_ci = boot_ci
    )
  }
  
  # Calculate variance-covariance matrix if requested
  vcov <- NULL
  if (effects > 0) {
    effect_cols <- 2:(effects + 1)
    vcov <- cov(boot_matrix[, effect_cols, drop = FALSE])
    rownames(vcov) <- paste0("Effect_", 1:effects)
    colnames(vcov) <- paste0("Effect_", 1:effects)
  }
  
  # Report failed replications
  if (failed_reps > 0 && verbose) {
    cat("\nWARNING: Restarted", failed_reps, 
        "bootstrap run(s) with no treated or control groups.\n")
  }
  
  return(list(
    estimates_with_se = estimates_with_se,
    H_combined = H_combined,
    boot_matrix = boot_matrix,
    vcov = vcov
  ))
}

#' Calculate standard errors from bootstrap
#'
#' @param estimates Original estimates matrix
#' @param boot_matrix Bootstrap results matrix
#' @param boot_ci Use empirical CI
#'
#' @return Matrix with estimates, SE, and CI
#' @keywords internal
calculate_se <- function(estimates, boot_matrix, boot_ci = FALSE) {
  
  n_estimates <- nrow(estimates)
  H_SE <- matrix(NA, nrow = n_estimates, ncol = 5)
  colnames(H_SE) <- c("Estimate", "SE", "LB_CI", "UB_CI", "Switchers")
  rownames(H_SE) <- rownames(estimates)
  
  # Copy estimates and switchers
  H_SE[, "Estimate"] <- estimates[, "Estimate"]
  H_SE[, "Switchers"] <- estimates[, "Switchers"]
  
  # Calculate SE and CI for each estimate
  for (j in 1:n_estimates) {
    boot_values <- boot_matrix[, j]
    boot_values <- boot_values[!is.na(boot_values)]
    
    if (length(boot_values) > 1) {
      # Standard error
      H_SE[j, "SE"] <- sd(boot_values)
      
      if (boot_ci) {
        # Empirical confidence intervals
        ci <- quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE)
        H_SE[j, "LB_CI"] <- ci[1]
        H_SE[j, "UB_CI"] <- ci[2]
      } else {
        # Normal approximation CI
        H_SE[j, "LB_CI"] <- H_SE[j, "Estimate"] - 1.96 * H_SE[j, "SE"]
        H_SE[j, "UB_CI"] <- H_SE[j, "Estimate"] + 1.96 * H_SE[j, "SE"]
      }
    }
  }
  
  return(H_SE)
}

#' Calculate standard errors for combinations
#'
#' @param H_cb Combination estimates matrix
#' @param boot_matrix_cb Bootstrap results for combinations
#' @param boot_ci Use empirical CI
#'
#' @return Matrix with estimates, SE, and CI
#' @keywords internal
calculate_se_combinations <- function(H_cb, boot_matrix_cb, boot_ci = FALSE) {
  
  if (is.null(H_cb)) return(NULL)
  
  n_combine <- nrow(H_cb)
  H_cb_SE <- matrix(NA, nrow = n_combine, ncol = 5)
  colnames(H_cb_SE) <- c("Estimate", "SE", "LB_CI", "UB_CI", "Switchers") 
  rownames(H_cb_SE) <- rownames(H_cb)
  
  # Copy estimates and switchers
  H_cb_SE[, "Estimate"] <- H_cb[, "Estimate"]
  H_cb_SE[, "Switchers"] <- H_cb[, "Switchers"]
  
  # Calculate SE and CI
  for (j in 1:n_combine) {
    boot_values <- boot_matrix_cb[, j]
    boot_values <- boot_values[!is.na(boot_values)]
    
    if (length(boot_values) > 1) {
      # Standard error
      H_cb_SE[j, "SE"] <- sd(boot_values)
      
      if (boot_ci) {
        # Empirical CI
        ci <- quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE)
        H_cb_SE[j, "LB_CI"] <- ci[1]
        H_cb_SE[j, "UB_CI"] <- ci[2]
      } else {
        # Normal approximation CI
        H_cb_SE[j, "LB_CI"] <- H_cb_SE[j, "Estimate"] - 1.96 * H_cb_SE[j, "SE"]
        H_cb_SE[j, "UB_CI"] <- H_cb_SE[j, "Estimate"] + 1.96 * H_cb_SE[j, "SE"]
      }
    }
  }
  
  return(H_cb_SE)
}