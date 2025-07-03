#' SDID Event Study Estimator
#'
#' @description
#' Implements Synthetic Difference-in-Differences (SDID) for event studies
#' with staggered adoption. Based on Arkhangelsky et al. (2021).
#'
#' This function requires the \code{synthdid} package to be installed.
#' Install it with: \code{devtools::install_github('synth-inference/synthdid')}
#'
#' @param Y Outcome variable (numeric vector)
#' @param G Unit identifier (factor or numeric)
#' @param T Time identifier (numeric)
#' @param D Treatment indicator (0/1)
#' @param effects Number of post-treatment effects to estimate (default: 0 = all)
#' @param placebo Number of pre-treatment placebo tests (default: NULL)
#' @param disag Show disaggregated cohort-level results (default: FALSE)
#' @param vce Variance calculation method: "bootstrap", "placebo", or "off" (default: "bootstrap")
#' @param brep Number of bootstrap replications (default: 50)
#' @param method Estimation method: "sdid", "did", or "sc" (default: "sdid")
#' @param covariates List with 'vars' (character vector) and 'method' ("optimized" or "projected")
#' @param boot_ci Use empirical bootstrap CI (default: FALSE)
#' @param combine Character vector of effect combinations to estimate (e.g., c("1 2", "3 4"))
#' @param verbose Show progress (default: TRUE)
#'
#' @return List containing:
#'   - estimates: Main treatment effects with standard errors
#'   - H_cohort: Cohort-specific results (if disag = TRUE)
#'   - H_combined: Combined effects (if combine specified)
#'   - vcov: Variance-covariance matrix (if requested)
#'   - boot_results: Bootstrap results (if requested)
#'
#' @import synthdid
#' @importFrom dplyr %>% group_by mutate ungroup
#' @export
sdid_event <- function(Y, G, T, D,
                       effects = 0,
                       placebo = NULL,
                       disag = FALSE,
                       vce = "bootstrap",
                       brep = 50,
                       method = "sdid",
                       covariates = NULL,
                       boot_ci = FALSE,
                       combine = NULL,
                       verbose = TRUE) {

  # Check that synthdid is installed
  if (!requireNamespace("synthdid", quietly = TRUE)) {
    stop("Package 'synthdid' is required but not installed.\n",
         "Please install it with:\n",
         "  devtools::install_github('synth-inference/synthdid')\n",
         "or\n",
         "  remotes::install_github('synth-inference/synthdid')")
  }

  # Input validation
  if (!method %in% c("sdid", "did", "sc")) {
    stop("Invalid method. Choose 'sdid', 'did', or 'sc'")
  }

  if (!vce %in% c("off", "bootstrap", "placebo")) {
    stop("Invalid vce. Choose 'off', 'bootstrap', or 'placebo'")
  }

  # Create data frame
  data <- data.frame(
    Y = Y,
    G = as.numeric(as.factor(G)),
    T = as.numeric(T),
    D = as.numeric(D)
  )

  # Remove any rows with missing Y values
  n_before <- nrow(data)
  data <- data[!is.na(data$Y), ]
  n_after <- nrow(data)
  if (n_before > n_after) {
    if (verbose) cat("Removed", n_before - n_after, "observations with missing outcome values.\n")
  }

  # Sort by G and T
  data <- data[order(data$G, data$T), ]

  # Create ever_treated indicator
  data <- data %>%
    dplyr::group_by(G) %>%
    dplyr::mutate(ever_treated = max(D)) %>%
    dplyr::ungroup()

  # Handle covariates if provided
  if (!is.null(covariates)) {
    if (verbose) cat("Adjusting for covariates...\n")

    cov_method <- ifelse(is.null(covariates$method), "optimized", covariates$method)

    # Residualize Y on covariates using control units
    control_data <- data[data$D == 0, ]

    # Create formula
    cov_formula <- as.formula(paste("Y ~", paste(covariates$vars, collapse = " + "),
                                    "+ factor(G) + factor(T)"))

    # Fit model on controls
    reg_model <- lm(cov_formula, data = control_data)

    # Get residuals for all data
    data$Y_res <- data$Y - predict(reg_model, newdata = data)
    data$Y <- data$Y_res
  }

  # Call core function
  if (verbose) cat(paste("Running", method, "estimation...\n"))

  core_results <- sdid_event_core(
    data = data,
    effects = effects,
    placebo = placebo,
    method = method,
    combine = combine
  )

  # Handle variance estimation
  if (vce != "off") {
    if (verbose) cat(paste("\nBootstrap replications (", brep, "), ", vce, " mode.\n", sep = ""))

    # Check for placebo vce validity
    if (vce == "placebo") {
      n_control <- length(unique(data$G[data$ever_treated == 0]))
      n_treated <- length(unique(data$G[data$ever_treated == 1]))

      if (n_control < n_treated) {
        stop("vce='placebo' requires more control units than treated units")
      }
    }

    # Run bootstrap
    boot_results <- bootstrap_sdid_event(
      data = data,
      core_results = core_results,
      effects = effects,
      placebo = placebo,
      method = method,
      vce = vce,
      brep = brep,
      combine = combine,
      boot_ci = boot_ci,
      verbose = verbose
    )

    # Add standard errors to results
    core_results$estimates <- boot_results$estimates_with_se

    if (!is.null(combine)) {
      core_results$H_combined <- boot_results$H_combined
    }

    if (!is.null(boot_results$vcov)) {
      core_results$vcov <- boot_results$vcov
    }

    core_results$boot_results <- boot_results$boot_matrix
  }

  # Print results
  if (verbose) {
    cat("\n", get_method_name(method), "\n", sep = "")
    print_results(core_results$estimates)

    if (!is.null(combine)) {
      cat("\nCombined ATTs\n")
      print_results(core_results$H_combined)
    }

    if (disag) {
      cat("\nDisaggregated ATTs - Cohort level\n")
      print(core_results$H_cohort)
    }
  }

  return(core_results)
}

#' Core SDID Event Study Function
#'
#' @param data Data frame with Y, G, T, D, ever_treated
#' @param effects Number of effects to estimate
#' @param placebo Number of placebo periods
#' @param method Estimation method
#' @param sampling Bootstrap sampling method (NULL, "bootstrap", or "placebo")
#' @param combine Effect combinations
#'
#' @return List with estimates and cohort results
#' @keywords internal
sdid_event_core <- function(data, effects, placebo, method,
                            sampling = NULL, combine = NULL) {

  # Handle sampling for bootstrap
  if (!is.null(sampling)) {
    if (sampling == "bootstrap") {
      # Cluster bootstrap by unit
      units <- unique(data$G)
      sampled_units <- sample(units, replace = TRUE)

      # Create new data with sampled units
      new_data <- list()
      for (i in seq_along(sampled_units)) {
        unit_data <- data[data$G == sampled_units[i], ]
        unit_data$G <- i  # Rename units to avoid duplicates
        new_data[[i]] <- unit_data
      }
      data <- do.call(rbind, new_data)

    } else if (sampling == "placebo") {
      # Placebo bootstrap: assign treatment to control units
      treated_pattern <- data[data$ever_treated == 1, c("G", "T", "D")]
      treated_pattern <- unique(treated_pattern[, c("T", "D")])

      # Keep only control units
      control_data <- data[data$ever_treated == 0, ]
      control_units <- unique(control_data$G)

      # Sample control units
      sampled_units <- sample(control_units, replace = TRUE)

      # Create new data
      new_data <- list()
      for (i in seq_along(sampled_units)) {
        unit_data <- control_data[control_data$G == sampled_units[i], ]
        unit_data$G <- i
        new_data[[i]] <- unit_data
      }
      data <- do.call(rbind, new_data)

      # Randomly assign treatment pattern to some units
      n_treat <- sum(treated_pattern$D[treated_pattern$T == min(treated_pattern$T[treated_pattern$D == 1])])
      treat_units <- sample(unique(data$G), n_treat)

      # Apply treatment pattern
      for (t_row in 1:nrow(treated_pattern)) {
        time_t <- treated_pattern$T[t_row]
        treat_t <- treated_pattern$D[t_row]
        data$D[data$T == time_t & data$G %in% treat_units] <- treat_t
      }

      # Update ever_treated
      data <- data %>%
        dplyr::group_by(G) %>%
        dplyr::mutate(ever_treated = max(D)) %>%
        dplyr::ungroup()
    }
  }

  # Calculate cohort (first treatment period)
  data <- data %>%
    dplyr::group_by(G) %>%
    dplyr::mutate(
      F_g_temp = D * T,
      C_temp = ifelse(any(D != 0), min(F_g_temp[D != 0]), 0),
      C = ifelse(ever_treated == 1, C_temp, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-F_g_temp, -C_temp)

  # Get time bounds
  T_max <- max(data$T)
  T_min <- min(data$T)

  # Calculate feasible effects
  data <- data %>%
    dplyr::mutate(T_g = ifelse(ever_treated == 1, T_max - C + 1, NA))

  L_g <- max(data$T_g, na.rm = TRUE)

  # Calculate feasible placebo periods
  L_pl_g <- max(data$C[data$ever_treated == 1]) - T_min

  # Set actual number of effects and placebos
  if (effects == 0 || effects > L_g) {
    effects <- L_g
  }

  if (!is.null(placebo)) {
    if (placebo == "all") {
      placebo <- L_pl_g
    } else if (placebo > L_pl_g) {
      placebo <- L_pl_g
    }
  } else {
    placebo <- 0
  }

  # Get unique cohorts
  cohorts <- sort(unique(data$C[data$C > 0]))
  n_cohorts <- length(cohorts)

  # Initialize results matrices
  res <- matrix(NA, nrow = 1 + L_g + placebo, ncol = n_cohorts)
  rownames(res) <- c("ATT_c",
                     paste0("Effect_c", 1:L_g),
                     if (placebo > 0) paste0("Placebo_c", 1:placebo) else NULL)
  colnames(res) <- paste0("c=", cohorts)

  c_weight <- numeric(n_cohorts)
  t_weight <- numeric(n_cohorts)

  # Loop over cohorts
  for (j in 1:n_cohorts) {
    cohort_c <- cohorts[j]

    # Get cohort data
    cohort_data <- data[data$C %in% c(0, cohort_c), ]
    N_Post_c <- T_max - cohort_c + 1
    N_Pre_c <- cohort_c - T_min
    N_Tr_c <- length(unique(cohort_data$G[cohort_data$C == cohort_c]))
    N_c <- length(unique(cohort_data$G))

    c_weight[j] <- N_Tr_c
    t_weight[j] <- N_Post_c * N_Tr_c

    # Estimate main ATT
    att_main <- estimate_att_cohort(
      data = cohort_data,
      cohort = cohort_c,
      periods = NULL,
      method = method
    )
    res[1, j] <- att_main

    # Estimate dynamic effects
    for (l in 1:N_Post_c) {
      if (l <= L_g) {
        att_l <- estimate_att_cohort(
          data = cohort_data,
          cohort = cohort_c,
          periods = cohort_c - 1 + l,
          method = method
        )
        res[1 + l, j] <- att_l
      }
    }

    # Estimate placebo effects
    if (placebo > 0 && N_Pre_c > 0) {
      for (l in 1:min(placebo, N_Pre_c)) {
        att_pl <- estimate_att_cohort_placebo(
          data = cohort_data,
          cohort = cohort_c,
          pl = l,
          method = method
        )
        res[1 + L_g + l, j] <- att_pl
      }
    }
  }

  # Aggregate results
  H <- aggregate_att(
    res = res,
    t_weight = t_weight,
    c_weight = c_weight,
    effects = effects,
    placebo = placebo,
    L_g = L_g
  )

  # Handle combinations
  H_combined <- NULL
  if (!is.null(combine)) {
    H_combined <- calculate_combinations(
      res = res,
      H = H,
      combine = combine,
      c_weight = c_weight
    )
  }

  return(list(
    estimates = H,
    H_cohort = if (exists("res")) res else NULL,
    H_combined = H_combined,
    cohorts = cohorts,
    c_weight = c_weight,
    t_weight = t_weight,
    effects = effects,
    placebo = placebo
  ))
}

#' Helper function to get method name
#' @keywords internal
get_method_name <- function(method) {
  switch(method,
         "sdid" = "Synthetic Difference-in-differences",
         "did" = "Difference-in-differences",
         "sc" = "Synthetic Control")
}

#' Helper function to print results table
#' @keywords internal
print_results <- function(H) {
  # Format the results nicely
  results_df <- as.data.frame(H)
  results_df$Estimate <- sprintf("%.5f", results_df$Estimate)
  if ("SE" %in% names(results_df)) {
    results_df$SE <- sprintf("%.5f", results_df$SE)
  }
  if ("LB_CI" %in% names(results_df)) {
    results_df$LB_CI <- sprintf("%.5f", results_df$LB_CI)
    results_df$UB_CI <- sprintf("%.5f", results_df$UB_CI)
  }
  print(results_df, row.names = TRUE)
}
