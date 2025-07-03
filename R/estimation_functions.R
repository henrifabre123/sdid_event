#' Estimate ATT for a cohort
#'
#' @param data Data frame with cohort data
#' @param cohort Cohort identifier
#' @param periods Specific periods to include (NULL for all post-treatment)
#' @param method Estimation method
#'
#' @return ATT estimate
#' @importFrom synthdid synthdid_estimate did_estimate sc_estimate
#' @keywords internal
estimate_att_cohort <- function(data, cohort, periods = NULL, method = "sdid") {
  # Create outcome matrix
  Y_wide <- data %>%
    dplyr::select(Y, G, T) %>%
    tidyr::pivot_wider(names_from = T, values_from = Y) %>%  # PAS de values_fill
    dplyr::select(-G) %>%
    as.matrix()

  # Create treatment indicator matrix
  W_wide <- data %>%
    dplyr::select(D, G, T) %>%
    tidyr::pivot_wider(names_from = T, values_from = D, values_fill = 0) %>%
    dplyr::select(-G) %>%
    as.matrix()

  # Get dimensions
  units <- sort(unique(data$G))  # IMPORTANT: Trier comme Stata
  times <- sort(unique(data$T))
  N <- length(units)
  TT <- length(times)

  # Find control and treated units
  control_units <- sort(unique(data$G[data$C == 0]))  # Trier
  treated_units <- sort(unique(data$G[data$C == cohort]))  # Trier
  N0 <- length(control_units)
  N1 <- length(treated_units)

  # Find pre and post periods - CORRIGER le calcul de T0
  pre_treatment_times <- times[times < cohort]
  post_treatment_times <- times[times >= cohort]
  T0 <- length(pre_treatment_times)
  T1 <- length(post_treatment_times)

  # Filter periods if specified
  if (!is.null(periods)) {
    period_cols <- which(times %in% periods)
    pre_cols <- which(times < cohort)
    keep_cols <- c(pre_cols, period_cols)
    Y_wide <- Y_wide[, keep_cols]
    W_wide <- W_wide[, keep_cols]
    times <- times[keep_cols]
    TT <- length(times)
    T1 <- length(which(times >= cohort))
  }

  # Reorder units: controls first, then treated
  unit_order <- c(which(units %in% control_units), which(units %in% treated_units))
  Y_wide <- Y_wide[unit_order, ]
  W_wide <- W_wide[unit_order, ]

  # Check dimensions
  if (nrow(Y_wide) != N || ncol(Y_wide) != TT) {
    stop(sprintf("Matrix dimension mismatch. Expected %d x %d, got %d x %d",
                 N, TT, nrow(Y_wide), ncol(Y_wide)))
  }

  # Check for sufficient pre/post periods and units
  if (T0 <= 0) {
    stop("No pre-treatment periods available for this cohort")
  }
  if (N0 <= 0) {
    stop("No control units available for this cohort")
  }

  # Vérifier si synthdid peut gérer les NA
  # Si synthdid R ne peut pas, on devra créer un wrapper

  if (method == "sdid") {
    # Si synthdid ne gère pas les NA, utiliser cette approche :
    if (any(is.na(Y_wide))) {
      # Créer un wrapper qui filtre comme Stata
      tau_hat <- synthdid_with_na_handling(Y_wide, N0, T0)
    } else {
      tau_hat <- synthdid::synthdid_estimate(Y_wide, N0, T0)
    }
    return(as.numeric(tau_hat))

  } else if (method == "did") {
    if (any(is.na(Y_wide))) {
      did_hat <- did_with_na_handling(Y_wide, N0, T0)
    } else {
      did_hat <- synthdid::did_estimate(Y_wide, N0, T0)
    }
    return(as.numeric(did_hat))

  } else if (method == "sc") {
    if (any(is.na(Y_wide))) {
      sc_hat <- sc_with_na_handling(Y_wide, N0, T0)
    } else {
      sc_hat <- synthdid::sc_estimate(Y_wide, N0, T0)
    }
    return(as.numeric(sc_hat))
  }
}

#' Estimate placebo ATT for a cohort
#'
#' @param data Data frame with cohort data
#' @param cohort Cohort identifier
#' @param pl Placebo period (periods before treatment)
#' @param method Estimation method
#'
#' @return Placebo ATT estimate
#' @keywords internal
estimate_att_cohort_placebo <- function(data, cohort, pl, method = "sdid") {

  # Keep only pre-treatment data
  pre_data <- data[data$T < cohort, ]

  # Check if we have enough pre-treatment periods
  if (length(unique(pre_data$T)) <= pl) {
    return(NA)  # Not enough pre-treatment periods for this placebo
  }

  # Create pseudo-treatment at cohort - pl
  pseudo_treat_time <- cohort - pl

  # Check if pseudo treatment time is valid
  if (pseudo_treat_time <= min(pre_data$T)) {
    return(NA)  # Pseudo treatment time too early
  }

  # Update the cohort indicator for the pseudo treatment
  pre_data$C_original <- pre_data$C
  pre_data$C <- ifelse(pre_data$C == cohort, pseudo_treat_time, 0)

  # Create pseudo treatment indicator
  pre_data$D_original <- pre_data$D
  pre_data$D <- ifelse(pre_data$T >= pseudo_treat_time & pre_data$C_original == cohort, 1, 0)

  # Estimate ATT for pseudo treatment
  att_placebo <- tryCatch({
    estimate_att_cohort(
      data = pre_data,
      cohort = pseudo_treat_time,
      periods = NULL,
      method = method
    )
  }, error = function(e) {
    warning(paste("Placebo estimation failed for period", pl, ":", e$message))
    return(NA)
  })

  return(att_placebo)
}

#' Aggregate ATT results across cohorts
#'
#' @param res Matrix of cohort-specific results
#' @param t_weight Time weights for aggregation
#' @param c_weight Cohort weights
#' @param effects Number of effects to report
#' @param placebo Number of placebo periods to report
#' @param L_g Maximum number of effects
#'
#' @return Aggregated results matrix
#' @keywords internal
aggregate_att <- function(res, t_weight, c_weight, effects, placebo, L_g) {

  # Prepare weights matrix
  n_cohorts <- ncol(res)
  W <- rbind(t_weight, c_weight, res)

  # Initialize results
  n_rows <- 1 + effects + placebo
  H <- matrix(NA, nrow = n_rows, ncol = 3)
  colnames(H) <- c("Estimate", "SE", "Switchers")

  row_names <- "ATT"
  if (effects > 0) {
    row_names <- c(row_names, paste0("Effect_", 1:effects))
  }
  if (placebo > 0) {
    row_names <- c(row_names, paste0("Placebo_", 1:placebo))
  }
  rownames(H) <- row_names

  # Calculate aggregated ATT
  valid_cols <- which(!is.na(W[3, ]))
  if (length(valid_cols) > 0) {
    weights <- W[1, valid_cols] / sum(W[1, valid_cols])
    H[1, 1] <- sum(weights * W[3, valid_cols])
    H[1, 3] <- sum(c_weight[valid_cols])
  }

  # Calculate dynamic effects
  for (j in 1:effects) {
    row_idx <- 3 + j
    valid_cols <- which(!is.na(W[row_idx, ]))
    if (length(valid_cols) > 0) {
      weights <- W[2, valid_cols] / sum(W[2, valid_cols])
      H[j + 1, 1] <- sum(weights * W[row_idx, valid_cols])
      H[j + 1, 3] <- sum(c_weight[valid_cols])
    }
  }

  # Calculate placebo effects
  if (placebo > 0) {
    for (j in 1:placebo) {
      row_idx <- 3 + L_g + j
      if (row_idx <= nrow(W)) {
        valid_cols <- which(!is.na(W[row_idx, ]))
        if (length(valid_cols) > 0) {
          weights <- W[2, valid_cols] / sum(W[2, valid_cols])
          H[1 + effects + j, 1] <- sum(weights * W[row_idx, valid_cols])
          H[1 + effects + j, 3] <- sum(c_weight[valid_cols])
        }
      }
    }
  }

  return(H)
}

#' Calculate combined effects
#'
#' @param res Matrix of cohort-specific results
#' @param H Main results matrix
#' @param combine Character vector of combinations
#' @param c_weight Cohort weights
#'
#' @return Combined effects matrix
#' @keywords internal
calculate_combinations <- function(res, H, combine, c_weight) {

  n_combine <- length(combine)
  H_cb <- matrix(NA, nrow = n_combine, ncol = 3)
  colnames(H_cb) <- c("Estimate", "SE", "Switchers")
  rownames(H_cb) <- paste0("Cmb_Effect_", 1:n_combine)

  for (k in 1:n_combine) {
    # Parse combination string
    effects_to_combine <- as.numeric(strsplit(combine[k], " ")[[1]])

    # Calculate weighted average
    total_weight <- 0
    weighted_sum <- 0

    for (s in effects_to_combine) {
      if (s > 0 && s < nrow(H)) {
        weight <- H[s + 1, 3]  # Number of switchers
        if (!is.na(weight) && !is.na(H[s + 1, 1])) {
          weighted_sum <- weighted_sum + H[s + 1, 1] * weight
          total_weight <- total_weight + weight
        }
      }
    }

    if (total_weight > 0) {
      H_cb[k, 1] <- weighted_sum / total_weight
      H_cb[k, 3] <- total_weight
    }
  }

  return(H_cb)
}
