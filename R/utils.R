#' Prepare data for SDID Event Study
#'
#' @description
#' Converts panel data into the format required for SDID event study estimation.
#' Handles unbalanced panels and creates necessary variables.
#'
#' @param data Data frame with panel data
#' @param unit_var Name of unit variable
#' @param time_var Name of time variable  
#' @param outcome_var Name of outcome variable
#' @param treatment_var Name of treatment indicator variable
#'
#' @return Data frame formatted for sdid_event
#' @export
prepare_sdid_data <- function(data, unit_var, time_var, 
                              outcome_var, treatment_var) {
  
  # Rename variables to standard names
  result <- data.frame(
    G = data[[unit_var]],
    T = data[[time_var]],
    Y = data[[outcome_var]],
    D = data[[treatment_var]]
  )
  
  # Ensure numeric types
  result$G <- as.numeric(as.factor(result$G))
  result$T <- as.numeric(result$T)
  result$Y <- as.numeric(result$Y)
  result$D <- as.numeric(result$D)
  
  # Check for validity
  if (any(result$D != 0 & result$D != 1)) {
    stop("Treatment variable must be binary (0/1)")
  }
  
  # Sort by unit and time
  result <- result[order(result$G, result$T), ]
  
  return(result)
}

#' Plot SDID Event Study Results
#'
#' @description
#' Creates an event study plot showing dynamic treatment effects over time.
#'
#' @param results Results object from sdid_event
#' @param type Type of plot: "effects" (default) or "coefficients"
#' @param include_placebo Include placebo tests in plot (default: TRUE)
#' @param confidence_level Confidence level for intervals (default: 0.95)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return ggplot2 object
#' @export
plot_sdid_event <- function(results, 
                            type = "effects",
                            include_placebo = TRUE,
                            confidence_level = 0.95,
                            ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  # Extract estimates
  estimates <- as.data.frame(results$estimates)
  estimates$period <- 1:nrow(estimates) - 1  # ATT is period 0
  
  # Separate effects and placebos
  n_effects <- results$effects
  n_placebo <- results$placebo
  
  # Label periods
  estimates$period_type <- "ATT"
  if (n_effects > 0) {
    effect_rows <- 2:(n_effects + 1)
    estimates$period_type[effect_rows] <- "Effect"
    estimates$period[effect_rows] <- 0:(n_effects - 1)
  }
  
  if (n_placebo > 0 && include_placebo) {
    placebo_rows <- (n_effects + 2):(n_effects + n_placebo + 1)
    estimates$period_type[placebo_rows] <- "Placebo"
    estimates$period[placebo_rows] <- -n_placebo:-1
  }
  
  # Filter data for plot
  plot_data <- estimates[estimates$period_type != "ATT", ]
  if (!include_placebo) {
    plot_data <- plot_data[plot_data$period_type != "Placebo", ]
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, 
                       ggplot2::aes(x = period, y = Estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "solid", alpha = 0.3) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Period Relative to Treatment",
      y = "Treatment Effect",
      title = "SDID Event Study"
    )
  
  # Add confidence intervals if available
  if ("SE" %in% names(plot_data) && !is.na(plot_data$SE[1])) {
    if ("LB_CI" %in% names(plot_data)) {
      # Use provided CI
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = LB_CI, ymax = UB_CI),
        width = 0.2
      )
    } else {
      # Calculate CI from SE
      z_score <- qnorm(1 - (1 - confidence_level) / 2)
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Estimate - z_score * SE,
          ymax = Estimate + z_score * SE
        ),
        width = 0.2
      )
    }
  }
  
  return(p)
}

#' Summary method for sdid_event results
#'
#' @param object Results from sdid_event
#' @param ... Additional arguments (ignored)
#'
#' @return Printed summary
#' @export
summary.sdid_event <- function(object, ...) {
  cat("\nSDID Event Study Results\n")
  cat("========================\n\n")
  
  # Main estimates
  cat("Treatment Effects:\n")
  print(object$estimates)
  
  # Combined effects if available
  if (!is.null(object$H_combined)) {
    cat("\nCombined Effects:\n")
    print(object$H_combined)
  }
  
  # Method info
  cat("\nEstimation details:\n")
  cat("- Number of cohorts:", length(object$cohorts), "\n")
  cat("- Total switchers:", sum(object$c_weight), "\n")
  cat("- Dynamic effects estimated:", object$effects, "\n")
  cat("- Placebo tests:", object$placebo, "\n")
  
  invisible(object)
}

#' Extract coefficients from sdid_event results
#'
#' @param object Results from sdid_event
#' @param ... Additional arguments (ignored)
#'
#' @return Named vector of coefficients
#' @export
coef.sdid_event <- function(object, ...) {
  coefs <- object$estimates[, "Estimate"]
  names(coefs) <- rownames(object$estimates)
  return(coefs)
}

#' Extract variance-covariance matrix
#'
#' @param object Results from sdid_event  
#' @param ... Additional arguments (ignored)
#'
#' @return Variance-covariance matrix
#' @export
vcov.sdid_event <- function(object, ...) {
  if (is.null(object$vcov)) {
    warning("Variance-covariance matrix not available. Run with vce='bootstrap'.")
    return(NULL)
  }
  return(object$vcov)
}

#' Create LaTeX table from results
#'
#' @param object Results from sdid_event
#' @param file Output file path (optional)
#' @param digits Number of decimal places
#' @param ... Additional arguments
#'
#' @return LaTeX table as character string
#' @export
latex_table <- function(object, file = NULL, digits = 3, ...) {
  
  # Format estimates
  est_df <- as.data.frame(object$estimates)
  
  # Round numeric columns
  numeric_cols <- sapply(est_df, is.numeric)
  est_df[numeric_cols] <- round(est_df[numeric_cols], digits)
  
  # Create table
  n_cols <- ncol(est_df)
  
  # Header
  latex <- "\\begin{table}[htbp]\n"
  latex <- paste0(latex, "\\centering\n")
  latex <- paste0(latex, "\\caption{SDID Event Study Results}\n")
  latex <- paste0(latex, "\\begin{tabular}{l", paste(rep("c", n_cols), collapse = ""), "}\n")
  latex <- paste0(latex, "\\hline\\hline\n")
  
  # Column names
  col_names <- colnames(est_df)
  col_names[col_names == "LB_CI"] <- "95\\% CI (Lower)"
  col_names[col_names == "UB_CI"] <- "95\\% CI (Upper)"
  
  latex <- paste0(latex, " & ", paste(col_names, collapse = " & "), " \\\\\n")
  latex <- paste0(latex, "\\hline\n")
  
  # Data rows
  for (i in 1:nrow(est_df)) {
    row_name <- rownames(est_df)[i]
    row_data <- est_df[i, ]
    
    # Format row
    row_str <- row_name
    for (j in 1:ncol(row_data)) {
      if (!is.na(row_data[, j])) {
        if (colnames(row_data)[j] == "Switchers") {
          row_str <- paste0(row_str, " & ", as.integer(row_data[, j]))
        } else {
          row_str <- paste0(row_str, " & ", sprintf(paste0("%.", digits, "f"), row_data[, j]))
        }
      } else {
        row_str <- paste0(row_str, " & ")
      }
    }
    
    latex <- paste0(latex, row_str, " \\\\\n")
  }
  
  # Footer
  latex <- paste0(latex, "\\hline\\hline\n")
  latex <- paste0(latex, "\\end{tabular}\n")
  latex <- paste0(latex, "\\end{table}\n")
  
  # Write to file if specified
  if (!is.null(file)) {
    writeLines(latex, file)
    cat("LaTeX table written to:", file, "\n")
  }
  
  return(invisible(latex))
}