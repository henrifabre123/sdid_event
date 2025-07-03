test_that("Different methods produce different results", {
  skip_if_not_installed("synthdid")
  
  data <- generate_test_data()# Unit tests for sdidevent package
  library(testthat)
  library(sdidevent)
  
  # Helper function to generate test data
  generate_test_data <- function(n_units = 30, n_periods = 20, 
                                 n_treated = 5, treatment_time = 10,
                                 effect_size = 5, seed = 123) {
    set.seed(seed)
    
    data <- expand.grid(
      unit = 1:n_units,
      time = 1:n_periods
    )
    
    # Assign treatment
    treated_units <- sample(1:n_units, n_treated)
    data$treated <- 0
    data$treated[data$unit %in% treated_units & data$time >= treatment_time] <- 1
    
    # Generate outcome
    data$unit_fe <- rep(rnorm(n_units), each = n_periods)
    data$time_fe <- rep(rnorm(n_periods), n_units)
    data$treatment_effect <- data$treated * effect_size
    data$outcome <- 10 + data$unit_fe + data$time_fe + 
      data$treatment_effect + rnorm(nrow(data))
    
    return(data)
  }
  
  test_that("synthdid package is required", {
    skip_if(requireNamespace("synthdid", quietly = TRUE), 
            "Test only relevant when synthdid not installed")
    
    data <- generate_test_data()
    
    expect_error(
      sdid_event(
        Y = data$outcome,
        G = data$unit,
        T = data$time,
        D = data$treated,
        vce = "off"
      ),
      "synthdid.*required"
    )
  })
  
  test_that("Basic sdid_event runs without errors", {
    skip_if_not_installed("synthdid")
    
    data <- generate_test_data()
    
    expect_no_error(
      results <- sdid_event(
        Y = data$outcome,
        G = data$unit,
        T = data$time,
        D = data$treated,
        vce = "off",
        verbose = FALSE
      )
    )
    
    expect_type(results, "list")
    expect_true("estimates" %in% names(results))
  })
  
  test_that("Different methods produce different results", {
    skip_if_not_installed("synthdid")
    
    data <- generate_test_data()
    
    results_sdid <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      method = "sdid", vce = "off", verbose = FALSE
    )
    
    results_did <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      method = "did", vce = "off", verbose = FALSE
    )
    
    results_sc <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      method = "sc", vce = "off", verbose = FALSE
    )
    
    # Methods should produce different estimates
    expect_false(
      all(results_sdid$estimates[, "Estimate"] == 
            results_did$estimates[, "Estimate"])
    )
    
    expect_false(
      all(results_sdid$estimates[, "Estimate"] == 
            results_sc$estimates[, "Estimate"])
    )
  })
  
  test_that("Effects parameter works correctly", {
    skip_if_not_installed("synthdid")
    
    data <- generate_test_data(n_periods = 30, treatment_time = 15)
    
    # Test with different numbers of effects
    results_5 <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 5, vce = "off", verbose = FALSE
    )
    
    results_10 <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 10, vce = "off", verbose = FALSE
    )
    
    # Should have correct number of rows (ATT + effects)
    expect_equal(nrow(results_5$estimates), 6)  # ATT + 5 effects
    expect_equal(nrow(results_10$estimates), 11) # ATT + 10 effects
  })
  
  test_that("Placebo tests work correctly", {
    data <- generate_test_data(n_periods = 30, treatment_time = 20)
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 5, placebo = 3,
      vce = "off", verbose = FALSE
    )
    
    # Should have ATT + 5 effects + 3 placebos = 9 rows
    expect_equal(nrow(results$estimates), 9)
    
    # Check row names
    expect_true("Placebo_1" %in% rownames(results$estimates))
    expect_true("Placebo_3" %in% rownames(results$estimates))
  })
  
  test_that("Bootstrap variance estimation works", {
    data <- generate_test_data()
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 3, vce = "bootstrap",
      brep = 20, verbose = FALSE
    )
    
    # Should have standard errors
    expect_true("SE" %in% colnames(results$estimates))
    expect_true(all(!is.na(results$estimates[, "SE"])))
    
    # Should have confidence intervals
    expect_true("LB_CI" %in% colnames(results$estimates))
    expect_true("UB_CI" %in% colnames(results$estimates))
  })
  
  test_that("Covariate adjustment works", {
    data <- generate_test_data()
    data$cov1 <- rnorm(nrow(data))
    data$cov2 <- runif(nrow(data))
    
    # Create data frame for covariate model
    cov_data <- data.frame(
      Y = data$outcome,
      G = data$unit,
      T = data$time,
      D = data$treated,
      cov1 = data$cov1,
      cov2 = data$cov2
    )
    
    expect_no_error(
      results <- sdid_event(
        Y = cov_data$Y, G = cov_data$G,
        T = cov_data$T, D = cov_data$D,
        covariates = list(vars = c("cov1", "cov2")),
        vce = "off", verbose = FALSE
      )
    )
  })
  
  test_that("prepare_sdid_data works correctly", {
    data <- generate_test_data()
    
    prepared <- prepare_sdid_data(
      data = data,
      unit_var = "unit",
      time_var = "time",
      outcome_var = "outcome",
      treatment_var = "treated"
    )
    
    expect_equal(ncol(prepared), 4)
    expect_equal(colnames(prepared), c("G", "T", "Y", "D"))
    expect_equal(nrow(prepared), nrow(data))
  })
  
  test_that("S3 methods work correctly", {
    data <- generate_test_data()
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 3, vce = "off", verbose = FALSE
    )
    
    # Test coef method
    coeffs <- coef(results)
    expect_type(coeffs, "double")
    expect_equal(length(coeffs), 4)  # ATT + 3 effects
    
    # Test summary method
    expect_output(summary(results), "SDID Event Study Results")
  })
  
  test_that("Multiple cohorts handled correctly", {
    # Generate data with two treatment cohorts
    set.seed(123)
    data <- expand.grid(unit = 1:40, time = 1:25)
    
    # First cohort treats at time 10, second at time 15
    cohort1_units <- 1:5
    cohort2_units <- 6:10
    
    data$treated <- 0
    data$treated[data$unit %in% cohort1_units & data$time >= 10] <- 1
    data$treated[data$unit %in% cohort2_units & data$time >= 15] <- 1
    
    # Generate outcome
    data$outcome <- rnorm(nrow(data)) + data$treated * 5
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 5, vce = "off", verbose = FALSE
    )
    
    expect_equal(length(results$cohorts), 2)
    expect_equal(results$cohorts, c(10, 15))
  })
  
  test_that("Error handling works correctly", {
    data <- generate_test_data()
    
    # Invalid method
    expect_error(
      sdid_event(
        Y = data$outcome, G = data$unit,
        T = data$time, D = data$treated,
        method = "invalid"
      ),
      "Invalid method"
    )
    
    # Invalid vce
    expect_error(
      sdid_event(
        Y = data$outcome, G = data$unit,
        T = data$time, D = data$treated,
        vce = "invalid"
      ),
      "Invalid vce"
    )
    
    # Non-binary treatment
    data$treated[1] <- 2
    expect_error(
      prepare_sdid_data(
        data = data,
        unit_var = "unit",
        time_var = "time",
        outcome_var = "outcome",
        treatment_var = "treated"
      ),
      "binary"
    )
  })
  
  test_that("Combinations work correctly", {
    data <- generate_test_data(n_periods = 30, treatment_time = 10)
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      effects = 6,
      combine = c("1 2 3", "4 5 6"),
      vce = "off", verbose = FALSE
    )
    
    expect_true(!is.null(results$H_combined))
    expect_equal(nrow(results$H_combined), 2)
    expect_equal(rownames(results$H_combined), c("Cmb_Effect_1", "Cmb_Effect_2"))
  })
  
  test_that("Disaggregated results work", {
    data <- generate_test_data()
    
    results <- sdid_event(
      Y = data$outcome, G = data$unit,
      T = data$time, D = data$treated,
      disag = TRUE, vce = "off", verbose = FALSE
    )
    
    expect_true(!is.null(results$H_cohort))
    expect_true(is.matrix(results$H_cohort))
  })
  
  # Run all tests
  test_dir("tests/testthat/")