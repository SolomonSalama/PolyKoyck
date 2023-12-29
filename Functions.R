adf_kpss_test <- function(time_series_variable) {
  # Compute partial autocorrelation function
  pacf_values <- pacf(time_series_variable, plot = FALSE)$acf

  # Calculate the critical value for the desired confidence level
  critical_value <- qnorm((1 + 0.95)/2) / sqrt(length(time_series_variable))

  # Initialize lag order
  lag_order <- 1

  # Sequentially test for significant PACF lags
  for (lag in 1:(length(pacf_values)-1)) {
    if (abs(pacf_values[lag]) > critical_value) {
      lag_order <- lag
    } else {
      break  # Break the loop when an insignificant lag is found
    }
  }

  # Printing PACF plot
  pacf(time_series_variable, main = paste(deparse(substitute(time_series_variable))))

  # Perform Augmented Dickey-Fuller test
  adf_result <- adf.test(time_series_variable, k = lag_order)

  # Perform Kwiatkowski-Phillips-Schmidt-Shin test with null = "Trend"
  kpss_result <- kpss.test(time_series_variable, null = "Trend")

  # Print ADF and KPSS test results
  cat("ADF Test Results:\n")
  print(adf_result)
  cat("\nKPSS Test Results:\n")
  print(kpss_result)

  # Print joint ADF and KPSS results
  if (adf_result$p.value < 0.05 & kpss_result$p.value >= 0.05) {
    cat("ADF test concludes stationarity and KPSS test concludes stationarity against trend stationarity null: stationary process\n")
  } else if (adf_result$p.value >= 0.05 & kpss_result$p.value >= 0.05) {
    cat("ADF test concludes non-stationarity and KPSS test concludes stationarity against trend stationarity null: trend-stationary process\n")
  } else if (adf_result$p.value < 0.05 & kpss_result$p.value < 0.05) {
    cat("ADF test concludes stationarity and KPSS test concludes non-stationarity against trend stationarity null: difference-stationary process\n")
  } else if (adf_result$p.value >= 0.05 & kpss_result$p.value < 0.05) {
    cat("ADF test concludes non-stationarity and KPSS test concludes non-stationarity against trend stationarity null: non-stationary process\n")
  }
}
auto_poly_koyck <- function(data, dependent_variable, independent_variables, max_lags = 10, bg_test_level = 0.10, L = "L", print_output = TRUE) {
  # Initial lag
  lag_number <- 1

  while (lag_number <= max_lags) {
    # Create lagged variables
    lagged_vars <- paste0(dependent_variable, "_lag", 1:lag_number)

    # Building formula without trend variable
    formula_str <- paste(c(dependent_variable, lagged_vars, independent_variables), collapse = " + ")
    formula_str <- gsub("trend", "", formula_str)  # Remove "trend" term

    # Ensuring that the formula is correctly formed
    formula <- as.formula(paste(dependent_variable, "~", formula_str))

    # Fit OLS/LM model
    ols_model <- lm(formula, data)

    # BG test for autocorrelation of order 1
    bg_test <- lmtest::bgtest(ols_model, order = 1)

    # Check p-value
    if (bg_test$p.value > bg_test_level) {
      # Stopping without printing
      break
    } else {
      lag_number <- lag_number + 1

      # Creating the next lag variable in the dataset
      next_lag_var <- paste0(dependent_variable, "_lag", lag_number)
      data[[next_lag_var]] <- lag(data[[dependent_variable]], lag_number)
    }
  }

  # Extracting coefficients from the final model for all lagged dependent variables
  poly <- round(coef(ols_model)[grep(paste0(dependent_variable, "_lag"), names(coef(ols_model)))],6)
  poly <- c(-1, poly)  # Fix the first coefficient at 1

  # Displaying the final model summary
  cat("\tDynamically Complete Koyck Lag Model\n")
  print(summary(ols_model))

  # Displaying the results of the final Breusch-Godfrey test
  print(bg_test)



  # Introduction to the function for polynomial root derivation and testing for unit roots
  cat("\tPolynomial Root Derivation and Testing for Unit Roots\n")
  cat("\n")

  # Print the Polynomial in question
  cat("Polynomial: ")
  cat(paste(ifelse(poly != 0, paste0(ifelse(poly > 0, "-", "+"), abs(poly), L, "^", seq_along(poly)-1), ""), collapse = " "), "\n\n")

  # Calculating polynomial roots and save them in a dataframe
  res <- data.frame(
    real = round(Re(polyroot(poly)), 6),
    complex = round(Im(polyroot(poly)), 6)
  )

  # Checking if roots outside the unit circle
  res <- data.frame(res, outside = c(sqrt(Re(polyroot(poly))^2 + Im(polyroot(poly))^2) > 1))

  # print results
  if (print_output) {
    print(res)
    cat("*Results are rounded to 6 digits.\n")

    # Check if all roots are outside the unit circle for polynomial decay function condition
    if (all(res$outside)) {
      cat("\n")
      cat("All roots are outside the unit circle: Koyck lag weights converge to zero\n")
    } else {
      cat("\n")
      cat("Some roots are inside the unit circle: Koyck lag weights do not converge to zero\n")
    }
  }

  invisible(res)
}

