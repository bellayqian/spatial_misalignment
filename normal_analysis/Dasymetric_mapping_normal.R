library(sp)
library(sf)
library(spdep)
library(MASS)
library(raster)
library(tidyverse)
library(ggplot2)

# Note: This script is adapted from Dasymetric_mapping_binomial.R for normal data

# Implement Dasymetric Mapping for Normal Data
dasymetric_mapping_normal <- function(misaligned_data) {
  # Extract components
  gridx <- misaligned_data$gridx
  gridy <- misaligned_data$gridy
  
  # Convert grids to centroids for gridx
  gridx_centroids <- st_centroid(gridx)
  
  # Find which y polygon contains each x centroid
  intersections <- st_within(gridx_centroids, gridy)
  
  # Create mapping table
  mapping_table <- data.frame(
    x_id = seq_len(nrow(gridx)),
    y_id = sapply(intersections, function(x) if(length(x) > 0) x[1] else NA)
  )
  
  # Function to aggregate normal data (using mean aggregation like Poisson version)
  aggregate_to_y_normal <- function(x_values, mapping) {
    y_values <- numeric(nrow(gridy))
    
    for(i in seq_along(y_values)) {
      x_indices <- mapping$x_id[mapping$y_id == i]
      if(length(x_indices) > 0) {
        y_values[i] <- mean(x_values[x_indices], na.rm = TRUE)
      } else {
        y_values[i] <- NA
      }
    }
    return(y_values)
  }
  
  # Transfer X variables to Y grid
  x_vars <- grep("covariate_x_", names(gridx), value = TRUE)
  result <- gridy
  
  for(var in x_vars) {
    # For normal data, we use mean aggregation
    result[[var]] <- aggregate_to_y_normal(gridx[[var]], mapping_table)
  }
  
  return(result)
}

# Function to analyze normal centroid mapping
analyze_centroid_mapping_normal <- function(original_data, mapped_data) {
  # Extract variables
  x_vars <- grep("covariate_x_", names(mapped_data), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data), value = TRUE)
  
  # Comparison summaries
  summaries <- list()
  
  # Show both original and mapped values for all variables
  for (var in c(x_vars, y_vars)) {
    # For x variables, get original from gridx
    if(grepl("covariate_x_", var)) {
      orig_values <- st_drop_geometry(original_data$gridx)[,var]
    } else {
      # For y variables, get original from gridy
      orig_values <- st_drop_geometry(original_data$gridy)[,var]
    }
    
    mapped_values <- st_drop_geometry(mapped_data)[,var]
    
    summaries[[var]] <- data.frame(
      Statistic = c("Mean", "Median", "SD", "IQR"),
      Original = c(
        mean(orig_values, na.rm = TRUE),
        median(orig_values, na.rm = TRUE),
        sd(orig_values, na.rm = TRUE),
        IQR(orig_values, na.rm = TRUE)
      ),
      Mapped = c(
        mean(mapped_values, na.rm = TRUE),
        median(mapped_values, na.rm = TRUE),
        sd(mapped_values, na.rm = TRUE),
        IQR(mapped_values, na.rm = TRUE)
      )
    )
  }
  
  # Calculate correlations for mapped variables (consistent with Poisson version)
  all_vars <- c(x_vars, y_vars)
  mapped_correlations <- mapped_data %>%
    st_drop_geometry() %>%
    dplyr::select(all_of(all_vars)) %>%
    cor(use = "complete.obs")
  
  return(list(
    variable_summaries = summaries,
    correlations = mapped_correlations
  ))
}

# Fit normal regression models to mapped data
fit_normal_models <- function(mapped_data) {
  # Remove geometry and convert to regular data frame
  mapped_data_df <- st_drop_geometry(mapped_data)
  
  # Get x and y variable names
  x_vars <- grep("covariate_x_", names(mapped_data_df), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data_df), value = TRUE)
  all_vars <- c(x_vars, y_vars)
  
  # Standardize variables for easier interpretation (same as Poisson version)
  for(var in all_vars) {
    mapped_data_df[[var]] <- scale(mapped_data_df[[var]])
  }
  
  # Remove any NA or infinite values
  for(var in c(all_vars, "y")) {
    mapped_data_df[[var]][is.na(mapped_data_df[[var]]) | is.infinite(mapped_data_df[[var]])] <- 0
  }
  
  # Fit normal linear model: y as continuous response
  formula <- as.formula(paste("y ~", paste(all_vars, collapse = " + ")))
  
  # Fit the model with error handling
  model <- lm(formula, data = mapped_data_df)
  
  # Calculate coefficients and CIs (consistent format with Poisson version)
  coef_table <- summary(model)$coefficients
  coefficients <- coef_table[,1]
  std_errors <- coef_table[,2]
  ci_lower <- coefficients - 1.96*std_errors
  ci_upper <- coefficients + 1.96*std_errors
  
  results <- data.frame(
    variable = names(coefficients),
    coefficient = coefficients,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  
  return(list(
    model = model,
    coefficients = results,
    r_squared = summary(model)$r.squared,
    adj_r_squared = summary(model)$adj.r.squared
  ))
}

# Compare true vs estimated for normal case
compare_true_vs_estimated_normal <- function(model_results, true_beta_x, true_beta_y) {
  if(is.null(model_results) || is.null(model_results$model)) {
    warning("Invalid model results provided")
    return(NULL)
  }
  
  # Get estimated coefficients and SEs from model results
  coef_table <- model_results$model$coefficients[-1] # Remove intercept
  se_table <- summary(model_results$model)$coefficients[-1, "Std. Error"]
  
  # Create data frame of true values
  true_values <- c(true_beta_x, true_beta_y)
  names(true_values) <- c(paste0("covariate_x_", 1:length(true_beta_x)),
                          paste0("covariate_y_", 1:length(true_beta_y)))
  
  # Ensure we have matching names (same approach as Poisson version)
  common_names <- intersect(names(coef_table), names(true_values))
  
  if(length(common_names) == 0) {
    warning("No matching variable names between model and true values")
    return(NULL)
  }
  
  # Subset to common variables
  coef_subset <- coef_table[common_names]
  se_subset <- se_table[common_names]
  true_subset <- true_values[common_names]
  
  # Combine into a data frame (consistent format with Poisson version)
  comparison <- data.frame(
    variable = common_names,
    true_beta = true_subset,
    estimated_beta = coef_subset,
    std_error = se_subset,
    bias = coef_subset - true_subset,
    relative_bias = ((coef_subset - true_subset) / true_subset) * 100
  )
  
  # Add 95% CI
  comparison$ci_lower <- comparison$estimated_beta - 1.96*comparison$std_error
  comparison$ci_upper <- comparison$estimated_beta + 1.96*comparison$std_error
  
  # Check if true value falls within CI
  comparison$within_ci <- comparison$true_beta >= comparison$ci_lower & 
    comparison$true_beta <= comparison$ci_upper
  
  # Remove row names
  rownames(comparison) <- NULL
  
  return(comparison)
}

# Visualization function for normal centroid mapping
plot_centroid_mapping_normal <- function(original_data, mapped_data, variable) {
  # For normal data, we plot the actual values
  if(grepl("covariate_x_", variable)) {
    # For X variables, use X grid data
    p1 <- ggplot() +
      geom_sf(data = original_data$gridx, aes_string(fill = variable)) +
      scale_fill_viridis_c(name = "Value") +
      ggtitle("Original Grid Values") +
      theme_minimal()
  } else {
    # For Y variables, use Y grid data
    p1 <- ggplot() +
      geom_sf(data = original_data$gridy, aes_string(fill = variable)) +
      scale_fill_viridis_c(name = "Value") +
      ggtitle("Original Grid Values") +
      theme_minimal()
  }
  
  p2 <- ggplot() +
    geom_sf(data = mapped_data, aes_string(fill = variable)) +
    scale_fill_viridis_c(name = "Value") +
    ggtitle("Mapped Values (Centroid Method)") +
    theme_minimal()
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# Visualize coefficients from normal model
plot_normal_coefficients <- function(model_results) {
  if(is.null(model_results) || is.null(model_results$coefficients)) {
    warning("Invalid model results for plotting")
    return(NULL)
  }
  
  coef_data <- model_results$coefficients
  
  # Remove intercept for plotting
  coef_data <- coef_data[coef_data$variable != "(Intercept)",]
  
  # Create plot for coefficients (similar to rate ratio plot in Poisson version)
  p1 <- ggplot(coef_data, aes(x = variable, y = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model Coefficients",
         x = "Variable", 
         y = "Coefficient (95% CI)") +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(p1)
}

#################### Example usage
# # Test with normal data
# test_data_normal <- simulate_misaligned_data(
#   res1 = c(5, 5),
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.6,
#   y_correlation = 0.8,
#   beta_x = c(0.5, -0.3, 0.8),
#   beta_y = c(0.4, 0.2, -0.6, 0.3)
# )
# 
# # Apply normal mapping and fit model
# mapped_result_normal <- dasymetric_mapping_normal(test_data_normal)
# model_results_normal <- fit_normal_models(mapped_result_normal)
# 
# # Compare true vs estimated values
# comparison_normal <- compare_true_vs_estimated_normal(
#   model_results_normal,
#   true_beta_x = c(0.5, -0.3, 0.8),
#   true_beta_y = c(0.4, 0.2, -0.6, 0.3)
# )
# print(comparison_normal)
# 
# # Analyze mapping results
# analysis_result_normal <- analyze_centroid_mapping_normal(test_data_normal, mapped_result_normal)
# 
# # Create visualizations
# plot_centroid_mapping_normal(test_data_normal, mapped_result_normal, "covariate_x_1")
# plot_normal_coefficients(model_results_normal)
