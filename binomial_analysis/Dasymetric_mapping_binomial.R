library(sp)
library(sf)
library(spdep)
library(MASS)
library(raster)
library(tidyverse)
library(ggplot2)

# Note: This script is adapted from Dasymetric_mapping.R for binomial data
# The key difference is that we now work with counts AND population sizes

# The data simulation functions remain the same since they're already in function_binomial_0508.R
# We just need the mapping and analysis functions

# Implement Dasymetric Mapping for Binomial Data
dasymetric_mapping_binomial <- function(misaligned_data) {
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
  
  # Function to aggregate binomial data (both counts and populations)
  aggregate_to_y_binomial <- function(x_counts, x_populations, mapping) {
    y_counts <- numeric(nrow(gridy))
    y_populations <- numeric(nrow(gridy))
    
    for(i in seq_along(y_counts)) {
      x_indices <- mapping$x_id[mapping$y_id == i]
      if(length(x_indices) > 0) {
        # Sum both counts and populations
        y_counts[i] <- sum(x_counts[x_indices], na.rm = TRUE)
        y_populations[i] <- sum(x_populations[x_indices], na.rm = TRUE)
      } else {
        y_counts[i] <- 0
        y_populations[i] <- 0
      }
    }
    return(list(counts = y_counts, populations = y_populations))
  }
  
  # Transfer X variables to Y grid
  x_vars <- grep("covariate_x_", names(gridx), value = TRUE)
  result <- gridy
  
  for(var in x_vars) {
    aggregated <- aggregate_to_y_binomial(
      gridx[[var]], 
      gridx$population, 
      mapping_table
    )
    result[[var]] <- aggregated$counts
    result[[paste0(var, "_pop")]] <- aggregated$populations
  }
  
  return(result)
}

# Function to analyze binomial centroid mapping
analyze_centroid_mapping_binomial <- function(original_data, mapped_data) {
  # Extract variables
  x_vars <- grep("covariate_x_", names(mapped_data), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data), value = TRUE)
  
  # Comparison summaries
  summaries <- list()
  
  # Show both original and mapped values for X variables
  for (var in x_vars) {
    # For x variables, get original from gridx
    orig_counts <- st_drop_geometry(original_data$gridx)[,var]
    orig_pops <- st_drop_geometry(original_data$gridx)[,"population"]
    orig_props <- orig_counts / orig_pops
    
    mapped_counts <- st_drop_geometry(mapped_data)[,var]
    mapped_pops <- st_drop_geometry(mapped_data)[,paste0(var, "_pop")]
    mapped_props <- mapped_counts / pmax(mapped_pops, 1)  # Avoid division by zero
    
    summaries[[var]] <- data.frame(
      Statistic = c("Mean_Count", "Mean_Pop", "Mean_Prop", "SD_Prop"),
      Original = c(
        mean(orig_counts, na.rm = TRUE),
        mean(orig_pops, na.rm = TRUE),
        mean(orig_props, na.rm = TRUE),
        sd(orig_props, na.rm = TRUE)
      ),
      Mapped = c(
        mean(mapped_counts, na.rm = TRUE),
        mean(mapped_pops, na.rm = TRUE),
        mean(mapped_props, na.rm = TRUE),
        sd(mapped_props, na.rm = TRUE)
      )
    )
  }
  
  # For Y variables, compare proportions
  for (var in y_vars) {
    orig_counts <- st_drop_geometry(original_data$gridy)[,var]
    orig_pops <- st_drop_geometry(original_data$gridy)[,"population"]
    orig_props <- orig_counts / orig_pops
    
    mapped_counts <- st_drop_geometry(mapped_data)[,var]
    mapped_pops <- st_drop_geometry(mapped_data)[,"population"]
    mapped_props <- mapped_counts / pmax(mapped_pops, 1)
    
    summaries[[var]] <- data.frame(
      Statistic = c("Mean_Count", "Mean_Pop", "Mean_Prop", "SD_Prop"),
      Original = c(
        mean(orig_counts, na.rm = TRUE),
        mean(orig_pops, na.rm = TRUE),
        mean(orig_props, na.rm = TRUE),
        sd(orig_props, na.rm = TRUE)
      ),
      Mapped = c(
        mean(mapped_counts, na.rm = TRUE),
        mean(mapped_pops, na.rm = TRUE),
        mean(mapped_props, na.rm = TRUE),
        sd(mapped_props, na.rm = TRUE)
      )
    )
  }
  
  # Calculate correlations for mapped proportions
  prop_vars <- c()
  for(var in c(x_vars, y_vars)) {
    if(var %in% x_vars) {
      counts <- mapped_data[[var]]
      pops <- mapped_data[[paste0(var, "_pop")]]
    } else {
      counts <- mapped_data[[var]]
      pops <- mapped_data[["population"]]
    }
    prop_var_name <- paste0(var, "_prop")
    mapped_data[[prop_var_name]] <- counts / pmax(pops, 1)
    prop_vars <- c(prop_vars, prop_var_name)
  }
  
  mapped_correlations <- mapped_data %>%
    st_drop_geometry() %>%
    dplyr::select(all_of(prop_vars)) %>%
    cor(use = "complete.obs")
  
  return(list(
    variable_summaries = summaries,
    correlations = mapped_correlations,
    proportion_data = mapped_data[prop_vars]
  ))
}

# Fit binomial regression models to mapped data
fit_binomial_models <- function(mapped_data) {
  # Remove geometry and convert to regular data frame
  mapped_data_df <- st_drop_geometry(mapped_data)
  
  # Get x and y variable names
  x_vars <- grep("covariate_x_", names(mapped_data_df), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data_df), value = TRUE)
  
  # Create proportion variables for X covariates
  for(var in x_vars) {
    pop_var <- paste0(var, "_pop")
    prop_var <- paste0(var, "_prop")
    if(pop_var %in% names(mapped_data_df)) {
      mapped_data_df[[prop_var]] <- mapped_data_df[[var]] / pmax(mapped_data_df[[pop_var]], 1)
    } else {
      warning(paste("Population variable", pop_var, "not found for", var))
      mapped_data_df[[prop_var]] <- mapped_data_df[[var]] / pmax(mapped_data_df$population, 1)
    }
  }
  
  # Create proportion variables for Y covariates
  for(var in y_vars) {
    prop_var <- paste0(var, "_prop")
    mapped_data_df[[prop_var]] <- mapped_data_df[[var]] / pmax(mapped_data_df$population, 1)
  }
  
  # Get proportion variable names
  x_prop_vars <- paste0(x_vars, "_prop")
  y_prop_vars <- paste0(y_vars, "_prop")
  all_prop_vars <- c(x_prop_vars, y_prop_vars)
  
  # Remove any NA or infinite values
  for(var in all_prop_vars) {
    mapped_data_df[[var]][is.na(mapped_data_df[[var]]) | is.infinite(mapped_data_df[[var]])] <- 0
    # Also bound proportions to [0,1]
    mapped_data_df[[var]] <- pmax(0, pmin(1, mapped_data_df[[var]]))
  }
  
  # Ensure y and population are valid for binomial
  mapped_data_df$y <- pmax(0, pmin(mapped_data_df$y, mapped_data_df$population))
  mapped_data_df$population <- pmax(1, mapped_data_df$population)  # Ensure population > 0
  
  # Fit binomial model: y as successes, population as trials
  formula <- as.formula(paste("cbind(y, population - y) ~", paste(all_prop_vars, collapse = " + ")))
  
  # Fit the model with error handling
  model <- tryCatch({
    glm(formula, family = binomial(link = "logit"), data = mapped_data_df)
  }, error = function(e) {
    warning("Binomial model fitting failed, trying with different approach")
    # Try with proportion as response and weights
    prop_formula <- as.formula(paste("I(y/population) ~", paste(all_prop_vars, collapse = " + ")))
    glm(prop_formula, family = binomial(link = "logit"), 
        data = mapped_data_df, weights = mapped_data_df$population)
  })
  
  # Calculate coefficients, odds ratios and CIs
  coef_table <- summary(model)$coefficients
  coefficients <- coef_table[,1]
  odds_ratios <- exp(coefficients)
  ci_lower_coef <- coefficients - 1.96*coef_table[,2]
  ci_upper_coef <- coefficients + 1.96*coef_table[,2]
  ci_lower_or <- exp(ci_lower_coef)
  ci_upper_or <- exp(ci_upper_coef)
  
  results <- data.frame(
    variable = names(coefficients),
    coefficient = coefficients,
    odds_ratio = odds_ratios,
    ci_lower_coef = ci_lower_coef,
    ci_upper_coef = ci_upper_coef,
    ci_lower_or = ci_lower_or, 
    ci_upper_or = ci_upper_or
  )
  
  return(list(
    model = model,
    coefficients = results,
    data_used = mapped_data_df
  ))
}

# Compare true vs estimated for binomial case
compare_true_vs_estimated_binomial <- function(model_results, true_beta_x, true_beta_y) {
  # Get estimated coefficients and SEs from model results
  coef_table <- model_results$model$coefficients[-1] # Remove intercept
  se_table <- summary(model_results$model)$coefficients[-1, "Std. Error"]
  
  # Create data frame of true values
  true_values <- c(true_beta_x, true_beta_y)
  names(true_values) <- c(paste0("covariate_x_", 1:length(true_beta_x), "_prop"),
                          paste0("covariate_y_", 1:length(true_beta_y), "_prop"))
  
  # Ensure we have matching names
  common_names <- intersect(names(coef_table), names(true_values))
  
  if(length(common_names) == 0) {
    warning("No matching variable names between model and true values")
    return(NULL)
  }
  
  # Subset to common variables
  coef_subset <- coef_table[common_names]
  se_subset <- se_table[common_names]
  true_subset <- true_values[common_names]
  
  # Combine into a data frame
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

# Visualization function for binomial centroid mapping
plot_centroid_mapping_binomial <- function(original_data, mapped_data, variable) {
  # Create proportion variables for plotting
  if(grepl("covariate_x_", variable)) {
    # For X variables, use X grid data
    orig_counts <- original_data$gridx[[variable]]
    orig_pops <- original_data$gridx$population
    orig_data <- original_data$gridx
    orig_data$proportion <- orig_counts / orig_pops
    
    mapped_counts <- mapped_data[[variable]]
    mapped_pops <- mapped_data[[paste0(variable, "_pop")]]
    mapped_data$proportion <- mapped_counts / pmax(mapped_pops, 1)
  } else {
    # For Y variables, use Y grid data
    orig_counts <- original_data$gridy[[variable]]
    orig_pops <- original_data$gridy$population
    orig_data <- original_data$gridy
    orig_data$proportion <- orig_counts / orig_pops
    
    mapped_counts <- mapped_data[[variable]]
    mapped_pops <- mapped_data$population
    mapped_data$proportion <- mapped_counts / pmax(mapped_pops, 1)
  }
  
  p1 <- ggplot() +
    geom_sf(data = orig_data, aes(fill = proportion)) +
    scale_fill_viridis_c(name = "Proportion") +
    ggtitle(paste("Original", variable, "Proportions")) +
    theme_minimal()
  
  p2 <- ggplot() +
    geom_sf(data = mapped_data, aes(fill = proportion)) +
    scale_fill_viridis_c(name = "Proportion") +
    ggtitle(paste("Mapped", variable, "Proportions")) +
    theme_minimal()
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# Visualize coefficients from binomial model
plot_binomial_coefficients <- function(model_results) {
  coef_data <- model_results$coefficients
  
  # Remove intercept for plotting
  coef_data <- coef_data[coef_data$variable != "(Intercept)",]
  
  # Create plot for coefficients
  p1 <- ggplot(coef_data, aes(x = variable, y = coefficient)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower_coef, ymax = ci_upper_coef), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Binomial Model Coefficients",
         x = "Variable", 
         y = "Coefficient (Log Odds)") +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  # Create plot for odds ratios
  p2 <- ggplot(coef_data, aes(x = variable, y = odds_ratio)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower_or, ymax = ci_upper_or), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Binomial Model Odds Ratios",
         x = "Variable", 
         y = "Odds Ratio") +
    geom_hline(yintercept = 1, linetype = "dashed")
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

#################### Example usage (commented out)
# # Test with binomial data
# test_data_binomial <- simulate_misaligned_binomial_data(
#   res1 = c(5, 5), 
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.6,
#   y_correlation = 0.8,
#   beta_x = c(0.03, -0.01, 0.06),
#   beta_y = c(0.04, 0.01, -0.05, 0.02),
#   pop_min = 100,
#   pop_max = 500
# )
# 
# # Apply binomial mapping and fit model
# mapped_result_binomial <- dasymetric_mapping_binomial(test_data_binomial)
# model_results_binomial <- fit_binomial_models(mapped_result_binomial)
# 
# # Compare true vs estimated values
# comparison_binomial <- compare_true_vs_estimated_binomial(
#   model_results_binomial,
#   true_beta_x = c(0.03, -0.01, 0.06),
#   true_beta_y = c(0.04, 0.01, -0.05, 0.02)
# )
# print(comparison_binomial)
# 
# # Analyze mapping results
# analysis_result_binomial <- analyze_centroid_mapping_binomial(test_data_binomial, mapped_result_binomial)
# 
# # Create visualizations
# plot_centroid_mapping_binomial(test_data_binomial, mapped_result_binomial, "covariate_x_1")
# plot_binomial_coefficients(model_results_binomial)