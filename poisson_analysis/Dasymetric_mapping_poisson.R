library(sp)
library(sf)
library(spdep)
library(MASS)
library(raster)
library(tidyverse)
library(nimble)
library(ggplot2)

# Data Simulation: No changes
gen_correlated_spat <- function(W, n_vars, rho = 0.6, var_spat = 1, correlation = 0.5, 
                                global_ints = NULL, verify = FALSE) {
  n <- nrow(W)
 
  # Create precision matrix for variables (Lambda)
  Sigma_vars <- matrix(correlation, n_vars, n_vars) # Correlation Matrix 
  diag(Sigma_vars) <- 1
  Lambda <- solve(Sigma_vars) 
  
  # Spatial precision matrix
  precision_matrix <- diag(colSums(W)) - rho * W
  
  # Create full precision matrix using Kronecker product
  full_prec <- kronecker(Lambda, precision_matrix)
  
  # Generate from MCAR distribution
  all_effects <- mvrnorm(1, mu = rep(0, n * n_vars), Sigma = solve(full_prec))
  spatial_effects <- matrix(all_effects, n, n_vars)
  
  if(verify) {
    # Check empirical correlation of spatial effects
    emp_cor <- cor(spatial_effects)
    cat("Theoretical correlation matrix:\n")
    print(Sigma_vars)
    cat("\nEmpirical correlation matrix of spatial effects:\n")
    print(emp_cor)
  }
  
  # Generate observations, put them back
  observed_values <- matrix(0, n, n_vars)
  for(j in 1:n_vars) {
    int <- ifelse(is.null(global_ints), 0.1, global_ints[j])
    observed_values[,j] <- rpois(n, lambda = exp(int + spatial_effects[,j]))
  }
  
  return(observed_values)
}

# Main simulation function
simulate_misaligned_data <- function(res1 = c(5, 5), res2 = c(10, 10), 
                                     seed = 2, n_covariates_x = 3, n_covariates_y = 4,
                                     x_correlation = 0.5, y_correlation = 0.5,
                                     beta_x = NULL, beta_y = NULL) {
  set.seed(seed)
  
  # Create the first spatial grid (coarser grid)
  grid1 <- GridTopology(cellcentre.offset = c(0.5, 0.5), 
                        cellsize = c(2, 2),
                        cells.dim = res1)
  sp_grid1 <- SpatialGrid(grid1)
  
  # Create the second spatial grid (finer grid)
  grid2 <- GridTopology(cellcentre.offset = c(0, 0), 
                        cellsize = c(1, 1),
                        cells.dim = res2)
  sp_grid2 <- SpatialGrid(grid2)
  
  # Convert grids to sf objects
  sp_grid2_poly <- st_as_sf(as(sp_grid2, 'SpatialPolygons'))
  sp_grid1_poly <- st_as_sf(as(sp_grid1, 'SpatialPolygons'))
  
  # Create non-nested misalignment
  pair_list <- list(c(2,3), c(4,5), c(6,7), c(8,9))
  union_ids <- matrix(0, res2[1], res2[2])
  id_num <- 1
  for (i in seq(1, res2[1], 2)) {
    union_ids[i, pair_list[[sample(1:4, size = 1)]]] <- id_num
    id_num <- id_num + 1
  }
  
  sp_grid2_poly$union_ids <- c(t(union_ids))
  
  # Union cells
  store_xunions <- NULL
  for (i in 1:max(sp_grid2_poly$union_ids)) {
    temp <- sp_grid2_poly[sp_grid2_poly$union_ids == i, ] %>% 
      st_union() %>% 
      st_sf()
    store_xunions <- rbind(store_xunions, temp)
  }
  
  # Merge remaining cells
  grid2_nounion <- subset(sp_grid2_poly, union_ids == 0)
  store_iunions <- NULL
  for (i in 1:nrow(sp_grid1_poly)) {
    temp <- grid2_nounion[c(st_contains(sp_grid1_poly[i,], grid2_nounion, sparse = F)) == T,] %>% 
      st_union() %>% 
      st_sf()
    store_iunions <- rbind(store_iunions, temp)
  }
  
  grid2_final <- rbind(store_iunions, store_xunions)
  
  # Rename grids
  gridx <- grid2_final
  gridy <- sp_grid1_poly
  gridy$ID <- 1:nrow(gridy)
  gridx$ID <- 1:nrow(gridx)
  
  # Generate X
  neighbors_x <- poly2nb(as(gridx, "Spatial"), queen = TRUE)
  W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)
  x_values <- gen_correlated_spat(W_x, n_covariates_x, 
                                  correlation = x_correlation,
                                  global_ints = rnorm(n_covariates_x, 3, 0.5))
  
  # Assign X values to grid
  for(i in 1:n_covariates_x) {
    gridx[[paste0("covariate_x_", i)]] <- x_values[,i]
  }
  
  # Generate correlated Y covariates
  neighbors_y <- poly2nb(as(gridy, "Spatial"), queen = TRUE)
  W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)
  y_values <- gen_correlated_spat(W_y, n_covariates_y, 
                                  correlation = y_correlation,
                                  global_ints = rnorm(n_covariates_y, 3, 0.5))
  
  # Assign Y covariate values
  for(i in 1:n_covariates_y) {
    gridy[[paste0("covariate_y_", i)]] <- y_values[,i]
  }
  
  # Create atoms for outcome generation
  atoms <- st_as_sf(raster::intersect(as(gridy, 'Spatial'), as(gridx, 'Spatial')))
  names(atoms)[which(names(atoms) == 'ID')] <- c("ID_y")
  names(atoms)[which(names(atoms) == 'ID.1')] <- c("ID_x")
  atoms$ID_atomorder <- 1:nrow(atoms)
  
  # Generate outcome with true associations
  spatial_effect_y <- gen_correlated_spat(W_y, n_vars = 1, rho = 0.6, 
                                          var_spat = 1, # Base spatial effect
                                          correlation = 1)[,1]
  
  # Calculate linear predictor incorporating covariate effects
  linear_pred_y <- rep(0.01, nrow(gridy))  # Smaller intercept on log scale
  
  # First aggregate X covariates to Y grid level using atoms
  for(i in 1:n_covariates_x) {
    # Get X covariate values
    x_vals <- gridx[[paste0("covariate_x_", i)]]
    
    # Use atoms to aggregate to Y grid
    for(j in 1:nrow(gridy)) {
      # Find atoms in this Y grid cell
      atom_indices <- which(atoms$ID_y == j)
      
      if(length(atom_indices) == 0) {
        warning(paste("No atoms found for Y grid cell", j))
        next
      }
      
      # Average X values over these atoms
      x_agg <- mean(x_vals[atoms$ID_x[atom_indices]])
      # Add effect to linear predictor
      linear_pred_y[j] <- linear_pred_y[j] + beta_x[i] * x_agg
    }
  }
  
  # Add Y covariate effects
  for(i in 1:n_covariates_y) {
    linear_pred_y <- linear_pred_y + 
      beta_y[i] * gridy[[paste0("covariate_y_", i)]]
  }
  
  # Generate final outcome
  gridy$y <- rpois(nrow(gridy), 
                   lambda = exp(linear_pred_y + spatial_effect_y))
  
  # Store true parameters
  true_params <- list(
    beta_x = beta_x,
    beta_y = beta_y,
    x_correlation = x_correlation,
    y_correlation = y_correlation
  )
  
  return(list(
    gridy = gridy, 
    gridx = gridx, 
    atoms = atoms, 
    true_params = true_params
  ))
}

# Implement Dasymetric Mapping 
dasymetric_mapping <- function(misaligned_data) {
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
  
  # Function to aggregate values
  aggregate_to_y <- function(x_values, mapping) {
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
    result[[var]] <- aggregate_to_y(gridx[[var]], mapping_table)
  }
  
  return(result)
}

# Function to analyze dasymetric mapping results
analyze_centroid_mapping <- function(original_data, mapped_data) {
  # Extract variables
  x_vars <- grep("covariate_x_", names(mapped_data), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data), value = TRUE)
  
  # Comparison summaries
  summaries <- list()
  
  # Show both original and mapped values
  for (var in c(x_vars, y_vars)) {
    # For x variables, get original from gridx
    if(grepl("covariate_x_", var)) {
      orig_vals <- st_drop_geometry(original_data$gridx)[,var]
    } else {
      # For y variables, get original from gridy
      orig_vals <- st_drop_geometry(original_data$gridy)[,var]
    }
    
    mapped_vals <- st_drop_geometry(mapped_data)[,var]
    
    summaries[[var]] <- data.frame(
      Statistic = c("Mean", "Median", "SD", "IQR"),
      Original = c(
        mean(orig_vals, na.rm = TRUE),
        median(orig_vals, na.rm = TRUE),
        sd(orig_vals, na.rm = TRUE),
        IQR(orig_vals, na.rm = TRUE)
      ),
      Mapped = c(
        mean(mapped_vals, na.rm = TRUE),
        median(mapped_vals, na.rm = TRUE),
        sd(mapped_vals, na.rm = TRUE),
        IQR(mapped_vals, na.rm = TRUE)
      )
    )
  }
  
  # Overall correlations in mapped data
  mapped_correlations <- mapped_data %>%
    st_drop_geometry() %>%
    dplyr::select(all_of(c(x_vars, y_vars))) %>%
    cor(use = "complete.obs")
  
  return(list(
    variable_summaries = summaries,
    correlations = mapped_correlations
  ))
}

# Fit Poisson regression models to mapped data
fit_poisson_models <- function(mapped_data) {
  # Remove geometry and convert to regular data frame
  mapped_data <- st_drop_geometry(mapped_data)
  
  # Get x and y variable names
  x_vars <- grep("covariate_x_", names(mapped_data), value = TRUE)
  y_vars <- grep("covariate_y_", names(mapped_data), value = TRUE)
  all_vars <- c(x_vars, y_vars)
  
  # Standardize variables for easier interpretation
  for(var in all_vars) {
    mapped_data[[var]] <- scale(mapped_data[[var]])
  }
  
  # Fit Poisson model: y as the count outcome
  formula <- as.formula(paste("y ~", paste(all_vars, collapse = " + ")))
  model <- glm(formula, family = poisson(link = "log"), data = mapped_data)
  
  # Calculate rate ratios and CIs
  coef_table <- summary(model)$coefficients
  rate_ratios <- exp(coef_table[,1])
  ci_lower <- exp(coef_table[,1] - 1.96*coef_table[,2])
  ci_upper <- exp(coef_table[,1] + 1.96*coef_table[,2])
  
  results <- data.frame(
    variable = names(rate_ratios),
    rate_ratio = rate_ratios,
    ci_lower = ci_lower, 
    ci_upper = ci_upper
  )
  
  return(list(
    model = model,
    rate_ratios = results
  ))
}

# Compared with true value
compare_true_vs_estimated <- function(model_results, true_beta_x, true_beta_y) {
  # Get estimated coefficients and SEs from model results
  coef_table <- model_results$model$coefficients[-1] # Remove intercept
  se_table <- summary(model_results$model)$coefficients[-1, "Std. Error"]
  
  # Create data frame of true values
  true_values <- c(true_beta_x, true_beta_y)
  names(true_values) <- c(paste0("covariate_x_", 1:length(true_beta_x)),
                          paste0("covariate_y_", 1:length(true_beta_y)))
  
  # Combine into a data frame
  comparison <- data.frame(
    variable = names(true_values),
    true_beta = true_values,
    estimated_beta = coef_table,
    std_error = se_table,
    bias = coef_table - true_values,
    relative_bias = ((coef_table - true_values) / true_values) * 100
  )
  
  # Add 95% CI
  comparison$ci_lower <- comparison$estimated_beta - 1.96*comparison$std_error
  comparison$ci_upper <- comparison$estimated_beta + 1.96*comparison$std_error
  
  # Check if true value falls within CI
  comparison$within_ci <- comparison$true_beta >= comparison$ci_lower & 
    comparison$true_beta <= comparison$ci_upper
  
  return(comparison)
}

### Visualization function for centroid mapping
plot_centroid_mapping <- function(original_data, mapped_data, variable) {
  p1 <- ggplot() +
    geom_sf(data = original_data$gridx, aes_string(fill = variable)) +
    scale_fill_viridis_c() +
    ggtitle("Original Grid Values") +
    theme_minimal()
  
  p2 <- ggplot() +
    geom_sf(data = mapped_data, aes_string(fill = variable)) +
    scale_fill_viridis_c() +
    ggtitle("Mapped Values (Centroid Method)") +
    theme_minimal()

  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# Visualize rate ratios
plot_rate_ratios <- function(model_results) {
  rr <- model_results$rate_ratios
  
  # Remove intercept
  rr <- rr[rr$variable != "(Intercept)",]
  
  # Create plot
  ggplot(rr, aes(x = variable, y = rate_ratio)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
    coord_flip() +
    theme_minimal() +
    xlab("Variable") + 
    ylab("Rate Ratio (95% CI)") +
    geom_hline(yintercept = 1, linetype = "dashed")
}

#################### Run code 
# test_data <- simulate_misaligned_data(
#   res1 = c(5, 5), 
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.6,
#   y_correlation = 0.8,
#   beta_x = c(0.03, -0.01, 0.06),
#   beta_y = c(0.04, 0.01, -0.05, 0.02)
# )
# 
# # Apply mapping and fit model
# mapped_result <- dasymetric_mapping(test_data)
# model_results <- fit_poisson_models(mapped_result)
# 
# # Compare true vs estimated values
# comparison <- compare_true_vs_estimated(
#   model_results,
#   true_beta_x = c(0.03, -0.01, 0.06),
#   true_beta_y = c(0.04, 0.01, -0.05, 0.02)
# )
# print(comparison)

#################### Analyze mapping results
# analysis_result <- analyze_centroid_mapping(test_data, mapped_result)
# 
# cat("Summary Statistics by Variable:\n")
# for(var in names(analysis_result$variable_summaries)) {
#   cat("\n", var, ":\n")
#   print(analysis_result$variable_summaries[[var]])
# }
# 
# cat("\nCorrelation Matrix of Mapped Variables:\n")
# print(round(analysis_result$correlations, 3))
# 
# # Create visualization for grids and mapping
# plot_centroid_mapping(test_data, mapped_result, "covariate_x_1")
# 
# # Comparison result
# ggplot(comparison, aes(x = variable)) +
#   geom_point(aes(y = estimated_beta), color = "blue") +
#   geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "blue") +
#   geom_point(aes(y = true_beta), color = "red", shape = 4, size = 3) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "Comparison of True vs Estimated Coefficients",
#        subtitle = "Blue points/bars = Estimated values with 95% CI\nRed X = True values",
#        x = "Variable",
#        y = "Coefficient Value")


