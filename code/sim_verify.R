################################################
source("function_0115.R")

# 1. Run multiple simulation study with different parameters
run_simulation_study <- function(sim_number, res1 = c(5, 5), res2 = c(10, 10),
                                 x_correlation, y_correlation,
                                 beta_x, beta_y, base_seed = 1) {
  
  set.seed(base_seed + sim_number)
  
  # Generate simulation data
  sim_data <- simulate_misaligned_data(
    res1 = res1,
    res2 = res2,
    n_covariates_x = length(beta_x),
    n_covariates_y = length(beta_y),
    x_correlation = x_correlation,
    y_correlation = y_correlation,
    beta_x = beta_x,
    beta_y = beta_y,
    seed = sim_number
  )
  
  # Run model and capture results
  model_output <- run_abrm(sim_data, paste0("sim_output_", sim_number, ".rds"))
  
  # Extract parameter estimates
  mcmc_summary <- model_output$summary$all.chains
  
  # Create results list
  results <- list(
    sim_number = sim_number,
    x_correlation = x_correlation,
    y_correlation = y_correlation,
    true_beta_x = beta_x,
    true_beta_y = beta_y,
    estimated_beta_x = mcmc_summary[paste0("beta_y[", 1:length(beta_x), "]"), "Mean"],
    estimated_beta_y = mcmc_summary[paste0("beta_y[", 
                                           (length(beta_x)+1):(length(beta_x)+length(beta_y)), "]"), "Mean"],
    beta_x_ci_lower = mcmc_summary[paste0("beta_y[", 1:length(beta_x), "]"), "95%CI_low"],
    beta_x_ci_upper = mcmc_summary[paste0("beta_y[", 1:length(beta_x), "]"), "95%CI_upp"],
    beta_y_ci_lower = mcmc_summary[paste0("beta_y[", 
                                          (length(beta_x)+1):(length(beta_x)+length(beta_y)), "]"), "95%CI_low"],
    beta_y_ci_upper = mcmc_summary[paste0("beta_y[", 
                                          (length(beta_x)+1):(length(beta_x)+length(beta_y)), "]"), "95%CI_upp"]
  )
  
  return(results)
}

# Function to run multiple simulations with different correlation strengths
run_correlation_study <- function(n_sims_per_setting = 10,
                                  correlation_grid = seq(0.2, 0.8, by = 0.2),
                                  beta_x = c(0.03, -0.05, 0.02),
                                  beta_y = c(0.04, 0.001, -0.08, 0.06)) {
  
  results_list <- list()
  counter <- 1
  
  # Loop through correlation settings
  for(x_cor in correlation_grid) {
    for(y_cor in correlation_grid) {
      cat(sprintf("\nRunning simulations for x_cor = %.2f, y_cor = %.2f\n", 
                  x_cor, y_cor))
      
      # Run multiple simulations for each correlation setting
      for(sim in 1:n_sims_per_setting) {
        cat(sprintf("  Simulation %d of %d\n", sim, n_sims_per_setting))
        
        # Changed from run_single_simulation to run_simulation_study
        results_list[[counter]] <- run_simulation_study(
          sim_number = counter,
          x_correlation = x_cor,
          y_correlation = y_cor,
          beta_x = beta_x,
          beta_y = beta_y
        )
        
        counter <- counter + 1
      }
    }
  }
  
  return(results_list)
}

# Function to analyze simulation results
analyze_simulation_results <- function(results_list) {
  # Convert list to data frame
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      sim_number = x$sim_number,
      x_correlation = x$x_correlation,
      y_correlation = x$y_correlation,
      parameter = c(paste0("beta_x", 1:length(x$true_beta_x)),
                    paste0("beta_y", 1:length(x$true_beta_y))),
      true_value = c(x$true_beta_x, x$true_beta_y),
      estimated_value = c(x$estimated_beta_x, x$estimated_beta_y),
      ci_lower = c(x$beta_x_ci_lower, x$beta_y_ci_lower),
      ci_upper = c(x$beta_x_ci_upper, x$beta_y_ci_upper)
    )
  }))
  
  # Calculate summary statistics
  summary_stats <- results_df %>%
    group_by(parameter, x_correlation, y_correlation) %>%
    summarise(
      true_value = first(true_value),
      mean_estimate = mean(estimated_value),
      sd_estimate = sd(estimated_value),
      bias = mean(estimated_value - true_value),
      rmse = sqrt(mean((estimated_value - true_value)^2)),
      coverage = mean(true_value >= ci_lower & true_value <= ci_upper),
      .groups = 'drop'
    )
  
  # Create visualization
  p <- ggplot(results_df, aes(x = true_value, y = estimated_value)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    facet_grid(x_correlation ~ y_correlation, 
               labeller = label_both) +
    labs(title = "Parameter Recovery by Correlation Structure",
         x = "True Value",
         y = "Estimated Value") +
    theme_minimal()
  
  return(list(
    summary_stats = summary_stats,
    plot = p,
    full_results = results_df
  ))
}

# Run the full simulation study
correlation_grid <- c(0.2, 0.4, 0.6, 0.8)
simulation_results <- run_correlation_study(
  n_sims_per_setting = 10,
  correlation_grid = correlation_grid,
  beta_x = c(0.03, -0.01, 0.06),
  beta_y = c(0.04, 0.01, -0.05, 0.02)
)

# Analyze and save results
results_analysis <- analyze_simulation_results(simulation_results)
print(results_analysis$summary_stats)
print(results_analysis$plot) # Display plot

# Save results
saveRDS(simulation_results, "full_simulation_results.rds")
saveRDS(results_analysis, "simulation_analysis.rds")


# Calculate overall performance metrics
overall_performance <- results_analysis$full_results %>%
  group_by(parameter) %>%
  summarise(
    bias = mean(estimated_value - true_value),
    rmse = sqrt(mean((estimated_value - true_value)^2)),
    coverage = mean(true_value >= ci_lower & true_value <= ci_upper)
  )

print("Overall Performance Metrics:")
print(overall_performance)

################################################
# 2. Begin to examine problems
source("sim2.R")

# a. Verify data generation. Debugging function to check simulated data structure
check_simulated_data <- function(data) {
  # Print dimensions
  cat("\nDimensions:\n")
  cat("X grid size:", nrow(data$gridx), "\n")
  cat("Y grid size:", nrow(data$gridy), "\n")
  cat("Number of atoms:", nrow(data$atoms), "\n")
  
  # Check X covariate correlations
  x_covs <- sapply(1:3, function(i) data$gridx[[paste0("covariate_x_", i)]])
  cat("\nX Covariate Correlations:\n")
  print(round(cor(x_covs), 3))
  
  # Check Y covariate correlations
  y_covs <- sapply(1:4, function(i) data$gridy[[paste0("covariate_y_", i)]])
  cat("\nY Covariate Correlations:\n")
  print(round(cor(y_covs), 3))
  
  # Basic statistics for outcome
  cat("\nOutcome (Y) Summary:\n")
  print(summary(data$gridy$y))
  
  # Check relationship between Y and covariates
  cat("\nCorrelations between Y and covariates:\n")
  # For X covariates (need to aggregate to Y grid first)
  x_agg <- matrix(0, nrow(data$gridy), 3)
  for(i in 1:3) {
    for(j in 1:nrow(data$gridy)) {
      atom_indices <- which(data$atoms$ID_y == j)
      x_agg[j,i] <- mean(x_covs[data$atoms$ID_x[atom_indices], i])
    }
  }
  cat("X covariates (aggregated):\n")
  print(round(cor(cbind(Y=data$gridy$y, x_agg))[1,-1], 3))
  
  cat("\nY covariates:\n")
  print(round(cor(cbind(Y=data$gridy$y, y_covs))[1,-1], 3))
  
  return(invisible(NULL))
}

# Run this on a single simulation
test_data <- simulate_misaligned_data(
  res1 = c(5, 5), res2 = c(10, 10),
  n_covariates_x = 3, n_covariates_y = 4,
  x_correlation = 0.6, y_correlation = 0.4,
  beta_x = c(0.5, 0.3, 0.7),
  beta_y = c(0.4, 0.6, 0.3, 0.5)
)
check_simulated_data(test_data)
## Outcome Y is extremely big, the linear predictor is exploding

# 2. Check MCMC model fit
check_model_fit <- function(model_output, true_params) {
  # Calculate ESS and Rhat from chains
  calculate_diagnostics <- function(chains) {
    n_chains <- length(chains)
    n_iter <- nrow(chains[[1]])
    
    # Combine chains for ESS
    combined <- do.call(rbind, chains)
    ess <- coda::effectiveSize(combined)
    
    # Calculate Rhat
    chain_means <- sapply(chains, colMeans)
    chain_vars <- sapply(chains, function(x) apply(x, 2, var))
    W <- mean(chain_vars)  # within-chain variance
    B <- n_iter * var(chain_means)  # between-chain variance
    Rhat <- sqrt((W + B/n_chains)/W)
    
    return(list(ESS = ess, Rhat = Rhat))
  }
  
  # Get diagnostics
  mcmc_diag <- calculate_diagnostics(model_output$samples)
  
  # Print diagnostics
  cat("\nMCMC Diagnostics:\n")
  cat("Effective Sample Sizes:\n")
  print(mcmc_diag$ESS)
  
  cat("\nGelman-Rubin Statistics:\n")
  print(mcmc_diag$Rhat)
  
  # Parameter recovery
  cat("\nParameter Recovery:\n")
  beta_names <- c(paste0("beta_x", 1:length(true_params$beta_x)),
                  paste0("beta_y", 1:length(true_params$beta_y)))
  true_betas <- c(true_params$beta_x, true_params$beta_y)
  
  estimates <- model_output$summary$all.chains[grep("beta_y\\[", rownames(model_output$summary$all.chains)),]
  
  recovery_df <- data.frame(
    Parameter = beta_names,
    True = true_betas,
    Estimate = estimates[,"Mean"],
    CI_Lower = estimates[,"95%CI_low"],
    CI_Upper = estimates[,"95%CI_upp"]
  )
  print(recovery_df)
  
  return(invisible(NULL))
}

# Run a single model fit with diagnostics
test_results <- run_abrm(test_data, "test_output.rds")
check_model_fit(test_results, test_data$true_params)

# 3. Nimble specification
# Modified initialization for better MCMC mixing
modify_nimble_inputs <- function(nimble_inputs) {
  # Adjust initial values
  nimble_inputs$inits$beta_0_y <- log(mean(nimble_inputs$data$y[,ncol(nimble_inputs$data$y)]))
  
  # Initialize betas closer to true values
  p_x <- nimble_inputs$constants$p_x
  p_y <- nimble_inputs$constants$p_y
  nimble_inputs$inits$beta_y <- rep(0.5, p_x + p_y)  # Start from reasonable values
  
  # Adjust spatial precision parameters
  nimble_inputs$inits$tau_x <- rep(1, p_x)
  nimble_inputs$inits$tau_y <- 1
  nimble_inputs$inits$tau_yx <- rep(1, p_y)
  
  return(nimble_inputs)
}

# Test with modified initialization
test_data <- simulate_misaligned_data(
  res1 = c(5, 5),
  res2 = c(10, 10),
  n_covariates_x = 3,
  n_covariates_y = 4,
  x_correlation = 0.6,
  y_correlation = 0.4,
  beta_x = c(0.5, 0.3, 0.7),
  beta_y = c(0.4, 0.6, 0.3, 0.5)
)

nimble_inputs <- prepare_nimble_inputs(
  bookkeeping = prepare_spatial_bookkeeping(test_data),
  adjacency = prepare_adjacency_matrices(
    prepare_spatial_bookkeeping(test_data)$gridy_yorder,
    prepare_spatial_bookkeeping(test_data)$gridx_xorder
  ),
  data = test_data
)

modified_inputs <- modify_nimble_inputs(nimble_inputs)
test_results <- run_nimble_model(
  modified_inputs$constants,
  modified_inputs$data,
  modified_inputs$inits
)

check_model_fit(test_results, test_data$true_params)



# Verify simulations
verify_simulation <- function(sim_data) {
  # Check X covariate correlations
  x_covs <- sapply(1:3, function(i) 
    sim_data$gridx[[paste0("covariate_x_", i)]])
  x_cor <- cor(x_covs)
  
  # Check Y covariate correlations
  y_covs <- sapply(1:4, function(i) 
    sim_data$gridy[[paste0("covariate_y_", i)]])
  y_cor <- cor(y_covs)
  
  # Check relationship between outcome and covariates
  y_outcome <- sim_data$gridy$y
  
  # Print results
  cat("X Covariate Correlations:\n")
  print(round(x_cor, 2))
  cat("\nY Covariate Correlations:\n")
  print(round(y_cor, 2))
  
  # Return correlation matrices for further analysis
  return(list(x_cor = x_cor, y_cor = y_cor))
}
verify_simulation(data)


covar_covid <- readRDS("../data/districts_covid19_icu_and_covariates.rds")
load("../data/final_misaligned_data.RData")


