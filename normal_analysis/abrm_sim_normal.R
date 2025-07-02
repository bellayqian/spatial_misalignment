# ABRM Simulation Script for Normal Distributions

source("function_normal_0610.R")

# Diagnostic function for normal ABRM
run_enhanced_diagnostics_normal <- function(results) {
  mcmc.out <- results$mcmc_output
  sim_data <- results$sim_data
  
  # Extract true parameters
  true_params <- sim_data$true_params
  
  # Calculate relative effective sample size
  calc_rel_ess <- function(ess, n_iterations) {
    return(ess / n_iterations)
  }
  
  # Calculate actual vs nominal effective sample sizes
  n_iterations <- length(mcmc.out$samples[[1]][,1])
  
  # Extract parameter estimates and diagnostics
  params <- rownames(mcmc.out$summary$all.chains)
  beta_params <- params[grep("^beta_y", params)]
  
  diagnostics_df <- data.frame(
    parameter = beta_params,
    true_value = c(true_params$beta_x, true_params$beta_y),
    estimate = mcmc.out[["convergence"]][["diagnostics"]][beta_params, "mean"],
    sd = mcmc.out[["convergence"]][["diagnostics"]][beta_params, "sd"],
    rhat = mcmc.out[["convergence"]][["diagnostics"]][beta_params, "rhat"],
    n_eff = mcmc.out[["convergence"]][["diagnostics"]][beta_params, "ess"],
    rel_ess = calc_rel_ess(mcmc.out[["convergence"]][["diagnostics"]][beta_params, "ess"], n_iterations)
  )
  
  # Add bias and coverage
  diagnostics_df$bias <- diagnostics_df$estimate - diagnostics_df$true_value
  diagnostics_df$rel_bias <- diagnostics_df$bias / diagnostics_df$true_value
  diagnostics_df$within_95ci <- 
    diagnostics_df$true_value >= mcmc.out$summary$all.chains[beta_params, "95%CI_low"] &
    diagnostics_df$true_value <= mcmc.out$summary$all.chains[beta_params, "95%CI_upp"]
  
  cat("\n=== Normal ABRM MCMC Diagnostic Summary ===\n")
  cat("\nConvergence Criteria:\n")
  cat("- Good: Rhat < 1.1, rel_ess > 0.05\n")
  cat("- Acceptable: Rhat < 1.2, rel_ess > 0.03\n")
  cat("- Poor: Rhat >= 1.2 or rel_ess <= 0.03\n\n")
  
  print(diagnostics_df)
  
  # Overall assessment
  cat("\nOverall Assessment:\n")
  cat("Proportion of parameters with good convergence: ", 
      mean(diagnostics_df$rhat < 1.1 & diagnostics_df$rel_ess > 0.05), "\n")
  cat("Coverage rate (% of parameters with true value in 95% CI): ", 
      mean(diagnostics_df$within_95ci) * 100, "%\n")
  
  return(diagnostics_df)
}

# Simulation function for normal ABRM
run_improved_simulation_normal <- function(
    sim_number,
    x_cor, y_cor,
    res1 = c(5, 5),
    res2 = c(10, 10),
    n_covariates_x = 3,
    n_covariates_y = 4,
    mcmc_params = list(
      niter = 30000,
      nburnin = 5000,
      thin = 10,
      nchains = 3),
    beta_x = c(0.5, -0.3, 0.8),
    beta_y = c(0.4, 0.2, -0.6, 0.3),
    seed = NULL, output_dir = NULL
) {
  
  # Generate simulation data for normal distributions
  sim_data <- simulate_misaligned_data(
    res1 = res1,
    res2 = res2,
    n_covariates_x = n_covariates_x,
    n_covariates_y = n_covariates_y,
    x_correlation = x_cor,
    y_correlation = y_cor,
    beta_x = beta_x,
    beta_y = beta_y,
    seed = ifelse(is.null(seed), sim_number, seed)
  )
  
  # Create metadata for this simulation
  sim_metadata <- list(
    sim_number = sim_number,
    x_correlation = x_cor,
    y_correlation = y_cor,
    output_dir = output_dir
  )
  
  # Preprocess data to ensure consistency for normal distributions
  sim_data <- preprocess_normal_data(sim_data)
  
  # Prepare spatial bookkeeping
  bookkeeping <- prepare_spatial_bookkeeping(sim_data)
  
  # Prepare adjacency matrices
  adjacency <- prepare_adjacency_matrices(
    bookkeeping$gridy_yorder,
    bookkeeping$gridx_xorder
  )
  
  # Prepare NIMBLE model inputs
  nimble_inputs <- prepare_nimble_inputs(bookkeeping, adjacency, sim_data)
  
  # Run NIMBLE model
  mcmc.out <- run_nimble_model(
    constants = nimble_inputs$constants,
    data = nimble_inputs$data,
    inits = nimble_inputs$inits,
    sim_metadata = sim_metadata,
    niter = mcmc_params$niter,
    nburnin = mcmc_params$nburnin,
    nchains = mcmc_params$nchains,
    save_plots = FALSE,
    output_dir = output_dir
  )
  
  return(list(
    mcmc_output = mcmc.out,
    sim_data = sim_data,
    sim_number = sim_number,
    x_correlation = x_cor,
    y_correlation = y_cor,
    output_dir = output_dir
  ))
}

# Run multiple simulations with different correlation strengths for normal ABRM
run_correlation_study_normal <- function(n_sims = 10,
                                         correlation_grid = c(0.2, 0.8),
                                         beta_x = c(0.5, -0.3, 0.8),
                                         beta_y = c(0.4, 0.2, -0.6, 0.3),
                                         base_seed = 1) {
  
  results_list <- list()
  counter <- 1
  
  # Create output directory for this run
  date_stamp <- format(Sys.Date(), "%Y%m%d")
  run_dir <- paste0("normal_sim_results_", date_stamp)
  dir.create(run_dir, showWarnings = FALSE)
  
  # Loop through correlation settings
  for(x_cor in correlation_grid) {
    for(y_cor in correlation_grid) {
      cat(sprintf("\nRunning normal ABRM simulations for x_cor = %.2f, y_cor = %.2f\n", x_cor, y_cor))
      
      # Run multiple simulations for each correlation setting
      for(sim in 1:n_sims) {
        cat(sprintf("  Simulation %d of %d\n", sim, n_sims))
        sim_seed <- base_seed + counter
        
        results <- run_improved_simulation_normal(
          sim_number = counter,
          x_cor = x_cor,
          y_cor = y_cor,
          beta_x = beta_x,
          beta_y = beta_y,
          seed = sim_seed,
          output_dir = run_dir
        )
        
        # Save individual simulation results
        saveRDS(results, file = file.path(run_dir, 
                                          sprintf("normal_sim_%d_xcor%.1f_ycor%.1f.rds", 
                                                  sim, x_cor, y_cor)))
        
        # Store in list and run diagnostics
        results_list[[counter]] <- results
        
        counter <- counter + 1
      }
    }
  }
  
  # Save complete results list
  saveRDS(results_list, file = file.path(run_dir, "complete_normal_results.rds"))
  
  # Set the output directory in the return value
  attr(results_list, "output_dir") <- run_dir
  
  return(results_list)
}

# Analyze simulation results for normal ABRM
analyze_simulation_results_normal <- function(results_list, output_prefix = NULL) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(gridExtra)
  require(reshape2)
  
  # Get the output directory from the results list
  output_dir <- attr(results_list, "output_dir")
  if(is.null(output_dir)) {
    output_dir <- "."
  }
  
  # If no output prefix provided, create one with timestamp
  if(is.null(output_prefix)) {
    output_prefix <- paste0("normal_sim_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Extract parameter estimation results
  results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
    x <- results_list[[i]]
    mcmc_summary <- x$mcmc_output$summary$all.chains
    true_params <- x$sim_data$true_params
    
    # Extract only the beta parameters
    beta_params <- grep("^beta_y\\[", rownames(mcmc_summary))
    
    # Create data frame with one row per parameter
    data.frame(
      sim_number = x$sim_number,
      x_correlation = x$x_correlation,
      y_correlation = x$y_correlation,
      parameter = c(paste0("beta_x", 1:length(true_params$beta_x)),
                    paste0("beta_y", 1:length(true_params$beta_y))),
      true_value = c(true_params$beta_x, true_params$beta_y),
      estimated_value = mcmc_summary[beta_params, "Mean"],
      ci_lower = mcmc_summary[beta_params, "95%CI_low"],
      ci_upper = mcmc_summary[beta_params, "95%CI_upp"],
      within_ci = mcmc_summary[beta_params, "95%CI_low"] <= c(true_params$beta_x, true_params$beta_y) & 
        mcmc_summary[beta_params, "95%CI_upp"] >= c(true_params$beta_x, true_params$beta_y),
      ci_width = mcmc_summary[beta_params, "95%CI_upp"] - mcmc_summary[beta_params, "95%CI_low"],
      bias = mcmc_summary[beta_params, "Mean"] - c(true_params$beta_x, true_params$beta_y),
      rel_bias = (mcmc_summary[beta_params, "Mean"] - c(true_params$beta_x, true_params$beta_y)) / 
        c(true_params$beta_x, true_params$beta_y) * 100
    )
  }))
  
  # Aggregate results by correlation combination and parameter
  parameter_summary <- results_df %>%
    group_by(x_correlation, y_correlation, parameter) %>%
    summarise(
      n_sims = n(),
      true_value = first(true_value),
      mean_estimate = mean(estimated_value),
      median_estimate = median(estimated_value),
      sd_estimate = sd(estimated_value),
      mean_bias = mean(bias),
      mean_rel_bias = mean(rel_bias),
      rmse = sqrt(mean(bias^2)),
      mean_ci_width = mean(ci_width),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Summary by correlation combination
  correlation_summary <- results_df %>%
    group_by(x_correlation, y_correlation) %>%
    summarise(
      n_params = n_distinct(parameter),
      n_sims = n_distinct(sim_number),
      overall_coverage = mean(within_ci) * 100,
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(rel_bias)),
      overall_rmse = sqrt(mean(bias^2)),
      avg_ci_width = mean(ci_width),
      .groups = 'drop'
    )
  
  # Save overall summary statistics
  write.csv(parameter_summary, 
            file = file.path(output_dir, paste0(output_prefix, "_parameter_summary.csv")), 
            row.names = FALSE)
  
  write.csv(correlation_summary, 
            file = file.path(output_dir, paste0(output_prefix, "_correlation_summary.csv")), 
            row.names = FALSE)
  
  return(list(
    parameter_summary = parameter_summary,
    correlation_summary = correlation_summary,
    full_results = results_df,
    output_dir = output_dir
  ))
}

# Analyze simulation results with simpler parametric approach for normal ABRM
analyze_simulation_results_simpler_parametric_normal <- function(results_list, output_prefix = NULL) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(gridExtra)
  require(reshape2)
  
  # Get the output directory from the results list
  output_dir <- attr(results_list, "output_dir")
  if(is.null(output_dir)) {
    output_dir <- "."  # Use current directory
  }
  
  # If no output prefix provided, create one with timestamp
  if(is.null(output_prefix)) {
    output_prefix <- paste0("normal_sim_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  # Extract parameter estimation results
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    mcmc_summary <- x$mcmc_output$summary$all.chains
    true_params <- x$sim_data$true_params
    
    data.frame(
      sim_number = x$sim_number,
      x_correlation = x$x_correlation,
      y_correlation = x$y_correlation,
      parameter = c(paste0("beta_x", 1:length(true_params$beta_x)),
                    paste0("beta_y", 1:length(true_params$beta_y))),
      true_value = c(true_params$beta_x, true_params$beta_y),
      estimated_value = mcmc_summary[grep("^beta_y\\[", rownames(mcmc_summary)), "Mean"],
      ci_lower = mcmc_summary[grep("^beta_y\\[", rownames(mcmc_summary)), "95%CI_low"],
      ci_upper = mcmc_summary[grep("^beta_y\\[", rownames(mcmc_summary)), "95%CI_upp"]
    )
  }))
  
  # Create comparison data separated by correlation combinations
  comparison_df <- results_df %>%
    group_by(x_correlation, y_correlation, parameter) %>%
    summarise(
      estimated_beta = mean(estimated_value),
      se_beta = sd(estimated_value) / sqrt(n()),
      true_beta = first(true_value),
      # Calculate CI using t-distribution
      ci_lower = estimated_beta - qt(0.975, n() - 1) * se_beta,
      ci_upper = estimated_beta + qt(0.975, n() - 1) * se_beta,
      .groups = 'drop'
    )
  
  # Create coefficient comparison plot for each correlation combination
  comparison_plots <- list()
  
  # Get unique correlation combinations
  corr_combos <- unique(comparison_df[, c("x_correlation", "y_correlation")])
  
  # Create a combined plot showing all combinations in facets
  p_combined <- ggplot(comparison_df, aes(x = parameter)) +
    geom_point(aes(y = estimated_beta), color = "blue") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "blue") +
    geom_point(aes(y = true_beta), color = "red", shape = 4, size = 3) +
    facet_grid(x_correlation ~ y_correlation, labeller = label_both) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal ABRM: True vs Estimated Coefficients by Correlation Structure",
         subtitle = "Blue points/bars = Estimated values with 95% CI\nRed X = True values",
         x = "Parameter",
         y = "Coefficient Value")
  
  # Create full file paths with the output directory
  pdf_file <- file.path(output_dir, paste0(output_prefix, "_plots.pdf"))
  csv_file <- file.path(output_dir, paste0(output_prefix, "_parameter_summary.csv"))
  
  # Save all results
  # Create single PDF with all plots
  pdf(pdf_file, width = 12, height = 8)
  print(p_combined)
  
  # Print MCMC diagnostics for each simulation
  for(i in seq_along(results_list)) {
    sim_data <- results_list[[i]]
    metadata <- list(
      sim_number = sim_data$sim_number,
      x_correlation = sim_data$x_correlation,
      y_correlation = sim_data$y_correlation
    )
    
    diag_plots <- create_diagnostic_plots(
      sim_data$mcmc_output$samples,
      metadata
    )
    
    print(diag_plots$trace)
    print(diag_plots$density)
  }
  
  dev.off()
  cat(sprintf("Normal ABRM plots saved to %s\n", pdf_file))
  
  # Save summary statistics
  write.csv(comparison_df, file = csv_file, row.names = FALSE)
  cat(sprintf("Normal ABRM summary statistics saved to %s\n", csv_file))
  
  return(list(
    summary_stats = comparison_df,
    comparison_plot_combined = p_combined,
    comparison_plots = comparison_plots,
    full_results = results_df,
    output_dir = output_dir,
    output_files = list(
      pdf = pdf_file,
      csv = csv_file
    )
  ))
}

# Define simulation parameters for normal distributions
DEFAULT_PARAMS_NORMAL <- list(
  res1 = c(5, 5),
  res2 = c(10, 10),
  n_covariates_x = 3,
  n_covariates_y = 4,
  mcmc_params = list(
    niter = 30000,
    nburnin = 5000,
    thin = 10,
    nchains = 3),
  beta_x = c(0.5, -0.3, 0.8),
  beta_y = c(0.4, 0.2, -0.6, 0.3)
)

# Run the full simulation study for normal ABRM
correlation_grid <- c(0.2, 0.6)
date_stamp <- format(Sys.Date(), "%Y%m%d")
output_prefix <- paste0("normal_corr", paste(correlation_grid, collapse="_"), "_", date_stamp)

# Run simulation
simulation_results_normal <- run_correlation_study_normal(
  n_sims = 3,
  correlation_grid = correlation_grid,
  beta_x = DEFAULT_PARAMS_NORMAL$beta_x,
  beta_y = DEFAULT_PARAMS_NORMAL$beta_y
)

# Analyze and save results
# results_analysis_normal <- analyze_simulation_results_normal(simulation_results_normal, output_prefix)
