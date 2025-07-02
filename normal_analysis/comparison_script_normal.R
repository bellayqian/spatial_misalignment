# Comparison of ABRM vs Dasymetric Mapping for Normal Data
# This script compares the results of two different approaches:
# 1. ABRM (Atom-based Regression Model) for normal distributions
# 2. Dasymetric mapping for normal distributions

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

source("function_normal_0610.R")
source("Dasymetric_mapping_normal.R")

set.seed(42)

# Function to run both methods and compare results
compare_methods_normal <- function(
    res1 = c(5, 5), 
    res2 = c(10, 10),
    n_covariates_x = 3,
    n_covariates_y = 4,
    x_correlation = 0.2,
    y_correlation = 0.6,
    beta_x = c(0.5, -0.3, 0.8),
    beta_y = c(0.4, 0.2, -0.6, 0.3),
    seed = 123,
    mcmc_params = list(
      niter = 30000,
      nburnin = 5000,
      thin = 10,
      nchains = 3)
) {
  cat("Generating simulated normal data...\n")
  
  # Create single simulated dataset for both methods
  sim_data <- simulate_misaligned_data(
    res1 = res1,
    res2 = res2,
    seed = seed,
    n_covariates_x = n_covariates_x,
    n_covariates_y = n_covariates_y,
    x_correlation = x_correlation,
    y_correlation = y_correlation,
    beta_x = beta_x,
    beta_y = beta_y
  )
  
  # Create output directory for results
  date_stamp <- format(Sys.Date(), "%Y%m%d")
  output_dir <- paste0("comparison_results_normal_", date_stamp)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Create metadata
  sim_metadata <- list(
    sim_number = 1,
    x_correlation = x_correlation,
    y_correlation = y_correlation,
    output_dir = output_dir
  )
  
  # Run ABRM method
  cat("\nRunning ABRM normal method...\n")
  
  # Preprocess data to ensure consistency
  sim_data <- preprocess_normal_data(sim_data)
  
  # Run the full ABRM pipeline
  abrm_results <- run_abrm(sim_data, 
                           file.path(output_dir, 'normal_abrm_output.rds'), 
                           sim_metadata)
  
  # Run Dasymetric mapping method
  cat("\nRunning Dasymetric mapping method (normal)...\n")
  mapped_data <- dasymetric_mapping_normal(sim_data)
  dasymetric_model <- fit_normal_models(mapped_data)
  
  # Get comparison with true values for Dasymetric approach
  dasymetric_comparison <- compare_true_vs_estimated_normal(
    dasymetric_model,
    true_beta_x = beta_x,
    true_beta_y = beta_y
  )
  
  # Extract ABRM parameter estimates
  abrm_parameters <- data.frame(
    variable = rownames(abrm_results$summary$all.chains),
    estimated_beta = abrm_results$summary$all.chains[, "Mean"],
    std_error = abrm_results$summary$all.chains[, "St.Dev."],
    ci_lower = abrm_results$summary$all.chains[, "95%CI_low"],
    ci_upper = abrm_results$summary$all.chains[, "95%CI_upp"]
  )
  
  # Filter beta parameters for covariates (beta_y contains effects of both X and Y covariates)
  beta_params <- grep("^beta_y\\[", abrm_parameters$variable)
  abrm_betas <- abrm_parameters[beta_params, ]
  
  # Rename to match with Dasymetric comparison
  true_values <- c(beta_x, beta_y)
  abrm_betas$true_beta <- true_values
  abrm_betas$variable <- c(paste0("covariate_x_", 1:length(beta_x)),
                           paste0("covariate_y_", 1:length(beta_y)))
  
  # Calculate bias metrics
  abrm_betas$bias <- abrm_betas$estimated_beta - abrm_betas$true_beta
  abrm_betas$relative_bias <- ((abrm_betas$estimated_beta - abrm_betas$true_beta) / abrm_betas$true_beta) * 100
  
  # Check if true value falls within CI
  abrm_betas$within_ci <- abrm_betas$true_beta >= abrm_betas$ci_lower & 
    abrm_betas$true_beta <= abrm_betas$ci_upper
  
  # Create combined comparison dataframe
  combined_comparison <- rbind(
    cbind(method = "ABRM", abrm_betas),
    cbind(method = "Dasymetric", dasymetric_comparison)
  )
  
  combined_comparison$x_correlation <- x_correlation
  combined_comparison$y_correlation <- y_correlation
  
  # Save combined comparison
  write.csv(combined_comparison, file = file.path(output_dir, "normal_method_comparison.csv"), row.names = FALSE)
  
  # Create comparison plots
  create_comparison_plots_normal(combined_comparison, output_dir)
  
  # Create summary statistics
  create_summary_statistics_normal(combined_comparison, output_dir)
  
  # Return results
  return(list(
    sim_data = sim_data,
    abrm_results = abrm_results,
    dasymetric_results = dasymetric_model,
    combined_comparison = combined_comparison
  ))
}

# Function to create comparison plots
create_comparison_plots_normal <- function(comparison_data, output_dir) {
  # Pivot the data for side by side comparison
  plot_data <- comparison_data %>%
    dplyr::select(method, variable, estimated_beta, ci_lower, ci_upper, true_beta) %>%
    mutate(variable = factor(variable, levels = unique(variable)))
  
  # Create coefficient plot
  coef_plot <- ggplot(plot_data, aes(x = variable, y = estimated_beta, color = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(aes(y = true_beta), color = "black", shape = 4, size = 3) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: Comparison of Coefficient Estimates by Method",
         subtitle = "Points = Estimated values with 95% CI\nX = True values",
         x = "Variable",
         y = "Coefficient Value")
  
  # Create bias plot
  bias_data <- comparison_data %>%
    dplyr::select(method, variable, relative_bias) %>%
    mutate(variable = factor(variable, levels = unique(variable)))
  
  bias_plot <- ggplot(bias_data, aes(x = variable, y = relative_bias, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: Relative Bias by Method (%)",
         x = "Variable",
         y = "Relative Bias (%)")
  
  # Create coverage plot (within_ci)
  coverage_data <- comparison_data %>%
    dplyr::select(method, variable, within_ci) %>%
    mutate(variable = factor(variable, levels = unique(variable)),
           covered = ifelse(within_ci, "Yes", "No"))
  
  coverage_plot <- ggplot(coverage_data, aes(x = variable, fill = covered)) +
    geom_bar(position = "dodge") +
    facet_wrap(~method) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: 95% CI Coverage by Method",
         subtitle = "Whether true parameter value falls within estimated 95% CI",
         x = "Variable",
         y = "Count")
  
  # Save plots
  pdf(file.path(output_dir, "normal_coefficient_comparison.pdf"), width = 10, height = 8)
  print(coef_plot)
  print(bias_plot)
  print(coverage_plot)
  dev.off()
  
  cat("Normal comparison plots saved to", file.path(output_dir, "normal_coefficient_comparison.pdf"), "\n")
}

# Function to create summary statistics
create_summary_statistics_normal <- function(comparison_data, output_dir) {
  # Calculate summary statistics by method
  summary_stats <- comparison_data %>%
    group_by(method) %>%
    summarize(
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Save summary
  write.csv(summary_stats, file = file.path(output_dir, "normal_method_summary_stats.csv"), row.names = FALSE)
  
  # Create a summary table for each parameter
  param_summary <- comparison_data %>%
    dplyr::select(method, variable, estimated_beta, true_beta, bias, relative_bias, within_ci) %>%
    pivot_wider(
      id_cols = variable,
      names_from = method,
      values_from = c(estimated_beta, bias, relative_bias, within_ci),
      names_sep = "_"
    )
  
  # Save parameter summary
  write.csv(param_summary, file = file.path(output_dir, "normal_parameter_summary.csv"), row.names = FALSE)
  
  # Print summary to console
  cat("\nNormal Method Comparison Summary:\n")
  print(summary_stats)
  
  cat("\nNormal Parameter-Specific Comparison:\n")
  print(param_summary)
  
  return(list(method_summary = summary_stats, param_summary = param_summary))
}

# Function to create summary plots for sensitivity analysis
create_sensitivity_summary_plots_normal <- function(combined_results, output_dir) {
  # Add correlation values to the combined results
  plot_data <- combined_results %>%
    group_by(method, variable, x_correlation, y_correlation) %>%
    summarize(
      mean_estimate = mean(estimated_beta),
      mean_lower = mean(ci_lower),
      mean_upper = mean(ci_upper),
      true_value = mean(true_beta),  # This should be the same for all
      mean_bias = mean(bias),
      mean_rel_bias = mean(relative_bias),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Create faceted plot by correlation values
  corr_plot <- ggplot(plot_data, aes(x = variable, y = mean_estimate, color = method)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = mean_lower, ymax = mean_upper), 
                  position = position_dodge(width = 0.5), width = 0.2) +
    geom_point(aes(y = true_value), color = "black", shape = 4, size = 3) +
    facet_grid(x_correlation ~ y_correlation, labeller = label_both) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: Parameter Estimates by Correlation Structure",
         subtitle = "Averaged across all simulations",
         x = "Variable",
         y = "Coefficient Value")
  
  # Create bias plot by correlation
  bias_plot <- ggplot(plot_data, aes(x = variable, y = mean_rel_bias, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(x_correlation ~ y_correlation, labeller = label_both) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: Relative Bias by Correlation Structure (%)",
         subtitle = "Averaged across all simulations",
         x = "Variable",
         y = "Relative Bias (%)")
  
  # Create coverage rate plot
  coverage_plot <- ggplot(plot_data, aes(x = variable, y = coverage_rate, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    facet_grid(x_correlation ~ y_correlation, labeller = label_both) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Normal Model: 95% CI Coverage Rate by Correlation Structure",
         subtitle = "Percentage of simulations where true value falls within CI",
         x = "Variable", 
         y = "Coverage Rate (%)")
  
  # Save plots
  pdf(file.path(output_dir, "normal_sensitivity_summary_plots.pdf"), width = 12, height = 10)
  print(corr_plot)
  print(bias_plot)
  print(coverage_plot)
  dev.off()
  
  cat("Normal sensitivity summary plots saved to", file.path(output_dir, "normal_sensitivity_summary_plots.pdf"), "\n")
}

# Function to run a parameter sensitivity analysis for normal data
run_sensitivity_analysis_normal <- function(
    correlation_grid = c(0.2, 0.6),
    n_sims = 10,
    base_seed = 123
) {
  # Create output directory
  date_stamp <- format(Sys.time(), "%Y%m%d")
  output_dir <- paste0("normal_sensitivity_analysis_", date_stamp)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Initialize results storage
  all_results <- list()
  counter <- 1
  
  # Run through correlation grid
  for(x_cor in correlation_grid) {
    for(y_cor in correlation_grid) {
      cat(sprintf("\nRunning normal analysis for x_cor = %.2f, y_cor = %.2f\n", x_cor, y_cor))
      
      # Generate multiple datasets upfront for this correlation combination
      sim_datasets <- list()
      for(sim in 1:n_sims) {
        sim_seed <- base_seed + ((sim-1) * 100) + counter
        
        cat(sprintf("  Generating dataset %d of %d (seed: %d)\n", sim, n_sims, sim_seed))
        
        # Generate the dataset
        sim_datasets[[sim]] <- simulate_misaligned_data(
          res1 = c(5, 5),
          res2 = c(10, 10),
          seed = sim_seed,
          n_covariates_x = 3,
          n_covariates_y = 4,
          x_correlation = x_cor,
          y_correlation = y_cor,
          beta_x = c(0.5, -0.3, 0.8),
          beta_y = c(0.4, 0.2, -0.6, 0.3)
        )
        
        # Store the true parameter values
        sim_datasets[[sim]]$true_params <- list(
          beta_x = c(0.5, -0.3, 0.8),
          beta_y = c(0.4, 0.2, -0.6, 0.3)
        )
      }
      
      # Now run both methods on each dataset
      for(sim in 1:n_sims) {
        sim_data <- sim_datasets[[sim]]
        cat(sprintf("  Analyzing dataset %d of %d\n", sim, n_sims))
        
        # Create sub-directory for this simulation
        sim_dir <- file.path(output_dir, sprintf("xcor%.1f_ycor%.1f_sim%d", x_cor, y_cor, sim))
        dir.create(sim_dir, showWarnings = FALSE)
        
        # Save the simulated dataset
        saveRDS(sim_data, file = file.path(sim_dir, "simulated_normal_data.rds"))
        
        # Run ABRM method
        cat("    Running ABRM normal method...\n")
        
        # Create metadata
        sim_metadata <- list(
          sim_number = sim,
          x_correlation = x_cor,
          y_correlation = y_cor,
          output_dir = sim_dir
        )
        
        # Preprocess data
        sim_data <- preprocess_normal_data(sim_data)
        
        # Run normal ABRM
        abrm_results <- run_abrm(sim_data, 
                                 file.path(sim_dir, 'normal_abrm_output.rds'), 
                                 sim_metadata)
        
        # Run Dasymetric mapping method
        cat("    Running Dasymetric mapping method (normal)...\n")
        mapped_data <- dasymetric_mapping_normal(sim_data)
        dasymetric_model <- fit_normal_models(mapped_data)
        
        # Get comparison with true values for Dasymetric approach
        true_beta_x <- sim_data$true_params$beta_x
        true_beta_y <- sim_data$true_params$beta_y
        
        dasymetric_comparison <- compare_true_vs_estimated_normal(
          dasymetric_model,
          true_beta_x = true_beta_x,
          true_beta_y = true_beta_y
        )
        
        # Extract ABRM parameter estimates
        abrm_parameters <- data.frame(
          variable = rownames(abrm_results$summary$all.chains),
          estimated_beta = abrm_results$summary$all.chains[, "Mean"],
          std_error = abrm_results$summary$all.chains[, "St.Dev."],
          ci_lower = abrm_results$summary$all.chains[, "95%CI_low"],
          ci_upper = abrm_results$summary$all.chains[, "95%CI_upp"]
        )
        
        # Filter beta parameters
        beta_params <- grep("^beta_y\\[", abrm_parameters$variable)
        abrm_betas <- abrm_parameters[beta_params, ]
        
        # Rename to match with Dasymetric comparison
        true_values <- c(true_beta_x, true_beta_y)
        abrm_betas$true_beta <- true_values
        abrm_betas$variable <- c(paste0("covariate_x_", 1:length(true_beta_x)),
                                 paste0("covariate_y_", 1:length(true_beta_y)))
        
        # Calculate bias metrics
        abrm_betas$bias <- abrm_betas$estimated_beta - abrm_betas$true_beta
        abrm_betas$relative_bias <- ((abrm_betas$estimated_beta - abrm_betas$true_beta) / abrm_betas$true_beta) * 100
        
        # Check if true value falls within CI
        abrm_betas$within_ci <- abrm_betas$true_beta >= abrm_betas$ci_lower & 
          abrm_betas$true_beta <= abrm_betas$ci_upper
        
        # Add simulation metadata
        abrm_betas$sim_number <- sim
        abrm_betas$x_correlation <- x_cor
        abrm_betas$y_correlation <- y_cor
        
        # Add method label
        abrm_betas$method <- "ABRM_Normal"
        
        # Add similar metadata to dasymetric results
        dasymetric_comparison$sim_number <- sim
        dasymetric_comparison$x_correlation <- x_cor
        dasymetric_comparison$y_correlation <- y_cor
        dasymetric_comparison$method <- "Dasymetric_Normal"
        
        # Create combined comparison for this simulation
        sim_comparison <- rbind(
          abrm_betas,
          dasymetric_comparison
        )
        
        # Save combined comparison for this simulation
        write.csv(sim_comparison, 
                  file = file.path(sim_dir, "normal_method_comparison.csv"), 
                  row.names = FALSE)
        
        # Create comparison plots for this simulation
        create_comparison_plots_normal(sim_comparison, sim_dir)
        
        # Add to overall results
        all_results[[counter]] <- sim_comparison
        counter <- counter + 1
      }
    }
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  
  # Create summary plots across all simulations
  create_sensitivity_summary_plots_normal(combined_results, output_dir)
  
  # Save combined results
  write.csv(combined_results, 
            file = file.path(output_dir, "normal_sensitivity_analysis_results.csv"), 
            row.names = FALSE)
  
  # Aggregate results by correlation, method, and variable
  summary_by_correlation_var <- combined_results %>%
    group_by(method, variable, x_correlation, y_correlation) %>%
    summarize(
      mean_estimate = mean(estimated_beta),
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Aggregate results by correlation and method
  summary_by_correlation <- combined_results %>%
    group_by(method, x_correlation, y_correlation) %>%
    summarize(
      mean_abs_bias = mean(abs(bias)),
      mean_rel_bias = mean(abs(relative_bias)),
      rmse = sqrt(mean(bias^2)),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Save summaries
  write.csv(summary_by_correlation_var, 
            file = file.path(output_dir, "normal_sensitivity_summary_by_variable.csv"), 
            row.names = FALSE)
  
  write.csv(summary_by_correlation, 
            file = file.path(output_dir, "normal_sensitivity_summary.csv"), 
            row.names = FALSE)
  
  # Create plot of method performance by correlation structure
  corr_plot <- ggplot(summary_by_correlation, 
                      aes(x = factor(x_correlation), y = factor(y_correlation), fill = rmse)) +
    geom_tile() +
    facet_wrap(~method) +
    scale_fill_viridis_c(direction = -1) +
    theme_minimal() +
    labs(title = "Normal Method Performance by Correlation Structure",
         subtitle = "Root Mean Square Error (RMSE)",
         x = "X Correlation",
         y = "Y Correlation")
  
  ggsave(file.path(output_dir, "normal_correlation_performance_heatmap.pdf"), 
         plot = corr_plot, 
         width = 10, height = 8)
  
  cat("\nNormal sensitivity analysis complete. Results saved to:", output_dir, "\n")
  
  return(list(
    combined_results = combined_results,
    summary_by_correlation = summary_by_correlation,
    summary_by_correlation_var = summary_by_correlation_var
  ))
}

# Single comparison
# single_comparison_normal <- compare_methods_normal(
#   res1 = c(5, 5),
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.3,
#   y_correlation = 0.7,
#   beta_x = c(0.5, -0.3, 0.8),
#   beta_y = c(0.4, 0.2, -0.6, 0.3),
#   seed = 123
# )

# Sensitivity analysis
sensitivity_results_normal <- run_sensitivity_analysis_normal(
  correlation_grid = c(0.2, 0.6),
  n_sims = 5,
  base_seed = 123
)