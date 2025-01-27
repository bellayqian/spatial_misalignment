# Load required packages
library(sp)
library(sf)
library(tidyverse)
library(ggplot2)
library(spdep)
library(MASS)
library(raster)

####### Data Generation ####### 
## All beta coefficients are close to 0 with confidence intervals centered around 0
## The precision matrices show minimal correlation between covariates
## The spatial parameters (tau) show some instability in estimation
# x_correlation = 0.6,  # Correlation between X covariates
# y_correlation = 0.4,  # Correlation between Y covariates
# beta_x = c(0.5, 0.3, 0.7),  # True effects of X covariates
# beta_y = c(0.4, 0.6, 0.3, 0.5)  # True effects of Y covariates

gen_correlated_spat <- function(W, n_vars, rho = 0.6, var_spat = 1, correlation = 0.5, 
                                global_ints = NULL, verify = TRUE) {
  n <- nrow(W)
  
  # Create precision matrix for variables (Lambda)
  Sigma_vars <- matrix(correlation, n_vars, n_vars) # Correlation Matrix 
  diag(Sigma_vars) <- 1
  Lambda <- solve(Sigma_vars) * (1/var_spat)
  
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
    int <- ifelse(is.null(global_ints), 3, global_ints[j])
    observed_values[,j] <- rpois(n, lambda = exp(int + spatial_effects[,j]))
  }
  
  return(observed_values)
}

# Main simulation function
simulate_misaligned_data <- function(res1 = c(5, 5), res2 = c(10, 10), 
                                     seed = 1, n_covariates_x = 3, n_covariates_y = 4,
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
  
  # 2-sim_misaligned_data.R
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
  
  # Generate outcome Y with true associations
  if(is.null(beta_x)) beta_x <- rnorm(n_covariates_x, 0.5, 0.1)
  if(is.null(beta_y)) beta_y <- rnorm(n_covariates_y, 0.5, 0.1) 

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
  linear_pred_y <- rep(3, nrow(gridy))  # Intercept
  
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

data <- simulate_misaligned_data(
  res1 = c(5, 5),
  res2 = c(10, 10),
  n_covariates_x = 3,
  n_covariates_y = 4,
  x_correlation = 0.6,  # Correlation between X covariates
  y_correlation = 0.4,  # Correlation between Y covariates
  beta_x = c(0.5, 0.3, 0.7),  # True effects of X covariates
  beta_y = c(0.4, 0.6, 0.3, 0.5)  # True effects of Y covariates
)
