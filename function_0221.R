# Load required libraries
library(sp)
library(sf)         # Simple feature for R
library(spdep)      # A package to calculate neighbors
library(MASS)
library(raster)
library(tidyverse)
library(nimble)     # A package for performing MCMC in R
library(ggplot2)
library(coda)

####
# Simulation part
gen_correlated_spat <- function(W, n_vars, rho = 0.6, var_spat = 1, correlation = 0.5, 
                                global_ints = NULL, verify = FALSE) {
  n <- nrow(W)
  
  # # BQ: Grid size adjustment
  # base_grid_size <- res1[1] * res1[2]
  # grid_factor <- sqrt(n/base_grid_size) # grid size scaling factor
  # # BQ: First analyze spatial structure
  # eigenvalues <- eigen(W)$values
  # max_rho <- 1/max(abs(eigenvalues))
  # # BQ: Adjust both rho and correlation
  # spatial_factor <- (mean(rowSums(W)) / n) * grid_factor
  # # BQ: Adjust rho to maintain valid spatial structure
  # adjusted_rho <- sign(rho) * min(abs(rho), 0.99 * max_rho)
  # # BQ: Adjust correlation to account for spatial effects
  # adjusted_correlation <- correlation / (1 + spatial_factor * abs(adjusted_rho))
  
  # Create precision matrix for variables (Lambda)
  Sigma_vars <- matrix(correlation, n_vars, n_vars) # Correlation Matrix 
  diag(Sigma_vars) <- 1
  Lambda <- solve(Sigma_vars) #* (1/var_spat)
  
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

# Spatial Bookkeeping part
create_adjacency_matrix <- function(grid) {
  neighbors <- spdep::poly2nb(as(grid, "Spatial"), queen = TRUE)
  W <- spdep::nb2mat(neighbors, style = "B", zero.policy = TRUE)
  return(W)
}

# Function to prepare spatial bookkeeping
prepare_spatial_bookkeeping <- function(data) {
  gridy <- data$gridy
  gridx <- data$gridx
  atoms <- data$atoms
  S_x <- nrow(data$gridx)
  S_y <- nrow(data$gridy)
  
  # Ensure IDs are present
  gridx$ID <- 1:nrow(gridx)
  gridy$ID <- 1:nrow(gridy)
  
  # Identify non-atom grids
  ## x grids that aren't atoms
  # x_vars <- c("x", grep("covariate_x_", names(gridx), value = TRUE))
  x_vars <- c(grep("covariate_x_", names(gridx), value = TRUE))
  x_nonatoms <- atoms$ID_x[which(duplicated(atoms$ID_x)==T)]
  x_nonatoms <- x_nonatoms[order(x_nonatoms)]
  J_x <- nrow(gridx) - length(x_nonatoms)
  
  ## y grids that aren't atoms
  y_vars <- c("y", grep("covariate_y_", names(gridy), value = TRUE))
  y_nonatoms <- atoms$ID_y[which(duplicated(atoms$ID_y)==T)]
  y_nonatoms <- y_nonatoms[order(y_nonatoms)]
  J_y <- nrow(gridy) - length(y_nonatoms)
  
  # Re-order x and y so that first J_x/J_y units are atom-equivalents 
  gridy_yorder <- gridy[c(gridy$ID[-y_nonatoms], y_nonatoms),]
  names(gridy_yorder)[which(names(gridy_yorder)=='ID')] <- 'ID_y'
  gridy_yorder$ID_yorder <- 1:nrow(gridy_yorder)
  
  gridx_xorder <- gridx[c(gridx$ID[-x_nonatoms], x_nonatoms),]
  names(gridx_xorder)[which(names(gridx_xorder)=='ID')] <- 'ID_x'
  gridx_xorder$ID_xorder <- 1:nrow(gridx_xorder)
  
  # Prepare ID crosswalk. Link the new IDs into the atoms dataset so we can find 
  # correspondence between new Y grid ids and atoms
  IDxwalk <- merge(st_drop_geometry(atoms), st_drop_geometry(gridy_yorder), by='ID_y')
  IDxwalk <- merge(IDxwalk, st_drop_geometry(gridx_xorder), by='ID_x')
  
  # Prepare mapping vectors, map atoms to grids
  ## get x/y_to_atom
  x_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder)]
  y_to_atom <- IDxwalk$ID_atomorder[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder)]
  
  ## get expand_x/y
  expand_x <- IDxwalk$ID_xorder[order(IDxwalk$ID_xorder)]
  expand_y <- IDxwalk$ID_yorder[order(IDxwalk$ID_yorder)]
  
  # Prepare latent indices
  ## get x/y_latentind
  IDxwalk <- IDxwalk[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder),]
  x_nonatoms <- IDxwalk$ID_xorder[which(duplicated(IDxwalk$ID_xorder)==T)]
  x_latentind <- matrix(0, length(x_nonatoms), 2)
  for (i in 1:length(x_nonatoms)) {
    x_latentind[i,] <- c(min(which(IDxwalk$ID_xorder==x_nonatoms[i])), max(which(IDxwalk$ID_xorder==x_nonatoms[i]))) - J_x
  }
  
  IDxwalk <- IDxwalk[order(IDxwalk$ID_yorder, IDxwalk$ID_atomorder),]
  y_nonatoms <- IDxwalk$ID_yorder[which(duplicated(IDxwalk$ID_yorder)==T)]
  y_latentind <- matrix(0, length(y_nonatoms), 2)
  for (i in 1:length(y_nonatoms)) {
    y_latentind[i,] <- c(min(which(IDxwalk$ID_yorder==y_nonatoms[i])), max(which(IDxwalk$ID_yorder==y_nonatoms[i]))) - J_y
  }
  
  # Prepare x_reorder
  IDxwalk <- IDxwalk[order(IDxwalk$ID_xorder, IDxwalk$ID_atomorder),]
  IDxwalk$xro <- 1:nrow(IDxwalk)
  IDxwalk <- IDxwalk[order(IDxwalk$ID_atomorder),]
  x_reorder <- IDxwalk$xro
  
  return(list(J_x = J_x, J_y = J_y, x_to_atom = x_to_atom, y_to_atom = y_to_atom,
              expand_x = expand_x, expand_y = expand_y, x_latentind = x_latentind,
              y_latentind = y_latentind, x_reorder = x_reorder,
              gridy_yorder = gridy_yorder, gridx_xorder = gridx_xorder,
              x_vars = x_vars, y_vars = y_vars))
}

# Function to prepare adjacency matrices
prepare_adjacency_matrices <- function(gridy_yorder, gridx_xorder) {
  ## Make spatial adjacency matrices for X and Y grids based on new orderings
  neighbors_x <- poly2nb(as(gridx_xorder, "Spatial"), queen = TRUE)
  W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)
  
  neighbors_y <- poly2nb(as(gridy_yorder, "Spatial"), queen = TRUE)
  W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)
  
  return(list(W_x = W_x, W_y = W_y))
}

# Function to prepare NIMBLE model inputs
prepare_nimble_inputs <- function(bookkeeping, adjacency, data) {
  # Basic dimensions 
  gridy_yorder <- bookkeeping$gridy_yorder
  gridx_xorder <- bookkeeping$gridx_xorder
  atoms <- data$atoms
  n_atoms <- if(is.null(nrow(atoms))) length(atoms) else nrow(atoms)
  
  # Identify variables
  x_vars <- c("x", grep("covariate_x_", names(gridx_xorder), value = TRUE))
  y_vars <- c("y", grep("covariate_y_", names(gridy_yorder), value = TRUE))
  
  # Create covariate matrices
  create_covariate_matrix <- function(grid_data, vars, atom_ids, grid_size) {
    if("ID_x" %in% names(grid_data)) {
      # For X grid
      covar_matrix <- matrix(NA, nrow = grid_size, ncol = length(vars))
      colnames(covar_matrix) <- vars
      
      # Fill in values
      for(i in 1:length(vars)) {
        covar_matrix[,i] <- as.numeric(grid_data[[vars[i]]])
      }
    } else {
      # For Y grid - p_y covariates first, then y
      p_y <- length(vars) - 1  # subtract 1 because 'y' is included in vars
      covar_matrix <- matrix(NA, nrow = grid_size, ncol = p_y + 1)
      
      # Get covariate names (excluding 'y')
      covar_names <- vars[vars != "y"]
      
      # Fill in covariates first
      for(i in 1:p_y) {
        covar_matrix[,i] <- as.numeric(grid_data[[covar_names[i]]])
      }
      
      # Fill in y as the last column
      covar_matrix[,p_y + 1] <- as.numeric(grid_data[["y"]])
      colnames(covar_matrix) <- c(covar_names, "y")
    }
    
    return(covar_matrix)
  }
  
  print("\nCreating and standardizing covariate matrices...")
  covar_x <- create_covariate_matrix(grid_data=gridx_xorder[order(gridx_xorder$ID_x), ], 
                                     vars=bookkeeping$x_vars, 
                                     atom_ids=unique(data$atoms$ID_x), 
                                     grid_size=nrow(gridx_xorder))
  
  covar_y <- create_covariate_matrix(gridy_yorder[order(gridy_yorder$ID_y), ], 
                                     bookkeeping$y_vars, 
                                     unique(data$atoms$ID_y), 
                                     nrow(gridy_yorder))
  
  # Create index matrices
  p_x <- ncol(covar_x)
  p_y <- ncol(covar_y) - 1  # subtract y column
  
  # Scale matrices for Wishart priors
  R_x <- diag(p_x)  # For X grid covariates 
  R_yx <- diag(p_y) # For Y grid covariates
  
  # Calculate reasonable initial values for lambda parameters
  calculate_initial_lambdas <- function(D, p_x, p_y) {
    # Use data means to initialize lambdas at reasonable values
    lambda_x_init <- matrix(0.01, nrow = D, ncol = p_x)
    lambda_yx_init <- matrix(0.01, nrow = D, ncol = p_y)
    lambda_y_init <- rep(0.01, D)
    
    return(list(
      lambda_x_init = lambda_x_init,
      lambda_yx_init = lambda_yx_init,
      lambda_y_init = lambda_y_init
    ))
  }
  
  # Create initial precision matrices that are guaranteed positive definite
  initialize_precision <- function(dim) {
    # Create a well-conditioned positive definite matrix
    scale <- diag(dim)
    df <- dim + 2  # Slightly more degrees of freedom for stability
    W <- rWishart(1, df, scale)[,,1]
    # Scale the matrix to have reasonable magnitudes
    W <- W / mean(diag(W))
    return(W)
  }
  
  # Initialize spatial random effects
  initialize_spatial_effects <- function(S, p) {
    # Start with small random values to break symmetry
    matrix(rnorm(S * p, 0, 0.1), nrow = S, ncol = p)
  }
  
  # Calculate dimensions
  D <- n_atoms
  S_x <- nrow(gridx_xorder)
  S_y <- nrow(gridy_yorder)
  initial_lambdas <- calculate_initial_lambdas(D, p_x, p_y)
  
  constants <- list(
    # Dimensions
    p_x = p_x,                  # Number of X covariates
    p_y = p_y,                  # Number of Y covariates
    D = as.integer(n_atoms),
    S_x = nrow(gridx_xorder),   # X grid size
    S_y = nrow(gridy_yorder),   # Y grid size
    J_x = bookkeeping$J_x,      # Number of atom-equivalent X grids
    J_y = bookkeeping$J_y,      # Number of atom-equivalent Y grids
    
    # Mapping between spaces
    x_to_atom = bookkeeping$x_to_atom,  # Maps X grids to atoms
    y_to_atom = bookkeeping$y_to_atom,  # Maps Y grids to atoms
    expand_x = bookkeeping$expand_x,     # Expands atoms to X grid
    expand_y = bookkeeping$expand_y,     # Expands atoms to Y grid
    x_reorder = bookkeeping$x_reorder,
    
    # Latent indices
    xlatent_ind = bookkeeping$x_latentind,
    ylatent_ind = bookkeeping$y_latentind,
    
    # CAR structures
    num_x = as.numeric(as.carAdjacency(adjacency$W_x)$num),
    weights_x = as.numeric(as.carAdjacency(adjacency$W_x)$weights),
    adj_x = as.numeric(as.carAdjacency(adjacency$W_x)$adj),
    num_y = as.numeric(as.carAdjacency(adjacency$W_y)$num),
    weights_y = as.numeric(as.carAdjacency(adjacency$W_y)$weights),
    adj_y = as.numeric(as.carAdjacency(adjacency$W_y)$adj),
    
    # Wishart prior matrices
    R_x = R_x,         # Add scale matrix for X covariates
    R_yx = R_yx,       # Add scale matrix for Y covariates
    df_x = p_x + 2,    # Degrees of freedom for Wishart
    df_yx = p_y + 2    # Degrees of freedom for Wishart
  )
  
  data <- list(
    x = covar_x,     # X covariates
    y = covar_y,     # Y covariates and response
    offs_x = rep(0, n_atoms),
    offs_y = rep(0, n_atoms),
    yx_obs = covar_y[, 1:p_y],  # Add explicit variables
    y_obs = covar_y[, p_y + 1]  # Y outcome
  )
  
  inits <- list(
    # Beta parameters with small non-zero values to avoid NA warnings
    beta_0_x = rep(0.01, p_x),
    beta_0_y = 0.01,
    beta_0_yx = rep(0.01, p_y),
    beta_y = rep(0.01, p_x + p_y),
    
    # Spatial random effects
    tau_x = rep(1, p_x),
    tau_y = 1,
    tau_yx = rep(1, p_y),
    
    # Matrices for spatial correlation
    Prec_x = initialize_precision(p_x),
    Prec_yx = initialize_precision(p_y),
    
    # Spatial random effects
    psi_x = matrix(rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x),
    psi_yx = matrix(rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y),
    phi_y = rnorm(S_y, 0, 0.1),
    phi_yx = matrix(rnorm(S_y * p_y, 0, 0.1), nrow = S_y, ncol = p_y),
    phi_x = matrix(rnorm(S_x * p_x, 0, 0.1), nrow = S_x, ncol = p_x),
    
    # Initialize temporary matrices
    temp_x = array(0, dim = c(n_atoms, p_x)),
    temp_yx = array(0, dim = c(n_atoms, p_y)),
    x_atomord = matrix(0.01, nrow = n_atoms, ncol = p_x),
    
    # Latent variables
    x_latent = matrix(1, nrow = D - constants$J_x, ncol = p_x),
    y_latent = rep(1, D - constants$J_y),
    yx_latent = matrix(1, nrow = D - constants$J_y, ncol = p_y),
    
    # Lambda values
    lambda_atom_x = initial_lambdas$lambda_x_init,
    lambda_atom_yx = initial_lambdas$lambda_yx_init,
    lambda_atom_y = initial_lambdas$lambda_y_init
  )
  
  # Now fill in the latent values with proper integer splits
  # For x_latent
  for(j in 1:p_x) {
    for(m in 1:(S_x-constants$J_x)) {
      total_count <- max(1L, data$x[constants$J_x + m, j])  # Ensure positive integer
      n_atoms <- diff(constants$xlatent_ind[m,]) + 1
      base_count <- as.integer(floor(total_count / n_atoms))
      remainder <- as.integer(total_count %% n_atoms)
      inits$x_latent[constants$xlatent_ind[m,1]:constants$xlatent_ind[m,2], j] <- base_count
      if(remainder > 0) {
        inits$x_latent[constants$xlatent_ind[m,1]:(constants$xlatent_ind[m,1]+remainder-1), j] <- base_count + 1L
      }
    }
  }
  
  # For y_latent
  for(m in 1:(S_y-constants$J_y)) {
    total_count <- max(1L, data$y[constants$J_y + m, p_y + 1])
    n_atoms <- diff(constants$ylatent_ind[m,]) + 1
    base_count <- as.integer(floor(total_count / n_atoms))
    remainder <- as.integer(total_count %% n_atoms)
    inits$y_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2]] <- base_count
    if(remainder > 0) {
      inits$y_latent[constants$ylatent_ind[m,1]:(constants$ylatent_ind[m,1]+remainder-1)] <- base_count + 1L
    }
  }
  
  # For yx_latent
  for(j in 1:p_y) {
    for(m in 1:(S_y-constants$J_y)) {
      total_count <- max(1L, data$y[constants$J_y + m, j])
      n_atoms <- diff(constants$ylatent_ind[m,]) + 1
      base_count <- as.integer(floor(total_count / n_atoms))
      remainder <- as.integer(total_count %% n_atoms)
      inits$yx_latent[constants$ylatent_ind[m,1]:constants$ylatent_ind[m,2], j] <- base_count
      if(remainder > 0) {
        inits$yx_latent[constants$ylatent_ind[m,1]:(constants$ylatent_ind[m,1]+remainder-1), j] <- base_count + 1L
      }
    }
  }
  
  return(list(
    constants = constants,
    data = data,
    inits = inits
  ))
}

# Function to run NIMBLE model
run_nimble_model <- function(constants, data, inits, sim_metadata, niter = 30000, 
                             nburnin = 5000, nchains = 3, save_plots = TRUE, output_dir = NULL) {
  source('nimble_abrm_0219.R')
  
  Rmodel <- nimbleModel(abrm, constants, data, inits, calculate = FALSE)
  compmod <- compileNimble(Rmodel)
  
  # Configure MCMC
  mcmcConf <- configureMCMC(compmod, enableWAIC = FALSE,
                            monitors = c('beta_0_y', 'beta_y',
                                         'tau_x', 'tau_yx', 'tau_y',
                                         'Prec_x', 'Prec_yx'),  
                            thin = 10,
                            useConjugacy = FALSE)
  
  # Add adaptive block sampling for betas
  beta_nodes <- c('beta_0_y', paste0('beta_y[', 1:(constants$p_x + constants$p_y), ']'))
  mcmcConf$addSampler(target = beta_nodes,
                      type = 'RW_block',
                      control = list(adaptive = TRUE,
                                     scale = 0.1,    # Smaller initial scale
                                     adaptInterval = 50)) # More conservative adaptation
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project = compmod)
  
  # Run MCMC
  cat("\nRunning MCMC with", nchains, "chains...\n")
  mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains, summary = TRUE)
  
  # Calculate convergence diagnostics
  cat("\nCalculating convergence diagnostics...\n")
  diagnostics <- check_mcmc_diagnostics(mcmc.out, sim_metadata)
  
  # Save diagnostic plots if requested
  if(save_plots && !is.null(diagnostics$plots)) {
    if(is.null(sim_metadata)) {
      plot_file <- "mcmc_diagnostics.pdf"
    } else {
      plot_file <- sprintf("mcmc_diagnostics_sim%d_xcor%.1f_ycor%.1f.pdf",
                           sim_metadata$sim_number,
                           sim_metadata$x_correlation,
                           sim_metadata$y_correlation)
    }
    
    if(!is.null(output_dir)) {
      plot_file <- file.path(output_dir, plot_file)
    }
    
    pdf(plot_file, width = 12, height = 8)
    print(diagnostics$plots$trace)
    print(diagnostics$plots$density)
    dev.off()
    cat("\nDiagnostic plots saved to", plot_file, "\n")
  }
  
  # Add diagnostics to output
  mcmc.out$convergence <- diagnostics
  
  return(mcmc.out)
}

# Main function to run the entire process
run_abrm <- function(data, output_file, sim_metadata) {
  
  # Prepare spatial bookkeeping
  bookkeeping <- prepare_spatial_bookkeeping(data)
  
  # Prepare adjacency matrices
  adjacency <- prepare_adjacency_matrices(bookkeeping$gridy_yorder, bookkeeping$gridx_xorder)
  
  # Prepare NIMBLE model inputs
  nimble_inputs <- prepare_nimble_inputs(bookkeeping, adjacency, data)
  
  # Run NIMBLE model
  mcmc.out <- run_nimble_model(nimble_inputs$constants, nimble_inputs$data, nimble_inputs$inits, sim_metadata)
  
  # Save output
  saveRDS(mcmc.out, file = output_file)
  
  return(mcmc.out)
}

# Function to calculate MCMC diagnostics
# Enhanced diagnostic checks with proper data structure
check_mcmc_diagnostics <- function(mcmc_output, sim_metadata = NULL) {
  # Validate chain structure
  chains_list <- mcmc_output$samples
  if (!is.list(chains_list) || length(chains_list) < 2) {
    warning("Insufficient number of chains for convergence diagnostics")
    return(list(converged = FALSE, error = "insufficient_chains"))
  }
  
  # Convert to mcmc objects
  chains_mcmc <- lapply(chains_list, function(x) {
    if (is.null(colnames(x))) colnames(x) <- paste0("param", 1:ncol(x))
    mcmc(x)
  })
  
  mcmc_chains <- mcmc.list(chains_mcmc)
  
  # Calculate Gelman-Rubin statistics
  gelman_diag <- try(gelman.diag(mcmc_chains, autoburnin = FALSE, multivariate = FALSE)$psrf, silent = TRUE)
  
  # Calculate effective sample sizes
  eff_size <- try(effectiveSize(mcmc_chains), silent = TRUE)
  
  # Parameter means and SDs
  param_names <- colnames(chains_list[[1]])
  param_means <- colMeans(do.call(rbind, chains_list))
  param_sds <- apply(do.call(rbind, chains_list), 2, sd)
  
  # Create diagnostics data frame
  diagnostics <- data.frame(
    parameter = param_names,
    mean = param_means,
    sd = param_sds,
    rhat = if(!inherits(gelman_diag, "try-error")) gelman_diag[,1] else rep(NA, length(param_names)),
    ess = if(!inherits(eff_size, "try-error")) eff_size else rep(NA, length(param_names)),
    stringsAsFactors = FALSE
  )
  
  # Add convergence assessment
  diagnostics$converged <- with(diagnostics, rhat < 1.1 & ess > 400)
  
  # Group parameters
  param_groups <- list(
    beta = grep("^beta", param_names, value = TRUE),
    tau = grep("^tau", param_names, value = TRUE),
    prec = grep("^Prec", param_names, value = TRUE)
  )
  
  # Check group convergence
  group_convergence <- sapply(param_groups, function(group_params) {
    group_diag <- diagnostics[diagnostics$parameter %in% group_params,]
    all(group_diag$converged)
  })
  
  # Create diagnostic plots
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if(is.null(sim_metadata)) {
      # Create default metadata if none provided
      sim_metadata <- list(
        sim_number = 1,
        x_correlation = NA,
        y_correlation = NA
      )
    }
    
    # Generate plots using updated function
    plots <- create_diagnostic_plots(chains_list, sim_metadata)
  } else {
    plots <- NULL
  }
  
  return(list(
    diagnostics = diagnostics,
    group_convergence = group_convergence,
    overall_converged = all(diagnostics$converged),
    plots = plots,
    chain_correlation = cor(do.call(rbind, chains_list))
  ))
}

# Updated print function to match new data structure
print_convergence_summary <- function(convergence_results) {
  cat("\nConvergence Summary:\n")
  cat("===================\n")
  
  # Overall convergence
  cat(sprintf("\nOverall convergence: %s\n", 
              ifelse(convergence_results$overall_converged, "Achieved", "Not achieved")))
  
  # Group convergence
  cat("\nParameter Group Convergence:\n")
  for(group in names(convergence_results$group_convergence)) {
    cat(sprintf("%s parameters: %s\n", 
                group,
                ifelse(convergence_results$group_convergence[[group]], 
                       "Converged", "Not converged")))
  }
  
  # Print detailed summary for problematic parameters
  problem_params <- subset(convergence_results$diagnostics, 
                           !converged)
  
  if(nrow(problem_params) > 0) {
    cat("\nParameters requiring attention:\n")
    print(problem_params[, c("parameter", "rhat", "ess")])
  }
  
  # Print general diagnostics summary
  cat("\nDiagnostic Summary Statistics:\n")
  diagnostics <- convergence_results$diagnostics
  cat(sprintf("Median ESS: %.1f\n", median(diagnostics$ess, na.rm = TRUE)))
  cat(sprintf("Max Rhat: %.3f\n", max(diagnostics$rhat, na.rm = TRUE)))
  
  # Check for any unusually high variances
  high_var_params <- subset(diagnostics, sd/abs(mean) > 0.5)
  if(nrow(high_var_params) > 0) {
    cat("\nParameters with high relative variance:\n")
    print(high_var_params[, c("parameter", "mean", "sd")])
  }
}

# Helper function for diagnostic plots
create_diagnostic_plots <- function(chains_list, sim_metadata) {
  require(ggplot2)
  require(reshape2)
  
  # Get parameters from chain data
  params <- colnames(chains_list[[1]])
  
  # Create title with metadata
  plot_title <- sprintf("MCMC Diagnostics (x_cor=%.1f, y_cor=%.1f, sim=%d)",
                        sim_metadata$x_correlation, 
                        sim_metadata$y_correlation, 
                        sim_metadata$sim_number)
  
  # Prepare data for plotting
  plot_data <- lapply(seq_along(chains_list), function(chain) {
    data <- as.data.frame(chains_list[[chain]][, params, drop = FALSE])  # Use params names here
    data$iteration <- 1:nrow(data)
    data$chain <- factor(chain)
    melt(data, id.vars = c("iteration", "chain"))
  })
  
  plot_data <- do.call(rbind, plot_data)
  
  # Create trace plots
  trace_plot <- ggplot(plot_data, aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~variable, scales = "free_y") +
    theme_minimal() +
    labs(title = plot_title,
         subtitle = "Trace Plots for All Parameters",
         x = "Iteration", y = "Parameter Value") +
    theme(strip.text = element_text(size = 8),
          plot.title = element_text(size = 11))
  
  # Create density plots
  density_plot <- ggplot(plot_data, aes(x = value, fill = chain)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~variable, scales = "free") +
    theme_minimal() +
    labs(title = plot_title,
         subtitle = "Density Plots for All Parameters",
         x = "Parameter Value", y = "Density") +
    theme(strip.text = element_text(size = 8),
          plot.title = element_text(size = 11))
  
  return(list(trace = trace_plot, density = density_plot))
}

###
# # Load in simulated data
# data <- simulate_misaligned_data(
#   res1 = c(5, 5),
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.3,  # Correlation between X covariates
#   y_correlation = 0.7,  # Correlation between Y covariates
#   beta_x = c(0.03, -0.01, 0.06),  # True effects of X covariates
#   beta_y = c(0.04, 0.01, -0.05, 0.02)  # True effects of Y covariates
# )
# 
# sim_metadata <- list(
#   sim_number = 1,
#   x_correlation = 0.3,
#   y_correlation = 0.7
# )
# 
# result <- run_abrm(data, 'model_output.rds', sim_metadata)
