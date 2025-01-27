# Load required libraries
library(sp)
library(spdep)
library(sf)
library(MASS)
library(raster)
library(tidyverse)
library(nimble)

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
  x_vars <- c("x", grep("covariate_x_", names(gridx), value = TRUE))
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
  # Extract basic components
  gridy_yorder <- bookkeeping$gridy_yorder
  gridx_xorder <- bookkeeping$gridx_xorder
  atoms <- data$atoms
  n_atoms <- if(is.null(nrow(atoms))) length(atoms) else nrow(atoms)
  
  # Identify variables
  x_vars <- c("x", grep("covariate_x_", names(gridx_xorder), value = TRUE))
  y_vars <- c("y", grep("covariate_y_", names(gridy_yorder), value = TRUE))
  
  # Get unique atom IDs
  x_atom_ids <- unique(atoms$ID_x)
  y_atom_ids <- unique(atoms$ID_y)
  
  # Order grids
  gridx_ordered <- gridx_xorder[order(gridx_xorder$ID_x), ]
  gridy_ordered <- gridy_yorder[order(gridy_yorder$ID_y), ]
  
  # Pre-compute indices for NIMBLE
  # This replaces the dynamic (1:p_x)[-j] indexing
  create_index_matrices <- function(p) {
    result <- matrix(0, nrow = p, ncol = p-1)
    for(j in 1:p) {
      idx <- setdiff(1:p, j)  # All indices except j
      result[j,] <- idx
    }
    return(result)
  }
  
  # Create expanded indices for multinomial calculations
  create_expanded_indices <- function(latent_ind, J, max_atoms_per_grid) {
    n_nonatom_grids <- nrow(latent_ind)
    indices <- matrix(0, nrow = n_nonatom_grids, ncol = max_atoms_per_grid)
    
    for(m in 1:n_nonatom_grids) {
      start_idx <- J + latent_ind[m,1]
      end_idx <- J + latent_ind[m,2]
      indices[m,1:(end_idx - start_idx + 1)] <- start_idx:end_idx
    }
    
    return(indices)
  }
  
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
  
  # Create initial lambda values
  create_initial_lambdas <- function(D, p_x, p_y) {
    lambda_x_init <- matrix(1, nrow = D, ncol = p_x)
    lambda_yx_init <- matrix(1, nrow = D, ncol = p_y)
    lambda_y_init <- rep(1, D)
    
    return(list(
      lambda_x_init = lambda_x_init,
      lambda_yx_init = lambda_yx_init,
      lambda_y_init = lambda_y_init
    ))
  }
  
  # Create X covariate matrix
  print("\nCreating covariate matrices...")
  covar_x <- create_covariate_matrix(gridx_ordered, x_vars, x_atom_ids, nrow(gridx_xorder))
  covar_y <- create_covariate_matrix(gridy_ordered, y_vars, y_atom_ids, nrow(gridy_yorder))
  
  # Create index matrices
  p_x <- ncol(covar_x)
  p_y <- ncol(covar_y) - 1  # subtract y column
  
  # Create index matrices
  x_index_mat <- create_index_matrices(p_x)
  y_index_mat <- create_index_matrices(p_y)
  
  # Calculate max atoms per grid
  max_atoms_x <- max(bookkeeping$x_latentind[,2] - bookkeeping$x_latentind[,1] + 1)
  max_atoms_y <- max(bookkeeping$y_latentind[,2] - bookkeeping$y_latentind[,1] + 1)
  
  # Create expanded indices
  x_expanded_indices <- create_expanded_indices(bookkeeping$x_latentind, 
                                                bookkeeping$J_x, max_atoms_x)
  y_expanded_indices <- create_expanded_indices(bookkeeping$y_latentind, 
                                                bookkeeping$J_y, max_atoms_y)
  
  # Create initial values
  initial_lambdas <- create_initial_lambdas(n_atoms, p_x, p_y)
  
  constants <- list(
    p_x = p_x,
    p_y = p_y,
    D = n_atoms,
    S_x = nrow(gridx_xorder),
    S_y = nrow(gridy_yorder),
    J_x = bookkeeping$J_x,
    J_y = bookkeeping$J_y,
    x_to_atom = bookkeeping$x_to_atom,
    y_to_atom = bookkeeping$y_to_atom,
    xlatent_ind = bookkeeping$x_latentind,
    ylatent_ind = bookkeeping$y_latentind,
    expand_x = bookkeeping$expand_x,
    expand_y = bookkeeping$expand_y,
    x_reorder = bookkeeping$x_reorder,
    num_x = as.numeric(as.carAdjacency(adjacency$W_x)$num),
    weights_x = as.numeric(as.carAdjacency(adjacency$W_x)$weights),
    adj_x = as.numeric(as.carAdjacency(adjacency$W_x)$adj),
    num_y = as.numeric(as.carAdjacency(adjacency$W_y)$num),
    weights_y = as.numeric(as.carAdjacency(adjacency$W_y)$weights),
    adj_y = as.numeric(as.carAdjacency(adjacency$W_y)$adj),
    x_index_mat = x_index_mat,
    y_index_mat = y_index_mat,
    x_expanded_indices = x_expanded_indices,
    y_expanded_indices = y_expanded_indices,
    max_atoms_x = max_atoms_x,
    max_atoms_y = max_atoms_y
  )
  
  data <- list(
    x = covar_x,
    y = covar_y,
    offs_x = rep(0, n_atoms),
    offs_y = rep(0, n_atoms)
  )
  
  inits <- list(
    beta_0_x = rep(0, p_x),
    beta_x = matrix(0, nrow = p_x, ncol = p_x - 1),
    beta_0_yx = rep(0, p_y),
    beta_yx = matrix(0, nrow = p_y, ncol = p_x + p_y - 1),
    beta_0_y = 0,
    beta_y = rep(0, p_x + p_y),
    tau_x = rep(1, p_x),
    tau_y = 1,
    tau_yx = rep(1, p_y),
    phi_x = matrix(0, nrow = constants$S_x, ncol = p_x),
    phi_y = rep(0, constants$S_y),
    phi_yx = matrix(0, nrow = constants$S_y, ncol = p_y),
    x_latent = matrix(1, nrow = n_atoms - constants$J_x, ncol = p_x),
    y_latent = rep(1, n_atoms - constants$J_y),
    yx_latent = matrix(1, nrow = n_atoms - constants$J_y, ncol = p_y),
    lambda_atom_x = initial_lambdas$lambda_x_init,
    lambda_atom_yx = initial_lambdas$lambda_yx_init,
    lambda_atom_y = initial_lambdas$lambda_y_init,
    lambda_atom_x_init = initial_lambdas$lambda_x_init,
    lambda_atom_yx_init = initial_lambdas$lambda_yx_init,
    lambda_atom_y_init = initial_lambdas$lambda_y_init,
    
    beta_sum_x_temp1 = matrix(0, nrow = n_atoms, ncol = p_x),
    beta_sum_x_temp2 = matrix(0, nrow = n_atoms, ncol = p_x),
    beta_sum_x = matrix(0, nrow = n_atoms, ncol = p_x),
    
    # Add intermediate sum variables for Y covariates
    beta_sum_yx_x_temp1 = matrix(0, nrow = n_atoms, ncol = p_y),
    beta_sum_yx_x_temp2 = matrix(0, nrow = n_atoms, ncol = p_y),
    beta_sum_yx_x_temp3 = matrix(0, nrow = n_atoms, ncol = p_y),
    beta_sum_yx_x = matrix(0, nrow = n_atoms, ncol = p_y),
    
    beta_sum_yx_y_temp1 = matrix(0, nrow = n_atoms, ncol = p_y),
    beta_sum_yx_y_temp2 = matrix(0, nrow = n_atoms, ncol = p_y),
    beta_sum_yx_y = matrix(0, nrow = n_atoms, ncol = p_y),
    
    # Add intermediate sum variables for final Y model
    beta_sum_y_x_temp1 = rep(0, n_atoms),
    beta_sum_y_x_temp2 = rep(0, n_atoms),
    beta_sum_y_x_temp3 = rep(0, n_atoms),
    beta_sum_y_x = rep(0, n_atoms),
    
    beta_sum_y_y_temp1 = rep(0, n_atoms),
    beta_sum_y_y_temp2 = rep(0, n_atoms),
    beta_sum_y_y_temp3 = rep(0, n_atoms),
    beta_sum_y_y = rep(0, n_atoms)
  )
  
  print("\n=== Model Component Summary ===")
  print(paste("X covariate matrix:", nrow(covar_x), "x", ncol(covar_x)))
  print(paste("Y covariate matrix:", nrow(covar_y), "x", ncol(covar_y)))
  print(paste("Number of atoms:", n_atoms))
  print(paste("X grid size:", constants$S_x))
  print(paste("Y grid size:", constants$S_y))
  print(paste("J_x:", constants$J_x))
  print(paste("J_y:", constants$J_y))
  print(paste("p_y =", p_y))
  print(paste("Dimensions of covar_y:", paste(dim(covar_y), collapse=" x ")))
  print("Column names of covar_y:")
  print(colnames(covar_y))
  
  return(list(
    constants = constants,
    data = data,
    inits = inits
  ))
}

# Function to run NIMBLE model
run_nimble_model <- function(constants, data, inits, niter = 10000, nburnin = 5000, nchains = 3) {
  source('nimble_abrm_1114success.R')
  Rmodel <- nimbleModel(abrm, constants, data, inits, calculate = FALSE)
  compmod <- compileNimble(Rmodel)
  
  # Configure MCMC
  mcmcConf <- configureMCMC(compmod, enableWAIC = FALSE,
                            monitors = c('beta_0_x', 'beta_x', 
                                         'beta_0_yx', 'beta_yx',
                                         'beta_0_y', 'beta_y',
                                         'tau_x', 'tau_yx', 'tau_y'), 
                            thin = 10,
                            useConjugacy = FALSE)
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project = compmod)
  
  # Run MCMC
  mcmc.out <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, 
                      nchains = nchains, summary = TRUE)
  
  return(mcmc.out)
}

# Main function to run the entire process
run_abrm <- function(data, output_file) {
  
  # Prepare spatial bookkeeping
  bookkeeping <- prepare_spatial_bookkeeping(data)
  
  # Prepare adjacency matrices
  adjacency <- prepare_adjacency_matrices(bookkeeping$gridy_yorder, bookkeeping$gridx_xorder)
  
  # Prepare NIMBLE model inputs
  nimble_inputs <- prepare_nimble_inputs(bookkeeping, adjacency, data)
  
  # Run NIMBLE model
  mcmc.out <- run_nimble_model(nimble_inputs$constants, nimble_inputs$data, nimble_inputs$inits)
  
  # Save output
  saveRDS(mcmc.out, file = output_file)
  
  return(mcmc.out)
}

################################################
# Load in simulated data
source("sim.R")
data <- simulate_misaligned_data(res1 = c(5, 5), res2 = c(10, 10), seed = 1,
                                 n_covariates_x = 3, n_covariates_y = 4)
# saveRDS(sim_data, file = 'sim_data.rds')
# sim_data <- readRDS('sim_data.rds')
# load('sim_data.RData')
result <- run_abrm(data, 'model_output.rds')