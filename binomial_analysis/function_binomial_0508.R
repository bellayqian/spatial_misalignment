library(sp)
library(sf)
library(spdep)
library(MASS)
library(raster)
library(tidyverse)
library(nimble)
library(ggplot2)
library(coda)

# Simulation part
gen_correlated_spat_binomial <- function(W, n_vars, rho = 0.6, var_spat = 1, correlation = 0.5, 
                                         global_ints = NULL, verify = FALSE, population_sizes = NULL) {
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
  
  # Generate observations using binomial
  observed_values <- matrix(0, n, n_vars)
  for(j in 1:n_vars) {
    int <- ifelse(is.null(global_ints), 0, global_ints[j])
    # Calculate probabilities using logistic function for binomial
    probabilities <- plogis(int + spatial_effects[,j])
    # Ensure probabilities aren't too extreme, to reduce chance of constraint violations
    probabilities <- pmin(pmax(probabilities, 0.1), 0.9)
    # Generate binomial data
    observed_values[,j] <- rbinom(n, size = population_sizes, prob = probabilities)
    
    # Add constraint check to ensure counts don't exceed population size
    observed_values[,j] <- pmin(observed_values[,j], population_sizes)
  }
  
  return(list(
    counts = observed_values,
    probabilities = probabilities,
    spatial_effects = spatial_effects
  ))
}

# Main simulation function
simulate_misaligned_binomial_data <- function(res1 = c(5, 5), res2 = c(10, 10), 
                                              seed = 2, n_covariates_x = 3, n_covariates_y = 4,
                                              x_correlation = 0.5, y_correlation = 0.5,
                                              beta_x = NULL, beta_y = NULL, 
                                              pop_min = 100, pop_max = 500) {
  set.seed(seed)
  
  # Create the spatial grids
  grid1 <- GridTopology(cellcentre.offset = c(0.5, 0.5), cellsize = c(2, 2), cells.dim = res1)
  sp_grid1 <- SpatialGrid(grid1)
  grid2 <- GridTopology(cellcentre.offset = c(0, 0), cellsize = c(1, 1), cells.dim = res2)
  sp_grid2 <- SpatialGrid(grid2)
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
  
  # Create atoms for outcome generation
  atoms <- st_as_sf(raster::intersect(as(gridy, 'Spatial'), as(gridx, 'Spatial')))
  names(atoms)[which(names(atoms) == 'ID')] <- c("ID_y")
  names(atoms)[which(names(atoms) == 'ID.1')] <- c("ID_x")
  atoms$ID_atomorder <- 1:nrow(atoms)
  
  # Generate atom populations directly
  atoms$population <- sample(pop_min:pop_max, nrow(atoms), replace = TRUE)
  
  # Sum up to get grid populations
  gridx$population <- sapply(1:nrow(gridx), function(i) {
    atom_indices <- which(atoms$ID_x == i)
    sum(atoms$population[atom_indices])
  })
  
  gridy$population <- sapply(1:nrow(gridy), function(i) {
    atom_indices <- which(atoms$ID_y == i)
    sum(atoms$population[atom_indices])
  })
  
  # Generate atom-level covariates and outcomes using CAR models
  # First create a spatial adjacency matrix for atoms
  atoms_sp <- as(atoms, "Spatial")
  neighbors_atoms <- poly2nb(atoms_sp, queen = TRUE)
  W_atoms <- nb2mat(neighbors_atoms, style = "B", zero.policy = TRUE)
  
  # Generate X covariates at atom level
  x_atom_results <- gen_correlated_spat_binomial(
    W_atoms, n_covariates_x, 
    correlation = x_correlation,
    global_ints = rnorm(n_covariates_x, 0, 0.5),
    population_sizes = atoms$population
  )
  
  # Assign X values to atoms
  for(i in 1:n_covariates_x) {
    atoms[[paste0("covariate_x_", i)]] <- x_atom_results$counts[,i]
  }
  
  # Generate Y covariates at atom level
  y_cov_atom_results <- gen_correlated_spat_binomial(
    W_atoms, n_covariates_y, 
    correlation = y_correlation,
    global_ints = rnorm(n_covariates_y, 0, 0.5),
    population_sizes = atoms$population
  )
  
  # Assign Y covariate values to atoms
  for(i in 1:n_covariates_y) {
    atoms[[paste0("covariate_y_", i)]] <- y_cov_atom_results$counts[,i]
  }
  
  # Generate atom-level outcome with true associations
  linear_pred_y_atom <- rep(0, nrow(atoms))
  
  # Add X covariate effects
  for(i in 1:n_covariates_x) {
    x_vals <- atoms[[paste0("covariate_x_", i)]]
    x_pops <- atoms$population
    x_props <- x_vals / x_pops
    # # Option 1: Log-odds
    # linear_pred_y_atom <- linear_pred_y_atom + beta_x[i] * log(x_props / (1 - x_props + 1e-10))
    # Raw Proportions
    linear_pred_y_atom <- linear_pred_y_atom + beta_x[i] * x_props
  }
  
  # Add Y covariate effects
  for(i in 1:n_covariates_y) {
    y_vals <- atoms[[paste0("covariate_y_", i)]]
    y_pops <- atoms$population
    y_props <- y_vals / y_pops
    # linear_pred_y_atom <- linear_pred_y_atom + beta_y[i] * log(y_props / (1 - y_props + 1e-10))
    linear_pred_y_atom <- linear_pred_y_atom + beta_y[i] * y_props
  }
  
  # Add spatial effect
  spatial_effect_y_atom <- gen_correlated_spat_binomial(
    W_atoms, n_vars = 1, rho = 0.6, var_spat = 1,
    correlation = 1, population_sizes = atoms$population
  )$spatial_effects[,1]
  
  # Generate final outcome (binomial) at atom level
  outcome_prob_atom <- plogis(linear_pred_y_atom + spatial_effect_y_atom)
  atoms$y <- rbinom(nrow(atoms), size = atoms$population, prob = outcome_prob_atom)
  atoms$true_prob <- outcome_prob_atom
  
  # Aggregate to grid level
  # Aggregate X covariates to X grid
  for(i in 1:n_covariates_x) {
    gridx[[paste0("covariate_x_", i)]] <- sapply(1:nrow(gridx), function(j) {
      atom_indices <- which(atoms$ID_x == j)
      sum(atoms[[paste0("covariate_x_", i)]][atom_indices])
    })
  }
  
  # Aggregate Y covariates to Y grid
  for(i in 1:n_covariates_y) {
    gridy[[paste0("covariate_y_", i)]] <- sapply(1:nrow(gridy), function(j) {
      atom_indices <- which(atoms$ID_y == j)
      sum(atoms[[paste0("covariate_y_", i)]][atom_indices])
    })
  }
  
  # Aggregate outcome to Y grid
  gridy$y <- sapply(1:nrow(gridy), function(j) {
    atom_indices <- which(atoms$ID_y == j)
    sum(atoms$y[atom_indices])
  })
  
  # Store true parameters
  true_params <- list(
    beta_0 = 0,
    beta_x = beta_x,
    beta_y = beta_y,
    x_correlation = x_correlation,
    y_correlation = y_correlation,
    spatial_effect_y = spatial_effect_y_atom
  )
  
  return(list(
    gridy = gridy, 
    gridx = gridx, 
    atoms = atoms, 
    true_params = true_params
  ))
}

# Pre-process data to ensure population summation
preprocess_data <- function(data) {
  atoms <- data$atoms
  gridx <- data$gridx
  gridy <- data$gridy
  
  # 1. Make sure population sizes exist
  if(!"population" %in% names(gridx)) {
    cat("Adding population sizes to gridx\n")
    gridx$population <- sample(100:500, nrow(gridx), replace = TRUE)
  }
  
  if(!"population" %in% names(gridy)) {
    cat("Adding population sizes to gridy\n")
    gridy$population <- sample(100:500, nrow(gridy), replace = TRUE)
  }
  
  if(!"population" %in% names(atoms)) {
    cat("Adding population sizes to atoms\n")
    atoms$area <- st_area(atoms)
    atoms$population <- round(atoms$area / sum(atoms$area) * sum(gridx$population))
  }
  
  # 2. Check each atom's population vs. the counts in corresponding grid cells
  atom_pop_sum_x <- numeric(nrow(gridx))
  atom_pop_sum_y <- numeric(nrow(gridy))
  
  # Calculate sum of atom populations for each grid
  for(i in 1:nrow(atoms)) {
    atom_pop_sum_x[atoms$ID_x[i]] <- atom_pop_sum_x[atoms$ID_x[i]] + atoms$population[i]
    atom_pop_sum_y[atoms$ID_y[i]] <- atom_pop_sum_y[atoms$ID_y[i]] + atoms$population[i]
  }
  
  # Check for zero or NA population sums and fix them
  if(any(atom_pop_sum_x <= 0) || any(is.na(atom_pop_sum_x))) {
    cat("Warning: Found zero, negative, or NA values in atom_pop_sum_x. Fixing...\n")
    problem_indices <- which(atom_pop_sum_x <= 0 | is.na(atom_pop_sum_x))
    atom_pop_sum_x[problem_indices] <- median(atom_pop_sum_x[atom_pop_sum_x > 0], na.rm=TRUE)
    if(is.na(atom_pop_sum_x[1])) atom_pop_sum_x <- rep(100, length(atom_pop_sum_x))
  }
  
  if(any(atom_pop_sum_y <= 0) || any(is.na(atom_pop_sum_y))) {
    cat("Warning: Found zero, negative, or NA values in atom_pop_sum_y. Fixing...\n")
    problem_indices <- which(atom_pop_sum_y <= 0 | is.na(atom_pop_sum_y))
    atom_pop_sum_y[problem_indices] <- median(atom_pop_sum_y[atom_pop_sum_y > 0], na.rm=TRUE)
    if(is.na(atom_pop_sum_y[1])) atom_pop_sum_y <- rep(100, length(atom_pop_sum_y))
  }
  
  # 3. Check and adjust X covariate counts to not exceed population
  cat("\nChecking X grid variables for population constraints:\n")
  for(i in 1:length(grep("covariate_x_", names(gridx), value=TRUE))) {
    var_name <- paste0("covariate_x_", i)
    violations <- 0
    for(j in 1:nrow(gridx)) {
      if(gridx[[var_name]][j] > atom_pop_sum_x[j]) {
        violations <- violations + 1
        cat("Adjusting", var_name, "for grid cell", j, "from", 
            gridx[[var_name]][j], "to", atom_pop_sum_x[j], "\n")
        gridx[[var_name]][j] <- atom_pop_sum_x[j]
      }
    }
    if(violations == 0) {
      cat("No violations found for", var_name, "\n")
    }
  }
  
  # 4. NEW: Check and adjust Y covariate counts and outcome to not exceed population
  cat("\nChecking Y grid variables for population constraints:\n")
  
  # Check Y covariates
  for(i in 1:length(grep("covariate_y_", names(gridy), value=TRUE))) {
    var_name <- paste0("covariate_y_", i)
    violations <- 0
    for(j in 1:nrow(gridy)) {
      if(gridy[[var_name]][j] > atom_pop_sum_y[j]) {
        violations <- violations + 1
        cat("Adjusting", var_name, "for grid cell", j, "from", 
            gridy[[var_name]][j], "to", atom_pop_sum_y[j], "\n")
        gridy[[var_name]][j] <- atom_pop_sum_y[j]
      }
    }
    if(violations == 0) {
      cat("No violations found for", var_name, "\n")
    }
  }
  
  # Check Y outcome
  if("y" %in% names(gridy)) {
    violations <- 0
    for(j in 1:nrow(gridy)) {
      if(gridy$y[j] > atom_pop_sum_y[j]) {
        violations <- violations + 1
        cat("Adjusting y (outcome) for grid cell", j, "from", 
            gridy$y[j], "to", atom_pop_sum_y[j], "\n")
        gridy$y[j] <- atom_pop_sum_y[j]
      }
    }
    if(violations == 0) {
      cat("No violations found for y (outcome)\n")
    }
  } else {
    cat("Warning: No 'y' column found in gridy\n")
  }
  
  # 5. NEW: Check if any atom has a population of 0
  if(any(atoms$population <= 0) || any(is.na(atoms$population))) {
    cat("\nWarning: Found", sum(atoms$population <= 0, na.rm=TRUE), "atoms with zero or negative population\n")
    cat("Warning: Found", sum(is.na(atoms$population)), "atoms with NA population\n")
    
    # Fix atom populations
    problem_atoms <- which(atoms$population <= 0 | is.na(atoms$population))
    median_pop <- median(atoms$population[atoms$population > 0], na.rm=TRUE)
    if(is.na(median_pop)) median_pop <- 100  # Fallback
    
    cat("Setting population of problematic atoms to", median_pop, "\n")
    atoms$population[problem_atoms] <- median_pop
  }
  
  # Return the adjusted data
  data$gridx <- gridx
  data$gridy <- gridy
  data$atoms <- atoms
  return(data)
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

# Function to prepare NIMBLE model inputs for binomial ABRM
prepare_nimble_inputs_binomial <- function(bookkeeping, adjacency, data) {
  # Basic dimensions - unchanged
  gridy_yorder <- bookkeeping$gridy_yorder
  gridx_xorder <- bookkeeping$gridx_xorder
  atoms <- data$atoms
  n_atoms <- if(is.null(nrow(atoms))) length(atoms) else nrow(atoms)
  
  # Identify variables - unchanged
  x_vars <- grep("covariate_x_", names(gridx_xorder), value = TRUE)
  y_vars <- c("y", grep("covariate_y_", names(gridy_yorder), value = TRUE))
  
  # Create covariate matrices - unchanged
  create_covariate_matrix <- function(grid_data, vars) {
    matrix_data <- matrix(NA, nrow = nrow(grid_data), ncol = length(vars))
    colnames(matrix_data) <- vars
    
    for(i in 1:length(vars)) {
      matrix_data[,i] <- as.numeric(grid_data[[vars[i]]])
    }
    
    return(matrix_data)
  }
  
  print("\nCreating covariate matrices...")
  
  # Create X covariate matrix
  covar_x <- create_covariate_matrix(
    gridx_xorder[order(gridx_xorder$ID_x), ], 
    bookkeeping$x_vars
  )
  
  # Create Y covariate matrix (includes response)
  covar_y <- matrix(NA, nrow = nrow(gridy_yorder), ncol = length(bookkeeping$y_vars))
  colnames(covar_y) <- bookkeeping$y_vars
  
  for(i in 1:length(bookkeeping$y_vars)) {
    covar_y[,i] <- as.numeric(gridy_yorder[order(gridy_yorder$ID_y), ][[bookkeeping$y_vars[i]]])
  }
  
  # Population sizes for binomial distribution
  n_x <- gridx_xorder$population[order(gridx_xorder$ID_x)]
  n_y <- gridy_yorder$population[order(gridy_yorder$ID_y)]
  
  # # Verify population sizes are positive
  # if(any(n_x <= 0) || any(is.na(n_x))) {
  #   cat("Warning: Found invalid X population sizes. Fixing...\n")
  #   n_x[n_x <= 0 | is.na(n_x)] <- median(n_x, na.rm=TRUE)
  #   if(is.na(n_x[1])) n_x <- rep(100, length(n_x))
  # }
  # 
  # if(any(n_y <= 0) || any(is.na(n_y))) {
  #   cat("Warning: Found invalid Y population sizes. Fixing...\n")
  #   n_y[n_y <= 0 | is.na(n_y)] <- median(n_y, na.rm=TRUE)
  #   if(is.na(n_y[1])) n_y <- rep(100, length(n_y))
  # }
  # 
  # # Check that observed counts don't exceed populations
  # for(j in 1:ncol(covar_x)) {
  #   for(i in 1:nrow(covar_x)) {
  #     if(covar_x[i,j] > n_x[i]) {
  #       cat("Warning: X count", covar_x[i,j], "exceeds population", n_x[i], 
  #           "at grid", i, "variable", j, ". Adjusting...\n")
  #       covar_x[i,j] <- n_x[i]
  #     }
  #   }
  # }
  # 
  # for(j in 1:ncol(covar_y)) {
  #   for(i in 1:nrow(covar_y)) {
  #     if(covar_y[i,j] > n_y[i]) {
  #       cat("Warning: Y count", covar_y[i,j], "exceeds population", n_y[i], 
  #           "at grid", i, "variable", j, ". Adjusting...\n")
  #       covar_y[i,j] <- n_y[i]
  #     }
  #   }
  # }
  
  # Calculate dimensions
  p_x <- ncol(covar_x)
  p_y <- ncol(covar_y) - 1  # subtract y column
  S_x <- nrow(gridx_xorder)
  S_y <- nrow(gridy_yorder)
  D <- n_atoms
  
  # Create index matrices - unchanged
  print("\nPreparing index matrices for hypergeometric distribution...")
  
  # Calculate total latent size
  x_latent_size <- D - bookkeeping$J_x
  y_latent_size <- D - bookkeeping$J_y
  
  # For X grid cells - create index matrix
  xlatent_ind <- matrix(0, nrow = S_x - bookkeeping$J_x, ncol = 2)
  x_non_atom_sizes <- numeric(S_x - bookkeeping$J_x)
  
  curr_idx <- 1
  for(m in 1:(S_x - bookkeeping$J_x)) {
    # Calculate sizes using the original latent indices
    start_idx <- bookkeeping$x_latentind[m, 1]
    end_idx <- bookkeeping$x_latentind[m, 2]
    cell_size <- end_idx - start_idx + 1
    
    # Store the start and end indices
    xlatent_ind[m, 1] <- curr_idx
    xlatent_ind[m, 2] <- curr_idx + cell_size - 1
    
    # Store sizes for future reference
    x_non_atom_sizes[m] <- cell_size
    
    # Update current index
    curr_idx <- curr_idx + cell_size
  }
  
  # For Y grid cells - create index matrix
  ylatent_ind <- matrix(0, nrow = S_y - bookkeeping$J_y, ncol = 2)
  y_non_atom_sizes <- numeric(S_y - bookkeeping$J_y)
  
  curr_idx <- 1
  for(m in 1:(S_y - bookkeeping$J_y)) {
    # Calculate sizes using the original latent indices
    start_idx <- bookkeeping$y_latentind[m, 1]
    end_idx <- bookkeeping$y_latentind[m, 2]
    cell_size <- end_idx - start_idx + 1
    
    # Store the start and end indices
    ylatent_ind[m, 1] <- curr_idx
    ylatent_ind[m, 2] <- curr_idx + cell_size - 1
    
    # Store sizes for future reference
    y_non_atom_sizes[m] <- cell_size
    
    # Update current index
    curr_idx <- curr_idx + cell_size
  }
  
  # Get population sizes directly from atoms
  n_x_latent <- numeric(x_latent_size)
  n_y_latent <- numeric(y_latent_size)
  
  # Map using x_to_atom and y_to_atom
  latent_x_count <- 1
  for(m in 1:(S_x - bookkeeping$J_x)) {
    start_idx <- bookkeeping$x_latentind[m, 1]
    end_idx <- bookkeeping$x_latentind[m, 2]
    
    for(i in start_idx:end_idx) {
      atom_idx <- bookkeeping$x_to_atom[bookkeeping$J_x + i]
      # IMPORTANT CHECK: Ensure atom_idx is valid
      if(is.na(atom_idx) || atom_idx < 1 || atom_idx > nrow(atoms)) {
        cat("Warning: Invalid atom index", atom_idx, "for X latent", latent_x_count, "\n")
        n_x_latent[latent_x_count] <- 100  # Default value
      } else {
        n_x_latent[latent_x_count] <- atoms$population[atom_idx]
      }
      latent_x_count <- latent_x_count + 1
    }
  }
  
  latent_y_count <- 1
  for(m in 1:(S_y - bookkeeping$J_y)) {
    start_idx <- bookkeeping$y_latentind[m, 1]
    end_idx <- bookkeeping$y_latentind[m, 2]
    
    for(i in start_idx:end_idx) {
      # Map to the correct atom index
      atom_idx <- bookkeeping$y_to_atom[bookkeeping$J_y + i]
      # IMPORTANT CHECK: Ensure atom_idx is valid
      if(is.na(atom_idx) || atom_idx < 1 || atom_idx > nrow(atoms)) {
        cat("Warning: Invalid atom index", atom_idx, "for Y latent", latent_y_count, "\n")
        n_y_latent[latent_y_count] <- 100  # Default value
      } else {
        n_y_latent[latent_y_count] <- atoms$population[atom_idx]
      }
      latent_y_count <- latent_y_count + 1
    }
  }
  
  # Check for zero or NA population sizes and fix them
  if(any(is.na(n_x_latent)) || any(n_x_latent <= 0)) {
    cat("Warning: Found", sum(is.na(n_x_latent)), "NA and", sum(n_x_latent <= 0, na.rm=TRUE), 
        "zero or negative values in n_x_latent\n")
    
    # Use median of non-NA, positive values as replacement
    median_val <- median(n_x_latent[!is.na(n_x_latent) & n_x_latent > 0])
    if(is.na(median_val)) median_val <- 100  # Fallback
    
    n_x_latent[is.na(n_x_latent) | n_x_latent <= 0] <- median_val
  }
  
  if(any(is.na(n_y_latent)) || any(n_y_latent <= 0)) {
    cat("Warning: Found", sum(is.na(n_y_latent)), "NA and", sum(n_y_latent <= 0, na.rm=TRUE), 
        "zero or negative values in n_y_latent\n")
    
    # Use median of non-NA, positive values as replacement
    median_val <- median(n_y_latent[!is.na(n_y_latent) & n_y_latent > 0])
    if(is.na(median_val)) median_val <- 100  # Fallback
    
    n_y_latent[is.na(n_y_latent) | n_y_latent <= 0] <- median_val
  }
  
  # Remaining code unchanged
  R_x <- diag(p_x)
  R_yx <- diag(p_y)
  offs_x <- rep(0, D)
  offs_y <- rep(0, D)
  
  # Prepare constants
  constants <- list(
    # Dimensions
    p_x = p_x,
    p_y = p_y,
    D = as.integer(n_atoms),
    S_x = S_x,
    S_y = S_y,
    J_x = bookkeeping$J_x,
    J_y = bookkeeping$J_y,
    
    # Index matrices
    xlatent_ind = xlatent_ind,
    ylatent_ind = ylatent_ind,
    
    # Mapping arrays
    x_to_atom = bookkeeping$x_to_atom,
    y_to_atom = bookkeeping$y_to_atom,
    expand_x = bookkeeping$expand_x,
    expand_y = bookkeeping$expand_y,
    x_reorder = bookkeeping$x_reorder,
    
    # Population sizes
    n_x = n_x,
    n_y = n_y,
    n_x_latent = n_x_latent,
    n_y_latent = n_y_latent,
    
    # Offsets
    offs_x = offs_x,
    offs_y = offs_y,
    
    # CAR structures
    num_x = as.numeric(as.carAdjacency(adjacency$W_x)$num),
    weights_x = as.numeric(as.carAdjacency(adjacency$W_x)$weights),
    adj_x = as.numeric(as.carAdjacency(adjacency$W_x)$adj),
    num_y = as.numeric(as.carAdjacency(adjacency$W_y)$num),
    weights_y = as.numeric(as.carAdjacency(adjacency$W_y)$weights),
    adj_y = as.numeric(as.carAdjacency(adjacency$W_y)$adj),
    
    # Wishart prior matrices
    R_x = R_x,
    R_yx = R_yx,
    df_x = p_x + 2,
    df_yx = p_y + 2
  )
  
  # Prepare data
  data_list <- list(
    x = covar_x,                # X covariates
    yx_obs = covar_y[, 1:p_y],  # Y covariates
    y_obs = covar_y[, p_y + 1],  # Y outcome
    Pop_atom_true = data$atoms$population[order(data$atoms$ID_atomorder)]
  )
  
  # Create initial values with more stable defaults
  inits <- list(
    # Beta parameters - use zeros for stability
    beta_0_x = rep(0, p_x),
    beta_0_y = 0,
    beta_0_yx = rep(0, p_y),
    beta_y = rep(0, p_x + p_y),
    
    # Precision parameters - use moderate values
    tau_x = rep(5, p_x),
    tau_y = 5,
    tau_yx = rep(5, p_y),
    
    # Matrices for spatial correlation - use identity matrices
    Prec_x = diag(p_x) + 0.1,
    Prec_yx = diag(p_y) + 0.1,
    
    # Spatial random effects - use smaller values for stability
    psi_x = matrix(rnorm(S_x * p_x, 0, 0.01), nrow = S_x, ncol = p_x),
    psi_yx = matrix(rnorm(S_y * p_y, 0, 0.01), nrow = S_y, ncol = p_y),
    phi_y = as.vector(rnorm(S_y, 0, 0.01)),
    phi_x = matrix(rnorm(S_x * p_x, 0, 0.01), nrow = S_x, ncol = p_x),
    phi_yx = matrix(rnorm(S_y * p_y, 0, 0.01), nrow = S_y, ncol = p_y),
    
    # Latent variables - will be filled by generate_valid_latent_inits
    x_latent = matrix(1, nrow = x_latent_size, ncol = p_x),
    yx_latent = matrix(1, nrow = y_latent_size, ncol = p_y),
    y_latent = rep(1, y_latent_size),
    
    # Temporary arrays - needed for derived quantities
    temp_x = matrix(0, nrow = D, ncol = p_x),
    temp_yx = matrix(0, nrow = D, ncol = p_y),
    x_atomord = matrix(0, nrow = D, ncol = p_x),
    
    # Linear predictors and transformed values - set to reasonable defaults
    linear_pred_x_base = matrix(0, nrow = D, ncol = p_x),
    linear_pred_yx_base = matrix(0, nrow = D, ncol = p_y),
    linear_pred_y_base = rep(0, D),
    # Odds values - set to 1 (which gives probability of 0.5)
    odds_atom_x = matrix(1, nrow = D, ncol = p_x),
    odds_atom_yx = matrix(1, nrow = D, ncol = p_y),
    odds_atom_y = rep(1, D),
    # Probability values - set to 0.5 for stability
    prob_atom_x = matrix(0.5, nrow = D, ncol = p_x),
    prob_atom_yx = matrix(0.5, nrow = D, ncol = p_y),
    prob_atom_y = rep(0.5, D),
    # Effect sums - initialize as zeros
    x_effect_sum = rep(0, D),
    y_effect_sum = rep(0, D)
  )
  
  # Improved function to generate valid initial values for latent variables
  generate_valid_latent_inits <- function(observed_totals, n_indices, population_sizes) {
    # Create output matrix with proper dimensions
    latent_values <- matrix(0, nrow = sum(n_indices[, 2] - n_indices[, 1] + 1), ncol = ncol(observed_totals))
    
    # Track positions in the output matrix
    curr_pos <- 1
    
    for(j in 1:ncol(observed_totals)) {
      for(m in 1:nrow(n_indices)) {
        start_idx <- n_indices[m, 1]
        end_idx <- n_indices[m, 2]
        cell_size <- end_idx - start_idx + 1
        
        # Get relevant population sizes
        if(start_idx <= length(population_sizes) && end_idx <= length(population_sizes)) {
          pop_sizes <- population_sizes[start_idx:end_idx]
        } else {
          cat("Warning: Index out of bounds for population sizes at cell", m, 
              "indices", start_idx, "-", end_idx, "vs max", length(population_sizes), "\n")
          pop_sizes <- rep(100, cell_size)  # Use default
        }
        
        # Get observed total for this cell
        if(m <= nrow(observed_totals)) {
          total <- observed_totals[m, j]
          if(is.na(total)) total <- 0
        } else {
          cat("Warning: Index out of bounds for observed totals at cell", m, 
              "vs max", nrow(observed_totals), "\n")
          total <- 0  # Use default
        }
        
        # Ensure pop_sizes are valid
        if(any(is.na(pop_sizes)) || any(pop_sizes <= 0)) {
          cat("Warning: Invalid population sizes for cell", m, "variable", j, "\n")
          pop_sizes[is.na(pop_sizes) | pop_sizes <= 0] <- 100  # Use default
        }
        
        # Ensure total doesn't exceed sum of population sizes
        total_pop <- sum(pop_sizes)
        if(total > total_pop) {
          cat("Warning: Total", total, "exceeds population sum", total_pop, 
              "for cell", m, "variable", j, "\n")
          total <- total_pop
        }
        
        # Calculate proportional allocation
        if(cell_size > 0 && total > 0) {
          props <- pop_sizes / total_pop
          values <- floor(props * total)
          
          # Distribute remainder to keep total exact
          remainder <- total - sum(values)
          if(remainder > 0) {
            # Sort indices by proportions
            sorted_idx <- order(props, decreasing = TRUE)
            for(i in 1:remainder) {
              idx <- sorted_idx[(i-1) %% cell_size + 1]
              if(values[idx] < pop_sizes[idx]) {
                values[idx] <- values[idx] + 1
              } else {
                # Find next available index
                for(k in 1:cell_size) {
                  next_idx <- sorted_idx[(i+k-1) %% cell_size + 1]
                  if(values[next_idx] < pop_sizes[next_idx]) {
                    values[next_idx] <- values[next_idx] + 1
                    break
                  }
                }
              }
            }
          }
          
          # Final check to ensure constraints are met
          if(sum(values) != total) {
            cat("Warning: Sum mismatch after allocation:", sum(values), "vs", total, 
                "for cell", m, "variable", j, "\n")
            
            # Find adjustable values
            under_capacity <- which(values < pop_sizes)
            over_capacity <- which(values > 0)
            
            # Adjust to match total
            diff <- total - sum(values)
            if(diff > 0 && length(under_capacity) > 0) {
              # Need to increase some values
              for(i in 1:diff) {
                idx <- under_capacity[(i-1) %% length(under_capacity) + 1]
                if(values[idx] < pop_sizes[idx]) {
                  values[idx] <- values[idx] + 1
                }
              }
            } else if(diff < 0 && length(over_capacity) > 0) {
              # Need to decrease some values
              for(i in 1:abs(diff)) {
                idx <- over_capacity[(i-1) %% length(over_capacity) + 1]
                if(values[idx] > 0) {
                  values[idx] <- values[idx] - 1
                }
              }
            }
          }
        } else {
          # For zero total or empty cell, all values are zero
          values <- rep(0, cell_size)
        }
        
        # Final check - all values must be non-negative and not exceed population
        for(i in 1:length(values)) {
          if(values[i] < 0) values[i] <- 0
          if(values[i] > pop_sizes[i]) values[i] <- pop_sizes[i]
        }
        
        # Store in output matrix
        if(start_idx <= nrow(latent_values) && end_idx <= nrow(latent_values)) {
          latent_values[start_idx:end_idx, j] <- values
        } else {
          cat("Warning: Cannot assign latent values for indices", start_idx, "-", end_idx, 
              "vs matrix size", nrow(latent_values), "\n")
        }
      }
    }
    
    return(latent_values)
  }
  
  # Initialize latent variables
  inits$x_latent <- generate_valid_latent_inits(
    covar_x[bookkeeping$J_x + 1:nrow(xlatent_ind), ],
    xlatent_ind,
    n_x_latent
  )
  
  inits$yx_latent <- generate_valid_latent_inits(
    covar_y[bookkeeping$J_y + 1:nrow(ylatent_ind), 1:p_y],
    ylatent_ind,
    n_y_latent
  )
  
  inits$y_latent <- generate_valid_latent_inits(
    matrix(covar_y[bookkeeping$J_y + 1:nrow(ylatent_ind), p_y + 1], ncol=1),
    ylatent_ind,
    n_y_latent
  )[,1]
  
  # # Verify latent variables
  # verify_latent_sums(
  #   inits$x_latent, 
  #   covar_x[bookkeeping$J_x + 1:nrow(xlatent_ind), ],
  #   xlatent_ind,
  #   "X"
  # )
  # 
  # verify_latent_sums(
  #   inits$yx_latent,
  #   covar_y[bookkeeping$J_y + 1:nrow(ylatent_ind), 1:p_y],
  #   ylatent_ind,
  #   "Y covariate"
  # )
  # 
  # verify_latent_sums(
  #   matrix(inits$y_latent, ncol=1),
  #   matrix(covar_y[bookkeeping$J_y + 1:nrow(ylatent_ind), p_y + 1], ncol=1),
  #   ylatent_ind,
  #   "Y outcome"
  # )
  
  # Return model inputs
  return(list(
    constants = constants,
    data = data_list,
    inits = inits
  ))
}

# Updated run_nimble_model_binomial function with debugging
run_nimble_model_binomial <- function(constants, data, inits, sim_metadata, niter = 30000, 
                                      nburnin = 5000, nchains = 3, save_plots = TRUE, output_dir = NULL) {
  source('nimble_abrm_binomial_0508.R')
  cat("\nBuilding NIMBLE model...\n")
  
  # First validate inputs
  # validate_nimble_inputs(constants, data, inits)
  
  # Build model
  Rmodel <- nimbleModel(abrm_binomial, constants, data, inits, calculate = FALSE)
  
  # # Debug model structure
  # debug_model_structure(Rmodel, constants, data, inits)
  # check_nimble_nodes(Rmodel)
  # check_dmfnchypg_nodes(Rmodel)
  # 
  # cat("\nApplying numerical stability fixes...\n")
  # Rmodel <- ensure_integer_values(Rmodel)
  # Rmodel <- verify_hypergeometric_constraints(Rmodel, constants)
  # 
  # Calculate initial log-likelihood
  cat("Initial log-likelihood: ", Rmodel$calculate(), "\n")
  
  # Check for NA log probabilities
  cat("\n=== Checking for NA log probabilities ===\n")
  
  # Check x
  x_logprob <- Rmodel[["logProb_x"]]
  if(!is.null(x_logprob) && any(is.na(x_logprob))) {
    cat("WARNING: NA values in logProb_x\n")
    na_indices <- which(is.na(x_logprob), arr.ind = TRUE)
    cat("NA indices:\n")
    print(head(na_indices, 10))
  }
  
  # Check x_latent
  x_latent_logprob <- Rmodel[["logProb_x_latent"]]
  if(!is.null(x_latent_logprob) && any(is.na(x_latent_logprob))) {
    cat("WARNING: NA values in logProb_x_latent\n")
    na_indices <- which(is.na(x_latent_logprob), arr.ind = TRUE)
    cat("NA indices:\n")
    print(head(na_indices, 10))
    
    # For each NA, check the corresponding node
    for(i in 1:min(5, nrow(na_indices))) {
      row <- na_indices[i, 1]
      col <- na_indices[i, 2]
      
      # Find which grid cell this corresponds to
      for(m in 1:nrow(constants$xlatent_ind)) {
        if(row >= constants$xlatent_ind[m,1] && row <= constants$xlatent_ind[m,2]) {
          cat(sprintf("\nNA at x_latent[%d, %d] corresponds to grid cell %d\n", row, col, m))
          
          # Check the inputs to dmfnchypg for this cell
          total <- Rmodel$x[constants$J_x + m, col]
          start_idx <- constants$xlatent_ind[m,1]
          end_idx <- constants$xlatent_ind[m,2]
          odds <- Rmodel$odds_atom_x[(constants$J_x + start_idx):(constants$J_x + end_idx), col]
          ni <- constants$n_x_latent[start_idx:end_idx]
          x_vals <- Rmodel$x_latent[start_idx:end_idx, col]
          
          cat("  Total:", total, "\n")
          cat("  Odds:", odds, "\n")
          cat("  ni:", ni, "\n")
          cat("  x values:", x_vals, "\n")
          cat("  Sum of x:", sum(x_vals), "\n")
          
          break
        }
      }
    }
  }
  
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
                                     scale = 0.1,
                                     adaptInterval = 50))
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc, project = compmod)
  
  # Set MCMC diagnostics
  nimbleOptions(MCMCprogressBar = TRUE, MCMCsaveHistory = TRUE)
  
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
run_abrm_binomial <- function(data, output_file, sim_metadata, niter = 30000, nburnin = 5000) {
  data=sim_data
  data <- preprocess_data(data)
  
  # Prepare spatial bookkeeping
  bookkeeping <- prepare_spatial_bookkeeping(data)
  
  # Prepare adjacency matrices
  adjacency <- prepare_adjacency_matrices(bookkeeping$gridy_yorder, bookkeeping$gridx_xorder)
  
  # Prepare NIMBLE model inputs
  nimble_inputs <- prepare_nimble_inputs_binomial(bookkeeping, adjacency, data)
  
  # Run NIMBLE model
  mcmc.out <- run_nimble_model_binomial(
    constants = nimble_inputs$constants, 
    data = nimble_inputs$data, 
    inits = nimble_inputs$inits, 
    sim_metadata
  )
  
  # Save output
  saveRDS(mcmc.out, file = output_file)
  
  return(mcmc.out)
}

# Calculate MCMC diagnostics
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

# # Simulate binomial data
# sim_data <- simulate_misaligned_binomial_data(
#   res1 = c(5, 5),
#   res2 = c(10, 10),
#   n_covariates_x = 3,
#   n_covariates_y = 4,
#   x_correlation = 0.3,
#   y_correlation = 0.7,
#   beta_x = c(0.03, -0.01, 0.06),
#   beta_y = c(0.04, 0.01, -0.05, 0.02),
#   pop_min = 100,
#   pop_max = 500
# )
# 
# # Set simulation metadata
# sim_metadata <- list(
#   sim_number = 1,
#   x_correlation = 0.3,
#   y_correlation = 0.7
# )
# 
# # Run binomial ABRM on simulated data
# source('nimble_abrm_binomial_0508.R')
# result <- run_abrm_binomial(
#   sim_data,
#   'binomial_model_output.rds',
#   sim_metadata,
#   niter = 500,
#   nburnin = 100
# )