# Simpler version
prepare_nimble_inputs_ <- function(bookkeeping, adjacency, data) {
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
  
  print("\nCreating covariate matrices...")
  covar_x <- create_covariate_matrix(gridx_xorder[order(gridx_xorder$ID_x), ], 
                                     bookkeeping$x_vars, 
                                     unique(data$atoms$ID_x), 
                                     nrow(gridx_xorder))
  
  covar_y <- create_covariate_matrix(gridy_yorder[order(gridy_yorder$ID_y), ], 
                                     bookkeeping$y_vars, 
                                     unique(data$atoms$ID_y), 
                                     nrow(gridy_yorder))
  
  # Create index matrices
  p_x <- ncol(covar_x)
  p_y <- ncol(covar_y) - 1  # subtract y column
  
  # Calculate max atoms per grid
  max_atoms_x <- max(bookkeeping$x_latentind[,2] - bookkeeping$x_latentind[,1] + 1)
  max_atoms_y <- max(bookkeeping$y_latentind[,2] - bookkeeping$y_latentind[,1] + 1)
  
  # Expanded indices for multinomial calculations
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
  
  x_expanded_indices <- create_expanded_indices(bookkeeping$x_latentind,
                                                bookkeeping$J_x, max_atoms_x)
  y_expanded_indices <- create_expanded_indices(bookkeeping$y_latentind,
                                                bookkeeping$J_y, max_atoms_y)
  
  # Scale matrices for Wishart priors
  R_x <- diag(p_x)  # For X grid covariates 
  R_yx <- diag(p_y) # For Y grid covariates
  
  constants <- list(
    # Dimensions
    p_x = p_x,
    p_y = p_y,
    D = as.integer(n_atoms),
    S_x = nrow(gridx_xorder),
    S_y = nrow(gridy_yorder),
    J_x = bookkeeping$J_x,
    J_y = bookkeeping$J_y,
    
    # Mapping between spaces
    x_to_atom = bookkeeping$x_to_atom,
    y_to_atom = bookkeeping$y_to_atom,
    expand_x = bookkeeping$expand_x,
    expand_y = bookkeeping$expand_y,
    x_reorder = bookkeeping$x_reorder,
    
    # Latent indices
    xlatent_ind = bookkeeping$x_latentind,
    ylatent_ind = bookkeeping$y_latentind,
    x_expanded_indices = x_expanded_indices,
    y_expanded_indices = y_expanded_indices,
    max_atoms_x = max_atoms_x,
    max_atoms_y = max_atoms_y,
    
    # CAR structures
    num_x = as.numeric(as.carAdjacency(adjacency$W_x)$num),
    weights_x = as.numeric(as.carAdjacency(adjacency$W_x)$weights),
    adj_x = as.numeric(as.carAdjacency(adjacency$W_x)$adj),
    num_y = as.numeric(as.carAdjacency(adjacency$W_y)$num),
    weights_y = as.numeric(as.carAdjacency(adjacency$W_y)$weights),
    adj_y = as.numeric(as.carAdjacency(adjacency$W_y)$adj),
    
    # Add explicit dimensions for matrix nodes
    temp_x = array(0, dim = c(n_atoms, p_x)),
    temp_yx = array(0, dim = c(n_atoms, p_y)),
    
    # Wishart prior matrices
    R_x = R_x,      # Add scale matrix for X covariates
    R_yx = R_yx,    # Add scale matrix for Y covariates
    df_x = p_x,     # Degrees of freedom for Wishart
    df_yx = p_y   # Degrees of freedom for Wishart
  )
  
  data <- list(
    x = covar_x,
    y = covar_y,
    offs_x = rep(0, n_atoms),
    offs_y = rep(0, n_atoms)
  )
  
  inits <- list(
    beta_0_y = 0,
    beta_y = rep(0, p_x + p_y),
    
    # Spatial random effects
    tau_x = rep(1, p_x),
    tau_y = 1,
    tau_yx = rep(1, p_y),
    
    # Initialize base spatial random effects (psi)
    psi_x = matrix(0, nrow = constants$S_x, ncol = p_x),
    psi_yx = matrix(0, nrow = constants$S_y, ncol = p_y),
    
    # Matrices for spatial correlation
    Prec_x = diag(p_x),  # Identity matrix
    Prec_yx = diag(p_y), # Identity matrix
    phi_x = matrix(0, nrow = constants$S_x, ncol = p_x),
    phi_y = rep(0, constants$S_y),
    phi_yx = matrix(0, nrow = constants$S_y, ncol = p_y),
    
    # Initialize matrices with proper dimensions
    temp_x = array(0, dim = c(n_atoms, p_x)),
    temp_yx = array(0, dim = c(n_atoms, p_y)),
    
    # Add initial values for latent counts
    x_latent = matrix(1, nrow = n_atoms - constants$J_x, ncol = p_x),
    y_latent = rep(1, n_atoms - constants$J_y),
    yx_latent = matrix(1, nrow = n_atoms - constants$J_y, ncol = p_y)
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

