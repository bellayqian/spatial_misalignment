### Nimble Code
## two spatial grids, X and Y, from which all variables arise
## nothing fully measured at the atom level
## all variables are poisson

abrm <- nimbleCode({
  
  ## constants
  ## p_x=number of x-grid predictors (does not include an intercept)
  ## p_y=number of y-grid predictors (does not include an intercept)
  ## D=number of atoms
  ## S_x= total number of X grids
  ## S_y= total number of Y grids
  ## J_x= first J_x x values are atom-equivalents
  ## J_y= first J_y y values are atom-equivalents
  ## x_to_atom= D-length vector of indices of atoms (in atom-order) corresponding to each X-grid
  ## x_latentind=(S_x-J_x) x 2 matrix with start and end indices of atoms for each non-atom-equivalent x-grid
  ## expand_x= D-length vector with x indices of atoms in x-order
  ## atom_to_x= D-length vector of x indices of atoms in atom-order
  ## y_to_atom= D-length vector of indices of atoms (in atom-order) corresponding to each y-grid
  ## y_latentind=(S_y-J_y) x 2 matrix with start and end indices of atoms for each non-atom-equivalent y-grid
  ## expand_y= D-length vector with y indices of atoms in y-order
  ## num = car.Adjacency$num
  ## weights = car.Adjacency$weights
  ## adj = car.Adjacency$adj
  
  ## data
  ## x is a S_x x p_x matrix ordered such that the atom-equivalents are first, then the non-atom-equivalents
  ## y is a S_y x (p_y+1) matrix with the first p_y columns the y-grid covariates and the last column y itself, ordered such that the atom-equivalents are first, then the non-atom-equivalents
  ## offs_x= X model offset
  ## offs_y= Y model offset
  
  ## betas for x-grid predictor imputation models
  for (j in 1:p_x){ 
    beta_0_x[j] ~ dnorm(0,1)
    for (k in 1:(p_x-1)){
      beta_x[j,k] ~ dnorm(0,1)
    }
  }
  
  ## betas for y-grid predictor imputation models
  for (j in 1:p_y){ 
    beta_0_yx[j] ~ dnorm(0,1)
    for (k in 1:(p_x+p_y-1)){
      beta_yx[j,k] ~ dnorm(0,1)
    }
  }
  
  ## betas for outcome model
  beta_0_y ~ dnorm(0,1)
  for (j in 1:(p_x+p_y)){
    beta_y[j] ~ dnorm(0,1)
  }
  
  ## spatial random effects
  tau_y ~ dgamma(0.001, 0.001)
  phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y, zero_mean = 1)
  
  for (j in 1:p_x){ 
    tau_x[j] ~ dgamma(0.001, 0.001)
    phi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j], zero_mean = 1)
  }
  
  for (j in 1:p_y){ 
    tau_yx[j] ~ dgamma(0.001, 0.001)
    phi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j], zero_mean = 1)
  }
  
  ################################
  ## MODELING X-GRID COVARIATES ##
  ################################
  
  ## initialize base lambda values
  for(j in 1:p_x) {
    for(d in 1:D) {
      linear_pred_x_base[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[x_to_atom[d]]
      lambda_atom_x_init[d,j] <- exp(linear_pred_x_base[d,j])
    }
  }
  
  ## model observed atom-level values
  ## observed atom level values (first J_x rows of S_x × p_x matrix)
  for (j in 1:p_x){ # p_x = 4
    for (i in 1:J_x){ # J_x = 15
      x[i,j] ~ dpois(lambda_atom_x_init[i,j]) 
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){ # S_x-J_x = 5
      # calculate numerators using initial lambdas
      for(k in 1:max_atoms_x) {
        numer_x[m,k,j] <- lambda_atom_x_init[x_expanded_indices[m,k],j]
      }
      denom_x[m,j] <- sum(numer_x[m,1:max_atoms_x,j])
      for(k in 1:max_atoms_x) {
        prob_vecx[m,k,j] <- numer_x[m,k,j]/denom_x[m,j]
      }
      x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ 
        dmulti(size=x[J_x+m,j], prob=prob_vecx[m,1:max_atoms_x,j])
    }
  }
  
  ## create atom-level matrices as derived quantities
  for (j in 1:p_x){ # p_x = 4
    # First fill in observed values
    for(i in 1:J_x) {
      temp_x[i,j] <- x[i,j]
    }
    # Then fill in latent values
    for(i in 1:(D-J_x)) {
      temp_x[J_x+i,j] <- x_latent[i,j]
    }
    # Create atom-ordered version
    for(d in 1:D) {
      x_atomord[d,j] <- temp_x[x_reorder[d],j]
    }
  }
  
  ## For X-grid covariates
  for (j in 1:p_x){
    for (d in 1:D){
      # Use intermediate sums
      beta_sum_x_temp1[d,j] <- temp_x[d,x_index_mat[j,1]] * beta_x[j,1]
      beta_sum_x_temp2[d,j] <- beta_sum_x_temp1[d,j] + 
        temp_x[d,x_index_mat[j,2]] * beta_x[j,2]
      beta_sum_x[d,j] <- beta_sum_x_temp2[d,j] + 
        temp_x[d,x_index_mat[j,3]] * beta_x[j,3]
      
      lambda_atom_x[d,j] <- exp(linear_pred_x_base[d,j] + beta_sum_x[d,j])
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## initialize base lambda values for Y covariates
  for(j in 1:p_y) {
    for(d in 1:D) {
      linear_pred_yx_base[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[y_to_atom[d]]
      lambda_atom_yx_init[d,j] <- exp(linear_pred_yx_base[d,j])
    }
  }
  
  ## observed atom level values (first J_y rows of S_y × p_y matrix)
  for (j in 1:p_y){ # p_y = 4
    for (i in 1:J_y){ 
      y[i,j] ~ dpois(lambda_atom_yx_init[i,j])
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_y){ # p_y = 4
    for (m in 1:(S_y-J_y)){ # S_y-J_y = 10
      # Calculate numerators using initial lambdas
      for(k in 1:max_atoms_y) {
        numer_yx[m,k,j] <- lambda_atom_yx_init[y_expanded_indices[m,k],j]
      }
      # Calculate denominator
      denom_yx[m,j] <- sum(numer_yx[m,1:max_atoms_y,j])
      # Normalize probabilities
      for(k in 1:max_atoms_y) {
        prob_vecyx[m,k,j] <- numer_yx[m,k,j]/denom_yx[m,j]
      }
      # Multinomial for latent values
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ 
        dmulti(size=y[J_y+m,j], prob=prob_vecyx[m,1:max_atoms_y,j])
    }
  }
  
  ## Create atom-level matrix as derived quantity
  for (j in 1:p_y){ # p_y = 4
    # Fill in observed values
    for(i in 1:J_y) {
      temp_yx[i,j] <- y[i,j]
    }
    # Fill in latent values
    for(i in 1:(D-J_y)) {
      temp_yx[J_y+i,j] <- yx_latent[i,j]
    }
  }
  
  ## For Y-grid covariates
  for (j in 1:p_y){
    for (d in 1:D){
      # X covariate effects - explicit steps
      beta_sum_yx_x_temp1[d,j] <- x_atomord[y_to_atom[d],1] * beta_yx[j,1]
      beta_sum_yx_x_temp2[d,j] <- beta_sum_yx_x_temp1[d,j] + 
        x_atomord[y_to_atom[d],2] * beta_yx[j,2]
      beta_sum_yx_x_temp3[d,j] <- beta_sum_yx_x_temp2[d,j] + 
        x_atomord[y_to_atom[d],3] * beta_yx[j,3]
      beta_sum_yx_x[d,j] <- beta_sum_yx_x_temp3[d,j] + 
        x_atomord[y_to_atom[d],4] * beta_yx[j,4]
      
      # Y covariate effects - explicit steps
      beta_sum_yx_y_temp1[d,j] <- temp_yx[d,y_index_mat[j,1]] * beta_yx[j,5]
      beta_sum_yx_y_temp2[d,j] <- beta_sum_yx_y_temp1[d,j] + 
        temp_yx[d,y_index_mat[j,2]] * beta_yx[j,6]
      beta_sum_yx_y[d,j] <- beta_sum_yx_y_temp2[d,j] + 
        temp_yx[d,y_index_mat[j,3]] * beta_yx[j,7]
      
      # Final lambda calculation
      lambda_atom_yx[d,j] <- exp(linear_pred_yx_base[d,j] + 
                                   beta_sum_yx_x[d,j] + 
                                   beta_sum_yx_y[d,j])
    }
  }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################
  
  ## initialize base lambda values for final outcome
  for(d in 1:D) {
    linear_pred_y_base[d] <- beta_0_y + phi_y[expand_y[d]] + offs_y[y_to_atom[d]]
    lambda_atom_y_init[d] <- exp(linear_pred_y_base[d])
  }
  
  ## observed atom level values (first J_y elements)
  for (i in 1:J_y){ # J_y = 25
    y[i,p_y+1] ~ dpois(lambda_atom_y_init[i])
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (m in 1:(S_y-J_y)){
    # Calculate numerators using initial lambdas
    for(k in 1:max_atoms_y) {
      numer_y[m,k] <- lambda_atom_y_init[y_expanded_indices[m,k]]
    }
    # Calculate denominator
    denom_y[m] <- sum(numer_y[m,1:max_atoms_y])
    # Normalize probabilities
    for(k in 1:max_atoms_y) {
      prob_vecy[m,k] <- numer_y[m,k]/denom_y[m]
    }
    # Multinomial for latent values
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ 
      dmulti(size=y[J_y+m,p_y+1], prob=prob_vecy[m,1:max_atoms_y])
  }
  
  ## calculate final lambda values
  for(d in 1:D) {
    # X predictor effects - explicit steps
    beta_sum_y_x_temp1[d] <- x_atomord[y_to_atom[d],1] * beta_y[1]
    beta_sum_y_x_temp2[d] <- beta_sum_y_x_temp1[d] + 
      x_atomord[y_to_atom[d],2] * beta_y[2]
    beta_sum_y_x_temp3[d] <- beta_sum_y_x_temp2[d] + 
      x_atomord[y_to_atom[d],3] * beta_y[3]
    beta_sum_y_x[d] <- beta_sum_y_x_temp3[d] + 
      x_atomord[y_to_atom[d],4] * beta_y[4]
    
    # Y predictor effects - explicit steps
    beta_sum_y_y_temp1[d] <- temp_yx[d,1] * beta_y[5]
    beta_sum_y_y_temp2[d] <- beta_sum_y_y_temp1[d] + 
      temp_yx[d,2] * beta_y[6]
    beta_sum_y_y_temp3[d] <- beta_sum_y_y_temp2[d] + 
      temp_yx[d,3] * beta_y[7]
    beta_sum_y_y[d] <- beta_sum_y_y_temp3[d] + 
      temp_yx[d,4] * beta_y[8]
    
    # Final lambda calculation
    lambda_atom_y[d] <- exp(linear_pred_y_base[d] + 
                              beta_sum_y_x[d] + 
                              beta_sum_y_y[d])
  }
})