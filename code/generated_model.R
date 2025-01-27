
abrm <- nimbleCode({
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
  
  ## Spatial random effects
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
  
  ## Initialize base lambda values
  for(j in 1:p_x) {
    for(d in 1:D) {
      linear_pred_x_base[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[x_to_atom[d]]
      lambda_atom_x_init[d,j] <- exp(linear_pred_x_base[d,j])
    }
  }
  
  ## model observed atom-level values
  for (j in 1:p_x){
    for (i in 1:J_x){
      x[i,j] ~ dpois(lambda_atom_x_init[i,j]) 
    }
  }
  
  ## multinomial probabilities and latent values
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){
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
  
  ## create temp matrices
  for (j in 1:p_x){
    for(i in 1:J_x) {
      temp_x[i,j] <- x[i,j]
    }
    for(i in 1:(D-J_x)) {
      temp_x[J_x+i,j] <- x_latent[i,j]
    }
    for(d in 1:D) {
      x_atomord[d,j] <- temp_x[x_reorder[d],j]
    }
  }
  
  ## X-grid final calculations
  for (j in 1:p_x){
    for (d in 1:D){
      beta_sum_x[d,j] <- 0
      lambda_atom_x[d,j] <- exp(linear_pred_x_base[d,j] + beta_sum_x[d,j])
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## Initialize Y base lambda values
  for(j in 1:p_y) {
    for(d in 1:D) {
      linear_pred_yx_base[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[y_to_atom[d]]
      lambda_atom_yx_init[d,j] <- exp(linear_pred_yx_base[d,j])
    }
  }
  
  ## Y observed values
  for (j in 1:p_y){
    for (i in 1:J_y){ 
      y[i,j] ~ dpois(lambda_atom_yx_init[i,j])
    }
  }
  
  ## Y multinomial probabilities
  for (j in 1:p_y){
    for (m in 1:(S_y-J_y)){
      for(k in 1:max_atoms_y) {
        numer_yx[m,k,j] <- lambda_atom_yx_init[y_expanded_indices[m,k],j]
      }
      denom_yx[m,j] <- sum(numer_yx[m,1:max_atoms_y,j])
      for(k in 1:max_atoms_y) {
        prob_vecyx[m,k,j] <- numer_yx[m,k,j]/denom_yx[m,j]
      }
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ 
        dmulti(size=y[J_y+m,j], prob=prob_vecyx[m,1:max_atoms_y,j])
    }
  }
  
  ## Y temp matrices
  for (j in 1:p_y){
    for(i in 1:J_y) {
      temp_yx[i,j] <- y[i,j]
    }
    for(i in 1:(D-J_y)) {
      temp_yx[J_y+i,j] <- yx_latent[i,j]
    }
  }
  
  ## Y-grid final calculations
  for (j in 1:p_y){
    for (d in 1:D){
      beta_sum_yx_x[d,j] <- x_atomord[y_to_atom[d],1] * beta_yx[j,1] + x_atomord[y_to_atom[d],2] * beta_yx[j,2] + x_atomord[y_to_atom[d],3] * beta_yx[j,3]
      beta_sum_yx_y[d,j] <- temp_yx[d,y_covar_indices[j,1]] * beta_yx[j,4] + temp_yx[d,y_covar_indices[j,2]] * beta_yx[j,5] + temp_yx[d,y_covar_indices[j,3]] * beta_yx[j,6] + temp_yx[d,y_covar_indices[j,4]] * beta_yx[j,7]
      lambda_atom_yx[d,j] <- exp(linear_pred_yx_base[d,j] + beta_sum_yx_x[d,j] + beta_sum_yx_y[d,j])
    }
  }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################
  
  ## Initialize final Y base lambda values
  for(d in 1:D) {
    linear_pred_y_base[d] <- beta_0_y + phi_y[expand_y[d]] + offs_y[y_to_atom[d]]
    lambda_atom_y_init[d] <- exp(linear_pred_y_base[d])
  }
  
  ## Final Y observed values
  for (i in 1:J_y){
    y[i,p_y+1] ~ dpois(lambda_atom_y_init[i])
  }
  
  ## Final Y multinomial probabilities
  for (m in 1:(S_y-J_y)){
    for(k in 1:max_atoms_y) {
      numer_y[m,k] <- lambda_atom_y_init[y_expanded_indices[m,k]]
    }
    denom_y[m] <- sum(numer_y[m,1:max_atoms_y])
    for(k in 1:max_atoms_y) {
      prob_vecy[m,k] <- numer_y[m,k]/denom_y[m]
    }
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ 
      dmulti(size=y[J_y+m,p_y+1], prob=prob_vecy[m,1:max_atoms_y])
  }
  
  ## Final Y calculations
  for(d in 1:D) {
    beta_sum_y_x[d] <- x_atomord[y_to_atom[d],1] * beta_y[1] + x_atomord[y_to_atom[d],2] * beta_y[2] + x_atomord[y_to_atom[d],3] * beta_y[3]
    beta_sum_y_y[d] <- temp_yx[d,1] * beta_y[4] + temp_yx[d,2] * beta_y[5] + temp_yx[d,3] * beta_y[6] + temp_yx[d,4] * beta_y[7] + temp_yx[d,5] * beta_y[8]
    lambda_atom_y[d] <- exp(linear_pred_y_base[d] + beta_sum_y_x[d] + beta_sum_y_y[d])
  }
})
