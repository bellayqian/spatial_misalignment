### Nimble Code
## two spatial grids, X and Y, from which all variables arise
## nothing fully measured at the atom level
## all variables are poisson

abrm <- nimbleCode({
  
  for(j in 1:p_x) {
    beta_0_x[j] ~ dnorm(0.5, sd = 1)  # Intercepts for X covariates
  }
  
  for(j in 1:p_y) {
    beta_0_yx[j] ~ dnorm(0.5, sd = 1) # Intercepts for Y covariates
  }

  beta_0_y ~ dnorm(3, sd = 1)  # Intercept for Y outcome
  
  for (j in 1:(p_x+p_y)){
    #beta_y[j] ~ dnorm(0.5, sd = 1) # Effects of X' and Y covariates on Y
    beta_y[j] ~ dnorm(0, sd = 1) # Effects of X' and Y covariates on Y
  }
  
  ## Spatial random effects for outcome Y
  tau_y ~ dgamma(2, 2)
  phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y)
  
  ## Spatial random effects for X-grid covariates with correlation
  for (j in 1:p_x){ 
    tau_x[j] ~ dgamma(2, 2)
    psi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j])
  }
  
  ## Correlation in the spatial random effects across x-grid covariates
  Prec_x[1:p_x,1:p_x] ~ dwish(R_x[1:p_x,1:p_x], df_x) # Wishart prior on precision matrix
  Achol_x[1:p_x,1:p_x] <- chol(inverse(Prec_x[1:p_x,1:p_x])) # Cholesky decomposition
  
  # Transform to correlated spatial effects
  for(i in 1:S_x){
    phi_x[i,1:p_x] <- Achol_x[1:p_x,1:p_x]%*%psi_x[i,1:p_x]
  }
  
  for (j in 1:p_y){ 
    tau_yx[j] ~ dgamma(2, 2)
    psi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j])
  }
  
  ## Correlation in the spatial random effects across y-grid covariates
  Prec_yx[1:p_y,1:p_y] ~ dwish(R_yx[1:p_y,1:p_y], df_yx)
  Achol_yx[1:p_y,1:p_y] <- chol(inverse(Prec_yx[1:p_y,1:p_y]))
  
  for(i in 1:S_y){
    phi_yx[i,1:p_y] <- Achol_yx[1:p_y,1:p_y]%*%psi_yx[i,1:p_y]
  }
  
  ################################
  ## MODELING X-GRID COVARIATES ##
  ################################
  ## For X-grid covariates (calculate lambda values first)
  for(j in 1:p_x) {
    for(d in 1:D) {
      linear_pred_x_base[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[x_to_atom[d]]
      lambda_atom_x[d,j] <- exp(linear_pred_x_base[d,j]) + 1e-10
    }
  }
  
  ## observed atom-level values (where we directly observe at atom level)
  for (j in 1:p_x){
    for (i in 1:J_x){
      x[i,j] ~ dpois(lambda_atom_x[i,j]) 
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){
      # Calculate probabilities directly from lambda values for atoms in this grid cell
      prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2],j] <- 
        lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j] / 
        sum(lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j])
      
      # Multinomial for latent values
      x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ 
        dmulti(size=x[J_x+m,j], 
               prob=prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2],j])
    }
  }
  
  ## create atom-level matrices as derived quantities
  for (j in 1:p_x){
    for(i in 1:J_x) {
      temp_x[i,j] <- x[i,j] # First fill in observed values
    }
    for(i in 1:(D-J_x)) {
      temp_x[J_x+i,j] <- x_latent[i,j] # Then fill in latent values
    }
  }
  
  # Single definition of x_atomord
  for (j in 1:p_x) {
    for(d in 1:D) {
      x_atomord[d,j] <- temp_x[x_reorder[d],j]
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## For Y-grid covariates (calculate lambda values first)
  for(j in 1:p_y) {
    for(d in 1:D) {
      linear_pred_yx_base[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[y_to_atom[d]]
      lambda_atom_yx[d,j] <- exp(linear_pred_yx_base[d,j]) + 1e-10
    }
  }
  
  ## observed atom level values (first J_y rows of S_y × p_y matrix)
  for (j in 1:p_y){
    for (i in 1:J_y){ 
      yx_obs[i,j] ~ dpois(lambda_atom_yx[i,j])
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_y){
    for (m in 1:(S_y-J_y)){
      # Calculate probabilities directly from lambda values for atoms in this grid cell
      prob_vecyx[ylatent_ind[m,1]:ylatent_ind[m,2],j] <- 
        lambda_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j] / 
        sum(lambda_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j])
      
      # Multinomial for latent values
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ 
        dmulti(size=yx_obs[J_y+m,j],
               prob=prob_vecyx[ylatent_ind[m,1]:ylatent_ind[m,2],j])
    }
  }
  
  ## Create atom-level Y covariates
  for (j in 1:p_y){
    for(i in 1:J_y) {
      temp_yx[i,j] <- yx_obs[i,j]
    }
    for(i in 1:(D-J_y)) {
      temp_yx[J_y+i,j] <- yx_latent[i,j]
    }
  }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################
  
  ## calculate final lambda values for Y
  for(d in 1:D) {
    linear_pred_y_base[d] <- beta_0_y + 
      sum(x_atomord[y_to_atom[d], 1:p_x] * beta_y[1:p_x]) +
      sum(temp_yx[d, 1:p_y] * beta_y[(p_x+1):(p_x+p_y)]) +
      phi_y[expand_y[d]] + offs_y[y_to_atom[d]]
    lambda_atom_y[d] <- exp(linear_pred_y_base[d]) + 1e-10
  }
  
  ## observed atom level values
  for (i in 1:J_y){
    y_obs[i] ~ dpois(lambda_atom_y[i])
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (m in 1:(S_y-J_y)){
    # Calculate raw probabilities for each set of atoms in this grid cell
    for(k in ylatent_ind[m,1]:ylatent_ind[m,2]) {
      raw_probs_y[k] <- lambda_atom_y[J_y + k]
    }
    
    # Calculate sum for normalization (for this specific m)
    sum_probs_y[m] <- sum(raw_probs_y[ylatent_ind[m,1]:ylatent_ind[m,2]])
    
    # Calculate normalized probabilities
    for(k in ylatent_ind[m,1]:ylatent_ind[m,2]) {
      prob_vecy[k] <- raw_probs_y[k]/sum_probs_y[m]
    }
    
    # Multinomial for latent values
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ 
      dmulti(size=y_obs[J_y+m], 
             prob=prob_vecy[ylatent_ind[m,1]:ylatent_ind[m,2]])
  }
})