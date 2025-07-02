### NIMBLE Code for Normal/Gaussian ABRM - FIXED VERSION
## Two spatial grids, X and Y, from which all variables arise
## Nothing fully measured at the atom level
## All variables are Normal/Gaussian distributed

abrm <- nimbleCode({
  
  for(j in 1:p_x) {
    beta_0_x[j] ~ dnorm(0, sd = 1)  # Intercepts for X covariates
  }
  
  for(j in 1:p_y) {
    beta_0_yx[j] ~ dnorm(0, sd = 1) # Intercepts for Y covariates
  }
  
  beta_0_y ~ dnorm(0, sd = 1)  # Intercept for Y outcome
  
  for (j in 1:(p_x+p_y)){
    beta_y[j] ~ dnorm(0, sd = 0.5) # Effects of X and Y covariates on Y
  }
  
  # Variance parameters for Normal distributions
  for(j in 1:p_x) {
    sigma2_x[j] ~ dinvgamma(2, 1)  # Variance for X covariates
  }
  
  for(j in 1:p_y) {
    sigma2_yx[j] ~ dinvgamma(2, 1) # Variance for Y covariates
  }
  
  sigma2_y ~ dinvgamma(2, 1)  # Variance for Y outcome
  
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
  
  ## Calculate mean parameters for atom-level X covariates
  for(j in 1:p_x) {
    for(d in 1:D) {
      mu_atom_x[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[x_to_atom[d]]
    }
  }
  
  ## Observed atom-level values (where we directly observe at atom level)
  for (j in 1:p_x){
    for (i in 1:J_x){
      x[i,j] ~ dnorm(mu_atom_x[i,j], sd = sqrt(sigma2_x[j])) 
    }
  }
  
  ## CONG ET AL ALGORITHM FOR X COVARIATES
  ## Sample latent atom-level values using constrained multivariate normal
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){
      nat_x <- length(xlatent_ind[m,1]:xlatent_ind[m,2])
      
      ## Step 1: Sample from unconstrained multivariate normal
      # w_x[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ dmnorm_chol(
      #   mu_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j], 
      #   cholesky = sqrt(sigma2_x[j])*diag(nat_x), 
      #   prec_param = FALSE
      # )
      
      # First define the covariance matrix as deterministic
      cov_x[xlatent_ind[m,1]:xlatent_ind[m,2], xlatent_ind[m,1]:xlatent_ind[m,2], j] <- sigma2_x[j] * diag(nat_x)
      
      # Then use it in the distribution
      w_x[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ dmnorm(
        mu_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j], 
        cov = cov_x[xlatent_ind[m,1]:xlatent_ind[m,2], xlatent_ind[m,1]:xlatent_ind[m,2], j]
      )
      
      ## Step 2: Project onto constraint hyperplane (Cong Algorithm 2)
      x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],j] <- 
        w_x[xlatent_ind[m,1]:xlatent_ind[m,2],j] + 
        (1/nat_x) * (x[J_x+m,j] - sum(w_x[xlatent_ind[m,1]:xlatent_ind[m,2],j])) * rep(1, nat_x)
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
  
  ## Calculate mean parameters for atom-level Y covariates
  for(j in 1:p_y) {
    for(d in 1:D) {
      mu_atom_yx[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[y_to_atom[d]]
    }
  }
  
  ## Observed atom level values (first J_y rows)
  for (j in 1:p_y){
    for (i in 1:J_y){ 
      yx_obs[i,j] ~ dnorm(mu_atom_yx[i,j], sd = sqrt(sigma2_yx[j]))
    }
  }
  
  ## CONG ET AL ALGORITHM FOR Y COVARIATES  
  ## Sample latent atom-level values using constrained multivariate normal
  for (j in 1:p_y){
    for (m in 1:(S_y-J_y)){
      nat_yx <- length(ylatent_ind[m,1]:ylatent_ind[m,2])
      
      ## Step 1: Sample from unconstrained multivariate normal
      # w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ dmnorm_chol(
      #   mu_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j], 
      #   cholesky = sqrt(sigma2_yx[j])*diag(nat_yx), 
      #   prec_param = FALSE
      # )
      
      # First define the covariance matrix as deterministic
      cov_yx[ylatent_ind[m,1]:ylatent_ind[m,2], ylatent_ind[m,1]:ylatent_ind[m,2], j] <- sigma2_yx[j] * diag(nat_yx)
      
      # Then use it in the distribution
      w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ dmnorm(
        mu_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j], 
        cov = cov_yx[ylatent_ind[m,1]:ylatent_ind[m,2], ylatent_ind[m,1]:ylatent_ind[m,2], j]
      )
      
      ## Step 2: Project onto constraint hyperplane (Cong Algorithm 2)
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] <- 
        w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j] + 
        (1/nat_yx) * (yx_obs[J_y+m,j] - sum(w_yx[ylatent_ind[m,1]:ylatent_ind[m,2],j])) * rep(1, nat_yx)
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
  
  ## Calculate mean parameters for atom-level Y outcome
  for(d in 1:D) {
    
    ## Convert covariates to rates/proportions by dividing by population
    for(k_x in 1:p_x) {
      x_mean[d, k_x] <- x_atomord[y_to_atom[d], k_x] / Pop_atom_true[y_to_atom[d]]
    }
    
    for(k_y in 1:p_y) {
      y_mean[d, k_y] <- temp_yx[d, k_y] / Pop_atom_true[y_to_atom[d]]
    }
    
    x_effect_sum[d] <- inprod(x_mean[d, 1:p_x], beta_y[1:p_x])
    y_effect_sum[d] <- inprod(y_mean[d, 1:p_y], beta_y[(p_x+1):(p_x+p_y)])
    
    mu_atom_y[d] <- beta_0_y + 
      x_effect_sum[d] + y_effect_sum[d] +
      phi_y[expand_y[d]] + offs_y[y_to_atom[d]]
  }
  
  ## Observed atom level values
  for (i in 1:J_y){
    y_obs[i] ~ dnorm(mu_atom_y[i], sd = sqrt(sigma2_y))
  }
  
  ## CONG ET AL ALGORITHM FOR Y OUTCOME
  ## Sample latent atom-level values using constrained multivariate normal
  for (m in 1:(S_y-J_y)){
    nat_y <- length(ylatent_ind[m,1]:ylatent_ind[m,2])
    
    ## Step 1: Sample from unconstrained multivariate normal
    # w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmnorm_chol(
    #   mu_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])], 
    #   cholesky = sqrt(sigma2_y)*diag(nat_y), 
    #   prec_param = FALSE
    # )
    
    # First define the covariance matrix as deterministic
    cov_y[ylatent_ind[m,1]:ylatent_ind[m,2], ylatent_ind[m,1]:ylatent_ind[m,2]] <- sigma2_y * diag(nat_y)
    
    # Then use it in the distribution
    w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmnorm(
      mu_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])], 
      cov = cov_y[ylatent_ind[m,1]:ylatent_ind[m,2], ylatent_ind[m,1]:ylatent_ind[m,2]]
    )
    
    ## Step 2: Project onto constraint hyperplane (Cong Algorithm 2)
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] <- 
      w_y[ylatent_ind[m,1]:ylatent_ind[m,2]] + 
      (1/nat_y) * (y_obs[J_y+m] - sum(w_y[ylatent_ind[m,1]:ylatent_ind[m,2]])) * rep(1, nat_y)
  }
})