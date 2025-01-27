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
  
  ## RCN: I DON'T THINK WE NEED THESE IMPUTATION MODEL BETAS 
  ## SINCE WE ARE NO LONGER INCLUDING THE OTHER COVARIATES IN THE IMPUTATION MODELS

  ## BQ: Absorb means, add intercepts for each process:
  for(j in 1:p_x) {
    # Center around 0.5 with reasonable variance based on simulation values (0.5, 0.3, 0.7)
    beta_0_x[j] ~ dnorm(0.5, sd = 0.2)  
  }
  
  for(j in 1:p_y) {
    # Center around 0.45 (mean of simulation values 0.4, 0.6, 0.3, 0.5)
    beta_0_yx[j] ~ dnorm(0.45, sd = 0.15)
  }

  ## Modify priors for better mixing
  # beta_0_y ~ dnorm(log(mean(y[,p_y+1])), sd = 1)  # More informative prior
  beta_0_y ~ dnorm(3, sd = 0.5)  # More informative prior based on simulation offset
  for (j in 1:(p_x+p_y)){
    # For X covariates (simulation uses 0.5, 0.3, 0.7)
    beta_y[j] ~ dnorm(0.5, sd = 0.2)
  }
  
  ## spatial random effects for outcome Y
  tau_y ~ dgamma(4, 2)
  phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y)
  
  ## RCN: SINCE WE CAN'T INCLUDE THE COVARIATES IN THE IMPUTATION MODELS, I THINK A REASONABLE COMPROMISE WOULD BE TO ALLOW FOR THE SPATIAL RANDOM EFFECTS TO BE CORRELATED ACROSS COVARIATES MEASURED ON THE SAME GRID
  ## BELOW I ADDED A FEW LINES TO ALLOW FOR THAT (FOLLOWING AN APPROACH SIMILAR TO THIS PAPER: https://doi.org/10.1016/j.sste.2023.100588)
  ## spatial random effects for X-grid covariates with correlation
  for (j in 1:p_x){ 
    tau_x[j] ~ dgamma(4, 2)
    psi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j])
  }
  
  ## this part allows for correlation in the spatial random effects across x-grid covariates
  Prec_x[1:p_x,1:p_x] ~ dwish(R_x[1:p_x,1:p_x], df_x) # Wishart prior on precision matrix
  Achol_x[1:p_x,1:p_x] <- chol(inverse(Prec_x[1:p_x,1:p_x])) # Cholesky decomposition
  
  # Transform to correlated spatial effects
  for(i in 1:S_x){
    phi_x[i,1:p_x] <- Achol_x[1:p_x,1:p_x]%*%psi_x[i,1:p_x]
  }
  
  for (j in 1:p_y){ 
    tau_yx[j] ~ dgamma(2, 1)
    # tau_yx[j] ~ dgamma(0.001, 0.001)
    psi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j])
  }
  
  ## this part allows for correlation in the spatial random effects across y-grid covariates
  Prec_yx[1:p_y,1:p_y] ~ dwish(R_yx[1:p_y,1:p_y], df_yx)
  Achol_yx[1:p_y,1:p_y] <- chol(inverse(Prec_yx[1:p_y,1:p_y]))
  
  for(i in 1:S_y){
    phi_yx[i,1:p_y] <- Achol_yx[1:p_y,1:p_y]%*%psi_yx[i,1:p_y]
  }
  
  ################################
  ## MODELING X-GRID COVARIATES ##
  ################################
  
  ## model observed atom-level values
  ## observed atom level values (first J_x rows of S_x × p_x matrix)
  for (j in 1:p_x){ # p_x = 4
    for (i in 1:J_x){ # J_x = 15
      x[i,j] ~ dpois(lambda_atom_x[i,j]) 
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){ # S_x-J_x = 5
      # calculate numerators using initial lambdas
      for(k in 1:max_atoms_x) {
        numer_x[m,k,j] <- lambda_atom_x[x_expanded_indices[m,k],j]
      }
      denom_x[m,j] <- sum(numer_x[m,1:max_atoms_x,j])
      for(k in 1:max_atoms_x) {
        prob_vecx[m,k,j] <- (numer_x[m,k,j] + 1e-10)/(denom_x[m,j] + 1e-10 * max_atoms_x)
      }
      ## RCN: I DON'T THINK IN GENERAL THE DIMENSIONS IN THIS CODE WILL WORK
      ## IN THE SIMULATED DATA, I THINK ALL OF THE NON-ATOM GRIDS HAVE THE SAME NUMBER OF ATOMS
      ## BUT IN REAL APPLICATIONS THAT WON'T GENERALLY BE THE CASE
      ## THEN THE LENGTH OF (xlatent_ind[m,1]:xlatent_ind[m,2]) WON'T BE THE SAME AS (1:max_atoms_x) SO I DON'T THINK THE CODE BELOW WILL WORK
      ## I THINK WE WILL NEED TO REVERT BACK TO THE WAY I HAD INITIALLY CODED THIS PART, WHICH SHOULD WORK IN THE MORE GENERAL SETTING
      ## SAME COMMENT APPLIES FOR THE Y-GRID COVARIATES AND THE OUTCOME
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
  }
  
  # Single definition of x_atomord
  for (j in 1:p_x) {
    for(d in 1:D) {
      x_atomord[d,j] <- temp_x[x_reorder[d],j]
    }
  }
  
  ## For X-grid covariates
  for(j in 1:p_x) {
    for(d in 1:D) {
      linear_pred_x_base[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[x_to_atom[d]]
      lambda_atom_x[d,j] <- exp(linear_pred_x_base[d,j])
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## observed atom level values (first J_y rows of S_y × p_y matrix)
  for (j in 1:p_y){ # p_y = 4
    for (i in 1:J_y){ 
      y[i,j] ~ dpois(lambda_atom_yx[i,j])
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_y){ # p_y = 4
    for (m in 1:(S_y-J_y)){ # S_y-J_y = 10
      # Calculate numerators using initial lambdas
      for(k in 1:max_atoms_y) {
        numer_yx[m,k,j] <- lambda_atom_yx[y_expanded_indices[m,k],j]
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
  for(j in 1:p_y) {
    for(d in 1:D) {
      linear_pred_yx_base[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[y_to_atom[d]]
      lambda_atom_yx[d,j] <- exp(linear_pred_yx_base[d,j])
    }
  }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################
  
  ## observed atom level values (first J_y elements)
  for (i in 1:J_y){ # J_y = 25
    y[i,p_y+1] ~ dpois(lambda_atom_y[i])
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (m in 1:(S_y-J_y)){
    # Calculate numerators using initial lambdas
    for(k in 1:max_atoms_y) {
      numer_y[m,k] <- lambda_atom_y[y_expanded_indices[m,k]]
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
  ## RCN: WE STILL NEED TO INCLUDE THE COVARIATES IN THE MODEL FOR THE OUTCOME (WHICH SHOULDN'T CREATE ANY CYCLES)
  ## THE OUTCOME MODEL COEFFICIENTS ARE THE KEY QUANTITIES WE WANT TO ESIMATE, SO CAN'T OMIT THOSE
  ## I EDITTED THE CODE BELOW TO ADD THOSE BACK IN
  for(d in 1:D) {
    linear_pred_y_base[d] <- beta_0_y + 
      sum(x_atomord[y_to_atom[d], 1:p_x] * beta_y[1:p_x]) +
      sum(temp_yx[d, 1:p_y] * beta_y[(p_x+1):(p_x+p_y)]) +
      phi_y[expand_y[d]] + offs_y[y_to_atom[d]]
    lambda_atom_y[d] <- exp(linear_pred_y_base[d]) + 1e-10
  }
})