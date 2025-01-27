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
    ## intercepts for each imputation model
    beta_0_x[j] ~ dnorm(0,1)
    for (k in 1:(p_x-1)){
      ## regression coefficients for each imputation model
      beta_x[j,k] ~ dnorm(0,1)
    }
  }
  
  ## betas for y-grid predictor imputation models
  for (j in 1:p_y){
    ## intercepts for each imputation model
    beta_0_yx[j] ~ dnorm(0,1)
    for (k in 1:(p_x+p_y-1)){
      ## regression coefficients for each imputation model
      beta_yx[j,k] ~ dnorm(0,1)
    }
  }
  
  ## betas for outcome model
  beta_0_y ~ dnorm(0,1)
  for (j in 1:(p_x+p_y)){
    beta_y[j] ~ dnorm(0,1)
  }
  
  
  ## outcome model spatial random effect precision parameter
  tau_y ~ dgamma(0.001, 0.001)
  ## vector of spatial random effects for the outcome model
  phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y, zero_mean = 1)
  
  ## allow different spatial random effects for each imputation model
  ## predictors measured on x-grid
  for (j in 1:p_x){
    tau_x[j] ~ dgamma(0.001, 0.001)
    ## vector of spatial random effects for x-grid predictor j's imputation model
    phi_x[1:S_x,j] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x[j], zero_mean = 1)
  }
  
  ## predictors measured on y-grid
  for (j in 1:p_y){
    tau_yx[j] ~ dgamma(0.001, 0.001)
    ## vector of spatial random effects for y-grid predictor j's imputation model
    phi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j], zero_mean = 1)
  }
  
  ################################
  ## MODELING X-GRID COVARIATES ##
  ################################
  
  ## observed atom level values are poisson ##
  for (j in 1:p_x){
    for (i in 1:J_x){
      x[i,j] ~ dpois(lambda_atom_x[i,j])
    }
  }
  
  ## missing (latent) atom level values are multinomial conditional on X grid sum
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){
      prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2],j]<-lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j]/sum(lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j])
      x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ dmulti(size=x[J_x+m,j],prob=prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2],j])
    }
  }
  
  ## make a matrix that holds the atom-level x values (observed and imputed) called temp_x
  ## also a re-ordered matrix that uses atom order to plug into the y model called x_atomord
  for (j in 1:p_x){
    temp_x[1:D,j]<-c(x[1:J_x,j],x_latent[1:(D-J_x),j])
    for (d in 1:D){
      x_atomord[d,j]<-temp_x[x_reorder[d],j]
    }
  }
  
  ## expected values for latent atom-level x-grid covariates (in x-order)
  for (j in 1:p_x){
    for (d in 1:D){
      lambda_atom_x[d,j]<-exp(beta_0_x[j]+sum(temp_x[d,(1:p_x)[-j]]*beta_x[j,1:(p_x-1)])+
                                phi_x[expand_x[d],j]+
                                offs_x[x_to_atom[d]])
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## observed atom level values are poisson ##
  for (j in 1:p_y){
    for (i in 1:J_y){
      y[i,j] ~ dpois(lambda_atom_yx[i,j])
    }
  }
  
  ## missing (latent) atom level values are multinomial conditional on Y grid sum
  for (j in 1:p_y){
    for (m in 1:(S_y-J_y)){
      prob_vecyx[ylatent_ind[m,1]:ylatent_ind[m,2],j]<-lambda_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j]/sum(lambda_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j])
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ dmulti(size=y[J_y+m,j],prob=prob_vecyx[ylatent_ind[m,1]:ylatent_ind[m,2],j])
    }
  }
  
  ## make a matrix that holds the atom-level y-grid covariate values (observed and imputed) called temp_yx
  for (j in 1:p_y){
    temp_yx[1:D,j]<-c(y[1:J_y,j],yx_latent[1:(D-J_y),j])
  }
  
  ## expected values for latent atom-level x-grid covariates (in x-order)
  for (j in 1:p_y){
    for (d in 1:D){
      lambda_atom_yx[d,j]<-exp(beta_0_yx[j]+sum(x_atomord[y_to_atom[d],(1:p_x)]*beta_yx[j,1:(p_x)])+
                                sum(temp_yx[d,(1:p_y)[-j]]*beta_yx[j,(p_x+1):(p_x+p_y-1)])+
                                phi_yx[expand_y[d],j]+
                                offs_y[y_to_atom[d]])
    }
  }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################
  
  ## expected values for latent atom-level outcomes
  for (d in 1:D){
    lambda_atom_y[d]<-exp(beta_0_y+
                            sum(x_atomord[y_to_atom[d]]*beta_y[1:p_x])+
                            sum(temp_yx[d,1:p_y]*beta_y[(p_x+1):(p_x+p_y)])+
                            phi_y[expand_y[d]]+
                            offs_y[y_to_atom[d]])
  }
  
  ## observed atom level values are poisson
  for (i in 1:J_y){
    y[i,p_y+1] ~ dpois(lambda_atom_y[i])
  }
  
  ## missing (latent) atom level values are multinomial conditional on Y grid sum
  for (m in 1:(S_y-J_y)){
    prob_vecy[ylatent_ind[m,1]:ylatent_ind[m,2]]<-lambda_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]/sum(lambda_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmulti(size=y[J_y+m,p_y+1],prob=prob_vecy[ylatent_ind[m,1]:ylatent_ind[m,2]])
  }
  
})
