## assumes X and Y have different spatial grids
## rest of predictors are at atom level
## X and Y are both poisson

abrm <- nimbleCode({
  
  ## constants
  ## p=number of atom-level predictors (excludes X)
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
  ## x is ordered such that the atom-equivalents are first, then the non-atom-equivalents
  ## y is ordered such that the atom-equivalents are first, then the non-atom-equivalents
  ## pred=Dxp matrix of atom-level predictors
  ## offs_x= X model offset
  ## offs_y= Y model offset

  
  ## betas for predictor model
  for (j in 1:p){
    beta_x[j] ~ dnorm(0,1)
  }
  
  ## betas for outcome model (there are p+1 b/c it's the number of predictors in the model for X plus 1 for X itself)
  for (j in 1:(p+1)){
    beta_y[j] ~ dnorm(0,1)
  }
  
  ## spatial random effect precision parameter
  tau_y ~ dgamma(0.001, 0.001)
  tau_x ~ dgamma(0.001, 0.001)
  
  ## vector of spatial random effects for the outcome model
  phi_y[1:S_y] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_y, zero_mean = 1)
  ## vector of spatial random effects for the predictor model
  phi_x[1:S_x] ~ dcar_normal(adj_x[], weights_x[], num_x[], tau_x, zero_mean = 1)
  
  ## expected values for latent atom-level covariate (in x-order)
  for (d in 1:D){
    lambda_atom_x[d]<-exp(sum(pred[x_to_atom[d],1:p]*beta_x[1:p])+phi_x[expand_x[d]]+offs_x[x_to_atom[d]])
  }
  
  
  ## observed atom level values are poisson ##
  for (i in 1:J_x){
    x[i] ~ dpois(lambda_atom_x[i])
  }
  
  ## missing (latent) atom level values are multinomial conditional on X grid sum
  for (m in 1:(S_x-J_x)){
    prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2]]<-lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])]/sum(lambda_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2])])
    x_latent[xlatent_ind[m,1]:xlatent_ind[m,2]] ~ dmulti(size=x[J_x+m],prob=prob_vecx[xlatent_ind[m,1]:xlatent_ind[m,2]])
  }
  
  ## put the atom-level x values (observed and latent) into atom order to plug into the y model
  temp_x[1:D]<-c(x[1:J_x],x_latent[1:(D-J_x)])
  for (d in 1:D){
    x_atomord[d]<-temp_x[x_reorder[d]]
  }
  ## note that in some cases we may want to divide this x by its denominator (exp(offset_x)) for use in the y model
  
  ## expected values for latent atom-level outcomes
  for (d in 1:D){
    lambda_atom_y[d]<-exp(sum(pred[y_to_atom[d],1:p]*beta_y[1:p])+x_atomord[y_to_atom[d]]*beta_y[p+1]+phi_y[expand_y[d]]+offs_y[y_to_atom[d]])
  }
  
  ## observed atom level values are poisson
  for (i in 1:J_y){
    y[i] ~ dpois(lambda_atom_y[i])
  }
  
  ## missing (latent) atom level values are multinomial conditional on Y grid sum
  for (m in 1:(S_y-J_y)){
    prob_vecy[ylatent_ind[m,1]:ylatent_ind[m,2]]<-lambda_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])]/sum(lambda_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])])
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ dmulti(size=y[J_y+m],prob=prob_vecy[ylatent_ind[m,1]:ylatent_ind[m,2]])
  }
  
})
