
county_cd_abrm <- nimbleCode({
  ## model coefficients for the outcome
  intercept_y ~ dnorm(0,50)
  beta_dwnom_y ~ dnorm(0, 1)
  beta_poverty ~ dnorm(0, 1)  
  beta_popdensity_y ~ dnorm(0, 1)  
  beta_black_y ~ dnorm(0, 1)  
  beta_asian_y ~ dnorm(0, 1)  
  beta_aian_y ~ dnorm(0, 1)  
  beta_hispanic_y ~ dnorm(0, 1)  
  
  ## model coefficients for the predictor
  intercept_x ~ dnorm(0,50)
  beta_dwnom_x ~ dnorm(0, 1)
  beta_popdensity_x ~ dnorm(0, 1)  
  beta_black_x ~ dnorm(0, 1)  
  beta_asian_x ~ dnorm(0, 1)  
  beta_aian_x ~ dnorm(0, 1)  
  beta_hispanic_x ~ dnorm(0, 1)  
  
  ## spatial random effect precision parameter
  tau_y ~ dgamma(0.001, 0.001)
  tau_x ~ dgamma(0.001, 0.001)
  
  ## vector of spatial random effects for the outcome model
  phi_y[1:B] ~ dcar_normal(adj[], weights[], num[], tau_y, zero_mean = 1)
  ## vector of spatial random effects for the predictor model
  phi_x[1:B] ~ dcar_normal(adj[], weights[], num[], tau_x, zero_mean = 1)
  
  ## expected values for latent atom-level covariate
  for (d in 1:D){
    lambda_atom_x[d]<-exp(intercept_x + 
                            beta_dwnom_x*dwnom[d]+
                            beta_popdensity_x * popdensity[d] +
                            beta_black_x * pct_black[d] + 
                            beta_asian_x * pct_asian[d]+
                            beta_aian_x * pct_aian[d]+
                            beta_hispanic_x * pct_hispanic[d]+
                            phi_x[out_ind[d]]+ 
                            offs_x[d])
  }
  
  ## observed atom level values are poisson ##
  for (j in 1:J){
    x[j] ~ dpois(lambda_atom_x[j])
  }
  
  ## latent atom-level values are multinomial conditional on the county level sum
  ## for each county that has multiple atoms ##
  prob_vec[1:K]<-createProbVec(vec=lambda_atom_x[1:D], a=J, M=M, mat = xlatent_ind[,], K=K)
  for (m in 1:M){
    x_latent[xlatent_ind[m,1]:xlatent_ind[m,2]] ~ dmulti(size=x[J+m],prob=prob_vec[xlatent_ind[m,1]:xlatent_ind[m,2]])
  }
  
  ## expected values for latent atom-level outcomes
  for (j in 1:J){
    lambda_atom_y[j]<-exp(intercept_y + 
                            beta_dwnom_y * dwnom[j] + 
                            beta_poverty*(100*x[j]/exp(offs_x[j])) + 
                            beta_popdensity_y * popdensity[j] +
                            beta_black_y * pct_black[j] + 
                            beta_asian_y * pct_asian[j]+
                            beta_aian_y * pct_aian[j]+
                            beta_hispanic_y * pct_hispanic[j]+
                            phi_y[j] + 
                            offs_y[j])
  }
  
  for (k in 1:K){
    lambda_atom_y[J+k]<-exp(intercept_y + 
                              beta_dwnom_y * dwnom[J+k] + 
                              beta_poverty*(100*x_latent[k]/exp(offs_x[J+k])) + 
                              beta_popdensity_y * popdensity[J+k] +
                              beta_black_y * pct_black[J+k] + 
                              beta_asian_y * pct_asian[J+k]+
                              beta_aian_y * pct_aian[J+k]+
                              beta_hispanic_y * pct_hispanic[J+k]+
                              phi_y[atom_ord[k]] + 
                              offs_y[J+k])
  }
  
  for (j in 1:J){
    lambda[j]<-lambda_atom_y[j]
    y[j] ~ dpois(lambda[j])
  }
  
  for (m in 1:M){
    lambda[J+m]<-sum(lambda_atom_y[(J+xlatent_ind[m,1]):(J+xlatent_ind[m,2])])
    y[J+m] ~ dpois(lambda[J+m])
  }
  
  ## save predicted values of the outcome for atoms in LA county ##
  la_preds[1:18]<-lambda_atom_y[3711:3728]
  
})
