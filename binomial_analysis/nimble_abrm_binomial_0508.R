### Nimble Code for Binomial ABRM
## two spatial grids, X and Y, from which all variables arise
## nothing fully measured at the atom level
## all variables are binomial

# Define the dmfnchypg function
dmfnchypg <- nimbleFunction(
  run = function(x = double(1), total = double(0), odds = double(1), ni = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Round all input values to handle floating point issues
    x_round <- round(x)
    total_round <- round(total)
    
    # Calculate sum of x using explicit loop
    sumX <- 0
    for(i in 1:length(x)) {
      sumX <- sumX + x[i]
    }
    
    # Check if sum(x) = total
    if(abs(sumX - total) > 0.1) {  # Allow for rounding differences
      if(log == 1) return(-Inf) else return(0)
    }
    
    # Check bounds constraints with explicit loops
    for(i in 1:length(x)) {
      if(x[i] < 0) {
        if(log == 1) return(-Inf) else return(0)
      }
      if(x_round[i] > round(ni[i])) {  # Round population sizes too
        if(log == 1) return(-Inf) else return(0)
      }
    }
    
    # Calculate kernel of log probability 
    logProb <- 0
    
    # Sum of log-factorials of ni
    for(i in 1:length(ni)) {
      logProb <- logProb + lfactorial(round(ni[i]))
    }
    
    # Subtract sum of log-factorials of x
    for(i in 1:length(x_round)) {
      logProb <- logProb - lfactorial(x_round[i])
    }
    
    # Subtract sum of log-factorials of (ni-x)
    for(i in 1:length(x_round)) {
      logProb <- logProb - lfactorial(round(ni[i]) - x_round[i])
    }
    
    # Add sum of x*log(odds)
    for(i in 1:length(x_round)) {
      # Ensure odds are positive
      safe_odds <- max(odds[i], 1e-10)
      logProb <- logProb + x_round[i] * log(safe_odds)
    }
    
    if(log == 1) return(logProb) else return(exp(logProb))
  }
)

# Ensure the BiasedUrn package is installed
# install.packages("BiasedUrn")
# Step 1: Define the R function that will wrap the BiasedUrn call
biasedUrn_rmfnc <- function(total, odds, ni) {
  # Validation and conversion
  total <- as.integer(round(total))
  ni <- as.integer(round(ni))
  
  # Ensure odds are positive
  odds[odds <= 0] <- 1e-10
  
  # Check total vs. available population
  sum_ni <- sum(ni)
  if(total > sum_ni) {
    warning(sprintf("Total (%d) exceeds sum of population sizes (%d). Adjusting total.", 
                    total, sum_ni))
    total <- sum_ni
  }
  
  # For single-atom cases, return the only possible solution
  if(length(ni) == 1) {
    return(total)
  }
  
  # BiasedUrn sampling with error handling
  # Double-check all inputs
  if(any(is.na(ni)) || any(is.na(odds)) || is.na(total)) {
    stop("NA values in inputs")
  }
  if(any(ni < 0) || any(odds < 0) || total < 0) {
    stop("Negative values in inputs")
  }
  
  result <- BiasedUrn::rMFNCHypergeo(1, ni, total, odds)
  return(as.numeric(result))
}

# Step 2: Create a nimble-compatible function using nimbleRcall
Rmfnchypg <- nimbleRcall(
  prototype = function(total = double(0), odds = double(1), ni = double(1)) {},
  returnType = double(1),
  Rfun = "biasedUrn_rmfnc"
)

# Step 3: Create the actual nimbleFunction wrapper that will be used in the model
rmfnchypg <- nimbleFunction(
  run = function(n = integer(0), total = double(0), odds = double(1), ni = double(1)) {
    returnType(double(1))
    # Simply call the nimbleRcall-created function
    return(Rmfnchypg(total, odds, ni))
  }
)

# Register the distribution
registerDistributions(list(
  dmfnchypg = list(
    BUGSdist = "dmfnchypg(total, odds, ni)",
    discrete = TRUE,
    types = c('value = double(1)', 'total = double(0)', 'odds = double(1)', 'ni = double(1)'),
    mixedSizes = TRUE,
    pqAvail = FALSE,
    range = c(0, Inf)
  )
))

abrm_binomial <- nimbleCode({
  
  for(j in 1:p_x) {
    beta_0_x[j] ~ dnorm(0, sd = 0.1)  # Intercepts for X covariates
  }
  
  for(j in 1:p_y) {
    beta_0_yx[j] ~ dnorm(0, sd = 0.1) # Intercepts for Y covariates
  }
  
  beta_0_y ~ dnorm(0, sd = 0.1)  # Intercept for Y outcome
  
  for (j in 1:(p_x+p_y)){
    beta_y[j] ~ dnorm(0, sd = 0.1) # Effects of X' and Y covariates on Y
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
    for(j in 1:p_x) {
      phi_x[i,j] <- inprod(Achol_x[j,1:p_x], psi_x[i,1:p_x])
    }
  }
  
  for (j in 1:p_y){ 
    tau_yx[j] ~ dgamma(2, 2)
    psi_yx[1:S_y,j] ~ dcar_normal(adj_y[], weights_y[], num_y[], tau_yx[j])
  }
  
  ## Correlation in the spatial random effects across y-grid covariates
  Prec_yx[1:p_y,1:p_y] ~ dwish(R_yx[1:p_y,1:p_y], df_yx)
  Achol_yx[1:p_y,1:p_y] <- chol(inverse(Prec_yx[1:p_y,1:p_y]))
  
  for(i in 1:S_y){
    for(j in 1:p_y) {
      phi_yx[i,j] <- inprod(Achol_yx[j,1:p_y], psi_yx[i,1:p_y])
    }
  }
  
  ################################
  ## MODELING X-GRID COVARIATES ##
  ################################
  
  ## For X-grid covariates (calculate linear predictor, odds, and probabilities)
  for(j in 1:p_x) {
    for(d in 1:D) {
      # Linear predictor (η_i in prompt)
      linear_pred_x_base[d,j] <- beta_0_x[j] + phi_x[expand_x[d],j] + offs_x[d]
      
      # Odds values (π_i/(1-π_i) in prompt)
      odds_atom_x[d,j] <- exp(min(max(linear_pred_x_base[d,j], -100), 100))
      
      # Probabilities (π_i in prompt) - for binomial observations
      prob_atom_x[d,j] <- min(0.9999, max(0.0001, odds_atom_x[d,j] / (1 + odds_atom_x[d,j])))
    }
  }
  
  ## observed atom-level values (where we directly observe at atom level)
  for (j in 1:p_x){
    for (i in 1:J_x){
      # Direct binomial observation - FIXED: Use correct atom index
      x[i,j] ~ dbinom(size=n_x[i], prob=prob_atom_x[x_to_atom[i],j])
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_x){
    for (m in 1:(S_x-J_x)){
      # Use multivariate non-central hypergeometric distribution
      x_latent[xlatent_ind[m,1]:xlatent_ind[m,2],j] ~ 
        dmfnchypg(total=x[J_x+m,j], 
                  odds=odds_atom_x[(J_x+xlatent_ind[m,1]):(J_x+xlatent_ind[m,2]),j],
                  ni=n_x_latent[xlatent_ind[m,1]:xlatent_ind[m,2]])
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
  
  # Create x_atomord for use in modeling Y
  for (j in 1:p_x) {
    for(d in 1:D) {
      x_atomord[d,j] <- temp_x[x_reorder[d],j]
    }
  }
  
  ################################
  ## MODELING Y-GRID COVARIATES ##
  ################################
  
  ## For Y-grid covariates (calculate linear predictor, odds, and probabilities)
  for(j in 1:p_y) {
    for(d in 1:D) {
      # Linear predictor (η_i in prompt)
      linear_pred_yx_base[d,j] <- beta_0_yx[j] + phi_yx[expand_y[d],j] + offs_y[d]
      
      # Odds values (π_i/(1-π_i) in prompt)
      odds_atom_yx[d,j] <- exp(min(max(linear_pred_yx_base[d,j], -100), 100))
      
      # Probabilities (π_i in prompt) - for binomial observations
      prob_atom_yx[d,j] <- min(0.9999, max(0.0001, odds_atom_yx[d,j] / (1 + odds_atom_yx[d,j])))
    }
  }
  
  ## observed atom level values (first J_y rows of S_y × p_y matrix)
  for (j in 1:p_y){
    for (i in 1:J_y){ 
      # Direct binomial observation - FIXED: Use correct atom index
      yx_obs[i,j] ~ dbinom(size=n_y[i], prob=prob_atom_yx[y_to_atom[i],j])
    }
  }
  
  ## multinomial probabilities and latent values for non-atom-equivalents
  for (j in 1:p_y){
    for (m in 1:(S_y-J_y)){
      # Use multivariate non-central hypergeometric distribution
      yx_latent[ylatent_ind[m,1]:ylatent_ind[m,2],j] ~ 
        dmfnchypg(total=yx_obs[J_y+m,j],
                  odds=odds_atom_yx[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2]),j],
                  ni=n_y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]])
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
  
  # #######################
  # ## MODELING Y ITSELF ##
  # #######################
  # 
  # ## calculate final odds and probability values for Y
  # for(d in 1:D) {
  #   # Calculate all components
  #   x_effect_sum[d] <- sum(x_atomord[y_to_atom[d], 1:p_x] * beta_y[1:p_x])
  #   y_effect_sum[d] <- sum(temp_yx[d, 1:p_y] * beta_y[(p_x+1):(p_x+p_y)])
  #   
  #   # Define linear predictor in a single statement
  #   linear_pred_y_base[d] <- beta_0_y + offs_y[d] + x_effect_sum[d] + y_effect_sum[d] + phi_y[expand_y[d]]
  #   
  #   # Calculate odds
  #   odds_atom_y[d] <- exp(min(max(linear_pred_y_base[d], -100), 100))
  #   
  #   # Calculate probability
  #   prob_atom_y[d] <- min(0.9999, max(0.0001, odds_atom_y[d] / (1 + odds_atom_y[d])))
  # }
  # 
  # ## observed atom level values
  # for (i in 1:J_y){
  #   # Direct binomial observation - FIXED: Use correct atom index
  #   y_obs[i] ~ dbinom(size=n_y[i], prob=prob_atom_y[y_to_atom[i]])
  # }
  # 
  # ## multinomial probabilities and latent values for non-atom-equivalents
  # for (m in 1:(S_y-J_y)){
  #   # Use multivariate non-central hypergeometric distribution
  #   y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~ 
  #     dmfnchypg(total=y_obs[J_y+m],
  #               odds=odds_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])],
  #               ni=n_y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]])
  # }
  
  #######################
  ## MODELING Y ITSELF ##
  #######################

  ## calculate final odds and probability values for Y
  for(d in 1:D) {
    # Calculate X covariate effects using proportions (vectorized)
    atom_idx <- y_to_atom[d]
    x_effect_sum[d] <- sum((x_atomord[atom_idx, 1:p_x] / Pop_atom_true[atom_idx]) * beta_y[1:p_x])
    
    # Calculate Y covariate effects using proportions (vectorized)  
    y_effect_sum[d] <- sum((temp_yx[d, 1:p_y] / Pop_atom_true[atom_idx]) * beta_y[(p_x+1):(p_x+p_y)])
    
    # Define linear predictor (same as before)
    linear_pred_y_base[d] <- beta_0_y + offs_y[d] + x_effect_sum[d] + y_effect_sum[d] + phi_y[expand_y[d]]
    
    # Calculate odds and probability (same as before)
    odds_atom_y[d] <- exp(min(max(linear_pred_y_base[d], -100), 100))
    prob_atom_y[d] <- min(0.9999, max(0.0001, odds_atom_y[d] / (1 + odds_atom_y[d])))
  }
  
  ## observed atom level values
  for (i in 1:J_y){
    # Direct binomial observation - FIXED: Use correct atom index
    y_obs[i] ~ dbinom(size=n_y[i], prob=prob_atom_y[y_to_atom[i]])
  }

  ## multinomial probabilities and latent values for non-atom-equivalents
  for (m in 1:(S_y-J_y)){
    # Use multivariate non-central hypergeometric distribution
    y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]] ~
      dmfnchypg(total=y_obs[J_y+m],
                odds=odds_atom_y[(J_y+ylatent_ind[m,1]):(J_y+ylatent_ind[m,2])],
                ni=n_y_latent[ylatent_ind[m,1]:ylatent_ind[m,2]])
  }
})