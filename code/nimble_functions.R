# Existing functions
nimbleCovarIndices <- nimbleFunction(
  run = function(j = integer(0), p = integer(0)) {
    returnType(integer(1))
    result <- numeric(p-1)
    current <- 1
    for(i in 1:p) {
      if(i != j) {
        result[current] <- i
        current <- current + 1
      }
    }
    return(result)
  }
)

calculateBetaSumX <- nimbleFunction(
  run = function(temp_x = double(2),    # D x p_x matrix
                 beta_x = double(2),     # p_x x (p_x-1) matrix
                 j = integer(0),         # current covariate index
                 p_x = integer(0)) {     # number of x covariates
    returnType(double(1))
    indices <- nimbleCovarIndices(j, p_x)
    beta_subset <- beta_x[j, 1:(p_x-1)]
    result <- numeric(nrow(temp_x))
    for(d in 1:nrow(temp_x)) {
      sum <- 0
      for(k in 1:(p_x-1)) {
        sum <- sum + temp_x[d, indices[k]] * beta_subset[k]
      }
      result[d] <- sum
    }
    return(result)
  }
)

calculateBetaSumY <- nimbleFunction(
  run = function(x_atomord = double(2),  # D x p_x matrix
                 temp_yx = double(2),     # D x p_y matrix
                 beta_yx = double(2),     # p_y x (p_x+p_y-1) matrix
                 j = integer(0),          # current covariate index
                 p_x = integer(0),        # number of x covariates
                 p_y = integer(0)) {      # number of y covariates
    returnType(double(1))
    result <- numeric(nrow(x_atomord))
    for(d in 1:nrow(x_atomord)) {
      sum_x <- 0
      for(k in 1:p_x) {
        sum_x <- sum_x + x_atomord[d,k] * beta_yx[j,k]
      }
      sum_y <- 0
      for(k in 1:(p_y-1)) {
        if(k < j) {
          sum_y <- sum_y + temp_yx[d,k] * beta_yx[j,p_x+k]
        } else if(k >= j) {
          sum_y <- sum_y + temp_yx[d,k+1] * beta_yx[j,p_x+k]
        }
      }
      result[d] <- sum_x + sum_y
    }
    return(result)
  }
)

# New functions
createCovariateSummer <- nimbleFunction(
  run = function(temp_mat = double(2),     # Input matrix (D x p)
                 beta_mat = double(2),      # Beta matrix
                 j = integer(0),            # Current index
                 p = integer(0),            # Number of covariates
                 type = character(0)) {     # Type of summation ("x", "y", or "yx")
    returnType(double(1))
    result <- numeric(nrow(temp_mat))
    
    if(type == "x") {
      indices <- nimbleCovarIndices(j, p)
      for(d in 1:nrow(temp_mat)) {
        sum <- 0
        for(k in 1:(p-1)) {
          sum <- sum + temp_mat[d, indices[k]] * beta_mat[j, k]
        }
        result[d] <- sum
      }
    } else if(type == "y") {
      for(d in 1:nrow(temp_mat)) {
        sum <- 0
        for(k in 1:p) {
          sum <- sum + temp_mat[d,k] * beta_mat[k]
        }
        result[d] <- sum
      }
    } else if(type == "yx") {
      p_x <- as.integer(ncol(temp_mat)/2)
      p_y <- p - p_x
      for(d in 1:nrow(temp_mat)) {
        sum_x <- 0
        for(k in 1:p_x) {
          sum_x <- sum_x + temp_mat[d,k] * beta_mat[j,k]
        }
        sum_y <- 0
        for(k in 1:(p_y-1)) {
          if(k < j) {
            sum_y <- sum_y + temp_mat[d,p_x+k] * beta_mat[j,p_x+k]
          } else if(k >= j) {
            sum_y <- sum_y + temp_mat[d,p_x+k+1] * beta_mat[j,p_x+k]
          }
        }
        result[d] <- sum_x + sum_y
      }
    }
    return(result)
  }
)

calculateMultinomProbs <- nimbleFunction(
  run = function(lambda = double(1),        # Vector of lambda values
                 latent_ind = double(2),    # Matrix of latent indices
                 m = integer(0)) {          # Current group index
    returnType(double(1))
    start_idx <- latent_ind[m,1]
    end_idx <- latent_ind[m,2]
    n <- end_idx - start_idx + 1
    
    probs <- numeric(n)
    sum_lambda <- 0
    
    # Calculate base probabilities and sum
    for(i in 1:n) {
      probs[i] <- lambda[start_idx + i - 1]
      sum_lambda <- sum_lambda + probs[i]
    }
    
    # Normalize probabilities
    for(i in 1:n) {
      probs[i] <- probs[i]/sum_lambda
    }
    
    return(probs)
  }
)

createTempMatrix <- nimbleFunction(
  run = function(obs_data = double(2),     # Observed data matrix
                 latent_data = double(2),   # Latent data matrix
                 J = integer(0)) {          # Number of atoms
    returnType(double(2))
    n_rows <- nrow(obs_data) + nrow(latent_data)
    n_cols <- ncol(obs_data)
    result <- matrix(0, nrow = n_rows, ncol = n_cols)
    
    # Fill observed data
    for(i in 1:J) {
      for(j in 1:n_cols) {
        result[i,j] <- obs_data[i,j]
      }
    }
    
    # Fill latent data
    for(i in 1:(n_rows-J)) {
      for(j in 1:n_cols) {
        result[J+i,j] <- latent_data[i,j]
      }
    }
    
    return(result)
  }
)