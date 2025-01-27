# Load required packages
library(sp)
library(sf)
library(tidyverse)
library(ggplot2)
library(spdep)
library(MASS)
library(raster)

####### Data Generation ####### 
## All beta coefficients are close to 0 with confidence intervals centered around 0
## The precision matrices show minimal correlation between covariates
## The spatial parameters (tau) show some instability in estimation
## All variables are generated independently, 
## so correlation structure between covariates is minimal ##

# Helper function for generating spatial data
gen_spat <- function(W, rho = 0.6, var_spat = 1, global_int = 3) {
  # Precision matrix for CAR model
  precision_matrix <- diag(colSums(W)) - rho * W  
  
  # Generate spatial effects
  spatial_effects <- mvrnorm(1, mu = rep(0, nrow(W)), Sigma = solve(precision_matrix) * var_spat)
  
  # Add observation noise
  observed_values <- rpois(nrow(W), lambda = exp(global_int + spatial_effects))
  
  return(observed_values)
}

# Main simulation function
simulate_misaligned_data <- function(res1 = c(5, 5), res2 = c(10, 10), seed = 1, n_covariates_x = 3, n_covariates_y = 4) {
  set.seed(seed)
  
  # Create the first spatial grid (coarser grid)
  grid1 <- GridTopology(cellcentre.offset = c(0.5, 0.5), 
                        cellsize = c(2, 2),
                        cells.dim = res1)
  sp_grid1 <- SpatialGrid(grid1)
  
  # Create the second spatial grid (finer grid)
  grid2 <- GridTopology(cellcentre.offset = c(0, 0), 
                        cellsize = c(1, 1),
                        cells.dim = res2)
  sp_grid2 <- SpatialGrid(grid2)
  
  # Convert grids to sf objects
  sp_grid2_poly <- st_as_sf(as(sp_grid2, 'SpatialPolygons'))
  sp_grid1_poly <- st_as_sf(as(sp_grid1, 'SpatialPolygons'))
  
  # Create non-nested misalignment
  pair_list <- list(c(2,3), c(4,5), c(6,7), c(8,9))
  union_ids <- matrix(0, res2[1], res2[2])
  id_num <- 1
  for (i in seq(1, res2[1], 2)) {
    union_ids[i, pair_list[[sample(1:4, size = 1)]]] <- id_num
    id_num <- id_num + 1
  }
  
  sp_grid2_poly$union_ids <- c(t(union_ids))
  
  # Union cells
  store_xunions <- NULL
  for (i in 1:max(sp_grid2_poly$union_ids)) {
    temp <- sp_grid2_poly[sp_grid2_poly$union_ids == i, ] %>% 
      st_union() %>% 
      st_sf()
    store_xunions <- rbind(store_xunions, temp)
  }
  
  # Merge remaining cells
  grid2_nounion <- subset(sp_grid2_poly, union_ids == 0)
  store_iunions <- NULL
  for (i in 1:nrow(sp_grid1_poly)) {
    temp <- grid2_nounion[c(st_contains(sp_grid1_poly[i,], grid2_nounion, sparse = F)) == T,] %>% 
      st_union() %>% 
      st_sf()
    store_iunions <- rbind(store_iunions, temp)
  }
  
  grid2_final <- rbind(store_iunions, store_xunions)
  
  # Rename grids
  gridx <- grid2_final
  gridy <- sp_grid1_poly
  
  # 2-sim_misaligned_data.R
  gridy$ID <- 1:nrow(gridy)
  gridx$ID <- 1:nrow(gridx)
  
  # Generate Y
  neighbors_y <- poly2nb(as(gridy, "Spatial"), queen = TRUE)
  W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)
  gridy$y <- gen_spat(W = W_y)
  
  for (i in 1:n_covariates_y) {
    gridy[[paste0("covariate_y_", i)]] <- gen_spat(W = W_y, rho = runif(1, 0.3, 0.8), var_spat = runif(1, 0.5, 1.5))
  }
  
  # Generate X
  neighbors_x <- poly2nb(as(gridx, "Spatial"), queen = TRUE)
  W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)
  gridx$x <- gen_spat(W = W_x)
  
  for (i in 1:n_covariates_x) {
    gridx[[paste0("covariate_x_", i)]] <- gen_spat(W = W_x, rho = runif(1, 0.3, 0.8), var_spat = runif(1, 0.5, 1.5))
  }
  
  # Generate atom-level predictor
  atoms <- st_as_sf(raster::intersect(as(gridy, 'Spatial'), as(gridx, 'Spatial')))
  names(atoms)[which(names(atoms) == 'ID')] <- c("ID_y")
  names(atoms)[which(names(atoms) == 'ID.1')] <- c("ID_x")
  atoms$ID_atomorder <- 1:nrow(atoms)
  
  # Generate a predictor over the atoms
  # Create atom neighbor list using Queen's case
  neighbors_pred <- poly2nb(as(atoms, "Spatial"), queen = TRUE)
  # Convert neighbor list to adjacency matrix
  W_pred <- nb2mat(neighbors_pred, style = "B", zero.policy = TRUE)
  pred <- gen_spat(W = W_pred)

  # Return simulated data
  return(list(gridy = gridy, gridx = gridx, atoms = atoms, pred = pred))
}

