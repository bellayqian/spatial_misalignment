library(sp)
library(spdep)
library(sf)
library(MASS)
library(raster)
library(tidyverse)
library(nimble)

## read in the misaligned grids
load('sim_grids.RData')

gridy$ID<-1:nrow(gridy)
gridx$ID<-1:nrow(gridx)

## for now i'll just simulate y (on y grid), x (on x grid), and one predictor (on atoms) randomly from spatial processes
## later will need to use more sophisticated DGP

## helper function for generating spatial data
gen_spat<-function(W,rho=.6,var_spat=1,global_int=3){
  # Precision matrix for CAR model
  precision_matrix <- diag(colSums(W)) - rho * W  
  
  # Generate spatial effects
  spatial_effects <- mvrnorm(1, mu = rep(0, nrow(W)), Sigma = solve(precision_matrix) * var_spat)
  
  # Add observation noise
  observed_values <- rpois(nrow(W), lambda=exp(global_int+spatial_effects))
  
  return(observed_values)
}

################
## generate Y ##
################

S_y<-nrow(gridy)

# Create Y grid neighbor list using Queen's case
neighbors_y <- poly2nb(as(gridy, "Spatial"), queen = TRUE)
# Convert neighbor list to adjacency matrix
W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)

## generate Y
gridy$y<-gen_spat(W=W_y)

################
## generate X ##
################

S_x<-nrow(gridx)

# Create Y grid neighbor list using Queen's case
neighbors_x <- poly2nb(as(gridx, "Spatial"), queen = TRUE)
# Convert neighbor list to adjacency matrix
W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)

## generate Y
gridx$x<-gen_spat(W=W_x)

###################################
## generate atom-level predictor ##
###################################

## get atoms
atoms<-st_as_sf(raster::intersect(as(gridy,'Spatial'),as(gridx,'Spatial')))
plot(st_geometry(atoms))
names(atoms)[which(names(atoms)=='ID')]<-c("ID_y")
names(atoms)[which(names(atoms)=='ID.1')]<-c("ID_x")
atoms$ID_atomorder<-1:nrow(atoms)

## now generate a predictor over the atoms
# Create atom neighbor list using Queen's case
neighbors_pred <- poly2nb(as(atoms, "Spatial"), queen = TRUE)
# Convert neighbor list to adjacency matrix
W_pred <- nb2mat(neighbors_pred, style = "B", zero.policy = TRUE)

## generate pred
pred<-gen_spat(W=W_pred)

## save simulated data
# save(gridy,gridx,atoms,pred,file='sim_data.RData')