library(sp)
library(spdep)
library(sf)
library(MASS)
library(raster)
library(tidyverse)
library(nimble)

## load in simulated data
load('sim_data.RData')

S_x<-nrow(gridx)
S_y<-nrow(gridy)

#####################
## prep for nimble ##
#####################

## get the spatial bookkeeping items ##
## x grids that aren't atoms
x_nonatoms<-atoms$ID_x[which(duplicated(atoms$ID_x)==T)]
x_nonatoms<-x_nonatoms[order(x_nonatoms)]
J_x<-nrow(gridx)-length(x_nonatoms)
## y grids that aren't atoms
y_nonatoms<-atoms$ID_y[which(duplicated(atoms$ID_y)==T)]
y_nonatoms<-y_nonatoms[order(y_nonatoms)]
J_y<-nrow(gridy)-length(y_nonatoms)

## re-order x and y so that first J_x/J_y units are atom-equivalents
gridy_yorder<-gridy[c(gridy$ID[-y_nonatoms],y_nonatoms),]
names(gridy_yorder)[which(names(gridy_yorder)=='ID')]<-'ID_y'
gridy_yorder$ID_yorder<-1:nrow(gridy_yorder)

gridx_xorder<-gridx[c(gridx$ID[-x_nonatoms],x_nonatoms),]
names(gridx_xorder)[which(names(gridx_xorder)=='ID')]<-'ID_x'
gridx_xorder$ID_xorder<-1:nrow(gridx_xorder)

## link the new IDs into the atoms dataset so we can find correspondence between new Y grid ids and atoms
IDxwalk<-merge(st_drop_geometry(atoms),st_drop_geometry(gridy_yorder),by='ID_y')
IDxwalk<-merge(IDxwalk,st_drop_geometry(gridx_xorder),by='ID_x')

## get x/y_to_atom
x_to_atom<-IDxwalk$ID_atomorder[order(IDxwalk$ID_xorder,IDxwalk$ID_atomorder)]
y_to_atom<-IDxwalk$ID_atomorder[order(IDxwalk$ID_yorder,IDxwalk$ID_atomorder)]
## get expand_x/y
expand_x<-IDxwalk$ID_xorder[order(IDxwalk$ID_xorder)]
expand_y<-IDxwalk$ID_yorder[order(IDxwalk$ID_yorder)]
## get x/y_latentind
IDxwalk<-IDxwalk[order(IDxwalk$ID_xorder,IDxwalk$ID_atomorder),]
x_nonatoms<-IDxwalk$ID_xorder[which(duplicated(IDxwalk$ID_xorder)==T)]
x_latentind<-matrix(0,length(x_nonatoms),2)
for (i in 1:length(x_nonatoms)){
  x_latentind[i,]<-c(min(which(IDxwalk$ID_xorder==x_nonatoms[i])),max(which(IDxwalk$ID_xorder==x_nonatoms[i])))-J_x
}

IDxwalk<-IDxwalk[order(IDxwalk$ID_yorder,IDxwalk$ID_atomorder),]
y_nonatoms<-IDxwalk$ID_yorder[which(duplicated(IDxwalk$ID_yorder)==T)]
y_latentind<-matrix(0,length(y_nonatoms),2)
for (i in 1:length(y_nonatoms)){
  y_latentind[i,]<-c(min(which(IDxwalk$ID_yorder==y_nonatoms[i])),max(which(IDxwalk$ID_yorder==y_nonatoms[i])))-J_y
}

## get x_reorder
IDxwalk<-IDxwalk[order(IDxwalk$ID_xorder,IDxwalk$ID_atomorder),]
IDxwalk$xro<-1:nrow(IDxwalk)
IDxwalk<-IDxwalk[order(IDxwalk$ID_atomorder),]
x_reorder<-IDxwalk$xro

## make spatial adjacency matrices for X and Y grids based on new orderings
neighbors_x <- poly2nb(as(gridx_xorder, "Spatial"), queen = TRUE)
W_x <- nb2mat(neighbors_x, style = "B", zero.policy = TRUE)

neighbors_y <- poly2nb(as(gridy_yorder, "Spatial"), queen = TRUE)
W_y <- nb2mat(neighbors_y, style = "B", zero.policy = TRUE)


## inputs to nimble
constants<-list(p=2,
                D=nrow(atoms),
                S_x= S_x,
                S_y= S_y,
                J_x= J_x,
                J_y= J_y,
                x_to_atom=x_to_atom,
                y_to_atom=y_to_atom,
                xlatent_ind=x_latentind,
                ylatent_ind=y_latentind,
                expand_x=expand_x,
                expand_y=expand_y,
                x_reorder=x_reorder,
                num_x = as.carAdjacency(W_x)$num,
                weights_x = as.carAdjacency(W_x)$weights,
                adj_x = as.carAdjacency(W_x)$adj,
                num_y = as.carAdjacency(W_y)$num,
                weights_y = as.carAdjacency(W_y)$weights,
                adj_y = as.carAdjacency(W_y)$adj)

data<-list(
  x=gridx_xorder$x,
  y=gridy_yorder$y,
  pred=cbind(1,pred),
  offs_x=rep(0,nrow(atoms)),
  offs_y=rep(0,nrow(atoms))
)

inits<-list(
  beta_x=rep(0,2),
  beta_y=rep(0,3),
  tau_x = 1, tau_y = 1,
  phi_x=rep(0,S_x),
  phi_y=rep(0,S_y)
)

######################
## run nimble model ##
######################

source('nimble_abrm.R')

## create nimble model ##
Rmodel <- nimbleModel(abrm, constants, data, inits,calculate=FALSE)

## compile nimble model ##
compmod <- compileNimble(Rmodel)

## configure mcmc ##
mcmcConf <- configureMCMC(compmod, enableWAIC = F,
                          monitors=c('beta_x',
                                     'beta_y'),thin=10,
                          useConjugacy=FALSE)

## build mcmc ##
Rmcmc <- buildMCMC(mcmcConf)

## compile mcmc ##
Cmcmc <- compileNimble(Rmcmc, project = compmod)

## run mcmc ##
mcmc.out <- runMCMC(Cmcmc, niter=10000,nburnin=5000,nchains=1,summary=T)

save(mcmc.out,file='model_output.RData')