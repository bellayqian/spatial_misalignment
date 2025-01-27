## CODE TO RUN ABRM ON COVID-19 AND DW-NOMINATE DATA ##

## version for fasrc batch submission ##
## this has to be run within singularity container ##

####################################################
## LOAD LIBRARIES AND FUNCTIONS AND SET CONSTANTS ##
####################################################

library(nimble)
library(sp)
library(spdep)
library(sf)

niter<-100000
nburn<-50000

## read in nimble functions ##
source('../nimble/nimble_functions.R')

## read in nimble models ##
source('../nimble/nimble_model.R')

#####################################################
## LOAD DATA AND CONSTRUCT NECESSARY NIMBLE INPUTS ##
#####################################################

## load in file with pre-processed county and atom level datasets
## both datasets must be sorted by
## (1) number of atoms
## (2) county geoid/fips

load('../data_processing/final_misaligned_data.RData')

## make the county spatial adjacency matrix ##
nb <- poly2nb(as_Spatial(county),row.names = as_Spatial(county)@data$geoid)
W=nb2mat(nb,style='B')
## apply the nimble function to put it in format nimble recognizes ##
adj_nimble<-as.carAdjacency(W)

## county indexes for sub-county atoms
atom_ord<-atom$out_ind[which(atom$num_atoms>1)]
## number of atoms for each county with multiple atoms
atom_count<-as.numeric(table(atom_ord))

## create matrix with start and end indices of atoms for each county with multiple atoms
ctr<-unique(atom_ord)
latx_indx<-matrix(0,nrow=length(atom_count),ncol=2)
for (i in 1:nrow(latx_indx)){
  latx_indx[i,]<-c(min(which(atom_ord==ctr[i])), max(which(atom_ord==ctr[i])))
}

## get reasonable init values for the latent poverty counts (to help convergence) ##
county_noshp<-county
st_geometry(county_noshp)<-NULL
initval_pov<-merge(atom[which(atom$num_atoms>1),c('county_geoid','pop_total','num_atoms')],
                   county_noshp[which(county$num_atoms>1),c('geoid','poverty_count_adjusted','pop_total')],
                   by.x='county_geoid',by.y='geoid')
initval_pov<-initval_pov[order(initval_pov$num_atoms,initval_pov$county_geoid),]
initval_pov$mult_probs<-round(initval_pov$pop_total.x/initval_pov$pop_total.y,2)

initval_pov$x_init<-NA
for (i in unique(initval_pov$county_geoid)){
  inds<-which(initval_pov$county_geoid==i)
  initval_pov$x_init[inds]<-rmultinom(n=1,size=initval_pov$poverty_count_adjusted[inds[1]],
                                      prob=initval_pov$mult_probs[inds])
}

############################################################################
## SET UP OBJECTS OF CONSTANTS, DATA, AND INIT VALUES TO PASS INTO NIMBLE ##
############################################################################

## list of constants to be passed into nimble ##
constants <- list(## number of counties ##
                  B = nrow(county), 
                  ## number of atoms ##
                  D = nrow(atom),
                  ## first J atoms in the dataset are county-equivalent atoms
                  J = max(which(atom$num_atoms==1)),
                  ## remaining K atoms are sub-county atoms
                  K = length(which(atom$num_atoms>1)),
                  ## number of counties split into multiple atoms
                  M=length(unique(atom_ord)),
                  ## unique county id for each atom
                  out_ind = atom$out_ind,
                  ## county ids for sub-county atoms only (atoms >J)
                  atom_ord = atom_ord,
                  ## number of atoms for counties with >1 atom
                  atom_count = atom_count,
                  ## matrix with start and end indices of atoms for each county with multiple atoms
                  xlatent_ind = latx_indx)

## list containing all data to be passed into nimble ##
data <- list(
  ## data at the county level
  y = county$deaths_round,
  x = county$poverty_count_adjusted,
  ## data at the atom level
  dwnom = c(scale(atom$nominate_dim1)),
  popdensity = c(scale(atom$popdensity)),
  pct_black= c(scale(atom$pct_black)),
  pct_asian= c(scale(atom$pct_asian)),
  pct_aian= c(scale(atom$pct_aian)),
  pct_hispanic=c(scale(atom$pct_hispanic)),
  offs_x = log(atom$pop_total),
  offs_y = log(atom$expected_count),
  ## objects needed for the spatial random effects ##
  num = adj_nimble$num, weights = adj_nimble$weights,adj = adj_nimble$adj)

## list of initial values of parameters to be passed into nimble ##
inits <- list(intercept_x = -2, intercept_y = 0, 
              beta_dwnom_y = 0, beta_poverty = 0,
              beta_popdensity_y =0, 
              beta_black_y =0, 
              beta_asian_y =0, 
              beta_aian_y =0, 
              beta_hispanic_y =0,
              beta_dwnom_x = 0,
              beta_popdensity_x =0, 
              beta_black_x =0, 
              beta_asian_x =0, 
              beta_aian_x =0, 
              beta_hispanic_x =0,
              tau_x = 1, tau_y = 1,
              phi_y = rep(0,nrow(county)), phi_x = rep(0,nrow(county)),
              x_latent=initval_pov$x_init)

#################################################
## CREATE AND RUN NIMBLE MODEL AND SAVE OUTPUT ##
#################################################

## create nimble model ##
Rmodel <- nimbleModel(county_cd_abrm, constants, data, inits,calculate=FALSE)

## compile nimble model ##
compmod <- compileNimble(Rmodel)

## configure mcmc ##
mcmcConf <- configureMCMC(compmod, enableWAIC = F,
                          monitors=c('intercept_y', 
                                     'beta_dwnom_y','beta_poverty',
                                     'beta_popdensity_y', 
                                     'beta_black_y', 
                                     'beta_asian_y', 
                                     'beta_aian_y', 
                                     'beta_hispanic_y',
                                     'intercept_x', 
                                     'beta_dwnom_x',
                                     'beta_popdensity_x', 
                                     'beta_black_x', 
                                     'beta_asian_x', 
                                     'beta_aian_x', 
                                     'beta_hispanic_x',
                                     'la_preds'),thin=10,
                          useConjugacy=FALSE)

## build mcmc ##
Rmcmc <- buildMCMC(mcmcConf)

## compile mcmc ##
Cmcmc <- compileNimble(Rmcmc, project = compmod)

## run mcmc ##
mcmc.out <- runMCMC(Cmcmc, niter=niter,nburnin=nburn,nchains=1,summary=T)

## save output ##
save(mcmc.out,file='results_abrm.RData')
