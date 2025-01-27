## master script for ABRM data pre-processing ##


###########################################
## EXTRACT AND PROCESS COUNTY-LEVEL DATA ##
###########################################

source('pull_county_data.R')

rm(list = ls())

###################################################
## CREATE ATOMS BY INTERSECTING COUNTIES AND CDS ##
###################################################

source('create_atoms.R')

##########################################################################
## EXTRACT AND PROCESS CENSUS BLOCK DATA AND AGGREGATE/MERGE WITH ATOMS ##
##########################################################################

source('pull_block_data.R')

######################################################################
## DO SOME FINAL DATA MANIPULATION STEPS AND SAVE THE ANALYTIC DATA ##
######################################################################

source('create_final_variables.R')
