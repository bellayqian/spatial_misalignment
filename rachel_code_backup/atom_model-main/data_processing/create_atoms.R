## get county and CD intersections ##
library(sf)
library(s2)
library(USAboundaries)
library(sp)
library(tidycensus)
library(raster)
library(dplyr)
library(data.table)
library(spdep)
library(tibble)
library(tidyr)
library(stringr)
library(magrittr)

###########################################################
## CREATE ATOMS BY INTERSECTING CD AND COUNTY SHAPEFILES ##
###########################################################

## get congressional district boundaries
cd<-us_congressional()
cd<-subset(cd,!(state_abbr%in% c('AK','HI','PR')))
#cd<-st_transform(cd,  3857)
cd2<-as(cd,'Spatial')

## get county boundaries
county<-us_counties()
county<-subset(county,!(state_abbr%in% c('AK','HI','PR')))
#county<-st_transform(county,  3857)
county2<-as(county,'Spatial')

## intersect
aa<-intersect(county2,cd2)
plot(aa)

## make back into an sf object for easier processing
ints<-st_as_sf(aa)

## read in CD-level Lancet data and merge in DW-NOMINATE since its atom-level values are known ##
lancet<-readRDS('districts_covid19_icu_and_covariates.rds')
lancet_nogeom<-lancet[,c('GEOID','nominate_dim1')]
st_geometry(lancet_nogeom)<-NULL
lancet_final<-as.data.frame(lancet_nogeom)

ints<-merge(ints,lancet_final,by.x='geoid.1',by.y='GEOID')
ints<-subset(ints,select=c('geoid.1','state_abbr.1','nominate_dim1','aland.1'))
names(ints)[1:4]<-c('county_geoid','state','nominate_dim1','county_area')

save(ints,file='county_cd_intersections.RData')

rm(list=ls())

#########################################
## GET NUMBER OF ATOMS FOR EACH COUNTY ##
#########################################

## read in the county level data
load('county_data.RData')

county<-county[order(county$geoid),]

county_4merge<-county
st_geometry(county_4merge)<-NULL

## read in atom data
load('county_cd_intersections.RData')
atom<-ints
atom<-merge(atom,county_4merge,
            by.x='county_geoid',by.y='geoid')

## order the county and atom datasets so that the counties that are atoms are first ##
## to do this, first count the number of times each county's geoid appears in the atom dataset (how many atoms it's split into) ##
countyXatom<-data.frame(table(atom$county_geoid))
names(countyXatom)<-c('geoid','num_atoms')

## merge the number of atoms for each county into the county-level dataset ##
county<-merge(county,countyXatom,by='geoid')
county<-county[order(county$num_atoms,county$geoid),]

## merge the number of atoms for each county into the atom-level dataset ##
atom<-merge(atom,countyXatom,by.x='county_geoid',by.y='geoid')
atom<-atom[order(atom$num_atoms,atom$county_geoid),]
atom$out_ind<-as.numeric(factor(atom$county_geoid,levels=county$geoid))

## identify/remove island counties (due to issues with spatial adjacency matrices) and update both datasets
county_sp1<-as_Spatial(county)
nb1 <- poly2nb(county_sp1)
islands<-which(unlist(lapply(nb1, function(x) identical(x,as.integer(0))))==T)
## remove islands from county data
county<-county[-islands,]
## remove islands from atom data ##
atom<-atom[which(atom$county_geoid %in% county$geoid),]

