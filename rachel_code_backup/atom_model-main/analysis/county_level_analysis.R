library(tidycensus)
library(ggplot2)
library(USAboundaries)
library(sf)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)
library(sp)
library(spdep)

## pull in previously created county data with poverty and age distn from ACS ##
load('../data_processing/county_data.RData')
names(county)[grep('deaths',names(county))]<-'covid_deaths'

## identify/remove islands and update dataset (for consistency with use of counties in ABRM)
county_sp1<-as_Spatial(county)
nb1 <- poly2nb(county_sp1)
islands<-which(unlist(lapply(nb1, function(x) identical(x,as.integer(0))))==T)
county<-county[-islands,]

## data pulled from 2020 census redistricting files ##
variables_dict <-
  tibble::tribble(
    ~var,          ~shortname,      ~desc,
    "P1_001N",     'pop_total',                      'total population estimates (race composition table)',
    "P1_003N",     'race_white_alone',               "white alone population estimate",
    "P1_004N",     'race_black_alone',               "black alone population estimate",
    "P1_005N",     'race_am_indian_ak_native_alone', "american indian alone population estimate",
    "P1_006N",     'race_asian_alone',               "asian alone population estimate",
    "P1_007N",     'race_native_hawaiian_opi_alone', "native hawaiian or other pacific islander alone",
    "P2_001N",     'hispanic_denom',                 "denominator population for hispanic estimate",
    "P2_002N",     'hispanic_count',                 "population estimate for hispanic",
  )

## most counties are atoms, so pull the county-level dataset for the whole US ##
dc2020_county<- get_decennial(geography = "county",
                              variables = variables_dict$var,year = 2020,sumfile = "pl")

# pivot to a wide format for renaming
dc2020_county <- dc2020_county %>% 
  pivot_wider(names_from = variable, values_from = value)

rename_vars <- setNames(variables_dict$var, variables_dict$shortname)
dc2020_county <- dc2020_county %>% rename(!!rename_vars)

# process racial composition data
dc2020_county <- dc2020_county %>%
  mutate(
    pct_white = 100*race_white_alone / pop_total,
    pct_black = 100*race_black_alone / pop_total,
    pct_aian = 100*race_am_indian_ak_native_alone / pop_total,
    pct_asian = 100*race_asian_alone / pop_total,
    pct_hispanic = 100*hispanic_count / hispanic_denom) %>%
  select(geoid=GEOID,pop_total,pct_white,pct_black,pct_aian,pct_asian,pct_hispanic)

dc2020_county<-data.frame(dc2020_county)


## merge acs and residtricting data ##
county<-merge(county,dc2020_county,by='geoid',all.x=T)

## compute pop density ##
county$popdensity<-100000*county$pop_total/county$aland

## compute expected counts ##
## christian's function for getting pop standards ##
source('../data_processing/pull_pop_std_rates.R')
std<-load_std_population_sizes()

## now use these to get age-standardized denominators ##
county<-county%>%mutate(
  expected_deaths=age0_24*std$rate[1]+age25_44*std$rate[2]+age45_64*std$rate[3]+age65_74*std$rate[4]+age75plus*std$rate[5]
)

## read in CD level file and keep only dw nominate and geoids ##
lancet<-readRDS('../data_processing/districts_covid19_icu_and_covariates.rds')
cd<-lancet[,c('GEOID','state_abb','district_name','nominate_dim1')]
st_geometry(cd)<-NULL

cdshp<-us_congressional()
cdshp<-subset(cdshp,!(state_abbr%in% c('AK','HI','PR')))
cdshp<-st_transform(cdshp, 3857)

county<-st_transform(county, 3857)

## assign each county to the cd that its centroid lies in and use the DW-nominate from that cd ##
county_cents<-st_centroid(county)
county_in_cd <- st_join(county_cents, cdshp, join = st_within)
st_geometry(county_in_cd)<-NULL
county_cd<-county_in_cd[,c('geoid.x','geoid.y')]
names(county_cd)<-c('geoid','cd_id')

## merge in DW nominate ##
county_cd<-merge(county_cd,cd[,c('GEOID','nominate_dim1')],by.x='cd_id',by.y='GEOID')

## merge CD IDs into county shapefile
county<-merge(county,county_cd,by='geoid')

## center and scale all predictors besides poverty ##
temp<-county[,c('nominate_dim1','popdensity','pct_black','pct_aian','pct_asian','pct_hispanic')]
st_geometry(temp)<-NULL
rescaled<-sapply(temp,scale)
rescaled<-data.frame(rescaled)
#names(rescaled)<-paste0(c('nominate_dim1','popdensity','pct_black','pct_aian','pct_asian','pct_hispanic'),'_rs')

#county<-cbind(county,rescaled)

county_fit<-glm(county$covid_deaths~county$poverty_pct+rescaled$nominate_dim1+ rescaled$popdensity+rescaled$pct_black+
                  rescaled$pct_aian+rescaled$pct_asian+rescaled$pct_hispanic+offset(log(county$expected_deaths)),
                family=poisson(link = 'log'))


summary(county_fit)
confint(county_fit)

save(county_fit,county,file='results_county_level.RData')
