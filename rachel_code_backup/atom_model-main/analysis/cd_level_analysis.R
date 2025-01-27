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

####################################
## PULL CD %POVERTY FROM 2020 ACS ##
####################################
variables_dict <-
  tibble::tribble(
    ~var,          ~shortname,      ~desc,
    "B23024_001",  'poverty_denominator', "20-64 year old total population for poverty estimates",
    "B23024_002",  'poverty_count',    "20-64 year old estimate of income below poverty level",
  )


acs_data <- get_acs(
  geography = "congressional district",
  #geometry = TRUE,
  year = 2020,
  variables = variables_dict$var)

# pivot to a wide format for renaming, dropping the margin of error data
acs_data <- acs_data %>% select(-moe) %>%
  pivot_wider(names_from = variable, values_from = estimate)

rename_vars <- setNames(variables_dict$var, variables_dict$shortname)
acs_data <- acs_data %>% rename(!!rename_vars)

# calculate percent poverty
acs_data <- acs_data %>%
  mutate(poverty_pct = 100*poverty_count  / poverty_denominator) %>% 
  select(geoid=GEOID,poverty_count,poverty_denominator,poverty_pct)

############################################
## PULL CD AGE DISTRIBUTION FROM 2020 ACS ##
############################################

# age stratified variable names
age_variables <- paste0("B01001_0", stringr::str_pad(1:49, 2, side='left', pad=0))

# pull stratified population count data
popsizes <- tidycensus::get_acs(
  geography = "congressional district",
  year = 2020,
  state = state.abb,
  variables = age_variables
)

# separate out total population estimate
total_population <- popsizes %>% filter(variable == 'B01001_001')
total_population %<>% rename(total_population = estimate)
total_population %<>% select(-c(variable, moe))

# clean variables
variables <- tidycensus::load_variables(
  year = 2020,
  dataset = 'acs5')

# get specific variable definitions
variables %<>% filter(name %in% age_variables)

# parse estimate labels
variables %<>% tidyr::separate(label, into = c('estimate_label', 'total_label', 'sex_gender', 'age_group'), sep = '!!')
variables %<>% select(-estimate_label, -total_label, -concept)
variables$sex_gender %<>% stringr::str_remove_all(":")

# join in the variable labels
popsizes %<>% left_join(variables, by = c('variable' = 'name'))

# aggregate sex/gender
popsizes %<>% group_by(GEOID, NAME, age_group) %>%
  summarize(
    estimate = sum(estimate),
    moe = tidycensus::moe_sum(estimate, moe))

# assign 10-year age bands
popsizes %<>% mutate(
  age_group = case_when(
    age_group == 'Under 5 years' ~ '0-24',
    age_group %in% c('5 to 9 years', '10 to 14 years') ~ '0-24',
    age_group %in% c('15 to 17 years', '18 and 19 years', '20 years', '21 years', '22 to 24 years') ~ '0-24',
    age_group %in% c('25 to 29 years', '30 to 34 years') ~ '25-44',
    age_group %in% c('35 to 39 years', '40 to 44 years') ~ '25-44',
    age_group %in% c('45 to 49 years', '50 to 54 years') ~ '45-64',
    age_group %in% c('55 to 59 years', '60 and 61 years', '62 to 64 years') ~ '45-64',
    age_group %in% c('65 and 66 years', '67 to 69 years', '70 to 74 years') ~ '65-74',
    age_group %in% c('75 to 79 years', '80 to 84 years') ~ '75+',
    age_group == '85 years and over' ~ '75+'
  ))

# sum up age groups within 10-year age bands
popsizes %<>% group_by(GEOID, age_group) %>%
  summarize(estimate = sum(estimate))

# make the age groups a factor with ordered levels
popsizes$age_group %<>% factor(levels = c('0-24', '25-44', '45-64', '65-74', '75+'))

## remove NAs (these are total population*2)
popsizes<-subset(popsizes,!is.na(age_group))

## put in wide format by age group
popsizes<-popsizes%>%
  pivot_wider(names_from = age_group, values_from = estimate)%>%
  rename(geoid=GEOID,age0_24=`0-24`,age25_44=`25-44`,age45_64=`45-64`,age65_74=`65-74`,age75plus=`75+`)

## merge popsize and poverty data
acs_cd<-merge(acs_data,popsizes,by='geoid')

############################################################################
## PULL CD TOTAL POP AND RACIAL COMPOSITION FROM 2020 CENSUS REDIST FILES ##
############################################################################

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

dc2020_cd<- get_decennial(geography = "congressional district",
                          variables = variables_dict$var,
                          year = 2020,sumfile = "pl")

# pivot to a wide format for renaming
dc2020_cd <- dc2020_cd %>% 
  pivot_wider(names_from = variable, values_from = value)

rename_vars <- setNames(variables_dict$var, variables_dict$shortname)
dc2020_cd <- dc2020_cd %>% rename(!!rename_vars)

# process racial composition data
dc2020_cd <- dc2020_cd %>%
  mutate(
    pct_white = 100*race_white_alone / pop_total,
    pct_black = 100*race_black_alone / pop_total,
    pct_aian = 100*race_am_indian_ak_native_alone / pop_total,
    pct_asian = 100*race_asian_alone / pop_total,
    pct_hispanic = 100*hispanic_count / hispanic_denom) %>%
  select(geoid=GEOID,pop_total,pct_white,pct_black,pct_aian,pct_asian,pct_hispanic)

dc2020_cd<-data.frame(dc2020_cd)


## merge ACS and census data ##
cd<-merge(acs_cd,dc2020_cd,by='geoid',all.x=T)

#################################################################
## COMPUTE CD EXPECTED COVID-19 DEATHS FOR AGE STANDARDIZATION ##
#################################################################

## pull christian's function for getting pop standards ##
source('../data_processing/pull_pop_std_rates.R')
std<-load_std_population_sizes()

## now use these to get age-standardized denominators ##
cd<-cd%>%mutate(
  expected_deaths=age0_24*std$rate[1]+age25_44*std$rate[2]+age45_64*std$rate[3]+age65_74*std$rate[4]+age75plus*std$rate[5]
)

#######################################################################################
## MERGE IN DW-NOM AND CD REAPPORTIONED COVID-19 MORTALITY DATA FROM LANCET RA PAPER ##
#######################################################################################

## read in CD level file and keep only the covid outcomes dw nominate and geoids ##
lancet<-readRDS('../data_processing/districts_covid19_icu_and_covariates.rds')
lancet<-lancet[,c('GEOID','ALAND','deaths_crude','nominate_dim1')]
names(lancet)[1]<-'geoid'

st_geometry(lancet)<-NULL

cd<-merge(cd,lancet,by='geoid')

######################################################
## FINAL DATA PROCESSING AND FIT POISSON REGRESSION ##
######################################################

## compute pop density ##
cd$popdensity<-100000*cd$pop_total/cd$ALAND

## center and scale all predictors besides poverty ##
rescaled<-sapply(cd[,c('nominate_dim1','popdensity','pct_black','pct_aian','pct_asian','pct_hispanic')],scale)
rescaled<-data.frame(rescaled)
names(rescaled)<-paste0(c('nominate_dim1','popdensity','pct_black','pct_aian','pct_asian','pct_hispanic'),'_rs')

cd<-cbind(cd,rescaled)

## fit poisson regression to these data ##
fit_cd<-glm(deaths_crude~nominate_dim1_rs+popdensity_rs+
                         poverty_pct+pct_black_rs+pct_aian_rs+
                         pct_asian_rs+pct_hispanic_rs+offset(log(expected_deaths)),
            data=cd,
            family = poisson(link = 'log'))

summary(fit_cd)
confint(fit_cd)

save(fit_cd,cd,file='results_cd_level.RData')
