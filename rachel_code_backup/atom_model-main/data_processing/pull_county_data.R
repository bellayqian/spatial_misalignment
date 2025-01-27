library(tidycensus)
library(ggplot2)
library(USAboundaries)
library(sf)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)

## grab the poverty data from 2020 ACS (this code is modified from Christian, hence the tidy syntax)
variables_dict <-
  tibble::tribble(
    ~var,          ~shortname,      ~desc,
    "B23024_001",  'poverty_denominator', "20-64 year old total population for poverty estimates",
    "B23024_002",  'poverty_count',    "20-64 year old estimate of income below poverty level",
  )


acs_data <- get_acs(
  geography = 'county',
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


## using christian's function to pull the age distribution data by county ##
## https://github.com/ctesta01/MidtermAccountability/blob/main/pkgs/MidtermAccountibilityCOVIDOutcomes/R/merge_in_county_popsizes.R
## https://github.com/ctesta01/MidtermAccountability/blob/main/pkgs/CongressionalDistricts2022/R/create_district_covid_outcomes.R


# specify variables
age_variables <- paste0("B01001_0", stringr::str_pad(1:49, 2, side='left', pad=0))

# pull population data
popsizes <- tidycensus::get_acs(
  geography = 'county',
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

## merge with the other county data ##
county<-merge(acs_data,popsizes,by='geoid')

## merge with county shapefile ##
county_shp<-us_counties()
county_shp<-subset(county_shp,!(state_abbr%in% c('AK','HI','PR')))
county_shp<-county_shp[,c('geoid','aland','geometry')]

county<-merge(county_shp,county,by='geoid')

## add in the county-level covid 19 death data (this comes from Lancet paper, before re-allocation to CDs) ##
covid<-readRDS('county_deaths_imputed.rds')

covid<-subset(covid,age_group=='all_ages' & period=="apr21_to_mar22")

county<-merge(county,covid[,c('county_fips','deaths')],by.x='geoid',by.y='county_fips')

## export data
save(county,file='county_data.RData')
