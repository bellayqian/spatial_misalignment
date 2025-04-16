# updated pull_county_data.R
library(tidycensus)
library(tigris)
library(ggplot2)
library(sf)
library(tibble)
library(tidyr)
library(dplyr)
library(stringr)
library(magrittr)

# Set tigris options for caching and simplified geometries
options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

# Define the data year
DATA_YEAR <- 2022

# Define variables for poverty data (ACS 2022)
poverty_vars <- c(
  poverty_denominator = "B23024_001",  # 20-64 year old total population for poverty estimates
  poverty_count = "B23024_002"         # 20-64 year old estimate of income below poverty level
)

# Pull poverty data from ACS
message("Fetching county poverty data from ACS...")
acs_data <- get_acs(
  geography = 'county',
  year = DATA_YEAR,
  variables = poverty_vars,
  output = "wide"
)

# Process poverty data
acs_data <- acs_data %>%
  rename(
    poverty_denominator = poverty_denominatorE,
    poverty_count = poverty_countE,
    poverty_denominator_moe = poverty_denominatorM,
    poverty_count_moe = poverty_countM
  ) %>%
  mutate(
    poverty_pct = 100 * poverty_count / poverty_denominator,
    geoid = GEOID
  ) %>%
  select(geoid, NAME, poverty_count, poverty_denominator, poverty_pct)

# Specify variables for age distribution (ACS 2022)
message("Fetching county age distribution data from ACS...")
age_variables <- paste0("B01001_0", stringr::str_pad(1:49, 2, side='left', pad=0))

# Pull age distribution data
popsizes <- get_acs(
  geography = 'county',
  year = DATA_YEAR,
  variables = age_variables
)

# Get variable definitions for mapping
variables <- load_variables(
  year = DATA_YEAR,
  dataset = 'acs5'
)

# Filter and clean variable definitions
variables <- variables %>% 
  filter(name %in% age_variables) %>%
  tidyr::separate(label, into = c('estimate_label', 'total_label', 'sex_gender', 'age_group'), sep = '!!', fill = "right") %>%
  select(-estimate_label, -total_label, -concept)

# Clean up sex/gender labels
variables$sex_gender <- stringr::str_remove_all(variables$sex_gender, ":")

# Join variable definitions to the population data
popsizes <- popsizes %>% 
  left_join(variables, by = c('variable' = 'name'))

# Aggregate by sex/gender
popsizes <- popsizes %>% 
  group_by(GEOID, NAME, age_group) %>%
  summarize(
    estimate = sum(estimate),
    moe = tidycensus::moe_sum(moe, estimate),
    .groups = "drop"
  )

# Assign 10-year age bands
popsizes <- popsizes %>% 
  mutate(
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
      age_group == '85 years and over' ~ '75+',
      TRUE ~ NA_character_
    )
  )

# Remove rows with NA age groups
popsizes <- popsizes %>% 
  filter(!is.na(age_group))

# Sum up age groups within 10-year age bands
popsizes <- popsizes %>% 
  group_by(GEOID, age_group) %>%
  summarize(estimate = sum(estimate), .groups = "drop")

# Make the age groups a factor with ordered levels
popsizes$age_group <- factor(popsizes$age_group, levels = c('0-24', '25-44', '45-64', '65-74', '75+'))

# Pivot to wide format by age group
popsizes <- popsizes %>%
  pivot_wider(names_from = age_group, values_from = estimate) %>%
  rename(
    geoid = GEOID,
    age0_24 = `0-24`,
    age25_44 = `25-44`,
    age45_64 = `45-64`,
    age65_74 = `65-74`,
    age75plus = `75+`
  )

# Get county geographic data from tigris instead of USAboundaries
message("Fetching county geographical data from tigris...")
county_sf <- tigris::counties(cb = TRUE, year = DATA_YEAR) %>%
  filter(!(STATEFP %in% c("02", "15", "72")))  # Filter out AK, HI, PR

# Standardize geoid format
county_sf <- county_sf %>%
  mutate(geoid = GEOID) %>%
  select(geoid, ALAND, geometry)

# Merge all county data
message("Merging county data...")
county <- county_sf %>%
  left_join(acs_data, by = "geoid") %>%
  left_join(popsizes, by = "geoid") %>%
  rename(aland = ALAND)

# Add in 2022 COVID-19 death data 
message("Adding COVID-19 death data...")
covid <- readRDS('./data/county_deaths_imputed.rds')

# Filter the COVID data as in the original script
# Update these filters based on your actual 2022 data structure
covid <- covid %>%
  filter(age_group == 'all_ages' & period == "apr21_to_mar22")

# Join with county data
county <- county %>%
  left_join(
    covid %>% select(county_fips, deaths),
    by = c("geoid" = "county_fips")
  )

# Calculate total population for each county
county <- county %>%
  rowwise() %>%
  mutate(
    pop_total = sum(c(age0_24, age25_44, age45_64, age65_74, age75plus), na.rm = TRUE)
  ) %>%
  ungroup()

# Export data
message("Saving county data...")
save(county, file = './data/county_data.RData')

message("County data processing complete!")