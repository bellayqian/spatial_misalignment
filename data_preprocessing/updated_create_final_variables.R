# Fixed create_final_variables.R
# Final data processing to create analysis-ready dataset
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)

message("Starting final variable creation...")

# Load demographic atom data (contains pop_total, racial percentages)
load('./data/block_data.RData')
atom_demographic <- atom  # Save demographic data
message("Loaded demographic atom data with", nrow(atom_demographic), "atoms.")

# Load geographic atom data (contains boundaries, county associations)
load("./data/atom_county_data.RData")
atom_geo <- atom  # Save geographic data
county_orig <- county  # Save original county data
message("Loaded geographic atom data with", nrow(atom_geo), "atoms.")

# Source population standardization functions
source('./R/pull_pop_std_rates.R')

# Spatial join using a simpler approach more similar to the original
message("Performing spatial join using simplified approach...")

# Prepare the atom_demographic (no need to match IDs - use spatial relationship)
# Transform to same projection
atom_demographic <- st_transform(atom_demographic, st_crs(atom_geo))

# We'll be selective about which columns we keep from each dataset
geo_cols <- c("atom_id", "county_geoid", "cd_geoid", "state", "county_name", "cd_name", "county_area", "nominate_dim1")
demo_cols <- c("pop_total", "pct_white", "pct_black", "pct_aian", "pct_asian", "pct_hispanic")

# Do the spatial join
atom <- st_join(
  atom_geo %>% select(all_of(geo_cols)),  
  atom_demographic %>% select(all_of(demo_cols))
)

##### Problem: When using st_join(), it joins rows when geometries intersect. 
# If a polygon in atom_geo intersects with multiple polygons in atom_demographic, 
# I'll get one row for each intersection.
##### 

# Calculate county totals using base R approach (similar to PI's code)
message("Computing county population totals...")
atom_noshp <- atom
st_geometry(atom_noshp) <- NULL

# Use aggregate() to sum populations by county (like the original script)
dc_county_pop <- aggregate(pop_total ~ county_geoid, data = atom_noshp, FUN = sum, na.rm = TRUE)
names(dc_county_pop) <- c('geoid', 'pop_total')

# Merge with county data (use merge like the original)
county <- merge(county_orig, dc_county_pop, by = 'geoid')

# Calculate adjusted poverty count (same as original)
message("Computing adjusted poverty counts...")
county$poverty_count_adjusted <- round((county$poverty_pct / 100) * county$pop_total)

message("Creating age-standardized denominators at atom level...")
# Create age-standardized denominators at atom level using county age distributions

# Calculate age distribution proportions (same as original)
atom_new <- atom
atom_new <- atom_new %>% 
  mutate(
    age_denom = age0_24 + age25_44 + age45_64 + age65_74 + age75plus
  ) %>% 
  mutate(
    # Add safety check for division by zero
    age0_24_p = ifelse(age_denom > 0, age0_24 / age_denom, 0),
    age25_44_p = ifelse(age_denom > 0, age25_44 / age_denom, 0),
    age45_64_p = ifelse(age_denom > 0, age45_64 / age_denom, 0),
    age65_74_p = ifelse(age_denom > 0, age65_74 / age_denom, 0),
    age75plus_p = ifelse(age_denom > 0, age75plus / age_denom, 0)
  ) %>% 
  select(-c(age0_24, age25_44, age45_64, age65_74, age75plus, age_denom))

# Calculate age-specific counts for each atom based on pop_total from block data
atom_new <- atom_new %>%
  mutate(
    age0_24 = pop_total * age0_24_p,
    age25_44 = pop_total * age25_44_p,
    age45_64 = pop_total * age45_64_p,
    age65_74 = pop_total * age65_74_p,
    age75plus = pop_total * age75plus_p
  )

# Get standard population rates for age standardization
message("Loading standard population rates...")
std <- load_std_population_sizes()

# Calculate age-standardized denominators
atom_new <- atom_new %>%
  mutate(
    expected_count = age0_24 * std$rate[1] +
      age25_44 * std$rate[2] +
      age45_64 * std$rate[3] +
      age65_74 * std$rate[4] +
      age75plus * std$rate[5]
  )

# Replace atom data with new calculations
atom <- atom_new

# Debug: Check specific counties before filtering
message("Before filtering, checking specific counties...")
jackson_al <- atom %>% filter(county_geoid == "01071")
orange_nc <- atom %>% filter(county_geoid == "37135")  
message(paste("Jackson County, AL has", nrow(jackson_al), "atoms before filtering"))
message(paste("Orange County, NC has", nrow(orange_nc), "atoms before filtering"))

# Special handling for counties with single atoms
# Get the count of atoms per county
atoms_per_county <- as.data.frame(table(atom$county_geoid))
names(atoms_per_county) <- c("county_geoid", "atom_count")

# Identify counties with only one atom
single_atom_counties <- atoms_per_county$county_geoid[atoms_per_county$atom_count == 1]
message(paste("Found", length(single_atom_counties), "counties with only one atom"))

# Debug: Check if our problem counties have single atoms
if("01071" %in% single_atom_counties) {
  message("Jackson County, AL has only one atom - will be preserved")
}
if("37135" %in% single_atom_counties) {
  message("Orange County, NC has only one atom - will be preserved")
}

# Remove atoms with zero population, BUT preserve one atom per county
message("Removing atoms with zero population while preserving counties with single atoms...")
atom_preserved <- atom %>%
  filter(county_geoid %in% single_atom_counties)  # Keep ALL atoms from counties with only one atom

atom_filtered <- atom %>%
  filter(!(county_geoid %in% single_atom_counties)) %>%  # Process only multi-atom counties
  filter(pop_total > 0)  # Remove zero-pop atoms from multi-atom counties

# Combine the preserved and filtered datasets
atom <- bind_rows(atom_preserved, atom_filtered)

message("After filtering, checking specific counties...")
jackson_al <- atom %>% filter(county_geoid == "01071")
orange_nc <- atom %>% filter(county_geoid == "37135")  
message(paste("Jackson County, AL has", nrow(jackson_al), "atoms after filtering"))
message(paste("Orange County, NC has", nrow(orange_nc), "atoms after filtering"))

# Recompute num_atoms since we removed some atoms
message("Recalculating atom counts...")
countyXatom <- as.data.frame(table(atom$county_geoid))
names(countyXatom) <- c('geoid', 'num_atoms')

# Update county-level atom counts (using merge like original)
county$num_atoms <- NULL
county <- merge(county, countyXatom, by = "geoid")
county <- county[order(county$num_atoms, county$geoid),]

# Update atom-level atom counts (using merge like original)
atom$num_atoms <- NULL
atom <- merge(atom, countyXatom, by.x = "county_geoid", by.y = "geoid")
atom <- atom[order(atom$num_atoms, atom$county_geoid),]

# Create out_ind index (same as original)
atom$out_ind <- as.numeric(factor(atom$county_geoid, levels = county$geoid))

message("Computing population density at atom level...")
# Compute population density at atom level
atom$atom_area <- as.numeric(st_area(atom))
atom$popdensity <- 100000 * atom$pop_total / atom$atom_area

# Add safety check for Inf values in population density
atom$popdensity[is.infinite(atom$popdensity) | is.na(atom$popdensity)] <- 0

# Round death counts for Poisson modeling
if("deaths" %in% names(county)) {
  county$deaths_round <- round(county$deaths)
}

message("Saving final misaligned data...")
# Save the final misaligned data
save(county, atom, file = './data/final_misaligned_data.RData')

message("Final variable creation complete!")