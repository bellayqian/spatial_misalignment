# Updated create_atoms.R
# Create atoms by intersecting county and congressional district shapefiles
library(sf)
library(tigris)
library(dplyr)
library(data.table)
library(tibble)
library(tidyr)
library(spdep)
library(stringr)
library(magrittr)

# 1. Execute Pull County Data from R/pull_county_data.R 
# to get data/county_cd_intersections.RData (intermediate)
# 2. Execute Create Atoms from R/updated_create_atoms.R
# to get data/atom_county_data.RData (atoms with county associations)
# 3. Generate Block Data if we do not have that
# 4. Execute Final Variables Creation to combine all and create the final analytical dataset.

# Set tigris options
options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

# Define the data year
DATA_YEAR <- 2022

message("Starting atom creation process...")

###########################################################
## CREATE ATOMS BY INTERSECTING CD AND COUNTY SHAPEFILES ##
###########################################################

# Get congressional district boundaries using tigris
message("Getting congressional district boundaries...")
cd <- tigris::congressional_districts(cb = TRUE, year = DATA_YEAR) %>%
  filter(!(STATEFP %in% c("02", "15", "72")))  # Filter out AK, HI, PR

# Get county boundaries using tigris
message("Getting county boundaries...")
county <- tigris::counties(cb = TRUE, year = DATA_YEAR) %>%
  filter(!(STATEFP %in% c("02", "15", "72")))  # Filter out AK, HI, PR

# Ensure both datasets use the same CRS
cd <- st_transform(cd, st_crs(county))

# Perform the intersection to create atoms
message("Performing county and CD intersection to create atoms...")
# Using st_intersection instead of sp::intersect
ints <- st_intersection(county, cd)

# Create a unique identifier for each atom
ints <- ints %>%
  mutate(
    atom_id = row_number(),
    county_geoid = GEOID,        # County GEOID
    cd_geoid = GEOID.1,          # CD GEOID
    state = STATEFP              # State FIPS from county
  )

# Debug output for problematic counties
message("Checking specific counties of interest...")
orange_nc <- ints %>% 
  filter(grepl("Orange", NAMELSAD) & state == "37")  # NC is FIPS 37
jackson_al <- ints %>% 
  filter(grepl("Jackson", NAMELSAD) & state == "01")  # AL is FIPS 01

message(paste("Orange County, NC has", nrow(orange_nc), "atoms after filtering"))
message(paste("Jackson County, AL has", nrow(jackson_al), "atoms after filtering"))

# Simplify the dataset to only include necessary columns
ints <- ints %>%
  select(atom_id, county_geoid, cd_geoid, state, NAMELSAD, NAMELSAD.1, ALAND, geometry) %>%
  rename(
    county_name = NAMELSAD,
    cd_name = NAMELSAD.1,
    county_area = ALAND
  )

# Try to get DW-NOMINATE data if available
message("Checking for DW-NOMINATE data...")
lancet <- readRDS('./data/districts_covid19_icu_and_covariates.rds')

# Process DW-NOMINATE data if available
if ("nominate_dim1" %in% names(lancet)) {
  lancet_nogeom <- lancet %>%
    select(GEOID, nominate_dim1)
  
  st_geometry(lancet_nogeom) <- NULL
  lancet_final <- as.data.frame(lancet_nogeom)
  
  # Merge DW-NOMINATE data with atoms
  ints <- ints %>%
    left_join(lancet_final, by = c("cd_geoid" = "GEOID"))
} else {
  message("DW-NOMINATE data not found in districts file.")
  # Add placeholder column
  ints$nominate_dim1 <- NA
}
  
# Save preliminary atom data
message("Saving preliminary atom data...")
save(ints, file = './data/county_cd_intersections.RData')

#########################################
## GET NUMBER OF ATOMS FOR EACH COUNTY ##
#########################################

# Load county data
message("Loading county data...")
load('./data/county_data.RData')

# Ensure county data is ordered by geoid
county <- county %>% arrange(geoid)

# Create a version without geometry for merging
county_4merge <- county
st_geometry(county_4merge) <- NULL

# Load atom data
message("Processing atom data to count atoms per county...")
atom <- ints

# Merge county data into atoms
atom <- atom %>%
  left_join(county_4merge, by = c("county_geoid" = "geoid"))

# Count number of atoms per county
countyXatom <- data.frame(table(atom$county_geoid))
names(countyXatom) <- c('geoid', 'num_atoms')

# Merge number of atoms into county dataset
county <- county %>%
  left_join(countyXatom, by = "geoid") %>%
  arrange(num_atoms, geoid)

# Fill in missing num_atoms values with 0
county$num_atoms[is.na(county$num_atoms)] <- 0

# Merge number of atoms into atom dataset
atom <- atom %>%
  left_join(countyXatom, by = c("county_geoid" = "geoid")) %>%
  arrange(num_atoms, county_geoid)

# Create out_ind index (position in county dataset)
county_index <- data.frame(
  county_geoid = county$geoid,
  out_ind = 1:nrow(county)
)

atom <- atom %>%
  left_join(county_index, by = "county_geoid")

# Identify island counties (due to issues with spatial adjacency matrices)
message("Identifying island counties...")
county_nb <- spdep::poly2nb(as(county, 'Spatial'))
islands <- which(sapply(county_nb, function(x) length(x) == 0))

if (length(islands) > 0) {
  message(paste("Removing", length(islands), "island counties..."))
  # Remove islands from county data
  county <- county[-islands, ]
  
  # Remove islands from atom data
  atom <- atom %>%
    filter(county_geoid %in% county$geoid)
}

# Calculate atom area
atom$atom_area <- as.numeric(st_area(atom))

# Save atom data
message("Saving atom data with county associations...")
save(atom, county, file = './data/atom_county_data.RData')

message("Atom creation process complete!")