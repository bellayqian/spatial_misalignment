# preprint maps - edited
library(sf)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(tigris)
library(leaflet)

# Set tigris options to use cached data
options(tigris_use_cache = TRUE)

# Load the data
load('./data/final_misaligned_data.RData')

# Get US state boundaries excluding AK, HI, and PR
stbounds <- tigris::states(cb = TRUE) %>%
  filter(!(STUSPS %in% c("AK", "HI", "PR")))

# Get county boundaries from tigris
county_tigris <- tigris::counties(cb = TRUE) %>%
  filter(!(STATEFP %in% c("02", "15", "72")))  # Filter out AK, HI, PR

# Get congressional districts
cd <- tigris::congressional_districts(cb = TRUE) %>%
  filter(!(STATEFP %in% c("02", "15", "72")))  # Filter out AK, HI, PR

# Make sure the CRS matches the loaded county data
if (exists("county") && is.data.frame(county) && "geometry" %in% names(county)) {
  county_tigris <- st_transform(county_tigris, st_crs(county))
  
  # Map CRS for consistency
  stbounds <- st_transform(stbounds, st_crs(county))
  cd <- st_transform(cd, st_crs(county))
  
  # Add atom data to county_tigris if needed
  if ("num_atoms" %in% names(county)) {
    # Perform a spatial join to associate counties with atom data
    county_tigris <- st_join(county_tigris, select(county, num_atoms), join = st_intersects, left = TRUE)
    
    # Replace NAs with 0 for counties that don't have atom data
    county_tigris$num_atoms[is.na(county_tigris$num_atoms)] <- 0
    
    # Create a category variable for atom counts
    county_tigris$num_atoms_cat <- cut(
      county_tigris$num_atoms,
      breaks = c(0, 1.5, 2.5, 5.5, 10.5, 20),
      labels = c('1', '2', '3-5', '6-10', '10+')
    )
  }
} else {
  # If original county data doesn't exist or isn't spatial, just use the tigris data
  message("Using tigris county data without atom information")
}

# Choose color palette
aa <- brewer.pal(n = 8, name = "Dark2")

## make map of counties intersected by atoms using ggplot2 ##
countybds <- ggplot() +
  geom_sf(data = county_tigris, fill = 'white', lwd = .2, color = aa[2]) +
  geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
  theme_void() +
  labs(title = "US Counties")

## Congressional district boundaries map
cdbds <- ggplot() +
  geom_sf(data = cd, fill = 'white', lwd = .2, color = aa[3]) +
  geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
  theme_void() +
  labs(title = "US Congressional Districts")

## Atom count map (if atom data is available)
if ("num_atoms_cat" %in% names(county_tigris)) {
  numat_fill <- ggplot() +
    geom_sf(data = county_tigris, aes(fill = num_atoms_cat), lwd = .2, color = aa[2]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme_void() +
    labs(title = "Number of Atoms per County", fill = 'Number of\nAtoms')
  
  # Create combined plot
  top_row <- plot_grid(countybds, cdbds, labels = c('A', 'B'), label_size = 12)
  combined_plot <- plot_grid(top_row, numat_fill, labels = c('', 'C'), label_size = 12, ncol = 1)
  
  # Display the combined plot
  print(combined_plot)
  
  # Uncomment to save the plot
  # pdf('./figures/maps_together.pdf', width=10)
  # print(combined_plot)
  # dev.off()
}

# Create Leaflet interactive maps
# Transform to WGS84 for Leaflet
county_leaflet <- st_transform(county_tigris, 4326)
cd_leaflet <- st_transform(cd, 4326)
states_leaflet <- st_transform(stbounds, 4326)

# County map with Leaflet
if ("num_atoms_cat" %in% names(county_leaflet)) {
  # Create color palette for atom counts
  pal <- colorFactor(
    palette = "YlGnBu",
    domain = county_leaflet$num_atoms_cat
  )
  
  # Interactive county map with atom counts
  county_interactive <- leaflet(county_leaflet) %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(
      fillColor = ~pal(num_atoms_cat),
      weight = 1,
      opacity = 1,
      color = "gray",
      dashArray = "",
      fillOpacity = 0.7,
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE
      ),
      label = ~paste0(NAME, ": ", num_atoms, " atoms"),
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto"
      )
    ) %>%
    addPolylines(
      data = states_leaflet,
      color = "black",
      weight = 2,
      opacity = 1
    ) %>%
    addLegend(
      position = "bottomright",
      pal = pal,
      values = ~num_atoms_cat,
      title = "Number of Atoms",
      opacity = 1
    )
} else {
  # Basic county map without atom data
  county_interactive <- leaflet(county_leaflet) %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    addPolygons(
      fillColor = "lightgreen",
      weight = 1,
      opacity = 1,
      color = "gray",
      dashArray = "",
      fillOpacity = 0.5,
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#666",
        dashArray = "",
        fillOpacity = 0.7,
        bringToFront = TRUE
      ),
      label = ~NAME,
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "15px",
        direction = "auto"
      )
    ) %>%
    addPolylines(
      data = states_leaflet,
      color = "black",
      weight = 2,
      opacity = 1
    )
}

# Congressional district map with Leaflet
cd_interactive <- leaflet(cd_leaflet) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = "lightblue",
    weight = 1,
    opacity = 1,
    color = "blue",
    dashArray = "",
    fillOpacity = 0.5,
    highlightOptions = highlightOptions(
      weight = 2,
      color = "#666",
      dashArray = "",
      fillOpacity = 0.7,
      bringToFront = TRUE
    ),
    label = ~paste0("District: ", NAMELSAD),
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"
    )
  ) %>%
  addPolylines(
    data = states_leaflet,
    color = "black",
    weight = 2,
    opacity = 1
  )

# Create a combined interactive map with both counties and congressional districts
combined_interactive <- leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron, group = "Base Map") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  
  # Add county layer
  addPolygons(
    data = county_leaflet,
    fillColor = if ("num_atoms_cat" %in% names(county_leaflet)) {
      ~pal(num_atoms_cat)
    } else {
      "lightgreen"
    },
    weight = 1,
    opacity = 1,
    color = "gray",
    dashArray = "",
    fillOpacity = 0.7,
    group = "Counties",
    highlightOptions = highlightOptions(
      weight = 2,
      color = "#666",
      dashArray = "",
      fillOpacity = 0.7,
      bringToFront = TRUE
    ),
    label = if ("num_atoms" %in% names(county_leaflet)) {
      ~paste0(NAME, ": ", num_atoms, " atoms")
    } else {
      ~NAME
    },
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"
    )
  ) %>%
  
  # Add CD layer
  addPolygons(
    data = cd_leaflet,
    fillColor = "transparent",
    weight = 2,
    opacity = 1,
    color = "blue",
    dashArray = "5",
    fillOpacity = 0,
    group = "Congressional Districts",
    highlightOptions = highlightOptions(
      weight = 3,
      color = "blue",
      dashArray = "",
      fillOpacity = 0.2,
      bringToFront = TRUE
    ),
    label = ~paste0("District: ", NAMELSAD),
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", padding = "3px 8px"),
      textsize = "15px",
      direction = "auto"
    )
  ) %>%
  
  # Add state boundaries
  addPolylines(
    data = states_leaflet,
    color = "black",
    weight = 2,
    opacity = 1,
    group = "State Boundaries"
  ) %>%
  
  # Add layer controls
  addLayersControl(
    baseGroups = c("Base Map", "Satellite"),
    overlayGroups = c("Counties", "Congressional Districts", "State Boundaries"),
    options = layersControlOptions(collapsed = FALSE)
  )

# Add legend if atom data is available
if ("num_atoms_cat" %in% names(county_leaflet)) {
  combined_interactive <- combined_interactive %>%
    addLegend(
      position = "bottomright",
      pal = pal,
      values = ~county_leaflet$num_atoms_cat,
      title = "Number of Atoms",
      opacity = 1
    )
}

# Display the interactive maps
# Note: In R Studio or R Markdown, these will display automatically
# In a Shiny app, you'll need to use the renderLeaflet function
print("Interactive maps created. Access them through the variables:")
print("county_interactive, cd_interactive, combined_interactive")