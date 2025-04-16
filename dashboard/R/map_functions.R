# R/map_functions.R
library(sf)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(leaflet)
library(tigris)
library(dplyr)

# Set tigris options to cache data downloads
options(tigris_use_cache = TRUE)

# Function to load geographical data with caching
load_geo_data <- function(state_filter = NULL) {
  # Get state boundaries using tigris
  stbounds <- states(cb = TRUE)
  stbounds <- subset(stbounds, !(STUSPS %in% c("AK", "HI", "PR", "AS", "VI")))
  
  # Get county boundaries using tigris
  county_data <- counties(cb = TRUE)
  county_data <- subset(county_data, !(STATEFP %in% c("02", "15", "72", "60", "78"))) 
  # Exclude AK, HI, PR, American Samoa, Virgin Islands
  
  # Get congressional district boundaries using tigris
  cd <- congressional_districts(cb = TRUE)
  cd <- subset(cd, !(STATEFP %in% c("02", "15", "72", "60", "78")))
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    # Get the state FIPS code from the abbreviation
    state_info <- filter(stbounds, STUSPS == state_filter)
    if(nrow(state_info) > 0) {
      state_fips <- state_info$STATEFP[1]
      
      # Filter county and CD data
      county_data <- filter(county_data, STATEFP == state_fips)
      cd <- filter(cd, STATEFP == state_fips)
      stbounds <- filter(stbounds, STUSPS == state_filter)
    }
  }
  
  return(list(
    states = stbounds,
    counties = county_data,
    congressional_districts = cd
  ))
}

# Function to generate county boundaries map with leaflet
generate_county_leaflet <- function(state_filter = NULL, county_data = NULL, stbounds = NULL) {
  # If data not provided, load it
  if(is.null(county_data) || is.null(stbounds)) {
    geo_data <- load_geo_data(state_filter)
    county_data <- geo_data$counties
    stbounds <- geo_data$states
  }
  
  # Create a leaflet map
  map <- leaflet() %>%
    # Add base map tiles
    addProviderTiles(providers$CartoDB.Positron) %>%
    # Add county boundaries
    addPolygons(
      data = county_data,
      fillColor = "white",
      fillOpacity = 0.5,
      color = "#FF8C00",
      weight = 1,
      opacity = 0.7,
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#FF4500",
        fillOpacity = 0.7,
        bringToFront = TRUE
      ),
      label = ~NAME,
      popup = ~paste("<strong>County:</strong>", NAME, 
                     "<br><strong>State:</strong>", STUSPS)
    ) %>%
    # Add state boundaries
    addPolylines(
      data = stbounds,
      color = "black",
      weight = 1.5,
      opacity = 0.8,
      label = ~NAME
    )
  
  # Set view to continental US if no state filter
  if(is.null(state_filter) || state_filter == "") {
    map <- map %>% setView(-96, 38, zoom = 4)
  } else {
    # If state filter is provided, zoom to that state
    state_bounds <- st_bbox(stbounds)
    map <- map %>% fitBounds(
      state_bounds[["xmin"]], state_bounds[["ymin"]],
      state_bounds[["xmax"]], state_bounds[["ymax"]]
    )
  }
  
  return(map)
}

# Function to generate congressional district boundaries map with leaflet
generate_cd_leaflet <- function(state_filter = NULL, cd = NULL, stbounds = NULL) {
  # If data not provided, load it
  if(is.null(cd) || is.null(stbounds)) {
    geo_data <- load_geo_data(state_filter)
    cd <- geo_data$congressional_districts
    stbounds <- geo_data$states
  }
  
  # Create a leaflet map
  map <- leaflet() %>%
    # Add base map tiles
    addProviderTiles(providers$CartoDB.Positron) %>%
    # Add congressional district boundaries
    addPolygons(
      data = cd,
      fillColor = "white",
      fillOpacity = 0.5,
      color = "#4169E1",
      weight = 1,
      opacity = 0.7,
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#0000FF",
        fillOpacity = 0.7,
        bringToFront = TRUE
      ),
      label = ~paste("District:", CD),
      popup = ~paste("<strong>Congressional District:</strong>", CD, 
                     "<br><strong>State:</strong>", STUSPS)
    ) %>%
    # Add state boundaries
    addPolylines(
      data = stbounds,
      color = "black",
      weight = 1.5,
      opacity = 0.8,
      label = ~NAME
    )
  
  # Set view to continental US if no state filter
  if(is.null(state_filter) || state_filter == "") {
    map <- map %>% setView(-96, 38, zoom = 4)
  } else {
    # If state filter is provided, zoom to that state
    state_bounds <- st_bbox(stbounds)
    map <- map %>% fitBounds(
      state_bounds[["xmin"]], state_bounds[["ymin"]],
      state_bounds[["xmax"]], state_bounds[["ymax"]]
    )
  }
  
  return(map)
}

# Function to generate atom count choropleth map with leaflet
generate_atom_leaflet <- function(state_filter = NULL, county_data = NULL, stbounds = NULL) {
  # If data not provided, load it
  if(is.null(county_data) || is.null(stbounds)) {
    geo_data <- load_geo_data(state_filter)
    county_data <- geo_data$counties
    stbounds <- geo_data$states
  }
  
  # Ensure atom data is joined to county data
  # We assumes county_data has a 'num_atoms' field from the loaded RData

  # Create a color palette for the atom counts
  pal <- colorBin(
    palette = "YlGnBu",
    domain = county_data$num_atoms,
    bins = c(0, 1, 2, 5, 10, 20, Inf),
    na.color = "#FFFFFF"
  )
  
  # Create a leaflet map
  map <- leaflet() %>%
    # Add base map tiles
    addProviderTiles(providers$CartoDB.Positron) %>%
    # Add county polygons colored by number of atoms
    addPolygons(
      data = county_data,
      fillColor = ~pal(num_atoms),
      fillOpacity = 0.7,
      color = "#FF8C00",
      weight = 0.5,
      opacity = 0.7,
      highlightOptions = highlightOptions(
        weight = 2,
        color = "#FF4500",
        fillOpacity = 0.9,
        bringToFront = TRUE
      ),
      label = ~paste(NAME, ": ", num_atoms, "atoms"),
      popup = ~paste("<strong>County:</strong>", NAME, 
                     "<br><strong>State:</strong>", STUSPS,
                     "<br><strong>Number of Atoms:</strong>", num_atoms)
    ) %>%
    # Add state boundaries
    addPolylines(
      data = stbounds,
      color = "black",
      weight = 1.5,
      opacity = 0.8,
      label = ~NAME
    ) %>%
    # Add a legend
    addLegend(
      position = "bottomright",
      pal = pal,
      values = county_data$num_atoms,
      title = "Number of Atoms",
      opacity = 0.7,
      labFormat = function(type, cuts, p) {
        paste0(c("1", "2", "3-5", "6-10", "10+"))
      }
    )
  
  # Set view to continental US if no state filter
  if(is.null(state_filter) || state_filter == "") {
    map <- map %>% setView(-96, 38, zoom = 4)
  } else {
    # If state filter is provided, zoom to that state
    state_bounds <- st_bbox(stbounds)
    map <- map %>% fitBounds(
      state_bounds[["xmin"]], state_bounds[["ymin"]],
      state_bounds[["xmax"]], state_bounds[["ymax"]]
    )
  }
  
  return(map)
}

# Function to generate county boundaries map using ggplot (for backward compatibility)
generate_county_map <- function(state_filter = NULL) {
  # Load data
  geo_data <- load_geo_data(state_filter)
  stbounds <- geo_data$states
  county_data <- geo_data$counties
  
  aa <- brewer.pal(n = 8, name = "Dark2")
  
  countybds <- ggplot() +
    geom_sf(data = county_data, fill = 'white', lwd = .2, color = aa[2]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    theme_void() +
    labs(title = "Continental US County Boundaries")
  
  return(countybds)
}

# Function to generate congressional district boundaries map using ggplot (for backward compatibility)
generate_cd_map <- function(state_filter = NULL) {
  # Load data
  geo_data <- load_geo_data(state_filter)
  stbounds <- geo_data$states
  cd <- geo_data$congressional_districts
  
  aa <- brewer.pal(n = 8, name = "Dark2")
  
  cdbds <- ggplot() +
    geom_sf(data = cd, fill = 'white', lwd = .2, color = aa[3]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    theme_void() +
    labs(title = "Congressional District Boundaries")
  
  return(cdbds)
}

# Function to generate atom count map using ggplot (for backward compatibility)
generate_atom_map <- function(state_filter = NULL, county_data = NULL) {
  # Load data if not provided
  if(is.null(county_data)) {
    geo_data <- load_geo_data(state_filter)
    stbounds <- geo_data$states
    county_data <- geo_data$counties
  } else {
    # Get state boundaries using tigris
    stbounds <- states(cb = TRUE)
    stbounds <- subset(stbounds, !(STUSPS %in% c("AK", "HI", "PR")))
    
    if (!is.null(state_filter) && state_filter != "") {
      stbounds <- filter(stbounds, STUSPS == state_filter)
    }
  }
  
  # Ensure county_data has num_atoms field
  # Create the category variable if it doesn't exist
  if(!"num_atoms_cat" %in% names(county_data)) {
    county_data$num_atoms_cat <- cut(
      county_data$num_atoms,
      breaks = c(0, 1.5, 2.5, 5.5, 10.5, 20),
      labels = c('1', '2', '3-5', '6-10', '10+')
    )
  }
  
  aa <- brewer.pal(n = 8, name = "Dark2")
  
  numat_fill <- ggplot() +
    geom_sf(data = county_data, aes(fill = num_atoms_cat), lwd = .2, color = aa[2]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme_void() +
    labs(
      title = "Number of Atoms per County",
      fill = "Number of\nAtoms"
    )
  
  return(numat_fill)
}

# Keep 3D visualization functions for backward compatibility
# Can be updated depending on needs
generate_3d_county_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Load data
  geo_data <- load_geo_data(state_filter)
  stbounds <- geo_data$states
  county_data <- geo_data$counties
  
  # Simplify geometries for better performance
  county_data <- st_simplify(county_data, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a 3D representation using plotly
  fig <- plot_ly()
  
  # Add county polygons
  for (i in 1:min(nrow(county_data), 500)) { # Limit to 500 counties for performance
    county_geom <- county_data[i,]
    coords <- st_coordinates(county_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(0.1, nrow(coords)),  # Fixed height for counties
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "orange",
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("County:", county_data$NAME[i]),
        showlegend = FALSE
      )
    }
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(0, nrow(coords)),  # Base level for states
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "black",
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("State:", stbounds$NAME[i]),
        showlegend = FALSE
      )
    }
  }
  
  # Configure the 3D scene
  fig <- fig %>% layout(
    title = "3D County Boundaries",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.2))
    )
  )
  
  return(fig)
}

generate_3d_cd_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Load data
  geo_data <- load_geo_data(state_filter)
  stbounds <- geo_data$states
  cd <- geo_data$congressional_districts
  
  # Simplify geometries for better performance
  cd <- st_simplify(cd, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a 3D representation
  fig <- plot_ly()
  
  # Add congressional district polygons
  for (i in 1:min(nrow(cd), 500)) { # Limit to 500 districts for performance
    cd_geom <- cd[i,]
    coords <- st_coordinates(cd_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(0.1, nrow(coords)),  # Fixed height for districts
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "blue",
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("District:", cd$CD[i]),
        showlegend = FALSE
      )
    }
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(0, nrow(coords)),  # Base level for states
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "black",
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("State:", stbounds$NAME[i]),
        showlegend = FALSE
      )
    }
  }
  
  # Configure the 3D scene
  fig <- fig %>% layout(
    title = "3D Congressional District Boundaries",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.2))
    )
  )
  
  return(fig)
}

generate_3d_atom_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Load data
  geo_data <- load_geo_data(state_filter)
  stbounds <- geo_data$states
  county_data <- geo_data$counties
  
  # Simplify geometries for better performance
  county_data <- st_simplify(county_data, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a color palette for atom counts
  atom_colors <- colorRampPalette(brewer.pal(5, "YlGnBu"))(5)
  
  # Create a 3D representation
  fig <- plot_ly()
  
  # Add county polygons with height based on atom count
  for (i in 1:min(nrow(county_data), 500)) { # Limit to 500 counties for performance
    county_geom <- county_data[i,]
    coords <- st_coordinates(county_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      # Calculate height based on atom count (assuming num_atoms exists)
      # If it doesn't exist, use a default height of 0.1
      atom_count <- if("num_atoms" %in% names(county_data)) county_data$num_atoms[i] else 1
      height <- min(0.5, atom_count / 10 * 0.5)  # Cap at 0.5 for visualization
      if(is.na(height)) height <- 0.1
      
      # Determine color based on atom count category
      color_idx <- min(5, ceiling(atom_count / 2))
      if(is.na(color_idx) || color_idx < 1) color_idx <- 1
      
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(height, nrow(coords)),
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = atom_colors[color_idx],
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("County:", county_data$NAME[i], 
                          "<br>Atoms:", atom_count),
        showlegend = FALSE
      )
    }
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
    if(nrow(coords) > 0) {
      fig <- fig %>% add_trace(
        x = coords[,1], 
        y = coords[,2],
        z = rep(0, nrow(coords)),  # Base level for states
        type = "scatter3d",
        mode = "lines",
        line = list(
          color = "black",
          width = 2
        ),
        hoverinfo = "text",
        hovertext = paste("State:", stbounds$NAME[i]),
        showlegend = FALSE
      )
    }
  }
  
  # Configure the 3D scene
  fig <- fig %>% layout(
    title = "3D Map of Atom Counts by County",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "Atom Count", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.6))
    )
  )
  
  return(fig)
}

# Function to create forest plots for method comparison
create_forest_plot <- function(data, x_corr = 0.2, y_corr = 0.2) {
  # Add correlation values to the combined results
  plot_data <- data %>%
    filter(x_correlation == x_corr, y_correlation == y_corr) %>%
    group_by(method, variable, x_correlation, y_correlation) %>%
    summarize(
      mean_estimate = mean(estimated_beta),
      mean_lower = mean(ci_lower),
      mean_upper = mean(ci_upper),
      true_value = mean(true_beta), 
      mean_bias = mean(bias),
      mean_rel_bias = mean(relative_bias),
      coverage_rate = mean(within_ci) * 100,
      .groups = 'drop'
    )
  
  # Reorder the variable levels
  plot_data$variable <- factor(plot_data$variable, 
                               levels = c("covariate_x_1", "covariate_x_2", "covariate_x_3", 
                                          "covariate_y_1", "covariate_y_2", "covariate_y_3", "covariate_y_4"))
  
  # Create the forest plot
  ggplot(plot_data, aes(x = mean_estimate, y = variable, color = method)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = mean_lower, xmax = mean_upper), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(x = true_value), shape = 4, size = 3, color = "black") +
    scale_color_manual(values = c("ABRM" = "#F08080", "Dasymetric" = "#20B2AA")) +
    labs(
      x = "Coefficient Value", 
      y = "") +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = "top"
    )
}


