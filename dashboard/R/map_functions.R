# R/map_functions.R
library(sf)
library(ggplot2)
library(USAboundaries)
library(RColorBrewer)
library(plotly)

# Function to generate county boundaries map
generate_county_map <- function(state_filter = NULL) {
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Filter by state if specified
  if (!is.null(state_filter)) {
    county <- subset(county, state_abbr == state_filter)
  }
  
  aa <- brewer.pal(n = 8, name = "Dark2")
  
  countybds <- ggplot() +
    geom_sf(data = county, fill = 'white', lwd = .2, color = aa[2]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    theme_void() +
    labs(title = "Continental US County Boundaries")
  
  return(countybds)
}

# Function to generate congressional district boundaries map
generate_cd_map <- function(state_filter = NULL) {
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  cd <- us_congressional()
  cd <- subset(cd, !(state_abbr %in% c('AK', 'HI', 'PR')))
  
  # Filter by state if specified
  if (!is.null(state_filter)) {
    cd <- subset(cd, state_abbr == state_filter)
  }
  
  aa <- brewer.pal(n = 8, name = "Dark2")
  
  cdbds <- ggplot() +
    geom_sf(data = cd, fill = 'white', lwd = .2, color = aa[3]) +
    geom_sf(data = stbounds, lwd = .2, color = 'black', fill = NA) +
    theme_void() +
    labs(title = "2010 Congressional District Boundaries")
  
  return(cdbds)
}

# Function to generate atom count map
generate_atom_map <- function(state_filter = NULL) {
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Create a copy of county data for this function
  county_data <- county
  
  # Filter by state if specified
  if (!is.null(state_filter)) {
    county_data <- subset(county_data, state_abbr == state_filter)
  }
  
  # Create the category variable
  county_data$num_atoms_cat <- cut(
    county_data$num_atoms,
    breaks = c(0, 1.5, 2.5, 5.5, 10.5, 20),
    labels = c('1', '2', '3-5', '6-10', '10+')
  )
  
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

# Improved fortify function for sf objects
fortify <- function(sf_obj) {
  if (nrow(sf_obj) == 0) {
    return(data.frame(long = numeric(0), lat = numeric(0), id = character(0)))
  }
  
  # Create a unique identifier if not present
  if (!"plot_id" %in% names(sf_obj)) {
    sf_obj$plot_id <- 1:nrow(sf_obj)
  }
  
  # Extract coordinates
  coords <- st_coordinates(sf_obj)
  
  # Handle multipolygons correctly
  if ("L2" %in% colnames(coords)) {
    # For multipolygons
    df <- data.frame(
      long = coords[, "X"],
      lat = coords[, "Y"],
      id = coords[, "L1"],  # Feature ID
      piece = coords[, "L2"],  # Polygon part
      group = paste(coords[, "L1"], coords[, "L2"], sep = ".")
    )
  } else {
    # For simple polygons
    df <- data.frame(
      long = coords[, "X"],
      lat = coords[, "Y"],
      id = coords[, "L1"],  # Feature ID
      group = coords[, "L1"]
    )
  }
  
  # Now we need to correctly join the attribute data
  # First get the attribute data without geometry
  attr_data <- st_drop_geometry(sf_obj)
  
  # Create a mapping from the original feature ID to the row in attr_data
  id_map <- data.frame(
    id = 1:nrow(sf_obj),
    row_id = 1:nrow(sf_obj)
  )
  
  # Merge the attribute data using this mapping
  # We use merge with id_map first to get the right row
  df_with_attrs <- merge(df, id_map, by = "id", all.x = TRUE)
  
  # Now we can safely merge with the attribute data
  result <- cbind(df_with_attrs, attr_data[df_with_attrs$row_id, , drop = FALSE])
  result$row_id <- NULL  # Remove the temporary mapping column
  
  return(result)
}

# Updated generate_3d_interactive_map function
generate_3d_interactive_map_ <- function(state_filter = NULL, detail_level = 0.01) {
  # Get state boundaries
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Get county and congressional district data
  county_data <- county
  cd <- us_congressional()
  cd <- subset(cd, !(state_abbr %in% c('AK', 'HI', 'PR')))
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    # Check if the column exists before filtering
    if ("state_abbr" %in% colnames(county_data)) {
      county_data <- subset(county_data, state_abbr == state_filter)
    } else if ("state_name" %in% colnames(county_data)) {
      # Try with state_name if state_abbr is not available
      county_data <- subset(county_data, state_name == state_filter)
    }
    
    if ("state_abbr" %in% colnames(cd)) {
      cd <- subset(cd, state_abbr == state_filter)
    }
    
    stbounds <- subset(stbounds, stusps == state_filter)
  }
  
  # Check if we have data to display
  if (nrow(county_data) == 0 || nrow(cd) == 0) {
    # Return empty plot with message
    return(plot_ly() %>% 
             layout(title = "No data available for selected state",
                    annotations = list(
                      x = 0.5,
                      y = 0.5,
                      text = "Please select a different state or 'All States'",
                      showarrow = FALSE
                    )))
  }
  
  # Simplify geometries for better performance
  states_sf <- st_simplify(stbounds, dTolerance = detail_level)
  county_sf <- st_simplify(county_data, dTolerance = detail_level)
  cd_sf <- st_simplify(cd, dTolerance = detail_level)
  
  # Add height attributes before fortifying
  states_sf$height <- 0
  county_sf$height <- 0.1
  cd_sf$height <- 0.2
  
  # Convert to data frames for plotly using our improved fortify function
  states_df <- fortify(states_sf)
  counties_df <- fortify(county_sf)
  cds_df <- fortify(cd_sf)
  
  # Create 3D interactive plot
  p <- plot_ly(source = "map_3d") %>%
    # Add state boundaries as base layer
    add_trace(
      data = states_df,
      x = ~long, y = ~lat, z = ~height,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 2),
      name = "State Boundaries",
      hoverinfo = "name"
    ) %>%
    # Add county layer
    add_trace(
      data = counties_df,
      x = ~long, y = ~lat, z = ~height,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "orange", width = 1.5),
      name = "County Boundaries",
      hoverinfo = "text",
      hovertext = ~paste("County ID:", id, "<br>Atoms:", num_atoms)
    ) %>%
    # Add congressional district layer
    add_trace(
      data = cds_df,
      x = ~long, y = ~lat, z = ~height,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "blue", width = 1.5),
      name = "Congressional Districts",
      hoverinfo = "text",
      hovertext = ~paste("District ID:", id)
    ) %>%
    layout(
      title = "3D Map of Counties and Congressional Districts",
      scene = list(
        aspectmode = "data",
        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.3))
      )
    )
  
  return(p)
}

# Updated generate_3d_chloropleth function
generate_3d_chloropleth_ <- function(state_filter = NULL, detail_level = 0.01) {
  # Get state boundaries
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Get county data
  county_data <- county
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    # Check if the column exists before filtering
    if ("state_abbr" %in% colnames(county_data)) {
      county_data <- subset(county_data, state_abbr == state_filter)
    } else if ("state_name" %in% colnames(county_data)) {
      # Try with state_name if state_abbr is not available
      county_data <- subset(county_data, state_name == state_filter)
    }
    
    stbounds <- subset(stbounds, stusps == state_filter)
  }
  
  # Check if we have data to display
  if (nrow(county_data) == 0) {
    # Return empty plot with message
    return(plot_ly() %>% 
             layout(title = "No data available for selected state",
                    annotations = list(
                      x = 0.5,
                      y = 0.5,
                      text = "Please select a different state or 'All States'",
                      showarrow = FALSE
                    )))
  }
  
  # Simplify geometries for better performance
  states_sf <- st_simplify(stbounds, dTolerance = detail_level)
  county_sf <- st_simplify(county_data, dTolerance = detail_level)
  
  # Create categories for atom counts and add height
  county_sf$num_atoms_cat <- cut(
    county_sf$num_atoms,
    breaks = c(0, 1.5, 2.5, 5.5, 10.5, 20),
    labels = c('1', '2', '3-5', '6-10', '10+')
  )
  
  # Create a numeric value for 3D height based on number of atoms
  county_sf$height <- (county_sf$num_atoms / max(county_sf$num_atoms, na.rm = TRUE)) * 0.5
  
  # Set height for state boundaries
  states_sf$height <- 0
  
  # Convert to data frames for plotly using our improved fortify function
  states_df <- fortify(states_sf)
  counties_df <- fortify(county_sf)
  
  # Create 3D chloropleth map
  p <- plot_ly(source = "map_3d") %>%
    # Add state boundaries as base layer
    add_trace(
      data = states_df,
      x = ~long, y = ~lat, z = ~height,
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 2),
      name = "State Boundaries",
      hoverinfo = "name"
    ) %>%
    # Add counties with height based on atom count
    add_trace(
      data = counties_df,
      x = ~long, y = ~lat, z = ~height,
      type = "scatter3d",
      mode = "lines",
      line = list(color = ~num_atoms, colorscale = "YlGnBu", width = 2),
      name = "Counties by Atom Count",
      hoverinfo = "text",
      hovertext = ~paste("County ID:", id, "<br>Atoms:", num_atoms)
    ) %>%
    layout(
      title = "3D Chloropleth Map of Atom Counts by County",
      scene = list(
        aspectmode = "data",
        xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        zaxis = list(title = "Atom Count", showgrid = TRUE, zeroline = FALSE, range = c(0, 0.6))
      )
    )
  
  return(p)
}

# Simple 3D County Map function with fixed camera
generate_3d_county_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Get state boundaries
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Get county data
  county_data <- county
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    if ("state_abbr" %in% colnames(county_data)) {
      county_data <- subset(county_data, state_abbr == state_filter)
    }
    stbounds <- subset(stbounds, stusps == state_filter)
  }
  
  # Simplify geometries for better performance
  county_data <- st_simplify(county_data, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a more direct 3D representation using plotly's geo capabilities
  fig <- plot_ly()
  
  # Add county polygons
  for (i in 1:nrow(county_data)) {
    county_geom <- county_data[i,]
    coords <- st_coordinates(county_geom)[,1:2]
    
    # Use number of atoms for coloring
    atom_count <- county_data$num_atoms[i]
    
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
      hovertext = paste("County:", county_data$name[i], 
                        "<br>Atoms:", atom_count),
      showlegend = FALSE
    )
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
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
      hovertext = paste("State:", stbounds$name[i]),
      showlegend = FALSE
    )
  }
  
  # Configure the 3D scene with fixed camera
  fig <- fig %>% layout(
    title = "3D County Boundaries",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),  # Camera positioned above the map
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",  # Allow panning but not rotation
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.2))
    )
  )
  
  return(fig)
}

# Simple 3D Congressional District Map function with fixed camera
generate_3d_cd_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Get state boundaries
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Get congressional district data
  cd <- us_congressional()
  cd <- subset(cd, !(state_abbr %in% c('AK', 'HI', 'PR')))
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    if ("state_abbr" %in% colnames(cd)) {
      cd <- subset(cd, state_abbr == state_filter)
    }
    stbounds <- subset(stbounds, stusps == state_filter)
  }
  
  # Simplify geometries for better performance
  cd <- st_simplify(cd, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a more direct 3D representation
  fig <- plot_ly()
  
  # Add congressional district polygons
  for (i in 1:nrow(cd)) {
    cd_geom <- cd[i,]
    coords <- st_coordinates(cd_geom)[,1:2]
    
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
      hovertext = paste("District:", cd$cd_name[i]),
      showlegend = FALSE
    )
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
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
      hovertext = paste("State:", stbounds$name[i]),
      showlegend = FALSE
    )
  }
  
  # Configure the 3D scene with fixed camera
  fig <- fig %>% layout(
    title = "3D Congressional District Boundaries",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),  # Camera positioned above the map
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",  # Allow panning but not rotation
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.2))
    )
  )
  
  return(fig)
}

# Function for 3D Atom Chloropleth with fixed camera
generate_3d_atom_map <- function(state_filter = NULL, detail_level = 0.01) {
  # Get state boundaries
  stbounds <- USAboundaries::us_states()
  stbounds <- subset(stbounds, !(stusps %in% c("AK", "HI", "PR")))
  
  # Get county data
  county_data <- county
  
  # Filter by state if specified
  if (!is.null(state_filter) && state_filter != "") {
    if ("state_abbr" %in% colnames(county_data)) {
      county_data <- subset(county_data, state_abbr == state_filter)
    }
    stbounds <- subset(stbounds, stusps == state_filter)
  }
  
  # Simplify geometries for better performance
  county_data <- st_simplify(county_data, dTolerance = detail_level)
  stbounds <- st_simplify(stbounds, dTolerance = detail_level)
  
  # Create a color palette for atom counts
  atom_colors <- colorRampPalette(brewer.pal(5, "YlGnBu"))(5)
  
  # Create a more direct 3D representation
  fig <- plot_ly()
  
  # Add county polygons with height based on atom count
  for (i in 1:nrow(county_data)) {
    county_geom <- county_data[i,]
    coords <- st_coordinates(county_geom)[,1:2]
    
    # Calculate height based on atom count
    atom_count <- county_data$num_atoms[i]
    height <- min(0.5, atom_count / 10 * 0.5)  # Cap at 0.5 for visualization
    
    # Determine color based on atom count category
    color_idx <- min(5, ceiling(atom_count / 2))
    if (color_idx < 1) color_idx <- 1
    
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
      hovertext = paste("County:", county_data$name[i], 
                        "<br>Atoms:", atom_count),
      showlegend = FALSE
    )
  }
  
  # Add state boundaries for context
  for (i in 1:nrow(stbounds)) {
    state_geom <- stbounds[i,]
    coords <- st_coordinates(state_geom)[,1:2]
    
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
      hovertext = paste("State:", stbounds$name[i]),
      showlegend = FALSE
    )
  }
  
  # Configure the 3D scene with fixed camera
  fig <- fig %>% layout(
    title = "3D Map of Atom Counts by County",
    scene = list(
      aspectmode = "data",
      camera = list(
        eye = list(x = 0, y = -0.1, z = 2),  # Camera positioned above the map
        center = list(x = 0, y = 0, z = 0)
      ),
      dragmode = "pan",  # Allow panning but not rotation
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "Atom Count", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, range = c(0, 0.6))
    )
  )
  
  return(fig)
}

