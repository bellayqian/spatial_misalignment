# app.R
library(shiny)
library(shinydashboard)
library(sf)
library(ggplot2)
library(DT)
library(plotly)
library(leaflet)
library(tigris)
library(RColorBrewer)
library(cowplot)
library(shinyjs)
library(dplyr)

# Set tigris options to cache data and prevent duplicate downloads
options(tigris_use_cache = TRUE)

# Load helper functions
source("R/map_functions.R")

# Load data
load("data/atom_county_data.RData")

# Get state list only once for UI dropdown
state_list <- states(cb = TRUE)
state_list <- subset(state_list, !(STUSPS %in% c("AK", "HI", "PR", "AS", "VI")))
state_choices <- c("All States" = "", sort(unique(state_list$STUSPS)))

# UI Definition
ui <- dashboardPage(
  dashboardHeader(
    title = "Atom-Based Regression Models Dashboard",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      id = "sidebarMenu",
      menuItem("Maps Explorer", tabName = "maps", icon = icon("map")),
      menuItem("Data Explorer", tabName = "data", icon = icon("table")),
      menuItem("Methods Comparison", tabName = "methods", icon = icon("chart-bar")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # State selection dropdown
    selectizeInput(
      "state_filter",
      "Select a State (optional):",
      choices = state_choices,
      selected = "",
      options = list(
        placeholder = 'Type to search...',
        onInitialize = I('function() { this.setValue(""); }')
      )
    ),
    
    # Add congressional district selector when a state is selected
    conditionalPanel(
      condition = "input.state_filter != ''",
      uiOutput("cd_selector")
    ),
    
    # Add county selector when a state is selected
    conditionalPanel(
      condition = "input.state_filter != ''",
      uiOutput("county_selector")
    )
  ),
  
  dashboardBody(
    # Include custom CSS
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    
    # Add shinyjs
    useShinyjs(),
    
    tabItems(
      # Maps Explorer tab (now the main tab)
      tabItem(tabName = "maps",
              fluidRow(
                box(
                  width = 12,
                  title = "Atom-Based Regression Models for Misaligned Spatial Data",
                  status = "primary",
                  solidHeader = TRUE,
                  p("This dashboard provides interactive visualizations of counties, congressional districts, and atoms for exploring spatial misalignment."),
                  uiOutput("dataset_summary")
                )
              ),
              fluidRow(
                tabBox(
                  width = 12,
                  title = "Interactive Maps",
                  id = "mapTabBox",
                  tabPanel("County Boundaries", 
                           leafletOutput("county_map_leaflet", height = "600px")),
                  tabPanel("Congressional Districts", 
                           leafletOutput("cd_map_leaflet", height = "600px")),
                  tabPanel("Atom Counts", 
                           leafletOutput("atom_map_leaflet", height = "600px"))
                )
              ),
              fluidRow(
                conditionalPanel(
                  condition = "input.state_filter != '' && input.mapTabBox == 'Atom Counts'",
                  box(
                    width = 12,
                    title = "Atom Count Distribution",
                    plotOutput("atom_histogram", height = "300px")
                  )
                )
              )
      ),
      
      # Data Explorer tab
      tabItem(tabName = "data",
              fluidRow(
                box(
                  width = 12,
                  title = "Data Explorer",
                  tabBox(
                    width = 12,
                    tabPanel("County Data", DT::dataTableOutput("county_table")),
                    tabPanel("Atom Data", DT::dataTableOutput("atom_table"))
                  )
                )
              )
      ),
      
      # Methods Comparison tab
      tabItem(tabName = "methods",
              fluidRow(
                box(
                  width = 12,
                  title = "Comparison: Dasymetric Mapping vs ABRM",
                  status = "primary",
                  solidHeader = TRUE,
                  p("This section compares traditional dasymetric mapping approaches with Atom-Based Regression Models."),
                  p("Select a distribution type to view method comparison results.")
                )
              ),
              fluidRow(
                tabBox(
                  width = 12,
                  title = "Comparison Results",
                  id = "methodsTabBox",
                  tabPanel("Poisson Distribution - Count data", 
                           plotOutput("poisson_forest_plot", height = "600px")),
                  tabPanel("Binomial Distribution - Rate data", 
                           plotOutput("binomial_forest_plot", height = "600px")),
                  tabPanel("Normal Distribution - Continuous data", 
                           plotOutput("normal_forest_plot", height = "600px"))
                )
              )
      ),
      
      # About tab
      tabItem(tabName = "about",
              fluidRow(
                box(
                  width = 12,
                  title = "About Atom-Based Regression Models",
                  status = "info",
                  solidHeader = TRUE,
                  h3("Introduction"),
                  p("Spatial misalignment—which occurs when data on multiple variables are collected using mismatched geographic boundary definitions—is a longstanding challenge in public health research. For instance, congressional districts can cut across multiple counties, and environmental hazard zones may cross census tract boundaries, in both cases creating intersecting areas that complicate efforts to study the relationships between health outcomes and their social, political, and environmental determinants."),
                  
                  h3("The ABRM Approach"),
                  p("Atom-based regression models (ABRM) offer a promising alternative by using atoms, the intersecting areas of all relevant units, as the fundamental units of analysis. By preserving the original spatial resolution of the data, ABRM account for uncertainty in statistical relationships while offering a robust method for handling misaligned data."),
                  
                  h3("Project Goal"),
                  p("To address these challenges, our work focuses on creating a comprehensive computational framework for ABRM, including developing open-source tools in R, to make ABRM more accessible to government agencies, health institutions, academic public health researchers, policy makers, and community-based organizations that produce, analyze, and/or use data from spatially misaligned geographic areas to inform their work for health equity."),
                  
                  h3("Funding and Project Information"),
                  p("This work was funded by the Robert Wood Johnson Foundation, Grant 81746. Project details are provided below."),
                  
                  h4("Project Title:"),
                  p("Aligning spatially misaligned data for health equity analysis, action, and accountability"),
                  
                  h4("Principal Investigators:"),
                  p("Dr. Nancy Krieger (PI) and Dr. Rachel Nethery (co-PI)"),
                  
                  h4("Start date:"),
                  p("July 2024"),
                  
                  h4("Project team and collaborators:"),
                  tags$ul(
                    tags$li("Yunzhe Qian (Bella), MS (Research Assistant, Dept of Biostatistics, HSPH)"),
                    tags$li("Rachel Nethery, PhD (Assistant Professor, Dept of Biostatistics, HSPH)"),
                    tags$li("Nancy Krieger, PhD (Professor, Department of Social and Behavioral Sciences (SBS))"),
                    tags$li("Nykesha Johnson, MPH (Statistical Data Analyst/Data Manager, SBS, HSPH)")
                  )
                )
              )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  output$dataset_summary <- renderUI({
    # Count the number of counties
    county_count <- nrow(county)
    
    # Count the number of unique congressional districts
    # This assumes there's a cd or cd_id column in your atom dataset
    if ("cd" %in% names(atom)) {
      cd_count <- length(unique(atom$cd))
    } else if ("cd_id" %in% names(atom)) {
      cd_count <- length(unique(atom$cd_id))
    } else if ("cd_name" %in% names(atom)) {
      cd_count <- length(unique(atom$cd_name))
    } else {
      cd_count <- "unknown"
    }
    
    # Count the number of atoms
    atom_count <- nrow(atom)
    
    p(paste0("The analytic dataset contains ", 
             format(county_count, big.mark=","), " counties, ", 
             format(cd_count, big.mark=","), " congressional districts, and ", 
             format(atom_count, big.mark=","), " atoms."))
  })
  
  # Generate Congressional District dropdown based on selected state
  output$cd_selector <- renderUI({
    req(input$state_filter)
    
    if (input$state_filter == "") {
      return(NULL)
    }
    
    # Get congressional districts using tigris
    cd_data <- congressional_districts(cb = TRUE)
    state_info <- filter(state_list, STUSPS == input$state_filter)
    
    if(nrow(state_info) > 0) {
      state_fips <- state_info$STATEFP[1]
      cd_data <- filter(cd_data, STATEFP == state_fips)
    }
    
    # Find the correct CD column
    cd_col <- names(cd_data)[grep("CD\\d+FP", names(cd_data))]
    if(length(cd_col) == 0) cd_col <- "CD"
    
    selectInput(
      "cd_filter",
      "Select a Congressional District:",
      choices = c("All Districts" = "", paste0("District ", sort(unique(cd_data[[cd_col]])))),
      selected = ""
    )
  })
  
  # Generate County dropdown based on selected state
  output$county_selector <- renderUI({
    req(input$state_filter)
    
    if (input$state_filter == "") {
      return(NULL)
    }
    
    # Get counties for the selected state using tigris
    county_state_data <- counties(state = input$state_filter, cb = TRUE)
    print(head(county_state_data))
    selectInput(
      "county_filter",
      "Select a County:",
      choices = c("All Counties" = "", sort(unique(county_state_data$NAME))),
      selected = ""
    )
  })
  
  # Leaflet map outputs
  output$county_map_leaflet <- renderLeaflet({
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    
    # Load county data from tigris
    if(is.null(state_filter)) {
      county_data <- counties(cb = TRUE)
      county_data <- subset(county_data, !(STATEFP %in% c("02", "15", "72", "60", "78")))
    } else {
      county_data <- counties(state = state_filter, cb = TRUE)
    }
    
    # Get state boundaries
    stbounds <- states(cb = TRUE)
    stbounds <- subset(stbounds, !(STUSPS %in% c("AK", "HI", "PR", "AS", "VI")))
    
    if(!is.null(state_filter)) {
      stbounds <- subset(stbounds, STUSPS == state_filter)
    }
    
    # Create leaflet map
    map <- leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
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
        popup = ~paste("<strong>County:</strong>", NAME)
      ) %>%
      addPolylines(
        data = stbounds,
        color = "black",
        weight = 1.5,
        opacity = 0.8,
        label = ~NAME
      )
    
    # Set view based on filter
    if(is.null(state_filter)) {
      map <- map %>% setView(-96, 38, zoom = 4)
    } else {
      state_bounds <- st_bbox(stbounds)
      map <- map %>% fitBounds(
        state_bounds[["xmin"]], state_bounds[["ymin"]],
        state_bounds[["xmax"]], state_bounds[["ymax"]]
      )
    }
    
    return(map)
  })
  
  output$cd_map_leaflet <- renderLeaflet({
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    
    # Load congressional district data from tigris
    if(is.null(state_filter)) {
      cd_data <- congressional_districts(cb = TRUE)
      cd_data <- subset(cd_data, !(STATEFP %in% c("02", "15", "72", "60", "78")))
    } else {
      cd_data <- congressional_districts(state = state_filter, cb = TRUE)
    }

    # Get state boundaries
    stbounds <- states(cb = TRUE)
    stbounds <- subset(stbounds, !(STUSPS %in% c("AK", "HI", "PR", "AS", "VI")))
    
    if(!is.null(state_filter)) {
      stbounds <- subset(stbounds, STUSPS == state_filter)
    }
    
    # Find the correct CD column
    cd_col <- names(cd_data)[grep("CD\\d+FP", names(cd_data))]
    if(length(cd_col) == 0) cd_col <- "CD" # Fallback if pattern not found
    
    # Create leaflet map
    map <- leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      addPolygons(
        data = cd_data,
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
        label = ~paste("District:", get(cd_col)),
        popup = ~paste("<strong>Congressional District:</strong>", get(cd_col))
      ) %>%
      addPolylines(
        data = stbounds,
        color = "black",
        weight = 1.5,
        opacity = 0.8,
        label = ~NAME
      )
    
    # Set view based on filter
    if(is.null(state_filter)) {
      map <- map %>% setView(-96, 38, zoom = 4)
    } else {
      state_bounds <- st_bbox(stbounds)
      map <- map %>% fitBounds(
        state_bounds[["xmin"]], state_bounds[["ymin"]],
        state_bounds[["xmax"]], state_bounds[["ymax"]]
      )
    }
    
    return(map)
  })
  
  output$atom_map_leaflet <- renderLeaflet({
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    
    # Load county data from tigris
    if(is.null(state_filter)) {
      county_data <- counties(cb = TRUE)
      county_data <- subset(county_data, !(STATEFP %in% c("02", "15", "72", "60", "78", "66", "69", "74", "GU")))
    } else {
      county_data <- counties(state = state_filter, cb = TRUE)
    }
    
    # Get state boundaries
    stbounds <- states(cb = TRUE)
    stbounds <- subset(stbounds, !(STUSPS %in% c("AK", "HI", "PR", "AS", "VI", "GU", "MP")))
    
    if(!is.null(state_filter)) {
      stbounds <- subset(stbounds, STUSPS == state_filter)
      
      # Specific debug for Jackson County, Alabama
      if(state_filter == "AL") {
        morgan_county <- county_data %>% filter(NAME == "Jackson")
        
        # Check if Jackson County exists in tigris data
        if(nrow(morgan_county) > 0) {
          print("Jackson County GEOID in tigris:")
          print(morgan_county$GEOID)
          
          # Check if this GEOID exists in your RData county file
          morgan_in_county <- county %>% filter(geoid == morgan_county$GEOID[1])
          print(paste("Jackson County in RData? Found:", nrow(morgan_in_county) > 0))
          
          if(nrow(morgan_in_county) > 0) {
            print("Jackson County num_atoms in RData:")
            print(morgan_in_county$num_atoms)
          }
        } else {
          print("Jackson County not found in tigris data for Alabama")
        }
      }
    }
    
    # Fix the CRS inconsistency warning
    county_data <- st_transform(county_data, 4326)  # Transform to WGS84
    stbounds <- st_transform(stbounds, 4326)  # Transform to WGS84
    
    # Join with atom data
    # Prepare county data for joining
    county_atoms <- st_drop_geometry(county) %>%
      dplyr::select(geoid, num_atoms)
    
    # Debug - show distribution of atom counts before join
    print("num_atoms distribution in RData:")
    print(summary(county_atoms$num_atoms))
    
    # Join by GEOID
    county_data$geoid <- county_data$GEOID
    county_data <- left_join(county_data, county_atoms, by = "geoid")
    
    county_data$num_atoms[is.na(county_data$num_atoms)] <- 0

    # Create color palette for atom counts
    pal <- colorBin(
      palette = "YlGnBu",
      domain = c(0, max(county_data$num_atoms)),
      bins = c(1, 2, 3, 6, 11, Inf),
      na.color = "#FFFFFF"
    )
    
    # Check for undefined or NA values after color palette creation
    print("Any NA values in num_atoms after handling?")
    print(any(is.na(county_data$num_atoms)))
    
    # Check for any potential JavaScript "undefined" issues (NaN, Inf, etc.)
    print("Any infinite values in num_atoms?")
    print(any(is.infinite(county_data$num_atoms)))
    
    print("Any NaN values in num_atoms?")
    print(any(is.nan(county_data$num_atoms)))
    
    # Create leaflet map
    map <- leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
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
        # Use NAME.x or NAME based on what's available after the join
        label = ~paste(NAME, ": ", ifelse(num_atoms == 0, "No atom data", as.character(num_atoms)), " atoms"),
        popup = ~paste("<strong>County:</strong>", NAME, 
                       "<br><strong>Number of Atoms:</strong>", 
                       ifelse(num_atoms == 0, "No atom data", as.character(num_atoms)))
      ) %>%
      addPolylines(
        data = stbounds,
        color = "black",
        weight = 1.5,
        opacity = 0.8,
        label = ~NAME
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = county_data$num_atoms,
        title = "Number of Atoms",
        opacity = 0.7,
        labFormat = function(type, cuts, p) {
          c("1", "2", "3-5", "6-10", "10+")
        }
      )
    
    # Set view based on filter
    if(is.null(state_filter)) {
      map <- map %>% setView(-96, 38, zoom = 4)
    } else {
      state_bounds <- st_bbox(stbounds)
      map <- map %>% fitBounds(
        state_bounds[["xmin"]], state_bounds[["ymin"]],
        state_bounds[["xmax"]], state_bounds[["ymax"]]
      )
    }
    
    return(map)
  })
  
  # Histogram of atom counts
  output$atom_histogram <- renderPlot({
    req(input$state_filter)
    
    filtered_county <- county
    if (input$state_filter != "") {
      if ("state_abbr" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state_abbr == input$state_filter)
      } else if ("state" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state == input$state_filter)
      }
    }
    
    ggplot(filtered_county, aes(x = num_atoms)) +
      geom_histogram(fill = "steelblue", bins = 20) +
      labs(
        title = "Distribution of Atom Counts per County",
        x = "Number of Atoms",
        y = "Count of Counties"
      ) +
      theme_minimal()
  })
  
  # County data table
  output$county_table <- DT::renderDataTable({
    filtered_county <- county
    
    # Apply state filter if selected
    if (input$state_filter != "") {
      if ("state_abbr" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state_abbr == input$state_filter)
      } else if ("state" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state == input$state_filter)
      }
    }
    
    # Apply county filter if selected
    if (!is.null(input$county_filter) && input$county_filter != "") {
      name_col <- if("name" %in% colnames(filtered_county)) "name" else "NAME"
      filtered_county <- subset(filtered_county, get(name_col) == input$county_filter)
    }
    
    # Convert sf object to dataframe and select only non-geometry columns for display
    county_df <- st_drop_geometry(filtered_county)
    
    DT::datatable(
      county_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })
  
  # Atom data table
  output$atom_table <- DT::renderDataTable({
    filtered_atom <- atom
    
    # Apply state filter if selected
    if (input$state_filter != "") {
      if ("state_abbr" %in% colnames(filtered_atom)) {
        filtered_atom <- subset(filtered_atom, state_abbr == input$state_filter)
      } else if ("state" %in% colnames(filtered_atom)) {
        filtered_atom <- subset(filtered_atom, state == input$state_filter)
      }
    }
    
    # Apply county filter if selected
    if (!is.null(input$county_filter) && input$county_filter != "") {
      county_col <- if("county" %in% colnames(filtered_atom)) "county" else "county_name"
      filtered_atom <- subset(filtered_atom, get(county_col) == input$county_filter)
    }
    
    # Apply CD filter if selected
    if (!is.null(input$cd_filter) && input$cd_filter != "") {
      cd_col <- if("cd_name" %in% colnames(filtered_atom)) "cd_name" else "cd"
      filtered_atom <- subset(filtered_atom, get(cd_col) == input$cd_filter)
    }
    
    # Convert sf object to dataframe and select only non-geometry columns for display
    atom_df <- st_drop_geometry(filtered_atom)
    
    DT::datatable(
      atom_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })
  
  # Forest plot for Poisson distribution comparison
  poisson_data <- read.csv("data/poisson_sensitivity_analysis_results.csv")
  binomial_data <- read.csv("data/binomial_sensitivity_analysis_results.csv")
  normal_data <- read.csv("data/normal_sensitivity_analysis_results.csv")
  
  output$poisson_forest_plot <- renderPlot({
    create_poisson_forest_plot(poisson_data, 0.2, 0.6)
  })
  output$binomial_forest_plot <- renderPlot({
    create_binomial_forest_plot(binomial_data, 0.2, 0.6)
  })
  output$normal_forest_plot <- renderPlot({
    create_normal_forest_plot(normal_data, 0.2, 0.6)
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)