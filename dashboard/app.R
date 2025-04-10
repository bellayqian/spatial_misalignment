# app.R
library(shiny)
library(shinydashboard)
library(sf)
library(ggplot2)
library(DT)
library(plotly)
library(USAboundaries)
library(RColorBrewer)
library(cowplot)
library(shinyjs)

# Load helper functions
source("R/map_functions.R")

# Load data
load("data/final_misaligned_data.RData")

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
      menuItem("Dashboard Home", tabName = "home", icon = icon("dashboard")),
      menuItem("Maps Explorer", tabName = "maps", icon = icon("map")),
      menuItem("Data Explorer", tabName = "data", icon = icon("table")),
      menuItem("Methods Comparison", tabName = "methods", icon = icon("chart-bar")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    
    # State selection dropdown - changed to selectizeInput for better search
    selectizeInput(
      "state_filter",
      "Select a State (optional):",
      choices = c("All States" = "", sort(unique(us_states()$stusps))),
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
    ),
    
    # Add detail level slider for 3D maps
    conditionalPanel(
      condition = "input.sidebarMenu == 'home'",
      sliderInput(
        "detail_level",
        "Map Detail Level:",
        min = 0.001,
        max = 0.1,
        value = 0.01,
        step = 0.001
      )
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
      # Dashboard Home tab
      tabItem(tabName = "home",
              fluidRow(
                box(
                  width = 12,
                  title = "Atom-Based Regression Models for Misaligned Spatial Data",
                  status = "primary",
                  solidHeader = TRUE,
                  p("This dashboard provides interactive visualizations of counties, congressional districts, and atoms for exploring spatial misalignment."),
                  p("The analytic dataset contains 3,104 counties, 432 congressional districts, and 3,728 atoms.")
                )
              ),
              fluidRow(
                box(
                  width = 12,
                  title = "3D Visualization",
                  status = "primary",
                  radioButtons("map_type", "Select Map Type:",
                               choices = list(
                                 "3D County Boundaries" = "counties", 
                                 "3D Congressional Districts" = "districts",
                                 "3D Atom Counts" = "atoms"
                               ),
                               selected = "counties",
                               inline = TRUE),
                  div(
                    class = "map-container",
                    plotlyOutput("map_3d", height = "600px"),
                    tags$div(id = "loading-spinner", class = "loading-spinner")
                  ),
                  p("Hover over different elements to view details. Use the buttons below to navigate to detailed views.")
                )
              ),
              fluidRow(
                column(width = 4,
                       actionButton("view_counties", "View Counties", 
                                    icon = icon("map"), 
                                    class = "action-button",
                                    style="color: #fff; background-color: #ff9800; border-color: #e68a00; width: 100%")
                ),
                column(width = 4,
                       actionButton("view_districts", "View Congressional Districts", 
                                    icon = icon("map-marker"), 
                                    class = "action-button",
                                    style="color: #fff; background-color: #2196F3; border-color: #0c7cd5; width: 100%")
                ),
                column(width = 4,
                       actionButton("view_atoms", "View Atom Counts", 
                                    icon = icon("cubes"), 
                                    class = "action-button",
                                    style="color: #fff; background-color: #4CAF50; border-color: #3d8b40; width: 100%")
                )
              )
      ),
      
      # Maps Explorer tab
      tabItem(tabName = "maps",
              fluidRow(
                tabBox(
                  width = 12,
                  title = "Interactive Maps",
                  id = "mapTabBox",
                  tabPanel("County Boundaries", 
                           plotlyOutput("county_map_interactive", height = "600px")),
                  tabPanel("Congressional Districts", 
                           plotlyOutput("cd_map_interactive", height = "600px")),
                  tabPanel("Atom Counts", 
                           plotlyOutput("atom_map_interactive", height = "600px"))
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
                  p("This section compares traditional dasymetric mapping approaches with Atom-Based Regression Models."),
                  p("Interactive comparisons will be available here.")
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
                  p("To address these challenges, our work focuses on creating a comprehensive computational framework for ABRM, including developing open-source tools in R, to make ABRM more accessible to government agencies, health institutions, academic public health researchers, policy makers, and community-based organizations that produce, analyze, and/or use data from spatially misaligned geographic areas to inform their work for health equity.")
                )
              )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Show/hide loading spinner
  show_spinner <- function() {
    shinyjs::show("loading-spinner")
  }
  
  hide_spinner <- function() {
    shinyjs::hide("loading-spinner")
  }
  
  # Hide loading spinner initially
  observe({
    hide_spinner()
  })
  
  # Generate Congressional District dropdown based on selected state
  output$cd_selector <- renderUI({
    req(input$state_filter)
    
    if (input$state_filter == "") {
      return(NULL)
    }
    
    cd_data <- us_congressional()
    cd_data <- subset(cd_data, !(state_abbr %in% c('AK', 'HI', 'PR')))
    
    if (input$state_filter != "") {
      cd_data <- subset(cd_data, state_abbr == input$state_filter)
    }
    
    selectInput(
      "cd_filter",
      "Select a Congressional District:",
      choices = c("All Districts" = "", sort(unique(cd_data$cd_name))),
      selected = ""
    )
  })
  
  # Generate County dropdown based on selected state
  output$county_selector <- renderUI({
    req(input$state_filter)
    
    if (input$state_filter == "") {
      return(NULL)
    }
    
    county_data <- county
    
    # Check if state_abbr exists in county data
    if ("state_abbr" %in% colnames(county_data)) {
      county_data <- subset(county_data, state_abbr == input$state_filter)
    }
    
    selectInput(
      "county_filter",
      "Select a County:",
      choices = c("All Counties" = "", sort(unique(county_data$name))),
      selected = ""
    )
  })
  
  # Interactive county map
  output$county_map_interactive <- renderPlotly({
    # Create a basic ggplot map first
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    p <- generate_county_map(state_filter)
    
    # Convert to plotly
    ggplotly(p) %>%
      layout(title = "County Boundaries")
  })
  
  # Interactive congressional district map
  output$cd_map_interactive <- renderPlotly({
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    p <- generate_cd_map(state_filter)
    
    ggplotly(p) %>%
      layout(title = "Congressional District Boundaries")
  })
  
  # Interactive atom count map
  output$atom_map_interactive <- renderPlotly({
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    p <- generate_atom_map(state_filter)
    
    ggplotly(p) %>%
      layout(title = "Number of Atoms per County")
  })
  
  # Histogram of atom counts
  output$atom_histogram <- renderPlot({
    req(input$state_filter)
    
    filtered_county <- county
    if (input$state_filter != "") {
      if ("state_abbr" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state_abbr == input$state_filter)
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
  
  # 3D map output - responds to the map type selector
  output$map_3d <- renderPlotly({
    show_spinner()
    
    state_filter <- if(input$state_filter == "") NULL else input$state_filter
    detail_level <- input$detail_level
    
    # Select the appropriate 3D map based on radio button selection
    plot <- if(input$map_type == "counties") {
      generate_3d_county_map(state_filter, detail_level)
    } else if(input$map_type == "districts") {
      generate_3d_cd_map(state_filter, detail_level)
    } else {
      generate_3d_atom_map(state_filter, detail_level)
    }
    
    hide_spinner()
    return(plot)
  })
  
  # Add observers for the action buttons
  observeEvent(input$view_counties, {
    updateTabItems(session, "sidebarMenu", "maps")
    updateTabsetPanel(session, "mapTabBox", selected = "County Boundaries")
  })
  
  observeEvent(input$view_districts, {
    updateTabItems(session, "sidebarMenu", "maps")
    updateTabsetPanel(session, "mapTabBox", selected = "Congressional Districts")
  })
  
  observeEvent(input$view_atoms, {
    updateTabItems(session, "sidebarMenu", "maps")
    updateTabsetPanel(session, "mapTabBox", selected = "Atom Counts")
  })
  
  # County data table
  output$county_table <- DT::renderDataTable({
    filtered_county <- county
    
    # Apply state filter if selected
    if (input$state_filter != "") {
      if ("state_abbr" %in% colnames(filtered_county)) {
        filtered_county <- subset(filtered_county, state_abbr == input$state_filter)
      }
    }
    
    # Apply county filter if selected
    if (!is.null(input$county_filter) && input$county_filter != "") {
      filtered_county <- subset(filtered_county, name == input$county_filter)
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
      }
    }
    
    # Apply county filter if selected
    if (!is.null(input$county_filter) && input$county_filter != "") {
      filtered_atom <- subset(filtered_atom, county == input$county_filter)
    }
    
    # Apply CD filter if selected
    if (!is.null(input$cd_filter) && input$cd_filter != "") {
      filtered_atom <- subset(filtered_atom, cd_name == input$cd_filter)
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
}

# Run the application
shinyApp(ui = ui, server = server)