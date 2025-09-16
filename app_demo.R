library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(stringr)
library(sf)
library(leaflet)
library(shinyjs)
library(mgcv)
library(terra)
library(raster)
library(rnaturalearth)
library(viridis)
library(mapview)
library(leafsync)
library(leaflet.extras2)
library(stars)

# Load your data
lung_cancer_sf = st_read("Counties/shapefiles/lungcancer.shp")
lung_cancer_sf = lung_cancer_sf |> filter(Year == 2022)

air_quality_data = read_csv("Data files/Air Quality data/air_quality_data_complete.csv")
air_quality_data = air_quality_data |> filter(Year == 2022 & `Air Pollutant` == "PM10") |> na.omit()
air_quality_sf = st_as_sf(air_quality_data, coords = c("Longitude", "Latitude"))
st_crs(air_quality_sf) = 4326

ui = fluidPage(
  useShinyjs(),
  
  # Sidebar and main content in flex container
  div(class = "sidebar",
      radioButtons( 
        inputId = "actual_predicted", 
        label = "View", 
        choices = list("Actual Data" = 1, "Model Predictions" = 2)
      ),
      selectInput("model", "Model", choices = c("GAM", "FDAPDE")),
      hr()
  ),
  
  # Main content area
  div(class = "main-content",
      uiOutput("map_ui")
  )
)

server = function(input, output, session) {
  observe({
    if (input$actual_predicted == 2){
      enable("model")
    } else {
      disable("model")
    }
  })
  
  output$map_ui <- renderUI({
    if(input$actual_predicted == 1) {
      # Create mapview objects
      m1 <- mapview(air_quality_sf, zcol = "value", map.types = "CartoDB.Positron")
      m2 <- mapview(lung_cancer_sf, zcol = "value", map.types = "CartoDB.Positron")
    } else {
      if(input$model == "GAM"){
        r = readRDS(file = paste0("Output Files/rastors/gam/PM10_2022_rastor.rds"))
        lung_cancer_pred_sf = st_read("Output Files/areal_outputs/gam/lungcancer_predictions.shp")
      } else if (input$model == "FDAPDE") {
        r = readRDS(file = paste0("Output Files/rastors/fdapde/PM10_2022_rastor.rds"))
        lung_cancer_pred_sf = st_read("Output Files/areal_outputs/fdapde/lungcancer_predictions.shp")
      }
      m1 <- mapview(r, map.types = "CartoDB.Positron")
      m2 <- mapview(lung_cancer_pred_sf, zcol = "preds", map.types = "CartoDB.Positron")
    }
    # Use leafsync with proper synchronization
    synced <- leafsync::sync(m1, m2, sync = "all")
    
    # Return the synced maps directly - leafsync handles the side-by-side layout
    synced
  })
}

shinyApp(ui, server)