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

mortality_data = st_read("Counties/shapefiles/shapefile.shp")
weather_data = read_csv("Weather Data/weather_data_long_clean.csv")
ireland = st_read("Counties/shapefiles/ireland_boundary.shp")
grid = readRDS(file = "ireland_grid.rds")

grid_df <- grid |>
  sf::st_coordinates() |>
  as.data.frame() |>
  setNames(c("Longitude", "Latitude"))


dsf = st_as_sf(weather_data, coords = c("Longitude", "Latitude"))
st_crs(dsf) = 4326

ui = page_sidebar(
  sidebar = sidebar(
    useShinyjs(),
    radioButtons( 
      inputId = "data_select", 
      label = "Data", 
      choices = list(
        "Mortality Data" = 1, 
        "Weather Data" = 2
      ) 
    ),
    radioButtons( 
      inputId = "actual_predicted", 
      label = "View", 
      choices = list(
        "Actual Data" = 1, 
        "Model Predictions" = 2
      ) 
    ),
    selectInput("model", "Model", choices = c("GAM", "FDAPDE", "INLA")),
    hr(),
    selectInput("year1", "Year", choices = unique(mortality_data$Year)),
    selectInput("cod", "Cause of Death", choices = unique(mortality_data$ICD10DG)),
    selectInput("year2", "Year", choices = unique(weather_data$year)),
    selectInput("parameter", "Parameters", choices = unique(weather_data$parameter)),
    sliderInput("slider", "Opacity", 
                min = 0, max = 1, value = 0.8)
  ),
  leafletOutput("map")
)

server = function(input, output, session) {
  observe({
    if (input$data_select == 1) {
      enable("year1")
      enable("cod")
      disable("year2")
      disable("parameter")
    } else {
      enable("year2")
      enable("parameter")
      disable("year1")
      disable("cod")
    }
    
    if (input$actual_predicted == 2){
      enable("model")
    } else {
      disable("model")
    }
  })
  
  subsetted1 = bindCache(
    reactive({
      req(input$year1, input$cod)
      mortality_data |> filter(Year == input$year1 & ICD10DG == input$cod)
    }),
    input$year1,
    input$cod
  )
  
  subsetted2 = bindCache(
    reactive({
      req(input$year2, input$parameter)
      dsf |> filter(year == input$year2 & parameter == input$parameter & !is.na(value))
    }),
    input$year2,
    input$parameter
  )
  
  subsetted3 = bindCache(
    reactive({
      req(input$actual_predicted == 2)
      req(input$year2, input$parameter, input$model)
      if(input$model == "GAM"){
        print("GAM")
        #grid_df = grid_df |> mutate(year = as.numeric(input$year2))
        #model = readRDS(file = paste0("models/gam/", input$parameter, "_model.rds"))
        #grid_df$value = mgcv::predict.gam(model, newdata = grid_df)
        #r <- rast(grid_df[, c("Longitude", "Latitude", "value")], type="xyz")
        #crs(r) <- "EPSG:4326"
        r = readRDS(file = paste0("rastors/gam/", input$parameter,"_", (input$year2),"_rastor.rds"))
        return(r)
      } else if(input$model == "FDAPDE") {
        print("FDAPDE")
        r = readRDS(file = paste0("rastors/fdapde/", input$parameter,"_", input$year2,"_rastor.rds"))
        return(r)
      } else if(input$model == "INLA") {
        print("INLA")
        r = readRDS(file = paste0("rastors/inla/", input$parameter,"_", input$year2,"_rastor.rds"))
        return(r)
      }

    }),
    input$year2,
    input$parameter,
    input$model
  )
  
  output$map = renderLeaflet({
    if (input$data_select == 1) {
      pal = colorNumeric(palette = "YlOrRd", domain = subsetted1()$VALUE)
      l = leaflet(subsetted1()) |> 
        addTiles()
      l = l |> addPolygons(
        color = "white",
        fillColor = ~pal(VALUE),
        fillOpacity = input$slider,
        label = ~paste(ArofRsd, ":", VALUE),
        weight = 2,
        highlightOptions = highlightOptions(bringToFront = TRUE, weight = 4, color = "white")
      )
      l = l |> addLegend(pal = pal, values = ~VALUE, opacity = input$slider)
      l
    }
    else{
      if (input$actual_predicted == 1){
        pal = colorNumeric(palette = "viridis", domain = abs(subsetted2()$value))
        l = leaflet(subsetted2()) |> addTiles() |>
          addCircles(lng = st_coordinates(subsetted2())[, 1],
                     lat = st_coordinates(subsetted2())[, 2],
                     radius = abs(subsetted2()$value)/max(abs(subsetted2()$value), na.rm = TRUE)*5000,
                     color = ~pal(abs(subsetted2()$value)), popup = ~paste(AreaOfResidence, ":", subsetted2()$value)) |>
          addLegend(pal = pal, values = ~value)
        l
      } else {
        rb <- raster::brick(subsetted3())
        contours <- rasterToContour(rb, levels = pretty(range(values(rb), na.rm = TRUE), n = 10))
        contours_sf <- st_as_sf(contours)
        
        pal <- colorNumeric("viridis", values(subsetted3()),
                            na.color = "transparent")
        leaflet() %>% addTiles() %>%
          addRasterImage(rb, colors = pal, opacity = 0.8) %>%
          #addPolylines(data = contours_sf, color = "black", weight = 1, opacity = 0.7) %>%
          addLegend(pal = pal, values = values(subsetted3()), title = input$parameter)
      }
    }

  })
}
shinyApp(ui, server)
