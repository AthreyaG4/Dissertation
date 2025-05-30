library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(stringr)
library(sf)
library(leaflet)
library(shinyjs)

data = st_read("Counties/shapefiles/shapefile.shp")

weather_data = read_csv("Weather Data/weather_data_long_clean.csv")

dsf = st_as_sf(weather_data, coords = c("Longitude", "Latitude"))
st_crs(dsf) = 4326

ui = page_sidebar(
  sidebar = sidebar(
    useShinyjs(),
    radioButtons( 
      inputId = "radio", 
      label = "Radio buttons", 
      choices = list(
        "Mortality Data" = 1, 
        "Weather Data" = 2
      ) 
    ),
    selectInput("year1", "Year", choices = unique(data$Year)),
    selectInput("cod", "Cause of Death", choices = unique(data$ICD10DG)),
    selectInput("year2", "Year", choices = unique(weather_data$year)),
    selectInput("parameter", "Parameters", choices = unique(weather_data$parameter)),
    sliderInput("slider", "Slider", 
                min = 0, max = 1, value = 0.8)
  ),
  leafletOutput("map")
)

server = function(input, output, session) {
  observe({
    if (input$radio == 1) {
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
  })
  
  subsetted1 = bindCache(
    reactive({
      req(input$year1, input$cod)
      data |> filter(Year == input$year1 & ICD10DG == input$cod)
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
  
  output$map = renderLeaflet({
    if (input$radio == 1) {
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
      pal = colorNumeric(palette = "viridis", domain = abs(subsetted2()$value))
      l = leaflet(subsetted2()) |> addTiles() |>
        addCircles(lng = st_coordinates(subsetted2())[, 1],
                   lat = st_coordinates(subsetted2())[, 2],
                   radius = abs(subsetted2()$value)/max(abs(subsetted2()$value), na.rm = TRUE)*5000,
                   color = ~pal(abs(subsetted2()$value)), popup = ~paste(AreaOfResidence, ":", subsetted2()$value)) |>
        addLegend(pal = pal, values = ~value)
      l
    }

  })
}
shinyApp(ui, server)
