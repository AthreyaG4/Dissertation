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

lung_cancer = st_read("Counties/shapefiles/lungcancer.shp")
air_quality_data = read_csv("Data files/Air Quality data/air_quality_data_complete.csv")

air_quality_sf = st_as_sf(air_quality_data, coords = c("Longitude", "Latitude"))
st_crs(air_quality_sf) = 4326

ui = page_sidebar(
  sidebar = sidebar(
    useShinyjs(),
    radioButtons( 
      inputId = "data_select", 
      label = "Data", 
      choices = list(
        "Lung Cancer" = 1, 
        "Air Quality Data" = 2
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
    selectInput("model", "Model", choices = c("GAM", "FDAPDE")),
    checkboxInput("use_covariates", "Use covariates", value = FALSE),
    hr(),
    selectInput("pollutant", "Pollutant", choices = unique(air_quality_data$`Air Pollutant`)),
    sliderInput("slider", "Opacity", min = 0, max = 1, value = 0.8)
  ),
  leafletOutput("map")
)

server <- function(input, output, session) {
  LUNG <- 1
  AIR  <- 2
  
  # --- UI toggles ---
  observe({
    toggleState("pollutant", condition = (input$data_select == AIR))
    
    toggleState("model", condition = (input$actual_predicted == 2))
    
    show_cov <- (input$data_select == LUNG && input$actual_predicted == 2)
    shinyjs::toggle(id = "use_covariates", condition = show_cov)
    toggleState("use_covariates", condition = show_cov)
    
    if (!show_cov && isTRUE(input$use_covariates)) {
      updateCheckboxInput(session, "use_covariates", value = FALSE)
    }
  })
  
  # --- Reactives ---
  lung_actual_sf <- bindCache(
    reactive({
      req(input$data_select == LUNG, input$actual_predicted == 1)
      
      dat <- dplyr::filter(lung_cancer, Year == 2022)
      
      if (is.na(sf::st_crs(dat))) sf::st_crs(dat) <- 2157
      dat <- dat |>
        sf::st_make_valid() |>
        sf::st_transform(4326) |>
        sf::st_zm(drop = TRUE)
      
      validate(need("CsNmbrs" %in% names(dat),
                    "Expected column 'CsNmbrs' in lung cancer data."))
      dat
    }),
    input$data_select, input$actual_predicted
  )
  
  # Lung cancer AREAL predictions (SF with 'preds')
  subsetted2 <- bindCache(
    reactive({
      req(input$data_select == LUNG, input$actual_predicted == 2, input$model)
      if(isTRUE(input$use_covariates)){
        shp <- file.path("Output Files", "areal_outputs", tolower(input$model), "lungcancer_predictions_covariate_2022.shp")
      } else {
        shp <- file.path("Output Files", "areal_outputs", tolower(input$model), "lungcancer_predictions_2022.shp")
      }
      validate(need(file.exists(shp), paste("Missing shapefile:", shp)))
      st_read(shp, quiet = TRUE)
    }),
    input$data_select, input$actual_predicted, input$model, input$use_covariates
  )
  
  # Air quality prediction raster
  subsetted3 <- bindCache(
    reactive({
      req(input$data_select == AIR, input$actual_predicted == 2, input$pollutant, input$model)
      r_path <- file.path("Output Files", "rastors", tolower(input$model),
                          sprintf("%s_%d_rastor.rds", input$pollutant, 2022))
      validate(need(file.exists(r_path), paste("Missing raster:", r_path)))
      r <- readRDS(r_path)
      if (inherits(r, "SpatRaster")) r <- raster::raster(r)
      r
    }),
    input$data_select, input$actual_predicted, input$pollutant, input$model
  )
  
  # Air quality actual points (sf)
  subsetted4 <- bindCache(
    reactive({
      req(input$pollutant)
      aq <- air_quality_sf |>
        dplyr::filter(Year == 2022, `Air Pollutant` == input$pollutant, !is.na(value))
      validate(need(nrow(aq) > 0, "No air-quality points for this pollutant (2022)."))
      aq
    }),
    input$pollutant
  )
  
  # --- helpers ---
  base_map <- function() {
    leaflet(options = leafletOptions(minZoom = 4)) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(lng = -8, lat = 53.2, zoom = 6)
  }
  
  # --- render ---
  output$map <- renderLeaflet({
    m <- base_map()
    
    if (input$data_select == LUNG) {
      if (input$actual_predicted == 1) {
        dat  <- lung_actual_sf()   # already WGS84
        vals <- dat$CsNmbrs
        rng  <- range(vals, na.rm = TRUE)
        
        # log-scale colors (continuous)
        pal_log <- colorNumeric("YlOrRd", domain = log1p(rng), na.color = "transparent")
        
        # pick a name column if present
        name_col <- intersect(c("ArofRsd", "County", "NAME"), names(dat))
        nm <- if (length(name_col)) dat[[name_col[1]]] else rep("", nrow(dat))
        
        m |>
          addPolygons(
            data = dat,
            weight = 1, color = "white",
            fillOpacity = input$slider,
            fillColor = ~pal_log(log1p(CsNmbrs)),
            label = ~sprintf("%s: %s", nm, scales::comma(CsNmbrs)),
            labelOptions = labelOptions(
              direction = "auto", textsize = "12px", opacity = 0.9,
              style = list(
                "background-color" = "rgba(255,255,255,0.85)",
                "padding"          = "3px 6px",
                "border-radius"    = "4px",
                "box-shadow"       = "0 1px 3px rgba(0,0,0,0.2)"
              )
            ),
            highlightOptions = highlightOptions(
              weight = 2, color = "#333", bringToFront = TRUE
            )
          ) |>
          addLegend(
            "topright",
            pal     = pal_log,
            values  = log1p(vals),
            title   = "Lung cancer (observed)",
            opacity = 1,
            labFormat = labelFormat(transform = expm1, digits = 0, big.mark = ",")
          )
      } else {
        pred <- subsetted2() |> sf::st_transform(4326) |> sf::st_zm(drop = TRUE)
        req(nrow(pred) > 0)
        validate(need("preds" %in% names(pred), "Expected 'preds' column missing."))
        
        vals <- pred$preds
        rng  <- range(vals, na.rm = TRUE)
        
        pal_log <- colorNumeric("YlOrRd", domain = log1p(rng), na.color = "transparent")
        
        name_col <- intersect(c("ArofRsd", "County", "NAME"), names(pred))
        nm <- if (length(name_col)) pred[[name_col[1]]] else rep("", nrow(pred))
        
        m |>
          addPolygons(
            data = pred,
            weight = 1, color = "white",
            fillOpacity = input$slider,
            fillColor = ~pal_log(log1p(preds)),
            label = ~sprintf("%s: %s", nm, scales::comma(preds)),
            labelOptions = labelOptions(
              direction = "auto", textsize = "12px", opacity = 0.9,
              style = list(
                "background-color" = "rgba(255,255,255,0.85)",
                "padding"          = "3px 6px",
                "border-radius"    = "4px",
                "box-shadow"       = "0 1px 3px rgba(0,0,0,0.2)"
              )
            ),
            highlightOptions = highlightOptions(
              weight = 2, color = "#333", bringToFront = TRUE
            )
          ) |>
          addLegend(
            "topright",
            pal     = pal_log,
            values  = log1p(vals),
            title   = "Lung cancer (observed)",
            opacity = 1,
            labFormat = labelFormat(transform = expm1, digits = 0, big.mark = ",")
          )
      }
      
    } else if (input$data_select == AIR) {
      if (input$actual_predicted == 1) {
        pts <- subsetted4()
        pal <- colorNumeric("viridis", pts$value, na.color = "transparent")
        m |>
          addCircleMarkers(
            data = pts, radius = 5, stroke = FALSE,
            fillOpacity = input$slider, fillColor = ~pal(value),
            label = ~sprintf("%s: %.2f", `Air Quality Station Name`, value)
          ) |>
          addLegend("topright", pal = pal, values = pts$value,
                    title = paste0(input$pollutant, " (2022)"), opacity = 1)
      } else {
        r <- subsetted3()
        rng <- range(values(r), na.rm = TRUE)
        pal <- colorNumeric("viridis", rng, na.color = "transparent")
        m |>
          addRasterImage(r, colors = pal, opacity = input$slider, project = TRUE) |>
          addLegend("topright", pal = pal, values = rng,
                    title = paste(input$model, input$pollutant, "2022"), opacity = 1)
      }
    }
    
  })
}
shinyApp(ui, server)
