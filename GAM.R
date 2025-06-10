library(ggplot2)
library(tidyverse)
library(stringr)
library(sf)
library(leaflet)
library(viridis)
library(mgcv)
library(terra)
library(raster)
library(rnaturalearth)
library(INLA)

data = read_csv("./Weather Data/weather_data_long_clean.csv")
grid = readRDS(file = "ireland_grid.rds")
years = seq(2001,2025)

parameters = unique(data$parameter)

if (!dir.exists("models")) dir.create("models")
if (!dir.exists("models/gam")) dir.create("models/gam")

for(param in parameters){
  gam_data = data
  gam_data = gam_data |>
    filter(parameter == param) |>
    dplyr::select(year, Latitude, Longitude, value)
  
  K = length(unique(paste(gam_data$Latitude, gam_data$Longitude)))
  
  model = mgcv::gam(value ~ s(year) + s(Latitude, Longitude, k = K), data = gam_data)
  saveRDS(model, file = paste0("models/gam/", param, "_model.rds"))
  
  for (yr in years) {
    if (!dir.exists("rastors")) dir.create("rastors")
    if (!dir.exists("rastors/gam")) dir.create("rastors/gam")
    grid_df <- grid |>
      sf::st_coordinates() |>
      as.data.frame() |>
      setNames(c("Longitude", "Latitude")) |>
      dplyr::mutate(year = yr)
    grid_df$value <- predict(model, newdata = grid_df)
    
    r <- rast(grid_df[, c("Longitude", "Latitude", "value")], type="xyz")
    crs(r) <- "EPSG:4326"
    saveRDS(r, file = paste0("rastors/gam/",param, "_", yr,"_rastor.rds"))
  }
}

parameter = "Precipitation Amount"
year = "2001"

rastor_file = readRDS(paste0("rastors/gam/", parameter, "_", year, "_rastor.rds"))

rb <- raster::brick(rastor_file)
contours <- rasterToContour(rb, levels = pretty(range(values(rb), na.rm = TRUE), n = 10))
contours_sf <- st_as_sf(contours)

pal <- colorNumeric("viridis", values(rastor_file),
                    na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(rb, colors = pal, opacity = 0.8) %>%
  addPolylines(data = contours_sf, color = "black", weight = 1, opacity = 0.7) %>%
  addLegend(pal = pal, values = values(rastor_file), title = "precipitation")