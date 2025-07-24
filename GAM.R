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

# Get the map (this was missing in your code)
map <- ne_countries(type = "countries",
                    country = "Ireland",
                    scale = "large", returnclass = "sf")

# Load data
data = read_csv("Air Quality data/air_quality_data_complete.csv")
grid = readRDS(file = "ireland_grid.rds")

# Set parameters
years = c(2022)  # Add the years you actually want to process
parameters = c("PM2.5")

# Create directories
if (!dir.exists("models")) dir.create("models")
if (!dir.exists("models/gam")) dir.create("models/gam")
if (!dir.exists("rastors")) dir.create("rastors")
if (!dir.exists("rastors/gam")) dir.create("rastors/gam")

for(param in parameters){
  for (yr in years) {  # Move year loop outside to process each year
    
    # Filter data for specific parameter and year
    gam_data <- data |>
      filter(`Air Pollutant` == param, Year == yr) |>
      dplyr::select(Latitude, Longitude, value) |>
      mutate(
        weight = (value - min(value, na.rm = TRUE)) / (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) + 1e-4
      ) |>
      na.omit()
    
    
    # Check if we have enough data
    if(nrow(gam_data) < 10) {
      cat("Not enough data for", param, "in year", yr, "\n")
      next
    }
    
    K = length(unique(paste(gam_data$Latitude, gam_data$Longitude)))
    
    # FIX: Corrected 'fmaily' to 'family'
    model = mgcv::gam(value ~ s(Latitude, Longitude, bs="tp", k = K), 
                      data = gam_data, 
                      weights = weight,
                      family = "gaussian",  # Fixed typo
                      method = "REML")
    
    saveRDS(model, file = paste0("models/gam/", param, "_", yr, "_model.rds"))
    
    # Create prediction grid
    grid <- terra::rast(map, nrows = 200, ncols = 200)
    xy <- terra::xyFromCell(grid, 1:ncell(grid))
    dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
                   crs = st_crs(map))
    
    # Filter points within map boundaries
    dp <- st_filter(dp, map)
    pred_coords = st_coordinates(dp)
    
    # Create prediction dataframe
    grid_df = as.data.frame(pred_coords)
    names(grid_df) = c("Longitude", "Latitude")  # This is correct
    
    # Make predictions
    grid_df$value <- predict(model, newdata = grid_df)
    
    # Create raster
    r <- rast(grid_df[, c("Longitude", "Latitude", "value")], type="xyz")
    crs(r) <- "EPSG:4326"
    saveRDS(r, file = paste0("rastors/gam/", param, "_", yr, "_rastor.rds"))
    
    cat("Completed", param, "for year", yr, "\n")
  }
}

# Visualization
parameter = "PM2.5"
year = "2022"  # Make sure this matches a year you processed
rastor_file = readRDS(paste0("rastors/gam/", parameter, "_", year, "_rastor.rds"))
rb <- raster::brick(rastor_file)

pal <- colorNumeric("viridis", values(rastor_file),
                    na.color = "transparent")

leaflet() %>% 
  addTiles() %>%
  addRasterImage(rb, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(rastor_file), title = parameter)  # Fixed title
