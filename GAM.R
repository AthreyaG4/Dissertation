library(ggplot2)
library(tidyverse)
library(sf)
library(mgcv)
library(terra)
library(rnaturalearth)
library(tmap)
library(gratia)
library(mapview)

map = ne_countries(type = "countries",
                    country = "Ireland",
                    scale = "large", returnclass = "sf")

data = read_csv(
  "Data files/Air Quality data/air_quality_data_complete.csv"
  )

# create dirs
if (!dir.exists("Output Files/rastors/")) {
  dir.create("Output Files/rastors/")
} 

if (!dir.exists("Output Files/rastors/gam/")) {
  dir.create("Output Files/rastors/gam/")
}

# helper function to fit the gam
fit_and_raster = function(param, yr, nrows = 1000, ncols = 1000) {
  df = data |>
    filter(`Air Pollutant` == param, Year == yr) |>
    dplyr::select(Latitude, Longitude, value) |>
    mutate(
      weight = (value - min(value, na.rm = TRUE)) /
        (max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) + 1e-4
    ) |>
    drop_na()
  
  if (nrow(df) < 10) {
    stop(sprintf("Not enough data for %s in %s", param, yr))
  }
  
  K = min(150, length(unique(paste(df$Latitude, df$Longitude))))
  
  mod = mgcv::gam(
    value ~ s(Longitude, Latitude, bs = "tp", k = K),
    data   = df,
    weights = df$weight,
    family  = gaussian(),
    method  = "REML"
  )
  
  print(summary(mod))
  print(appraise(mod))
  
  y_obs  = df$value
  y_pred = predict(mod, type = "response")  # fitted counts
  
  rmse = sqrt(mean((y_obs - y_pred)^2))
  print(rmse)
  
  saveRDS(mod, paste0("models/gam/", param,".rds"))
  
  grid = terra::rast(map, nrows = 1000, ncols = 1000)
  xy = terra::xyFromCell(grid, 1:ncell(grid))
  dp = st_as_sf(as.data.frame(xy), coords = c("x", "y"),
                 crs = st_crs(map))
  
  # Filter points within map boundaries
  dp = st_filter(dp, map)
  pred_coords = st_coordinates(dp)
  
  # Create prediction dataframe
  grid_df = as.data.frame(pred_coords)
  names(grid_df) = c("Longitude", "Latitude")  # This is correct
  
  # Make predictions
  grid_df$value = predict(mod, newdata = grid_df)
  
  # Create raster
  r = rast(
    grid_df[, c("Longitude", "Latitude", "value")], type="xyz"
  )
  crs(r) = "EPSG:4326"
  
  # save and return
  path = paste0(
    "Output Files/rastors/gam/", param, "_", yr, "_rastor.rds"
    )
  saveRDS(r, path)
  r
}

# build rasters for both parameters
r_pm25 = fit_and_raster("PM2.5", 2022)
r_pm10 = fit_and_raster("PM10",  2022)

# If you prefer to load from disk instead of refitting:
r_pm25 = readRDS("Output Files/rastors/gam/PM2.5_2022_rastor.rds")
r_pm10 = readRDS("Output Files/rastors/gam/PM10_2022_rastor.rds")

m_pm25 = tm_shape(r_pm25) +
  tm_raster(style = "cont", palette = "viridis", title = "PM2.5") +
  tm_layout(legend.outside = TRUE, frame = FALSE)

m_pm10 = tm_shape(r_pm10) +
  tm_raster(style = "cont", palette = "viridis", title = "PM10") +
  tm_layout(legend.outside = TRUE, frame = FALSE)

tmap_arrange(m_pm25, m_pm10, ncol = 2)