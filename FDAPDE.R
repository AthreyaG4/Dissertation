library(tidyverse)
library(fdaPDE)
library(dplyr)
library(terra)
library(raster)
library(rnaturalearth)
library(INLA)
library(sf)
library(ggplot2)
library(leaflet)
library(viridis)

#get the map
map <- ne_countries(type = "countries",
                    country = "Ireland",
                    scale = "large", returnclass = "sf")

#plot the map
plot(sf::st_geometry(map))


#get the data
data = read_csv("./Weather Data/weather_data_long_clean.csv")

years = seq(2001, 2025)
param = "Mean Temperature"
for (yr in years) {
  #filter for particular year
  data_filtered = data |>
    filter(parameter == param & year == yr) |>
    dplyr::select(-AreaOfResidence, -year, -parameter)
  
  #data_filtered$value <- log1p(data_filtered$value)
  #view the points on the map
  points(data_filtered[,c("Longitude","Latitude")], col = "red")
  data_filtered |> summarise(n = n(), var = var(value), mean = mean(value))
  
  #make sf object
  data_sf = st_as_sf(data_filtered, coords = c("Longitude","Latitude"), crs=sf::st_crs(map))
  
  #filter out data points which lie outside the map boundaries
  data_sf = st_filter(data_sf, map)
  
  #z = scale(data_filtered$value)
  
  #make grid for mesh
  sq100 <- sf::st_make_grid(map, n = c(50,25)) # rect grid on bbox
  tri_pts <- sf::st_coordinates(sf::st_centroid(sq100[sf::st_union(map)])) # Triangulation points
  outer_pts <- sf::st_coordinates(sf::st_union(sq100[sf::st_union(map)])) # boundary points
  plot(outer_pts)
  
  tri_pts = rbind(outer_pts[, 1:2], tri_pts) # combine boundary points and inner points
  tri_pts = tri_pts[!duplicated(tri_pts),]
  tri_seg <- cbind(1:(nrow(outer_pts)-1), c(2:(nrow(outer_pts)-1),1))
  tri_mesh  = create.mesh.2D(nodes = tri_pts, segments = tri_seg)
  tri_meshr = refine.mesh.2D(tri_mesh,minimum_angle = 30,maximum_area = 0.03)
  plot(tri_meshr)
  points(data_sf, col = "red")
  FEM_basis <- fdaPDE::create.FEM.basis(tri_mesh)
  FEM_fit <- fdaPDE::smooth.FEM(as.matrix(data_filtered[,c(2,1)]), data_filtered[,3], FEM_basis, lambda = 10^-3,
                                DOF.evaluation = 'exact', 
                                lambda.selection.lossfunction = 'GCV')
  
  #prediction grid. 200 x 200 points. filtered out points which lie outside the map.
  grid <- terra::rast(map, nrows = 200, ncols = 200)
  xy <- terra::xyFromCell(grid, 1:ncell(grid))
  dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
                 crs = st_crs(map))
  indicespointswithin <- which(st_intersects(dp, map,
                                             sparse = FALSE))
  dp <- st_filter(dp, map)
  pred_coords = st_coordinates(dp) #predict for these points. around 20000 data points which lie inside the map.
  
  FEM_zfit <- fdaPDE::eval.FEM(FEM_fit$fit.FEM, pred_coords)
  
  #mean_z <- attr(z, "scaled:center")
  #sd_z <- attr(z, "scaled:scale")
  
  #pred_values <- FEM_zfit * sd_z + mean_z
  
  output = cbind(pred_coords, FEM_zfit)
  
  #make rastor object for plotting
  r <- terra::rast(output, type = "xyz")
  crs(r) <- "EPSG:4326"
  saveRDS(r, paste0("rastors/fdapde/", param, "_", yr, "_rastor.rds"))
}

rb <- raster::brick(r)
contours <- rasterToContour(rb, levels = pretty(range(values(rb), na.rm = TRUE), n = 10))
contours_sf <- st_as_sf(contours)
pal <- colorNumeric("viridis", values(r),
                    na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(rb, colors = pal, opacity = 0.8) %>%
  #addPolylines(data = contours_sf, color = "black", weight = 1, opacity = 0.7) %>%
  addLegend(pal = pal, values = values(r), title = "input$parameter")














