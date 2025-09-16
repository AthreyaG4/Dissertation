library(tidyverse)
library(fdaPDE)
library(dplyr)
library(terra)
library(raster)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(leaflet)
library(viridis)
library(tmap)

#get the map
map = ne_countries(type = "countries",
                    country = "Ireland",
                    scale = "large", returnclass = "sf")

map = st_transform(map, crs = 4326)

#plot the map
plot(sf::st_geometry(map))


map_m = st_transform(map, 2157)  
ire    = st_make_valid(st_union(map_m))


ire_parts = st_cast(ire, "POLYGON")


ire_main = ire_parts[which.max(st_area(ire_parts))]

# Simplify and densify boundary
ire = st_simplify(ire, dTolerance = 50)
ire_dense = st_segmentize(ire, dfMaxLength = 5000)

rings = ire_dense |>
  st_boundary() |>
  st_cast("MULTILINESTRING") |>
  st_cast("LINESTRING")

all_nodes = list()
all_segs  = list()
offset = 0

for (i in seq_len(length(rings))) {
  crd = st_coordinates(rings[i])[, 1:2, drop = FALSE]
  
  # drop duplicated closing vertex if present
  if (nrow(crd) > 1 && isTRUE(all(crd[1, ] == crd[nrow(crd), ]))) {
    crd = crd[-nrow(crd), , drop = FALSE]
  }
  
  n = nrow(crd)
  if (n < 2) next
  
  all_nodes[[i]] = crd
  all_segs[[i]]  = cbind(offset + 1: n, offset + c(2:n, 1))
  offset = offset + n
}

boundary_nodes = do.call(rbind, all_nodes)
segments       = do.call(rbind, all_segs)

grid  = st_make_grid(ire, n = c(50, 50))
cent  = st_centroid(grid)
inside = st_within(cent, ire, sparse = FALSE)[,1]
interior_nodes = st_coordinates(cent[inside])

nodes = rbind(boundary_nodes, interior_nodes)
nodes = nodes[!duplicated(nodes),]

mesh = create.mesh.2D(
  nodes    = nodes,
  segments = segments
)

mesh_sf = st_sf(
  geometry = st_sfc(
    lapply(1:nrow(mesh$triangles), function(i) {
      tri = mesh$nodes[mesh$triangles[i, ], ]
      st_polygon(list(rbind(tri, tri[1,])))
    }),
    crs = 2157
  )
)

tm_shape(map) + 
  tm_borders(lwd = 2, col = "black") +
  tm_shape(mesh_sf) +
  tm_borders(col = "blue", alpha = 0.5) +
  tm_layout(frame = FALSE)

n_triangles = nrow(mesh$triangles)
n_nodes     = nrow(mesh$nodes)
cat("Mesh has", n_triangles, "triangles and", n_nodes, "nodes.\n")

FEM_basis = fdaPDE::create.FEM.basis(mesh)

#get the data
data = read_csv("./Data files/Weather Data/weather_data_long_clean.csv")
data = read_csv("./Data files/Air Quality data/air_quality_data_complete.csv")

fit_fem = function(param, yr) {
  data_filtered = data |>
    filter(`Air Pollutant` == param & Year == yr) |>
    dplyr::select(-`Air Quality Station Name`, -Year, -`Air Pollutant`) |>
    na.omit()
  
  #view the points on the map
  points(data_filtered[,c("Longitude","Latitude")], col = "red")
  
  #make sf object
  data_sf = st_as_sf(data_filtered, 
                     coords = c("Longitude","Latitude"), 
                     crs=sf::st_crs(map)
                     )
  
  #filter out data points which lie outside the map boundaries
  data_sf = st_filter(data_sf, map)
  
  # Project data points to mesh CRS (2157)
  data_sf = st_transform(data_sf, 2157)

  # Extract coordinates for fitting
  data_coords = st_coordinates(data_sf)
  
  FEM_fit = fdaPDE::smooth.FEM(data_coords, 
                               data_sf$value, 
                               FEM_basis, 
                               lambda = 10^seq(-10,10),
                               DOF.evaluation = 'exact', 
                               lambda.selection.lossfunction = 'GCV'
                               )
  
  # Build grid in EPSG:2157 to match mesh
  map2157 = st_transform(map, 2157)
  
  grid = terra::rast(map2157, nrows = 400, ncols = 400)
  xy = terra::xyFromCell(grid, 1:ncell(grid))
  dp = st_as_sf(as.data.frame(xy), 
                coords = c("x", "y"), 
                crs = st_crs(map2157)
                )
  
  # keep only points inside Ireland
  dp = st_filter(dp, ire_main)
  pred_coords = st_coordinates(dp)
  z = fdaPDE::eval.FEM(FEM_fit$fit.FEM, pred_coords)
  print(z)
  
  y_obs = as.numeric(data_sf$value)
  y_hat = as.numeric(FEM_fit$solution$z_hat)
  
  resid = y_obs - y_hat
  n     = length(resid)
  
  RMSE  = sqrt(mean(resid^2, na.rm = TRUE))
  MAE   = mean(abs(resid), na.rm = TRUE)
  R2    = 1 - sum(resid^2, na.rm = TRUE) / 
    sum((y_obs - mean(y_obs))^2, na.rm = TRUE)
  GCV   = FEM_fit$optimization$GCV
  
  out = list(RMSE = RMSE, MAE = MAE, R2 = R2, GCV = GCV)
  print(out)
  
  r_2157 = terra::rast(cbind(pred_coords, z), type = "xyz")
  terra::crs(r_2157) = st_crs(map_m)$wkt
  r_ll = terra::project(r_2157, "EPSG:4326")
  r_ll = terra::mask(r_ll, terra::vect(map))
  r_ll = terra::trim(r_ll)
  
  path = paste0(
    "Output Files/rastors/fdapde/", param, "_", yr, "_rastor.rds"
  )
  saveRDS(r_ll, path)
  
  r_ll
}

r_pm25 = fit_fem("PM2.5", 2022)
r_pm10 = fit_fem("PM10",  2022)

rng  = range(c(values(r_pm25), values(r_pm10)), na.rm=TRUE)
brks = pretty(rng, n=7)

# r_pm25 = readRDS("Output Files/rastors/gam/PM2.5_2022_rastor.rds")
# r_pm10 = readRDS("Output Files/rastors/gam/PM10_2022_rastor.rds")

tmap_mode("plot")

m_pm25 = tm_shape(r_pm25) +
  tm_raster(style = "cont", 
            at = brks, 
            palette = "viridis", 
            title = "PM2.5"
            ) +
  tm_layout(legend.outside = TRUE, frame = FALSE)

m_pm10 = tm_shape(r_pm10) +
  tm_raster(style = "cont", 
            at = brks, 
            palette = "viridis", 
            title = "PM10"
            ) +
  tm_layout(legend.outside = TRUE, frame = FALSE)

tmap_arrange(m_pm25, m_pm10, ncol = 2)