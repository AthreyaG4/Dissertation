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
library(fmesher)

#read data
data = read_csv("./Weather Data/weather_data_long_clean.csv")
df = data |>
  rename(station.name = AreaOfResidence) |>
  mutate(station.id = as.numeric(as.factor(station.name)))
parameters = unique(data$parameter)
parameters = parameters[-4]

#get the map
map <- ne_countries(type = "countries",
                    country = "Ireland",
                    scale = "medium", returnclass = "sf")
#make boundary
map <- st_transform(map, crs = 4326)
map_sp <- as_Spatial(map)
bnd <- inla.sp2segment(map_sp)

#prediction grid
grid <- terra::rast(map, nrows = 200, ncols = 200)
# coordinates of all cells
xy <- terra::xyFromCell(grid, 1:ncell(grid))

dp <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
               crs = st_crs(map))

# indices points within the map
indicespointswithin <- which(st_intersects(dp, map,
                                           sparse = FALSE))

# points within the map
dp <- st_filter(dp, map)

#coordinates
pred_coords = st_coordinates(dp)

for (param in parameters) {
  #filter and add station id
  df_filtered = df |>
    filter(parameter == param)
  
  #number of years in the dataset
  n_year = length(unique(df_filtered$year))
  expanded_coords <- pred_coords[rep(1:nrow(pred_coords), times=n_year), ]
  
  # Find station names or IDs present in every year
  stations_in_all_years <- df_filtered %>%
    distinct(station_name = station.name, year) %>%   # or use station_id if you have it
    group_by(station_name) %>%
    summarise(n_years = n(), .groups = "drop") %>%
    filter(n_years == n_year) %>%
    pull(station_name)
  n_stations = length(stations_in_all_years)
  
  # Filter the original dataframe to keep only these stations
  df_filtered <- df_filtered %>%
    filter(station.name %in% stations_in_all_years) |>
    dplyr::select(-parameter, -station.name)
  
  #make sf object
  df_sf = st_as_sf(df_filtered, coords = c("Longitude", "Latitude"))
  st_crs(df_sf) <- "EPSG:4326" 
  
  #create mesh
  mesh = fmesher::fm_mesh_2d(loc = unique(st_coordinates(df_sf)), boundary = bnd, max.edge=c(40,100), min.angle=c(21,21), cutoff = 0.1)
  #plot mesh and the data points
  plot(mesh)
  points(st_coordinates(df_sf), col = "red")
  
  #number of nodes in the mesh
  mesh$n
  
  spde = inla.spde2.matern(mesh=mesh)
  
  df_filtered$time = rep(1:n_year, each = n_stations)
  
  #building estimation stack
  A.est = inla.spde.make.A(mesh,
                           loc=as.matrix(df_filtered[,c("Longitude","Latitude")]),
                           group=df_filtered$time,
                           n.group=n_year)
  
  field.indices = inla.spde.make.index("field",n.spde=spde$n.spde,n.group=n_year)
  
  stack.est <- inla.stack(tag = "est",
                          data = list(value = df_filtered$value), A = list(1, A.est),
                          effects = list(data.frame(Intercept = rep(1, nrow(A.est))),
                                         s = field.indices))
  
  #prediction for single year
  #A.pred = inla.spde.make.A(mesh, loc=as.matrix(pred_coords),
  #                          group=1, n.group=n_year)
  
  #stack.pred <- inla.stack(tag = "pred",
  #                         data = list(value = NA), A = list(1, A.pred),
  #                         effects = list(data.frame(Intercept = rep(1, nrow(A.pred))),
  #                                        s = field.indices))
  
  #building prediction stack
  A.pred = inla.spde.make.A(mesh, loc=as.matrix(expanded_coords),
                            group=rep(1:n_year, each=nrow(pred_coords)), n.group=n_year)
  
  stack.pred <- inla.stack(tag = "pred",
                           data = list(value = NA), A = list(1, A.pred),
                           effects = list(data.frame(Intercept = rep(1, nrow(A.pred))),
                                          s = field.indices))
  
  stack = inla.stack(stack.est, stack.pred)
  
  formula <- (value ~ -1 + Intercept + f(field, model=spde,
                                         group=field.group, control.group=list(model="ar1")))
  mod.mode = inla(formula,
                  data=inla.stack.data(stack.est, spde=spde),
                  family="gaussian",
                  control.predictor=list(A=inla.stack.A(stack.est), compute=FALSE))
  
  mod = inla(formula,
             data=inla.stack.data(stack, spde=spde),
             family="gaussian",
             control.predictor=list(A=inla.stack.A(stack), compute=TRUE),
             control.mode=list(theta=mod.mode$mode$theta, restart=FALSE))
  mod$summary.fixed
  
  index <- inla.stack.index(stack = stack, tag = "pred")$data
  pred_mean <- mod$summary.fitted.values[index, "mean"]
  
  
  output = cbind(expanded_coords, output = pred_mean)
  r = terra::rast(output, type = "xyz")
  crs(r) = "EPSG:4326"
  if(!dir.exists("rastors")) dir.create("rastors")
  if(!dir.exists("rastors/inla")) dir.create("rastors/inla")
  for (yr in 1:n_year){
    upper_bound = (yr - 1) * nrow(pred_coords) + nrow(pred_coords)
    lower_bound = (yr - 1) * nrow(pred_coords) + 1
    y = yr + 2000
    r <- terra::rast(output[lower_bound:upper_bound,], type = "xyz")
    crs(r) <- "EPSG:4326"
    saveRDS(r, paste0("rastors/inla/",param,"_",y,"_rastor.rds"))
  }
}


parameter = "Precipitation Amount"
year = "2005"

r = readRDS(paste0("rastors/inla/",parameter,"_",year,"_rastor.rds"))
  
rb <- raster::brick(r)
contours <- rasterToContour(rb, levels = pretty(range(values(rb), na.rm = TRUE), n = 10))
contours_sf <- st_as_sf(contours)
pal <- colorNumeric("viridis", values(r),
                    na.color = "transparent")
leaflet() %>% addTiles() %>%
  addRasterImage(rb, colors = pal, opacity = 0.8) %>%
  #addPolylines(data = contours_sf, color = "black", weight = 1, opacity = 0.7) %>%
  addLegend(pal = pal, values = values(r), title = "input$parameter")











