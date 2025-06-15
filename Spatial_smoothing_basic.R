library(INLA)
library(fmesher)

# simple data example
grid = seq(0, 1, len = 50)
xy <- expand.grid(x = grid, y = grid) # coordinates
z <- exp(-(xy$x-xy$y)^2/(0.125)) # observations (RBF example)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(z, 50)) |> plotly::add_surface()
# add noise 
z_err <- z + rnorm(50^2, 0, 1/16);
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(z_err, 50)) |> plotly::add_surface()
# add holes 
holes <- sample.int(50^2, 200)
z_holes <- z_err
z_holes[holes] <- NA
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(z_holes, 50)) |> plotly::add_surface()


# NOTE: In non-regular surfaces you can use sf package to extract boundaries.

# FDAPDE
mesh_grid <- seq(0, 1, len = 20) #creates mesh with fewer points that observations locations
tri_mesh <- fdaPDE::create.mesh.2D(expand.grid(mesh_grid, mesh_grid), order = 1)
plot(tri_mesh) # can use more mesh points for better fits

FEM_basis <- fdaPDE::create.FEM.basis(tri_mesh)
FEM_fit <- fdaPDE::smooth.FEM(as.matrix(xy[-holes, ]), z_err[-holes], FEM_basis, lambda = 10^(-6:6),
                   DOF.evaluation = 'exact', 
                   lambda.selection.lossfunction = 'GCV')
FEM_zfit <- fdaPDE::eval.FEM(FEM_fit$fit.FEM, xy)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(FEM_zfit, 50)) |> plotly::add_surface()
Metrics::rmse(FEM_zfit, z)


# GAM
xyz = cbind(xy[-holes,], z = z_err[-holes])
GAM_fit <- mgcv::gam(z ~ s(x,y,bs = 'tp'), data = xyz,
                     family = 'gaussian', method = 'REML')
GAM_zfit <- mgcv::predict.gam(GAM_fit, xy)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(GAM_zfit, 50)) |> plotly::add_surface()
Metrics::rmse(GAM_zfit, z)

# INLA (using inlabru)

# boundary
bnd <- inlabru::spoly(
  data.frame(
    x = c(0,1,1,0),
    y = c(0,0,1,1)
  ), 
  format = 'sf'
)
inla_mesh <- fm_mesh_2d_inla(boundary = bnd, max.edge = 0.2)
# inla mesh
#inla_mesh <- fmesher::fm_mesh_2d(
#  loc = fmesher::fm_hexagon_lattice(bnd, edge_len = 1),
#  boundary = bnd,
#  max.edge = 0.05

#)

plot(inla_mesh) # adjust max.edge to control number of nodes

matern <- INLA::inla.spde2.matern(inla_mesh)
xyz_inla <- sf::st_as_sf(xyz, coords = c('x', 'y')) #has to be sf object
cmp <- z ~ field(geometry,  model = matern) + Intercept(1)
INLA_fit <- inlabru::bru(cmp, xyz_inla, family = "gaussian")
INLA_zfit <- predict(INLA_fit, 
                     fmesher::fm_pixels(inla_mesh, dims = c(50, 50)), 
                     ~ field + Intercept)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(INLA_zfit$mean, 50)) |> plotly::add_surface()
Metrics::rmse(INLA_zfit$mean, z)








### Easy plot

quick_leaflet_plot <- function(sf_data, point_data, palette = "viridis"){
  pl <- leaflet::leaflet(sf_data) |> 
    leaflet::addTiles() |> 
    leaflet::addPolygons(data=sf_data,
                         color = "#444444",
                         weight = 1,
                         smoothFactor = 0.5,
                         opacity = 1.0,
                         fillOpacity = 0.5,
                         highlightOptions = leaflet::highlightOptions(color = "white", weight = 2,bringToFront = TRUE)
    ) 
  if (ncol(point_data) > 2){
    pal <- leaflet::colorNumeric(palette = palette, domain = point_data[,3])
    pl |> leaflet::addCircles(color = ~pal(point_data[,3]), opacity = 1,
                              lng=point_data[,1], lat=point_data[,2])
  }else{
    pl |> leaflet::addCircles(color = '#ff000088', opacity = 1,
                              lng=point_data[,1], lat=point_data[,2])
  }
  
}


#### Ireland 
IRL <- rnaturalearth::ne_states(geounit = 'ireland')
plot(sf::st_geometry(IRL))
IRL_sf <- sf::st_as_sf(IRL)

## data simulation
grid = seq(0, 1, len = 50)
xy <- expand.grid(x = grid, y = grid) # coordinates
z <- exp(-(xy$x-xy$y)^2/(0.125)) # observations (RBF example)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(z, 50)) |> plotly::add_surface()

IRL_bb <- sf::st_bbox(IRL_sf) # This just gives the values at the corners of a box over the map see plot below

# transform the xy values to range of lon + lat of ireland
IRL_xyz <- as.data.frame(cbind('Lon' = xy$x * (IRL_bb$xmax - IRL_bb$xmin) + IRL_bb$xmin , 
                 'Lat' = xy$y * (IRL_bb$ymax - IRL_bb$ymin) + IRL_bb$ymin, 
                 'obs' = z))
IRL_xyz_sf <- sf::st_as_sf(IRL_xyz, coords = c("Lon", "Lat"), crs=sf::st_crs(IRL_sf))
# get all the simulated points within the map
in_map <- sf::st_intersects(IRL_xyz_sf, IRL_sf)
in_map <- sapply(1:50^2, function(i) ifelse(length(in_map[[i]]) == 0, NA, in_map[[i]]))
IRL_xyz <- IRL_xyz[!is.na(in_map), ]
quick_leaflet_plot(IRL_sf, IRL_xyz)

# Add noise and holes
# add noise 
IRL_xyz$obs_err <- IRL_xyz$obs + rnorm(nrow(IRL_xyz), 0, 1/8)
quick_leaflet_plot(IRL_sf, IRL_xyz[, c(1,2,4)])
# add holes 
holes <- sample.int(nrow(IRL_xyz), 500)
z_holes <- IRL_xyz$obs_err
z_holes[holes] <- NA
IRL_xyz$obs_err_holes <- z_holes
quick_leaflet_plot(IRL_sf, IRL_xyz[, c(1,2,5)])


# Fit with FDAPDE
sq100 <- sf::st_make_grid(IRL_sf, n = c(50, 25)) # rect grid on bbox
plot(sq100)
plot(sf::st_geometry(IRL_sf), add = TRUE)
plot(sq100[sf::st_union(IRL_sf)], col = '#ff000088', add = TRUE)
tri_pts <- sf::st_coordinates(sf::st_centroid(sq100[sf::st_union(IRL)])) # Triangulation points
outer_pts <- sf::st_coordinates(sf::st_union(sq100[sf::st_union(IRL_sf)])) # boundary points
plot(outer_pts, type = 'l')
tri_pts = rbind(outer_pts[, 1:2], tri_pts) # combine boundary points and inner points
tri_pts = tri_pts[!duplicated(tri_pts),]
tri_seg <- cbind(1:(nrow(outer_pts)-1), c(2:(nrow(outer_pts)-1),1))
tri_mesh <- fdaPDE::create.mesh.2D(tri_pts, order = 1, segments = tri_seg)
plot(tri_mesh)
FEM_basis <- fdaPDE::create.FEM.basis(tri_mesh)
FEM_fit <- fdaPDE::smooth.FEM(as.matrix(IRL_xyz[-holes, 1:2]), IRL_xyz$obs_err_holes[-holes], FEM_basis, lambda = 10^(-6:6),
                              DOF.evaluation = 'exact', 
                              lambda.selection.lossfunction = 'GCV')
FEM_zfit <- fdaPDE::eval.FEM(FEM_fit$fit.FEM, IRL_xyz[, 1:2])
Metrics::rmse(FEM_zfit, IRL_xyz$obs)
IRL_xyz$FEMfit <- FEM_zfit
quick_leaflet_plot(IRL_sf, IRL_xyz[, c(1,2,6)])

# Fit INLA
bnd <- sf::st_boundary(sf::st_union(sq100[sf::st_union(IRL_sf)]))
inla_mesh <- fmesher::fm_mesh_2d(
  loc = fmesher::fm_hexagon_lattice(bnd, edge_len = 1),
  boundary = sf::st_boundary(bnd),
  max.edge = 0.5
  
)
inla_mesh <- fm_mesh_2d_inla(boundary = bnd, max.edge = 0.2)
plot(inla_mesh)


#### FRENCH
#### Ireland 
FRA <- rnaturalearth::ne_states(geounit = 'france')
plot(sf::st_geometry(FRA))
# For simplicity take out Corse island
FRA <- subset(FRA, region != 'Corse')
plot(sf::st_geometry(FRA))
FRA_sf <- sf::st_as_sf(FRA)

## data simulation
grid = seq(0, 1, len = 50)
xy <- expand.grid(x = grid, y = grid) # coordinates
z <- exp(-(xy$x-xy$y)^2/(0.125)) # observations (RBF example)
FRA_bb <- sf::st_bbox(FRA_sf) # This just gives the values at the corners of a box over the map see plot below

# transform the xy values to range of lon + lat of france
FRA_xyz <- as.data.frame(cbind('Lon' = xy$x * (FRA_bb$xmax - FRA_bb$xmin) + FRA_bb$xmin , 
                               'Lat' = xy$y * (FRA_bb$ymax - FRA_bb$ymin) + FRA_bb$ymin, 
                               'obs' = z))
FRA_xyz_sf <- sf::st_as_sf(FRA_xyz, coords = c("Lon", "Lat"), crs=sf::st_crs(FRA_sf))
# get all the simulated points within the map
in_map <- sf::st_intersects(FRA_xyz_sf, FRA_sf)
in_map <- sapply(1:50^2, function(i) ifelse(length(in_map[[i]]) == 0, NA, in_map[[i]]))
FRA_xyz <- FRA_xyz[!is.na(in_map), ]
quick_leaflet_plot(FRA_sf, FRA_xyz)

# Add noise and holes
# add noise 
FRA_xyz$obs_err <- FRA_xyz$obs + rnorm(nrow(FRA_xyz), 0, 1/8)
quick_leaflet_plot(FRA_sf, FRA_xyz[, c(1,2,4)])
# add holes 
holes <- sample.int(nrow(FRA_xyz), 500)
z_holes <- FRA_xyz$obs_err
z_holes[holes] <- NA
FRA_xyz$obs_err_holes <- z_holes
quick_leaflet_plot(FRA_sf, FRA_xyz[, c(1,2,5)])


# Fit with FDAPDE
sq100 <- sf::st_make_grid(FRA_sf, n = c(50, 25)) # rect grid on bbox
plot(sq100)
plot(sf::st_geometry(FRA_sf), add = TRUE)
plot(sq100[sf::st_union(FRA_sf)], col = '#ff000088', add = TRUE)
tri_pts <- sf::st_coordinates(sf::st_centroid(sq100[sf::st_union(FRA)])) # Triangulation points
outer_pts <- sf::st_coordinates(sf::st_union(sq100[sf::st_union(FRA_sf)])) # boundary points
plot(outer_pts, type = 'l')
tri_pts = rbind(outer_pts[, 1:2], tri_pts) # combine boundary points and inner points
tri_pts = tri_pts[!duplicated(tri_pts),]
tri_seg <- cbind(1:(nrow(outer_pts)-1), c(2:(nrow(outer_pts)-1),1))
tri_mesh <- fdaPDE::create.mesh.2D(tri_pts, order = 1, segments = tri_seg)
plot(tri_mesh)
FEM_basis <- fdaPDE::create.FEM.basis(tri_mesh)
FEM_fit <- fdaPDE::smooth.FEM(as.matrix(FRA_xyz[-holes, 1:2]), FRA_xyz$obs_err_holes[-holes], FEM_basis, lambda = 10^(-6:6),
                              DOF.evaluation = 'exact', 
                              lambda.selection.lossfunction = 'GCV')
FEM_zfit <- fdaPDE::eval.FEM(FEM_fit$fit.FEM, FRA_xyz[, 1:2])
Metrics::rmse(FEM_zfit, FRA_xyz$obs)
FRA_xyz$FEMfit <- FEM_zfit
quick_leaflet_plot(FRA_sf, FRA_xyz[, c(1,2,6)])

# Fit INLA
bnd <- sf::st_boundary(sf::st_union(sq100[sf::st_union(FRA_sf)]))
inla_mesh <- fmesher::fm_mesh_2d(
  loc = fmesher::fm_hexagon_lattice(bnd, edge_len = 1),
  boundary = sf::st_boundary(bnd),
  max.edge = 0.5
  
)
plot(inla_mesh)
