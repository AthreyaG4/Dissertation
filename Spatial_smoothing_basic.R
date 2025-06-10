library(INLA)
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

# inla mesh

inla_mesh <- fmesher::fm_mesh_2d_inla(boundary = bnd, max.edge = 0.1)
#inla_mesh <- fmesher::fm_mesh_2d(
#  loc = fmesher::fm_hexagon_lattice(bnd, edge_len = 1),
#  boundary = bnd,
#  max.edge = 0.05
#)

plot(inla_mesh) # adjust max.edge to control number of nodes

xyz = cbind(xy,z)

matern <- INLA::inla.spde2.matern(inla_mesh)
xyz_inla <- sf::st_as_sf(xyz, coords = c('x', 'y')) #has to be sf object
cmp <- z ~ field(geometry,  model = matern) + Intercept(1)
INLA_fit <- inlabru::bru(cmp, xyz_inla, family = "gaussian")
INLA_zfit <- predict(INLA_fit, 
                     fmesher::fm_pixels(inla_mesh, dims = c(50, 50)), 
                     ~ field + Intercept)
plotly::plot_ly(x = 1:50, y = 1:50, z = matrix(INLA_zfit$mean, 50)) |> plotly::add_surface()
Metrics::rmse(INLA_zfit$mean, z)



