library(sf)
library(dplyr)
library(ggplot2)
library(mapview)
library(fdaPDE)
library(rmapshaper)
library(tidyverse)
library(patchwork)
library(tmap)
library(scales)

data_sf = st_read("Counties/shapefiles/lungcancer.shp")

d <- data_sf |>
  filter(Year == 2022) |>
  st_make_valid()

d

u    <- st_union(d)
poly <- st_cast(u, "POLYGON")
mainland <- poly[which.max(st_area(poly))]

d_no_islands <- st_intersection(d, mainland) |>
  st_make_valid() |>
  st_collection_extract("POLYGON")

d_no_islands <- d_no_islands |>
  filter(st_area(geometry) > units::set_units(1, km^2))

d_simple <- ms_simplify(d_no_islands, keep = 0.001, keep_shapes = TRUE)

if (is.na(st_crs(d_simple)) || st_crs(d_simple)$epsg != 2157) {
  d_simple <- st_transform(d_simple, 2157)
}

touches_list <- st_touches(d_simple) 

cent_pts <- st_point_on_surface(d_simple)
cent_xy  <- st_coordinates(cent_pts)
y_cent   <- cent_xy[,2]
x_cent   <- cent_xy[,1]

county_col <- "Region"
start_idx <- which(grepl("\\bcork\\b", tolower(d_simple[[county_col]])))
stopifnot(length(start_idx) == 1)

n <- nrow(d_simple)
visited <- rep(FALSE, n)
visit_id <- rep(NA_integer_, n)

queue <- start_idx
rank <- 1L

while (length(queue) > 0) {
  i <- queue[1]
  queue <- queue[-1]
  if (visited[i]) next
  visited[i] <- TRUE
  visit_id[i] <- rank; rank <- rank + 1L
  
  nb <- touches_list[[i]]
  nb <- nb[!visited[nb]]
  
  if (length(nb) > 1) {
    nb <- nb[order(-y_cent[nb], x_cent[nb])]
  }
  
  queue <- c(queue, nb)
}

while (any(!visited)) {
  seed <- which(!visited)[which.max(y_cent[!visited])]
  queue <- seed
  while (length(queue) > 0) {
    i <- queue[1]
    queue <- queue[-1]
    if (visited[i]) next
    visited[i] <- TRUE
    visit_id[i] <- rank; rank <- rank + 1L
    nb <- touches_list[[i]]
    nb <- nb[!visited[nb]]
    if (length(nb) > 1) {
      nb <- nb[order(-y_cent[nb], x_cent[nb])]
    }
    queue <- c(queue, nb)
  }
}

d_simple$visit_id <- visit_id
d$visit_id <- visit_id

d_simple = d_simple |> arrange(visit_id)
d = d |> arrange(visit_id)

make_mesh <- function(poly) {
  county <- st_simplify(poly)
  county_dense <- st_segmentize(county, dfMaxLength = 5000)

  rings <- county_dense |>
    st_boundary() |>
    st_cast("MULTILINESTRING") |>
    st_cast("LINESTRING")

  all_nodes <- list()
  all_segs  <- list()
  offset <- 0
  
  for (i in seq_len(length(rings))) {
    crd <- st_coordinates(rings[i])[, 1:2, drop = FALSE]

    if (nrow(crd) > 1 && isTRUE(all(crd[1, ] == crd[nrow(crd), ]))) {
      crd <- crd[-nrow(crd), , drop = FALSE]
    }
    
    n <- nrow(crd)
    if (n < 2) next
    
    all_nodes[[i]] <- crd
    all_segs[[i]]  <- cbind(offset + 1: n, offset + c(2:n, 1))
    offset <- offset + n
  }
  
  boundary_nodes <- do.call(rbind, all_nodes)
  segments       <- do.call(rbind, all_segs)

  grid  <- st_make_grid(county, n = c(20,20))
  cent  <- st_centroid(grid)
  inside <- st_within(cent, county, sparse = FALSE)[,1]
  interior_nodes <- st_coordinates(cent[inside])

  nodes_all <- rbind(boundary_nodes, interior_nodes)

  nodes_round <- round(nodes_all, 6)
  
  nodes_unique <- unique(nodes_round)
  
  row_map <- match(
    data.frame(t(nodes_round)),
    data.frame(t(nodes_unique))
  )
  
  segments_unique <- matrix(row_map[as.vector(segments)], ncol = 2)
  
  valid_seg <- segments_unique[,1] != segments_unique[,2]
  segments_unique <- segments_unique[valid_seg,]
  
  nodes <- as.matrix(nodes_unique)
  
  nodes_sf <- st_as_sf(data.frame(x = nodes[,1], y = nodes[,2]),
                       coords = c("x","y"), crs = st_crs(county))
  
  mesh <- create.mesh.2D(
    nodes    = nodes,
    segments = segments_unique,
    order    = 1
  )
  mesh
}

all_meshes <- list()

for (i in seq_len(nrow(d_simple))) {
  county_name <- d_simple$Region[i]
  poly <- d_simple[i,]
  cat("Meshing:", county_name, "\n")
  mesh_sf <- make_mesh(poly)
  all_meshes[[i]] <- mesh_sf
}

merge_two_meshes <- function(meshA, meshB, tol = 6) {
  nA <- nrow(meshA$nodes)
  
  nodes_all <- rbind(meshA$nodes, meshB$nodes)
  
  tris_all <- rbind(meshA$triangles, meshB$triangles + nA)
  segs_all <- rbind(meshA$segments,  meshB$segments  + nA)
  
  nodes_round  <- round(nodes_all, tol)
  nodes_unique <- unique(nodes_round)
  
  row_map <- match(
    data.frame(t(nodes_round)),
    data.frame(t(nodes_unique))
  )
  
  tri_comb <- matrix(row_map[as.vector(tris_all)], ncol = 3)
  seg_comb <- matrix(row_map[as.vector(segs_all)], ncol = 2)
  
  tri_comb <- tri_comb[apply(tri_comb, 1, function(r) length(unique(r)) == 3), , drop=FALSE]
  seg_comb <- seg_comb[seg_comb[,1] != seg_comb[,2], , drop=FALSE]
  
  create.mesh.2D(
    nodes     = as.matrix(nodes_unique),
    triangles = tri_comb,
    segments  = seg_comb,
    order     = 1
  )
}

mesh_global <- all_meshes[[1]]
for (i in 2:length(all_meshes)) {
  cat("Merging county", i, "of", length(all_meshes), "\n")
  mesh_global <- merge_two_meshes(mesh_global, all_meshes[[i]])
}

FEMbasis = create.FEM.basis(mesh_global)

mesh_to_sf <- function(mesh, crs) {
  tris <- mesh$triangles
  nxy  <- mesh$nodes
  polys <- lapply(seq_len(nrow(tris)), function(i) {
    xy <- nxy[tris[i, ], , drop=FALSE]
    st_polygon(list(rbind(xy, xy[1, , drop=FALSE])))
  })
  st_as_sf(st_sfc(polys, crs = crs))
}

mesh_global_sf <- mesh_to_sf(mesh_global, st_crs(d_simple))

nrow(mesh_global$nodes)
nrow(mesh_global$triangles)

tmap_mode("plot")   # static
tm_shape(d_simple) +
  tm_borders(col = "black", lwd = 2) +
  tm_shape(mesh_global_sf) +
  tm_borders(col = "black", lwd = 0.2) +
  tm_layout(frame = FALSE)

tri_centroids <- st_centroid(mesh_global_sf)

cent_join <- st_join(tri_centroids, d_simple["Region"], left = TRUE)

county_names <- d_simple$Region
county_index <- match(cent_join$Region, county_names)

n_counties <- length(county_names)
n_tri      <- nrow(mesh_global_sf)

i <- county_index
j <- seq_len(n_tri)

valid <- !is.na(i)
i <- i[valid]; j <- j[valid]

incidence_matrix <- sparseMatrix(i = i, j = j, x = 1,
                                 dims = c(n_counties, n_tri),
                                 dimnames = list(county_names, paste0("T", j)))

solution <- smooth.FEM(observations = d_simple$CsNmbrs,
                       incidence_matrix = incidence_matrix, 
                       FEMbasis = FEMbasis, 
                       DOF.evaluation = "exact",
                       lambda = 10^seq(0,10),
                       lambda.selection.lossfunction = "GCV",
                       areal.data.avg = FALSE
)

solution$optimization$lambda_vector
solution$optimization$lambda_solution
solution$solution$z_hat
solution$solution$rmse
solution$optimization$GCV_vector
solution$optimization$dof[solution$optimization$lambda_position]

y_obs <- as.numeric(d$CsNmbrs)
y_hat <- as.numeric(solution$solution$z_hat) 

resid <- y_obs - y_hat
n     <- length(resid)

RMSE  <- sqrt(mean(resid^2, na.rm = TRUE))
R2    <- 1 - sum(resid^2, na.rm = TRUE) / sum( (y_obs - mean(y_obs))^2, na.rm = TRUE)
GCV   <- solution$optimization$GCV

out <- list(RMSE = RMSE, R2 = R2, GCV = GCV)
print(out)
print(solution$optimization$lambda_solution)

d_plot <- d %>%
  mutate(observed = CsNmbrs, 
         preds   = solution$solution$z_hat) |>
  st_transform(4326)

st_write(d_plot, "Output Files/areal_outputs/fdapde/lungcancer_predictions_2022.shp", delete_layer = TRUE)

rng   <- range(d_plot$preds, na.rm = TRUE)
rng[1] <- floor(rng[1])
rng[2] <- round(rng[2])
ticks <- expm1(seq(log1p(rng[1]), log1p(rng[2]), length.out = 6))
ticks <- round(ticks)

gg <- ggplot(d_plot) +
  geom_sf(aes(fill = preds), color = "black", linewidth = 0.5, alpha = 0.75) +
  scale_fill_distiller(
    palette   = "YlOrRd",
    direction = 1,
    trans     = log1p_trans(),
    limits    = rng,
    breaks    = ticks,
    labels    = label_comma(),
    na.value  = NA,
    name      = "",
    guide     = guide_colorbar(
      title.position = "top",
      barheight = unit(90, "pt"),
      ticks = TRUE
    )
  ) +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

gg
