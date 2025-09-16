library(sf)
library(spdep)
library(mgcv)
library(ggplot2)
library(scales)
library(grid)
library(tidyverse)
library(mapview)
library(gratia)

# read & ensure projected in meters
data_sf <- st_read("Counties/shapefiles/lungcancer.shp")
if (is.na(st_crs(data_sf))) st_crs(data_sf) <- 2157 else data_sf <- st_transform(data_sf, 2157)
data_sf = data_sf |> filter(Year == 2022)

data_sf <- st_make_valid(data_sf)

nb <- poly2nb(data_sf)
data_sf$region_id <- factor(seq_len(nrow(data_sf)))
names(nb) <- levels(data_sf$region_id)

mrf_model <- gam(CsNmbrs ~ s(region_id, bs = "mrf", xt = list(nb = nb), k = 10),
                 data = data_sf, family = poisson(), method = "REML")

summary(mrf_model)
appraise(mrf_model)

y_obs  <- data_sf$CsNmbrs
y_pred <- predict(mrf_model, type = "response")

rmse <- sqrt(mean((y_obs - y_pred)^2))
rmse

# predictions per area
data_sf$preds <- as.numeric(predict(mrf_model, type = "response"))
data_sf$residuals <- data_sf$CsNmbrs - data_sf$preds
d <- st_transform(data_sf, 4326) |>
  dplyr::select(geometry, preds, residuals)

st_write(d, "Output Files/areal_outputs/gam/lungcancer_predictions_2022.shp", delete_layer = TRUE)

rng   <- range(d$preds, na.rm = TRUE)
rng[1] <- floor(rng[1])
rng[2] <- ceiling(rng[2])
ticks <- expm1(seq(log1p(rng[1]), log1p(rng[2]), length.out = 6))
ticks <- round(ticks)

gg1 <- ggplot(d) +
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

gg1

rng <- range(d$residuals, na.rm = TRUE)
rng <- c(floor(rng[1]), ceiling(rng[2]))
ticks <- seq(rng[1], rng[2], length.out = 6)

gg2 <- ggplot(d) +
  geom_sf(aes(fill = residuals), color = "black", linewidth = 0.5, alpha = 0.75) +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = 1,
    limits    = c(-80,80),
    labels    = label_comma(),
    breaks    = ticks,
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

gg2
