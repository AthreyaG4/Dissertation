library(ggplot2)
library(tidyverse)
library(sf)
library(mgcv)
library(tmap)
library(gratia)
library(mapview)
library(spdep)
library(scales)
library(grid)

pm25_model = read_rds(file = "models/gam/PM2.5.rds")
pm10_model = read_rds(file = "models/gam/PM10.rds")

data_sf = st_read("Counties/shapefiles/lungcancer.shp")

if (is.na(st_crs(data_sf))) {
  st_crs(data_sf) = 4326
} else {
  data_sf = st_transform(data_sf, 4326)
}

data_sf = data_sf |> filter(Year == 2022)
data_sf = st_make_valid(data_sf)

centroids = st_centroid(data_sf)
centroids_df = data.frame(Longitude = st_coordinates(centroids)[,"X"],
                          Latitude = st_coordinates(centroids)[,"Y"])

pm25_preds = predict(pm25_model, newdata = centroids_df)
pm10_preds = predict(pm10_model, newdata = centroids_df)

data_sf$pm25_preds = pm25_preds
data_sf$pm10_preds = pm10_preds

data_sf = st_transform(data_sf, 2157)

nb = poly2nb(data_sf)
data_sf$region_id = factor(seq_len(nrow(data_sf)))
names(nb) = levels(data_sf$region_id)

mrf_model = gam(
  CsNmbrs ~ s(region_id, bs = "mrf", xt = list(nb = nb), k = 10),
  data = data_sf, family = poisson(), method = "REML"
  )

mrf_cov_model_pm25 = gam(
  CsNmbrs ~ 1 +
    s(region_id, bs = "mrf", xt = list(nb = nb), k = 10) + 
    pm25_preds,
  data = data_sf,
  family = poisson(), 
  method = "REML")

mrf_cov_model_pm10 = gam(
  CsNmbrs ~ 1 +
    s(region_id, bs = "mrf", xt = list(nb = nb), k = 10) + 
    pm10_preds,
  data = data_sf,
  family = poisson(), 
  method = "REML")

# compare the two with anova
anova(mrf_model, mrf_cov_model_pm25, test = "Chisq")

summary(mrf_model)
summary(mrf_cov_model_pm25)
summary(mrf_cov_model_pm10)

appraise(mrf_model)
appraise(mrf_cov_model_pm25)
appraise(mrf_cov_model_pm10)

y_obs  = data_sf$CsNmbrs
y_pred = predict(mrf_cov_model_pm10, type = "response")  # fitted counts

rmse = sqrt(mean((y_obs - y_pred)^2))
rmse

# predictions per area
data_sf$preds = as.numeric(
  predict(mrf_cov_model_pm25, type = "response")
  ) 

#calculate residuals
data_sf$residuals = data_sf$CsNmbrs - data_sf$preds
d = st_transform(data_sf, 4326)

d = d |> dplyr::select(geometry, preds, residuals)
st_write(d, "Output Files/areal_outputs/gam/lungcancer_predictions_covariate_2022.shp", delete_layer = TRUE)

#plot predictions
rng   = range(d$preds, na.rm = TRUE)
rng[1] = floor(rng[1])
rng[2] = ceiling(rng[2])
ticks = expm1(seq(log1p(rng[1]), log1p(rng[2]), length.out = 6))
ticks = round(ticks)

gg1 = ggplot(d) +
  geom_sf(aes(fill = preds), color = "black", 
          linewidth = 0.5, alpha = 0.75) +
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

gg2 = ggplot(d) +
  geom_sf(aes(fill = residuals), color = "black", 
          linewidth = 0.5, alpha = 0.75) +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = 1,
    limits    = c(-80,80),
    labels    = label_comma(),
    breaks    = ticks,
    na.value  = NA,
    guide     = guide_colorbar(
      title.position = "top",
      barheight = unit(90, "pt"),
      ticks = TRUE
    )
  ) +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

gg2
