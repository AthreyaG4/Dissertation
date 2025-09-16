library(tidyverse)
library(mapview)
library(sf)
library(viridis)
library(leafsync)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(tmap)

data <- readr::read_csv("Data files/Air Quality data/air_quality_data_complete.csv")
data_pm10 = data |> filter(`Air Pollutant` == "PM10", Year == 2022) |> drop_na()
data_pm25 = data |> filter(`Air Pollutant` == "PM2.5", Year == 2022) |> drop_na()

ireland <- ne_countries(country = "Ireland", scale = "large", returnclass = "sf")

# PM10 plot
p1 <- ggplot() +
  geom_sf(data = ireland, fill = "grey95", color = "black") +
  geom_point(data = data_pm10, aes(x = Longitude, y = Latitude, color = value), size=2) +
  scale_color_viridis_c(name = "PM10") +
  coord_sf(xlim=c(-11, -5.5), ylim=c(51.3, 55.5)) +
  theme_void()

# PM2.5 plot
p2 <- ggplot() +
  geom_sf(data = ireland, fill = "grey95", color = "black") +
  geom_point(data = data_pm25, aes(x = Longitude, y = Latitude, color = value), size=2) +
  scale_color_viridis_c(name = "PM2.5") +
  coord_sf(xlim=c(-11, -5.5), ylim=c(51.3, 55.5)) +
  theme_void()

# side by side
p2 + p1

data_sf <- st_read("Counties/shapefiles/lungcancer.shp")
if (is.na(st_crs(data_sf))) st_crs(data_sf) <- 2157 else data_sf <- st_transform(data_sf, 2157)
d = data_sf |> filter(Year == 2022)

rng   <- range(d$CsNmbrs, na.rm = TRUE)
rng[1] <- floor(rng[1])
rng[2] <- round(rng[2])
ticks <- expm1(seq(log1p(rng[1]), log1p(rng[2]), length.out = 6))
ticks <- round(ticks)  # clean integers

gg <- ggplot(d) +
  geom_sf(aes(fill = CsNmbrs), color = "black", linewidth = 0.5, alpha = 0.75) +
  scale_fill_distiller(
    palette   = "YlOrRd",      # ColorBrewer
    direction = 1,             # 1 = red high; -1 = yellow high
    trans     = log1p_trans(), # continuous log colours
    limits    = rng,
    breaks    = ticks,
    labels    = label_comma(), # 1,000 separators
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
