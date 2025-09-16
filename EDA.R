library(tidyverse)
library(ggplot2)
library(qqplotr)

#-----------Lung Cancer-----------------------

lung_cancer <- read_delim("Data files/Cancer Data/Cancers.csv", delim = "\t", locale = locale(encoding = "UTF-16LE"))

lung_cancer = lung_cancer |> filter(Sex == "Both") |> dplyr::select(Region, Year, `Case Numbers`)

lung_cancer$Region = toupper(lung_cancer$Region)

lung_cancer = lung_cancer |> filter(! (Region %in% c("DUBLIN NORTH", "DUBLIN SOUTH", "HSE DUBIN SOUTH EAST", "HSE DUBLIN MID LEINSTER", "HSE DUBLIN NORTH EAST", "HSE SOUTH WEST", "HSE WEST", "HSE WEST NORTH WEST", "IRELAND", "TIPPERARY NORTH", "TIPPERARY SOUTH")))

lung_cancer$Region[lung_cancer$Region == "DUBLIN (TOTAL)"] = "DUBLIN"

lung_cancer$Region[lung_cancer$Region == "TIPPERARY (TOTAL)"] = "TIPPERARY"

lung_cancer$`Case Numbers` = as.numeric(lung_cancer$`Case Numbers`)

lung_cancer_filtered = lung_cancer |> filter(Region != "DUBLIN" & Region != "CORK")

ggplot(lung_cancer, aes(x = Year, y = `Case Numbers`, color = Region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.8) +
  facet_wrap(~Region, scales = "free_y") +
  labs(
    x = "Year",
    y = "Number of cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",       # removes legend
    axis.text.x = element_text(angle = 45, hjust = 1) # rotates x labels if crowded
  ) +
  scale_x_continuous(breaks = seq(2000, 2022, by = 10))  # ticks every 5 years

lung_cancer_cummulative = lung_cancer |>
  group_by(Region) |>
  summarise(
    Region = first(Region),
    `Case Numbers` = sum(`Case Numbers`)
    )

ggplot(lung_cancer_cummulative, aes(x = Region, y = `Case Numbers`, fill = Region)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(`Case Numbers`, 2)), 
            vjust = -0.5, size = 3) +
  labs(
    title = "",
    x = "",
    y = "Cases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",       # removes legend
    axis.text.x = element_text(angle = 45, hjust = 1) # rotates x labels if crowded
  )

#------------------Pollution Data---------------------------------
data_pm25 = read_csv("Data files/Air Quality data/PM2.5.csv")
data_pm25 = data_pm25 |> filter(Year == 2022) |> dplyr::select(`Air Pollution Level`, Longitude, Latitude, Altitude, `City Population`)

data_pm10 = read_csv("Data files/Air Quality data/PM10.csv")
data_pm10 = data_pm10 |> filter(Year == 2022) |> dplyr::select(`Air Pollution Level`, Longitude, Latitude, Altitude, `City Population`)

# Q-Q plot
qqnorm(data_pm10$`Air Pollution Level`,
       main = "Q-Q Plot of PM2.5 Levels (2022)")
qqline(data_pm10$`Air Pollution Level`, col = "red", lwd = 2)


df1 <- data_pm25 %>%
  rename(value = `Air Pollution Level`) %>%
  mutate(value = as.numeric(value), pollutant = "PM2.5") %>%
  tidyr::drop_na(value)

df2 <- data_pm10 %>%
  rename(value = `Air Pollution Level`) %>%
  mutate(value = as.numeric(value), pollutant = "PM10") %>%
  tidyr::drop_na(value)

df <- rbind(df1, df2)

gg <- ggplot(data = df, mapping = aes(sample = pm25)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
gg

gg <- ggplot(data = df, mapping = aes(sample = value, color = pollutant, fill = pollutant)) +
  stat_qq_band(alpha=0.5) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~ pollutant) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
gg












