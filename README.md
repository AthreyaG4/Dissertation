# Spatial Epidemiology & Air Quality Analysis

A comprehensive R-based project for analyzing the relationship between air quality and health outcomes (lung cancer) in Ireland using advanced spatial statistical methods.

## Project Overview

This dissertation project develops and compares spatial smoothing techniques to model air quality pollution levels and their association with health outcomes across Ireland. The project integrates epidemiological data, environmental monitoring, and meteorological information to create predictive models and interactive visualizations.

## Key Features

- **Spatial Modeling**: Implementation of Generalized Additive Models (GAM) and Functional Data Analysis on Partial Differential Equations (FDAPDE) for spatial smoothing
- **Air Quality Analysis**: Processing and modeling of multiple air pollutants (PM10, PM2.5, and others)
- **Health Data Integration**: Lung cancer incidence data mapped against pollution levels at the county level
- **Interactive Visualization**: Shiny app for exploring actual data and model predictions
- **Automated Data Scraping**: Web scraping of Met Éireann historical weather data
- **Data Processing Pipeline**: Comprehensive Quarto documents for data cleaning, exploratory analysis, and integration

## Project Structure

```
├── app.R                          # Main Shiny application for visualization
├── app_demo.R                     # Demo version of the Shiny app
├── GAM.R                          # Generalized Additive Model implementation
├── FDAPDE.R                       # FDAPDE spatial smoothing models
├── GAM_areal_covariate.R         # GAM with areal covariate integration
├── FDAPDE_areal_covariate.R      # FDAPDE with areal covariate integration
├── EDA.R                          # Exploratory data analysis
├── view_data.R                    # Data visualization utilities
├── Data Processing.qmd           # Quarto document for data processing workflow
├── Data files/
│   ├── Air Quality data/         # Processed air quality measurements
│   ├── Cancer Data/              # Lung cancer incidence data
│   ├── Mortality Rates/          # Mortality statistics
│   ├── Weather Data/             # Meteorological data from Met Éireann
│   └── Population/               # Population statistics
├── Counties/
│   ├── shapefiles/               # Spatial boundary data for Irish counties
│   └── *.gpkg                    # GeoPackage administrative boundary files
├── models/
│   ├── gam/                      # Saved GAM model outputs
│   └── fdapde/                   # Saved FDAPDE model outputs
├── Output Files/
│   ├── areal_outputs/            # Areal analysis results
│   └── rastors/                  # Raster outputs for model predictions
├── scraping/
│   ├── scraping.py               # Python script for Met Éireann data scraping
│   └── *.csv                     # Scraping logs and tracking
└── archive/                      # Previous versions and experimental analyses
```

## Data Sources

- **Air Quality**: Point-based air quality measurements from monitoring stations across Ireland
- **Lung Cancer**: County-level incidence data
- **Weather Data**: Automated weather station data from Met Éireann (Irish meteorological service)
- **Geographic Data**: OSi (Ordnance Survey Ireland) administrative boundary data
- **Population**: Central Statistics Office population data

## Methodology

### Spatial Modeling Approaches

1. **GAM (Generalized Additive Models)**
   - Thin-plate spline basis functions for spatial smoothing
   - Weighted estimation based on data variability
   - Automatic smoothing parameter selection via REML

2. **FDAPDE (Functional Data Analysis on PDE)**
   - PDE-based spatial regularization
   - Handles areal and point-referenced data
   - Incorporates geographic covariates

### Analysis Features

- Comparison of model performance across spatial methods
- Incorporation of areal covariates (e.g., demographic factors)
- Raster prediction outputs for visualization and further analysis
- Cross-validation and model diagnostics

## Main Scripts

| Script | Purpose |
|--------|---------|
| `app.R` | Interactive Shiny dashboard for exploring predictions |
| `GAM.R` | Fit GAM models and generate raster predictions |
| `FDAPDE.R` | Fit FDAPDE models with spatial regularization |
| `EDA.R` | Exploratory analysis of health and environmental data |
| `Data Processing.qmd` | Data cleaning and integration workflow |
| `scraping/scraping.py` | Automated weather data collection |

## Requirements

### R Packages
- `shiny`, `bslib` - Interactive web application
- `tidyverse`, `dplyr` - Data manipulation
- `sf`, `terra`, `raster` - Spatial data handling
- `mgcv` - GAM implementation
- `ggplot2` - Visualization
- `INLA` - Integrated Nested Laplace Approximation (alternative method)
- `leaflet` - Interactive mapping

### Python Packages
- `selenium` - Web scraping automation
- `pandas` - Data manipulation

## Usage

### Running the Interactive App
```r
source("app.R")
```

### Fitting Models
```r
# Fit GAM models
source("GAM.R")

# Fit FDAPDE models
source("FDAPDE.R")
```

### Data Processing
Open and run the Quarto document:
```
Data Processing.qmd
```

### Scraping Weather Data
```bash
python scraping/scraping.py
```

## Output Files

- **Models**: Saved model objects in `models/gam/` and `models/fdapde/`
- **Predictions**: Raster files in `Output Files/rastors/` for visualization
- **Areal Analysis**: Results in `Output Files/areal_outputs/`

## Notes

- All spatial analyses are conducted at the county level for Ireland
- Model outputs are saved as `.rds` files for reproducibility
- Raster predictions enable continuous spatial prediction surfaces
- The Shiny app allows interactive exploration of both actual data and model predictions

## License & Attribution

This is a dissertation project. References to data sources and methodologies should be included in academic work utilizing this codebase.
