#' ---
#' title: "Automatic interpolation of spatial points data with Kriging"
#' #' abstract: "This script demonstrates how to use the autoKrigMap function for automated spatial interpolation of sensor point data (e.g. Geophilus-Multi-Sensor-Platform)."
#' author: "Ingmar Schröter, Sandra Post"
#' date: "`r format(Sys.Date(), '%Y-%m-%d')`"
#' warning: false
#' message: false
#' lang: eng
#' format: 
#'   typst:
#'    toc: true
#'    number-sections: true
#' ---

#' # Load packages 
# LOAD PACKAGES ----
library(cowplot)
library(ggplot2)
library(sf)
library(gstat)
library(terra)
library(dplyr)
library(concaveman)

#' # Load helper functions 
# LOAD HELPER FUNCTIONS ----
source("R/createRast.R")
source("R/autofitVariogram3.R")
source("R/autoKrigMap.R")

#' # Input data 
# INPUT DATA ----

# Farm info 
farm_name <- "LPP"
field_name <- "1392"

## Field boundary ----
field_sf <- st_read("data/1392/PP_1392_Schlaggrenze_Neu.shp", quiet = T)

## Sensor point data ----
sensor_sf <- st_read("data/1392/LPP_1392_EPSG_4326_Geophilus_2017-09-27.shp", quiet = T) 
#sensor_sf <- st_transform(x = data_sf, crs = 25833)
#sensor_sp <- as(sensor_sf, "Spatial")

# Spatial interpolation (Kriging) ----
ok_data <- autoKrigMap(x = sensor_sf[c(1, 2, 8)], boundary = field_sf) 

#| label: fig-interpolated-maps
#| fig-cap: "Results of the interpolation process for all selected sensor point data."
plot(ok_data$maps)

#| label: fig-variograms
#| fig-cap: "Experimental and fitted theoretical variograms for all selected sensor point data."
#| fig-width: 5.83
#| fig-height: 8.27 
ok_data$varioPlots

# Create base raster for kriging  ----
#base_r <- createRast(boundary = field_sf, res = c(10, 10), buffer = 10, value = -9999)
#base_df <- as.data.frame(base_r, xy=T)
#base_grid <- st_as_sf(x = base_df, coords = c("x", "y"), crs = 25833)

# Save Output as GeoTiff ---- 
r <- ok_data$maps
crs_info <- crs(r, describe = TRUE)
crs_info$name
crs_info$authority
crs_info$code
#crs(r)  <- "epsg:25832"

for (i in seq_len(nlyr(r))) {
  terra::writeRaster(r[[i]], paste0("output/", farm_name, "_", field_name, "_OK_",names(r)[i],".tif"), overwrite=TRUE)
}

# Generate report ----

# PDF
#quarto::quarto_render(input = "Vignette_autoKrigMap.R")

# HTML
#quarto::quarto_render(input = "Vignette_autoKrigMap.R", output_format = c("html"))



