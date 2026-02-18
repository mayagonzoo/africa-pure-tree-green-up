# ============================================================
# Author: Maya Gonzalez
# Date: July 2025
#
# Description:
# - Merges phenological metrics with environmental covariates
#   (soil, topography, rainfall)
# - Aggregates phenological metrics across years per tree
# - Performs stratified random spatial subsampling by country
#
# Inputs:
# - Pheno_rain_annual_100000_trees.csv
# - Full_pheno_params_cleaned.csv
# - soil_pca_0_60_data.csv
# - DEMData_100000_trees.csv
#
# Outputs:
# - Env_data_annual.csv
# - Final_8500_pts.csv
# - subset_dNDVI_pts.csv
#
# Notes:
# - Downsampling is necessary because subsequent spatial models are
#   computationally intensive and slow to converge at larger sample sizes
# - Stratified random sampling balances coverage across the latitudinal gradient
# ============================================================

library(dplyr)
library(sf)
library(terra)
library(purrr)
library(rnaturalearth)
library(RColorBrewer)

# ---- File paths ----
# Update these paths to match your local directory structure
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

#############################################################
### Data import
#############################################################
# Annual phenology and rainfall metrics
green_up <- read.csv(file.path(data_dir, "Pheno_rain_annual_100000_trees.csv"))
pheno_params <- read.csv(file.path(data_dir, "Full_pheno_params_cleaned.csv"))
# Environmental covariates
soil_pca <- read.csv(file.path(data_dir, "soil_pca_0_60_data.csv"))
# DEM-derived topography (retain relevant variables only)
slope_data <- read.csv(file.path(data_dir, "DEMData_100000_trees.csv")) %>% 
  dplyr::select(tree_ID, slope, elevation, curvature)

# Join phenological metrics with environmental covariates
data_annual = dat_GUD_SRS %>% 
  inner_join(gup_data, by = c("tree_ID", "rainyear", "GUD")) %>%
  left_join(slope_data, by = "tree_ID") %>% 
  left_join(soil_pca, by = "tree_ID") %>% 
  dplyr::filter(ant_rain_sum < 2000,                  # Remove outliers
                curvature < 50 & curvature > -50) %>%
  drop_na(ant_rain_sum, PC1_texture_shallow)

### Assign countries
# Load global country boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert data to spatial object and assign country labels
dat_sf <- data_annual %>%
  st_as_sf(coords = c("longitude", "latitude"),
           crs = 4326, remove = FALSE) %>%
  st_join(world %>% dplyr::select(admin)) %>%
  dplyr::filter(!admin %in% c("Democratic Republic of the Congo",
                              "Namibia", "Mozambique"))

# Drop geometry
data_annual <- dat_sf %>% st_drop_geometry()

# Save full annual dataset
write.csv(data_annual,
          file = file.path(out_dir, "Env_data_annual.csv"),
          row.names = FALSE)

#############################################################
### Aggregate phenology across years
#############################################################
dat.ag <- data_annual %>% 
  group_by(tree_ID, latitude, longitude) %>% 
  dplyr::summarise(GUD = mean(GUD),
                   SRS = mean(SRS),
                   GUD_DOY = mean(GUD_DOY),
                   SRS_DOY = mean(SRS_DOY),
                   pheno_lag = mean(pheno_lag),
                   dNDVI = mean(dNDVI),
                   NDVIseas = mean(NDVIseas),
                   ratio = mean(ratio),
                   gup_rate = mean(gup_rate),
                   rain_sum = mean(rain_sum), # based on 2019-2024 rain
                   MAP = mean(ant_rain_sum),  # based on 2018-2023 rain
                   soil_texture = mean(PC1_texture_shallow),
                   slope = mean(slope),
                   elevation = mean(elevation),
                   curvature = mean(curvature),
                   admin = first(admin),
                   .groups = "drop")


#############################################################
### Stratified spatial subsampling
#############################################################
# Sample up to 1000 trees per country, then fill to 8500 total
set.seed(123)
min_sample <- 1000
by_country <- dat.ag %>%
  st_drop_geometry() %>%
  dplyr::distinct(tree_ID, admin) %>%
  group_split(admin) %>%
  map_dfr(~ slice_sample(.x, n = min(nrow(.x), min_sample)))

remaining_n <- 8500 - nrow(by_country) # Remaining needed to hit 8500

# Draw additional samples from countries with surplus data
extra_pool <- dat.ag %>%
  distinct(tree_ID, admin) %>%
  filter(!tree_ID %in% by_country$tree_ID) %>%
  group_by(admin) %>%
  filter(n() > min_sample) %>%
  ungroup() %>%
  slice_sample(n = remaining_n)

final_ids <- bind_rows(by_country, extra_pool)

dat.ag <- dat.ag %>%
  filter(tree_ID %in% final_ids$tree_ID)

# Keep trees with early green-up for subsequent models of ∆NDVI
dat.ag.1 <- dat.ag %>%
  drop_na(dNDVI) %>% 
  filter(pheno_lag < 0, dNDVI > 0)

write.csv(dat.ag,
          file = file.path(data_dir, "Final_8500_pts.csv"),
          row.names = FALSE)

write.csv(dat.ag.1,
          file = file.path(out_dir, "subset_dNDVI_pts.csv"),
          row.names = FALSE)
