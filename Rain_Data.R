# ============================================================
# Authors: Maya Gonzalez
# Date: January 2025
#
# Description:
# - Extracts daily CHIRPS rainfall for tree locations
# - Cleans and aggregates rainfall by tree and rainyear
# - Computes annual rainfall totals, variance, and 1-year lagged rainfall
# - Identifies Start of the Rainy Season (SRS) using a threshold-based rule
# - Merges rainfall metrics with green-up data (GUD)
#
# SRS Definition:
# - First 10-day window with ≥ 25 mm rainfall
# - Followed by ≥ 20 mm rainfall in the subsequent 20 days
#
# Inputs:
# - rain_dat*_raw.csv (CHIRPS daily rainfall extracted in batches)
# - GUP.csv (tree-level green-up phenology from NDVI processing)
#
# Outputs:
# - Pheno_rain_annual_100000_trees.csv
#
# Notes:
# - Rainyear defined as July–June
# - CHIRPS data processed in several batches due to memory constraints
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(chirps)

############################################################
## Functions
############################################################

#-----------------------------------------------------------
# calculate_srs()
#
# Calculates the Start of the Rainy Season (SRS) for a single
# tree_ID × rainyear time series.
#
# Criteria:
#   1) First 10-day window has ≥ 25 mm rainfall
#   2) Following 20-day window has > 20 mm rainfall
#
# Returns:
#   - Date corresponding to the first day of the qualifying
#     10-day window
#   - NA if criteria are never met or time series is too short
#
# IMPORTANT:
#   - rain_data must be ordered by date BEFORE calling
#   - chirps must be daily rainfall (mm)
#-----------------------------------------------------------
calculate_srs <- function(rain_data) {
  n <- nrow(rain_data)  # Number of rows in the group
  
  # Return NA if the group has fewer than 30 rows
  if (n < 30) {
    return(NA)
  }
  
  for (i in 1:(n - 30)) {
    # Subset for the first 10 days
    window_10 <- rain_data$chirps[i:(i + 9)]
    total_rain_10 <- sum(window_10)
    
    # If the first 10-day criterion is met
    if (total_rain_10 >= 25) {
      # Subset for the subsequent 20 days
      window_20 <- rain_data$chirps[(i + 10):(i + 29)]
      total_rain_20 <- sum(window_20)
      
      # Check the second criterion
      if (total_rain_20 > 20) {
        srs <- rain_data$date[i]  # Store the start date of the 10-day window
        return(srs)  # Return the SRS once found
      }
    }
  }
  
  return(NA)  # Return NA if no SRS is found
}

############################################################
## Rainfall batch processing function
############################################################

#-----------------------------------------------------------
# Processes a single CHIRPS rainfall batch and returns
# SRS and annual rainfall metrics by tree_ID × rainyear.
#
# Inputs:
#   rain_dat           - CHIRPS rainfall data for a batch
#   coords_to_treeID   - lookup table mapping coordinates
#                        to tree_ID
#
# Output:
#   dat_srs            - one row per tree_ID × rainyear
#-----------------------------------------------------------
process_rain <- function(rain_dat, coords_to_treeID) {
  
  #---------------------------------------------------------
  # Clean rainfall data and attach tree_IDs
  #---------------------------------------------------------
  
  # Assign the same tree_IDs as GUP dat, modify date
  rain_clean <- rain_dat %>% 
    rename("longitude" = "lon",
           "latitude"  = "lat") %>% 
    dplyr::select(-1) %>%                   # Drop auto-generated index column
    mutate(latitude  = round(latitude, 5),  # Rounding is critical to ensure 
           longitude = round(longitude, 5)  # successful coordinate joins
           ) %>% 
    merge(coords_to_treeID,
          by = c("longitude", "latitude"),
          all.x = TRUE) %>%
    mutate(
      year = year(date),
      month = month(date),
      # Rain year starts in June
      rainyear = ifelse(month >= 6, year, year - 1),
      # Continuous day index (used for time-series operations)
      day = yday(date) + (lubridate::year(date) - 2018) * 365
      ) %>% 
    dplyr::filter(chirps >= 0)     # Remove invalid rainfall values
  
  #---------------------------------------------------------
  # Annual rainfall summaries (per tree_ID × rainyear)
  #---------------------------------------------------------
  rain_summary <- rain_clean %>%
    group_by(tree_ID, rainyear) %>%
    summarise(
      rain_sum = sum(chirps),
      rain_var = var(chirps)
      ) %>% 
    # Lagged (antecedent) rainfall is calculated *after*
    # summarisation and within tree_ID
    mutate(
      ant_rain_sum = lag(rain_sum,
                         order_by = rainyear,
                         default = NA)
      )
  
  #---------------------------------------------------------
  # Calculate SRS for each tree_ID × rainyear
  #---------------------------------------------------------
  dat_srs <- rain_clean %>%
    group_by(tree_ID, rainyear) %>%
    arrange(tree_ID, date) %>%     # Important: order dates before finding SRS
    mutate(SRS_DOY = calculate_srs(cur_data())) %>% 
    ungroup() %>% 
    group_by(tree_ID, rainyear, SRS_DOY, latitude, longitude) %>% 
    summarise() %>% 
    # Attach annual rainfall metrics
    left_join(rain_summary, 
              by = c("tree_ID", "rainyear"))
  
  return(dat_srs)
}

############################################################
## File paths
############################################################
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

############################################################
## Import green-up (phenology) data
############################################################
dat_Gup <- read.csv(file.path(data_dir, "GUP.csv")) %>% 
  mutate(latitude = round(latitude, 5),
         longitude = round(longitude, 5))

############################################################
## Prepare coordinate batches for CHIRPS extraction
############################################################
long_lat <- dat_Gup %>% 
  ungroup() %>% 
  distinct(longitude, latitude) %>% 
  rename("lon" = "longitude",
         "lat" = "latitude")

long_lat_1 <- long_lat[1:3000,]
long_lat_2 <- long_lat[3001:6000,]
long_lat_3 <- long_lat[6001:9000,]
long_lat_4 <- long_lat[9001:12000,]
long_lat_5 <- long_lat[12001:15000,]
long_lat_6 <- long_lat[15001:18000,]
long_lat_7 <- long_lat[18001:19521,]


############## Get CHIRPS Rainfall Data #################
############### Warning: SLOW ###########################
dates = c("2018-07-01", "2024-06-30")
rain_dat_1 = get_chirps(long_lat_1, dates, server = "CHC", as.matrix = FALSE)
rain_dat_2 = get_chirps(long_lat_2, dates, server = "CHC", as.matrix = FALSE)
rain_dat_3 = get_chirps(long_lat_3, dates, server = "CHC", as.matrix = FALSE)
rain_dat_4 = get_chirps(long_lat_4, dates, server = "CHC", as.matrix = FALSE)
rain_dat_5 = get_chirps(long_lat_5, dates, server = "CHC", as.matrix = FALSE)
rain_dat_6 = get_chirps(long_lat_6, dates, server = "CHC", as.matrix = FALSE)
rain_dat_7 = get_chirps(long_lat_7, dates, server = "CHC", as.matrix = FALSE)

# Save to avoid downloading data again
write.csv(rain_dat_1, file = file.path(out_dir, "rain_dat1_raw.csv"), row.names = FALSE)
write.csv(rain_dat_2, file = file.path(out_dir, "rain_dat2_raw.csv"), row.names = FALSE)
write.csv(rain_dat_3, file = file.path(out_dir, "rain_dat3_raw.csv"), row.names = FALSE)
write.csv(rain_dat_4, file = file.path(out_dir, "rain_dat4_raw.csv"), row.names = FALSE)
write.csv(rain_dat_5, file = file.path(out_dir, "rain_dat5_raw.csv"), row.names = FALSE)
write.csv(rain_dat_6, file = file.path(out_dir, "rain_dat6_raw.csv"), row.names = FALSE)
write.csv(rain_dat_7, file = file.path(out_dir, "rain_dat7_raw.csv"), row.names = FALSE)

############################################################
## Import CHIRPS Rain Data from WD
############################################################
rain_dat_1 = read.csv(file.path(data_dir, "rain_dat1_raw.csv"))
rain_dat_2 = read.csv(file.path(data_dir, "rain_dat2_raw.csv"))
rain_dat_3 = read.csv(file.path(data_dir, "rain_dat3_raw.csv"))
rain_dat_4 = read.csv(file.path(data_dir, "rain_dat4_raw.csv"))
rain_dat_5 = read.csv(file.path(data_dir, "rain_dat5_raw.csv"))
rain_dat_6 = read.csv(file.path(data_dir, "rain_dat6_raw.csv"))
rain_dat_7 = read.csv(file.path(data_dir, "rain_dat7_raw.csv"))

############################################################
## Coordinate → tree_ID lookup table
############################################################
coords_to_treeID = dat_Gup %>%
  dplyr::select(longitude, latitude, tree_ID) %>%
  distinct()

############################################################
## Process rainfall batches
############################################################
dat_srs_1 <- process_rain(rain_dat_1, coords_to_treeID)
dat_srs_2 <- process_rain(rain_dat_2, coords_to_treeID)
dat_srs_3 <- process_rain(rain_dat_3, coords_to_treeID)
dat_srs_4 <- process_rain(rain_dat_4, coords_to_treeID)
dat_srs_5 <- process_rain(rain_dat_5, coords_to_treeID)
dat_srs_6 <- process_rain(rain_dat_6, coords_to_treeID)
dat_srs_7 <- process_rain(rain_dat_7, coords_to_treeID)

############################################################
## Join rainfall metrics with green-up phenology
############################################################
# Join rain and green-up summaries based on tree_ID and 
# remove the trees that were omitted in the rain data cleaning
join_gud_srs <- function(dat_srs) {
  dat_Gup %>%
    dplyr::filter(tree_ID %in% unique(dat_srs$tree_ID)) %>%
    left_join(dat_srs,
              by = c("tree_ID", "rainyear", "latitude", "longitude")) %>%
    mutate(
      SRS_DOY = yday(SRS_DOY),
      # Circular transformation to align with phenological year
      SRS = ((SRS_DOY - 181 + 365) %% 365) + 1
    ) %>%
    drop_na(GUD, SRS)
}

dat_GUD_SRS_1 <- join_gud_srs(dat_srs_1)
dat_GUD_SRS_2 <- join_gud_srs(dat_srs_2)
dat_GUD_SRS_3 <- join_gud_srs(dat_srs_3)
dat_GUD_SRS_4 <- join_gud_srs(dat_srs_4)
dat_GUD_SRS_5 <- join_gud_srs(dat_srs_5)
dat_GUD_SRS_6 <- join_gud_srs(dat_srs_6)
dat_GUD_SRS_7 <- join_gud_srs(dat_srs_7)

############################################################
## Merge all batches into final dataset
############################################################
dat_GUD_SRS <- rbind(
  dat_GUD_SRS_1, dat_GUD_SRS_2, dat_GUD_SRS_3,
  dat_GUD_SRS_4, dat_GUD_SRS_5, dat_GUD_SRS_6,
  dat_GUD_SRS_7) %>%
  dplyr::select(
    tree_ID, latitude, longitude, rainyear,
    GUD, GUD_DOY,
    SRS, SRS_DOY,
    rain_sum, rain_var, ant_rain_sum
    )

############################################################
## Export final dataset
############################################################
write.csv(
  dat_GUD_SRS,
  file = file.path(out_dir, "Pheno_rain_annual_100000_trees.csv"),
  row.names = FALSE
)
