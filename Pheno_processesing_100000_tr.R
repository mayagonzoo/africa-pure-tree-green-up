# ============================================================
# Title: NDVI Time-Series Processing & Green-Up Detection
# Authors: Rico Holdo & Maya Gonzalez 
# Date: February 2025
#
# Description:
# - Merges pentad NDVI time series from multiple GEE exports
# - Applies cloud filtering and LOESS imputation
# - Smooths NDVI using Savitzky–Golay filtering
# - Computes green-up day (GUD) following Archibald & Scholes (2007)
#
# Inputs:
# - TREE_BATCH_*_NDVI.csv (exported from Google Earth Engine)
#
# Outputs:
# - hundred_thousand_NDVI_time_series.csv
# - GUP.csv
#
# Notes:
# - NDVI was extracted in batches due to server memory limits in GEE
# - Pentads are 5-day periods starting Jan 1
# - Rainyear defined as July–June
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(gsignal)   # Savitzky–Golay filtering for NDVI smoothing
library(TTR)       # Moving averages used in green-up calculation
library(ggplot2)

### Functions ########################################

# Function to identify time series without long gaps of missing data
# Series are rejected if they contain runs of >3 consecutive NA values,
# which are considered too long for reliable interpolation
NAseqs <- function(x){
  runs <- rle(is.na(x))
  NAruns <- runs$lengths[runs$values]
  keep <- TRUE
  if (length(NAruns) > 0){
    if (max(NAruns) > 3)
      keep <- FALSE
  }
  return(keep)
}

#### Loess impute function
# Fills short gaps in NDVI time series using LOESS interpolation.
# Intended only for short gaps (< ~15 days); long gaps are filtered earlier.
#
# Source adapted from:
# https://gis.stackexchange.com/questions/279354/ndvi-time-series-with-missing-values
#
# y            Numeric vector with NA values
# x.length     Length of resulting vector (defaults to length(y))
# s            LOESS span parameter (controls smoothness)
# smooth.data  If TRUE, smooths the entire series instead of only imputing NAs
impute.loess <- function(y, x.length = NULL, s = 0.2, 
                         smooth.data = FALSE, ...) {
  if(is.null(x.length)) { x.length = length(y) }
  options(warn=-1) # Suppress LOESS warnings for edge cases (e.g., sparse data)
  x <- 1:x.length
  p <- loess(y ~ x, span = s, data.frame(x=x, y=y))
  if(smooth.data == TRUE) {
    y <- predict(p, x)
  } else {
    na.idx <- which( is.na(y) )
    if( length(na.idx) > 1 ) {
      y[na.idx] <- predict(p, data.frame(x=na.idx))
    }
  }   
  return(y)
}

# Archibald green-up index function
# Green-up is defined as the first positive NDVI anomaly
# following the seasonal minimum in the smoothed NDVI time series
Gup.fn <- function(NDVIsg, yday){
  yday <- yday[!is.na(NDVIsg)]
  NDVIsg <- NDVIsg[!is.na(NDVIsg)]
  # Subtract a short-term moving average to identify rapid NDVI increases
  NDVIma <- TTR::SMA(NDVIsg, n = 4)
  Gup = NDVIsg - lag(NDVIma)
  # Identify seasonal minimum, assumed to precede green-up
  minidx <- which.min(NDVIsg)
  # Make greenup 0 up to this point (needs to come after)
  Gup[1:(minidx - 1)] <- 0
  # Find first positive Gup value after minidx
  Day.Gup <- yday[min(which(Gup > 0))]
  return(Day.Gup)
}

##################################################
# Change to user working directory
# ---- File paths ----
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

# Read and merge all TREE_BATCH NDVI exports from Google Earth Engine
T1 = read.csv(file.path(data_dir, "TREE_BATCH_1_NDVI.csv"))
T2 = read.csv(file.path(data_dir, "TREE_BATCH_2_NDVI.csv"))
T3 = read.csv(file.path(data_dir, "TREE_BATCH_3_NDVI.csv"))
T4 = read.csv(file.path(data_dir, "TREE_BATCH_4_NDVI.csv"))
T5 = read.csv(file.path(data_dir, "TREE_BATCH_5_NDVI.csv"))
dat = rbind(T1, T2, T3, T4, T5)
rm(T1, T2, T3, T4, T5) # clear space
# Preserve spatial metadata (lat/lon) separately to avoid duplication
# during time-series expansion and filtering
meta <- dat
meta <- dplyr::filter(dat, !duplicated(tree_ID))
meta <- meta %>% dplyr::select(tree_ID, latitude, longitude)
dat <- dat %>% dplyr::select(tree_ID, year, pentad, NDVI, cloud_prob)
# Pad dat with NAs for missing dates
dat <- dat %>% complete(tree_ID, year, pentad)
# Convert pentad index to calendar date
# Pentads are defined as 5-day periods starting on January 1
dat$date = as.Date((dat$pentad - 1) * 5 + 1, origin = paste0(dat$year, "-01-01"))
dat$month <- month(dat$date)
# Define rainyear (July–June) to align with regional wet–dry seasonality
dat <- dat %>% mutate(rainyear = ifelse(month > 6, year, year - 1))
dat$yday <- yday(dat$date)
dat$Day <- dat$yday + (year(dat$date) - 2017) * 365
dat <- dat %>% dplyr::filter(rainyear > 2016 & rainyear < 2024)
# Mask NDVI values with high cloud probability (>10%)
# Cloud threshold chosen to balance data availability and quality
dat <- dat %>% mutate(NDVI = ifelse(cloud_prob < 10, NDVI, NA))
# Remove negative NDVI values
dat <- dat %>% mutate(NDVI = ifelse(NDVI < 0, NA, NDVI))

# Find dry season pentad sequences that don't have NA runs
# that are more than 3 long (hard to interpolate if sequences are long)
dat <- dat %>% mutate(seas = ifelse(month > 6 & month < 12, 'Dry', 'Wet'))
dat.ag <- dat %>% dplyr::filter(seas == 'Dry') %>% 
  group_by(tree_ID, rainyear) %>% 
  summarise(NAseqs = NAseqs(NDVI))
# Remove locations/years with long NA runs
dat.ag <- dat.ag %>% dplyr::filter(NAseqs)
# We only want to keep these cases in dat
dat.ag <- dat.ag %>% mutate(IDyear = paste(tree_ID, rainyear, sep = '_'))
dat <- dat %>% mutate(IDyear = paste(tree_ID, rainyear, sep = '_'))
dat <- dat %>% dplyr::filter(IDyear %in% dat.ag$IDyear)
dat <- dat %>% dplyr::select(-IDyear)

# NDVI processing pipeline:
# 1) Impute short gaps using LOESS
# 2) Smooth full series using Savitzky–Golay filtering
# Before any smoothing, missing pentads need to be imputed
dat <- dat %>% group_by(tree_ID, rainyear) %>% 
  mutate(NDVIint = impute.loess(NDVI))
# Now smooth data
dat <- dat %>% group_by(tree_ID, rainyear) %>% 
  mutate(NDVIsg = sgolayfilt(NDVIint))

# Retain only tree-years with a clear seasonal NDVI signal:
# - Minimum NDVI < 0.3 during dry season
# - Maximum NDVI > 0.7 during wet season
dat.ag2 <- dat %>% 
  group_by(tree_ID, rainyear) %>%
  summarise(NDVImin = min(NDVIsg, na.rm = TRUE),
            NDVImax = max(NDVIsg, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(passes = NDVImin < 0.3 & NDVImax > 0.7)

# Retain only trees for which *all observed years* exhibit a 
# clear seasonal NDVI cycle
dat.ag2 <- dat.ag2 %>%
  group_by(tree_ID) %>%
  summarise(all_pass = all(passes), .groups = "drop") %>%
  dplyr::filter(all_pass) %>%
  pull(tree_ID)

dat <- dat %>% dplyr::filter(tree_ID %in% dat.ag2)

# Redefine day-of-year relative to start of rainyear (July 1 = day 1)
dat <- dat %>%
  mutate(rainday = if_else(yday(date) >= 182, yday(date) - 181, 
                           yday(date) + 184))

# Remove tree-years with any negative NDVIsg values
dat <- dat %>%
  group_by(tree_ID, rainyear) %>%
  dplyr::filter(all(NDVIsg >= 0 | is.na(NDVIsg))) %>%
  ungroup()

# Convert date column to Date
dat$date = as.Date(dat$date)

#export NDVI time series data
write.csv(dat, file = file.path(out_dir, "hundred_thousand_NDVI_time_series.csv"),
          row.names = FALSE)

dat = read.csv(file.path(data_dir, "hundred_thousand_NDVI_time_series.csv"))
# Compute green-up day in both rainyear day units and calendar DOY
dat.Gup <- dat %>% 
  group_by(tree_ID, rainyear) %>% 
  summarize(GUD = Gup.fn(NDVIsg, rainday),
            GUD_DOY = Gup.fn(NDVIsg, yday))

# Reconnect with metadata to create an annual summary of green-up
dat.Gup <- left_join(dat.Gup, meta)
dat.Gup$tree_ID <- factor(dat.Gup$tree_ID)
dat.Gup$tree_ID <- as.factor(dat.Gup$tree_ID)

# Export
write.csv(dat.Gup, file = file.path(out_dir, "GUP.csv"), row.names = FALSE)

