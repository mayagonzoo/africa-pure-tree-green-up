# ============================================================
# Title: Time series Screen and Additional Phenology Metrics
# Author: Maya Gonzalez
# Date: July 2025
#
# Description:
# - Screens tree-level NDVI time series for bimodal seasonal behavior
# - Merges cleaned NDVI data with rainfall-derived phenology metrics
# - Computes rates and magnitudes of NDVI change around green-up (GUD)
#   and start of the rainy season (SRS)
#
# Inputs:
# - hundred_thousand_NDVI_time_series.csv
#   (pentad NDVI time series derived from 
#   Pheno_processsing_10000_tr.R)
# - Pheno_rain_annual_100000_trees.csv
#
# Outputs:
# - Full_pheno_params_cleaned.csv
#   (tree-year phenological rate and magnitude metrics)
#
# ============================================================

library(dplyr)
library(ggplot2)
library(zoo)

# ---- File paths ----
# Users should update these paths to their local environment
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

#############################################################
### Functions
#############################################################

# Detect bimodal NDVI seasonal behavior
#
# This function identifies time series with ≥2 growth peaks separated
# by sufficiently deep senescence valleys. Valleys occurring during
# the wet season are only considered valid if NDVI drops below a
# conservative threshold, preventing false bimodality due to noise.
#
# Args:
# - ndvi: numeric vector of smoothed NDVI values
# - seas: character vector indicating season ("Wet" or "Dry")
# - wet_valley_thresh: NDVI threshold below which wet-season valleys
#   are considered true senescence events
#
# Returns:
# - TRUE if a valid bimodal growth cycle is detected
# - FALSE otherwise
detect_bimodal <- function(ndvi, seas, wet_valley_thresh = 0.5) {
  
  # Indices of local NDVI maxima (growth peaks)
  grw <- which(diff(sign(diff(ndvi))) == -2) + 1
  
  # Indices of local NDVI minima (candidate senescence valleys)
  sns_all <- which(diff(sign(diff(ndvi))) == 2) + 1
  
  # Retain dry-season valleys and sufficiently deep wet-season valleys
  sns <- sns_all[
    seas[sns_all] == "Dry" |
      (seas[sns_all] == "Wet" & ndvi[sns_all] < wet_valley_thresh)
  ]
  
  # Require at least two growth peaks and three senescence points
  if (length(grw) < 2 || length(sns) < 3) return(FALSE)
  
  valid_cycles <- 0  # Counter for valid bimodal cycles
  
  # Evaluate peak–valley–peak structure
  for (i in 2:length(grw)) {
    grw1 <- ndvi[grw[i - 1]]
    grw2 <- ndvi[grw[i]]
    
    sns_range <- sns[sns > grw[i - 1] & sns < grw[i]]
    if (length(sns_range) == 0) next
    
    sns_idx <- sns_range[which.min(ndvi[sns_range])]
    sns_value <- ndvi[sns_idx]
    
    # Require sufficient amplitude on both sides of the valley
    if ((grw1 - sns_value >= 0.4) &&
        (grw2 - sns_value >= 0.4)) {
      valid_cycles <- valid_cycles + 1
    }
  }
  
  # Enforce seasonal consistency
  if (!any(seas[sns] == "Wet")) return(FALSE)
  if (!any(seas[grw] == "Wet")) return(FALSE)
  
  return(valid_cycles >= 1)
}

#############################################################
### Data import
#############################################################
dat = read.csv(file.path(data_dir, "hundred_thousand_NDVI_time_series.csv")) %>% 
  mutate(date = as.Date(date))

dat.1 = dat # Backup copy to avoid re-importing

# Annual rainfall and SRS metrics
dat_GUD_SRS = read.csv(file.path(data_dir, "Pheno_rain_annual_100000_trees.csv"))

#############################################################
### NDVI smoothing and bimodality screening
#############################################################

# Apply a short rolling mean to reduce high-frequency noise
dat_roll <- dat %>%
  arrange(tree_ID, rainyear, date) %>%
  group_by(tree_ID, rainyear) %>%
  mutate(NDVIsg = rollmean(NDVIsg, k = 3, fill = NA, align = "center")) %>%
  ungroup()

# Identify tree-year series exhibiting bimodal growth patterns
bimodal_flags <- dat_roll %>%
  group_by(tree_ID, rainyear) %>%
  summarise(bimodal = detect_bimodal(NDVIsg, seas),
            .groups = "drop") %>%
  dplyr::filter(bimodal == TRUE)

# Visual inspection of a random subset of flagged series
bimodal_flags %>%
  slice_sample(n = 12) %>%
  inner_join(dat, by = c("tree_ID", "rainyear")) %>%
  ggplot(aes(x = date, y = NDVIsg)) +
  geom_line(group = 1) +
  facet_wrap(~ tree_ID + rainyear, scales = "free_x") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  labs(x = "Date", y = "NDVI") +
  theme_minimal()

# Remove all trees exhibiting bimodal behavior
dat <- dat %>%
  anti_join(bimodal_flags %>% 
              dplyr::select(tree_ID) %>%
              distinct(),
            by = "tree_ID")

#############################################################
### Phenological rate and magnitude metrics
#############################################################

dat_rate <- dat %>%
  inner_join(dat_GUD_SRS, by = c("tree_ID", "rainyear")) %>% 
  group_by(tree_ID, rainyear) %>%
  mutate(
    NDVIsg_max = max(NDVIsg, na.rm = TRUE),
    max_day = rainday[which.max(NDVIsg)]
  ) %>%
  filter(!is.na(GUD_DOY)) %>%
  mutate(
    # NDVI near green-up (±2 days)
    GUD_NDVI = {
      idx <- which(abs(yday - GUD_DOY) <= 2)
      if (length(idx) > 0) NDVIsg[idx[1]] else NA_real_
    },
    # NDVI near start of rainy season (±2 days)
    SRS_NDVI = {
      idx <- which(abs(yday - SRS_DOY) <= 2)
      if (length(idx) > 0) NDVIsg[idx[1]] else NA_real_
    },
    pheno_lag = GUD - SRS,              # days between SRS and GUD
    dNDVI = SRS_NDVI - GUD_NDVI,        # pre-rain NDVI increase
    NDVIseas = NDVIsg_max - GUD_NDVI,   # total seasonal increase
    ratio = dNDVI / NDVIseas,           # relative pre-rain contribution
    gup_rate = max_day - GUD            # days from GUD to NDVI peak
  ) %>%
  filter(rainday >= GUD & rainday <= max_day) %>%
  ungroup()

#############################################################
### Summarize and export
#############################################################

gup_data <- dat_rate %>% 
  distinct(tree_ID, rainyear, GUD, GUD_NDVI, SRS_NDVI, NDVIsg_max,
           max_day, pheno_lag, dNDVI, NDVIseas, ratio, gup_rate) %>% 
  filter(gup_rate > 0)

write.csv(gup_data,
          file = file.path(out_dir, "Full_pheno_params_cleaned.csv"),
          row.names = FALSE)