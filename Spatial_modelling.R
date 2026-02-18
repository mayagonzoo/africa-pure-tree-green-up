# ============================================================
# Author: Maya Gonzalez
# Date: August 2025
#
# Description:
# - Fits linear and spatial autoregressive (SAR) error models to explain
#   variation in GUD and ∆NDVI
# - Compares candidate models using AIC
#
# Inputs:
# - Final_8500_pts.csv
#
# Outputs:
# - Serialized SAR model objects (.rds)
#
# Notes:
# - Data were downsampled to 8,500 points in the previous script to improve 
#   convergence and run time of SAR models
# - All SAR models use spatial error structure (errorsarlm)
# - Distance threshold (350 km) chosen to ensure connectivity
# ============================================================

# ---- Required packages ----
library(dplyr)
library(sf)
library(spdep)
library(spatialreg)
library(gstat)

# Change to user working directory
# ---- File paths ----
data_dir <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"
out_dir  <- "/Volumes/Holdo_Lab/Projects/Africa tree grass phenology/Public Datasets"

# ============================================================
# Import data
# ============================================================
dat.ag = read.csv(file.path(data_dir, "Final_8500_pts.csv"))
dat.ag.1 = read.csv(file.path(data_dir, "subset_dNDVI_pts.csv"))

# ============================================================
# OLS models with spatial trend
# ============================================================
mod0 = lm(GUD ~ latitude + longitude, data = dat.ag)

mod1.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2), 
                data = dat.ag)

mod2.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP, data = dat.ag)

mod3.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP+ elevation, data = dat.ag)

mod4.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP * elevation, data = dat.ag)

mod5.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP + soil_texture, data = dat.ag)

mod6.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP * soil_texture, data = dat.ag)

mod7.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP* slope, data = dat.ag)

mod8.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP * curvature, data = dat.ag)

mod9.trend = lm(GUD ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                  MAP * slope + elevation, data = dat.ag)

aic_table = AIC(mod0, mod1.trend, mod2.trend, mod3.trend, mod4.trend,
                mod5.trend, mod6.trend, mod7.trend, mod8.trend, mod9.trend)
aic_table$delta = aic_table$AIC - min(aic_table$AIC)
aic_table_sorted = aic_table[order(aic_table$delta), ]
aic_table_sorted


# ============================================================
# Spatial autoregressive (SAR) error models for GUD
# Warning: May be SLOW 
# ============================================================
coords <- cbind(dat.ag$longitude, dat.ag$latitude) 

# Distance-based neighbors (km, longlat = TRUE)
nb_dist <- dnearneigh(coords, d1 = 0, d2 = 350, longlat = TRUE) 
lw_dist <- nb2listw(nb_dist, style = "W", zero.policy = TRUE) # Create spatial weights

mod1.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2),
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod1.sar, file = file.path(out_dir, "mod1.sar"))

mod2.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod2.sar, file = file.path(out_dir, "mod2.sar"))

mod3.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP + elevation,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod3.sar, file = file.path(out_dir, "mod3.sar"))

mod4.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * elevation,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod4.sar, file = file.path(out_dir, "mod4.sar"))

mod5.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP + soil_texture,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod5.sar, file = file.path(out_dir, "mod5.sar"))

mod6.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * soil_texture,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod6.sar, file = file.path(out_dir, "mod6.sar"))

mod7.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * slope,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod7.sar, file = file.path(out_dir, "mod7.sar"))

mod8.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * curvature,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod8.sar, file = file.path(out_dir, "mod8.sar"))

mod9.sar <- errorsarlm(GUD ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * slope + elevation,
                       data = dat.ag, listw = lw_dist, method = "eigen")
saveRDS(mod9.sar, file = file.path(out_dir, "mod9.sar"))

aic_table <- AIC(mod1.sar, mod2.sar, mod3.sar, mod4.sar,
                 mod5.sar, mod6.sar, mod7.sar, mod8.sar, mod9.sar)
aic_table$delta <- aic_table$AIC - min(aic_table$AIC)
aic_table_sorted <- aic_table[order(aic_table$delta), ]
aic_table_sorted


# ============================================================
# dNDVI models (subset: early green-up only)
# ============================================================

mod9.5.trend = lm(dNDVI ~ longitude + latitude, data = dat.ag.1)
mod10.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2),
                 data = dat.ag.1)
mod11.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP, data = dat.ag.1)

mod12.trend = lm(dNDVI~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP + elevation, data = dat.ag.1)

mod13.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP * elevation, data = dat.ag.1)

mod14.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP + soil_texture, data = dat.ag.1)

mod15.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP * soil_texture, data = dat.ag.1)

mod16.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP * slope, data = dat.ag.1)

mod17.trend = lm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2) +
                   MAP * curvature, data = dat.ag.1)

aic_table = AIC(mod9.5.trend, mod10.trend, mod11.trend, mod12.trend, mod13.trend, mod14.trend,
                mod15.trend, mod16.trend, mod17.trend)
aic_table$delta = aic_table$AIC - min(aic_table$AIC)
aic_table_sorted = aic_table[order(aic_table$delta), ]
aic_table_sorted

# ---- SAR models for dNDVI ----
coords_1 <- cbind(dat.ag.1$longitude, dat.ag.1$latitude) 

# Distance-based neighbors (km, longlat = TRUE)
nb_dist_1 <- dnearneigh(coords_1, d1 = 0, d2 = 350, longlat = TRUE) 
lw_dist_1 <- nb2listw(nb_dist_1, style = "W", zero.policy = TRUE) # Create spatial weights

mod10.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) + latitude + I(latitude^2), 
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod10.sar, file = file.path(out_dir,"mod10.sar"))

mod11.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) + 
                         latitude + I(latitude^2) + MAP, 
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod11.sar, file = file.path(out_dir,"mod11.sar"))

mod12.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) + 
                         latitude + I(latitude^2) + MAP + elevation, 
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod12.sar, file = file.path(out_dir,"mod12.sar"))

mod13.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * elevation,
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod13.sar, file = file.path(out_dir,"mod13.sar"))

mod14.sar = errorsarlm(dNDVI  ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP + soil_texture,
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod14.sar, file = file.path(out_dir,"mod14.sar"))

mod15.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * soil_texture,
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod15.sar, file = file.path(out_dir,"mod15.sar"))

mod16.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * slope,
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod16.sar, file = file.path(out_dir,"mod16.sar"))

mod17.sar = errorsarlm(dNDVI ~ longitude + I(longitude^2) +
                         latitude + I(latitude^2) + MAP * curvature,
                       data = dat.ag.1, listw = lw_dist_1, method = "eigen")
saveRDS(mod17.sar, file = file.path(out_dir,"mod17.sar"))


aic_table = AIC(mod10.sar, mod11.sar, mod12.sar, mod13.sar,
                mod14.sar, mod15.sar, mod16.sar, mod17.sar)
aic_table$delta = aic_table$AIC - min(aic_table$AIC)
aic_table_sorted = aic_table[order(aic_table$delta), ]
aic_table_sorted


# ============================================================
# Standardize coefficients of the best-fitting models
# ============================================================
dat.ag$MAP_s  <- scale(dat.ag$MAP)
dat.ag$elevation_s <- scale(dat.ag$elevation)
dat.ag$slope_s <- scale(dat.ag$slope)

mod9.std <- errorsarlm(GUD ~ longitude + I(longitude^2) + latitude + 
                         I(latitude^2) + MAP_s * slope_s + elevation_s,
                       data = dat.ag, listw = lw_dist, method = "eigen")
summary(mod9.std)


dat.ag.1$MAP_s  <- scale(dat.ag.1$MAP)
dat.ag.1$soil_texture_s <- scale(dat.ag.1$soil_texture)

mod15.std <- errorsarlm(dNDVI ~ longitude + I(longitude^2) +
                          latitude + I(latitude^2) + MAP_s * soil_texture_s,
                        data = dat.ag.1, listw = lw_dist_1, method = "eigen")
summary(mod15.std)

