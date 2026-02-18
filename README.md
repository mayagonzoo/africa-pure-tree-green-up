# Africa pure tree green-up

Code and datasets for the manuscript titled "Pre-Rain Green-Up in African woodlands and savannas is driven by hydrology rather than photoperiod".


Rasters of pure tree cover (resampled to 10m) were requested directly from the authors of [More than one quarter of Africa’s tree cover is found outside areas previously classified as forest](https://doi.org/10.1038/s41467-023-37880-4).
All other datasets used in this analysis were derived from the following publicly available sources:
- [ESA WorldCover](esa-worldcover.org/en)
- [MODIS Land Cover Type Yearly Global (MCD12Q1.061)](www.earthdata.nasa.gov/data/catalog/lpcloud-mcd12q1-061)
- [Harmonized Sentinel-2 Level 2-A (Surface Reflectance)](https://dataspace.copernicus.eu/)
- [Climate Hazards Center InfraRed Precipitation with Station data (CHIRPS v2)](data.chc.ucsb.edu/products/CHIRPS-2.0)
- [NASA DEM Global 1 Arc Second V001](www.earthdata.nasa.gov/data/catalog/lpcloud-nasadem-hgt-001)
- [WoSIS SoilGrids 2.0](www.isric.org/explore/soilgrids)

## Workflow

Scripts must be run in the following order. Outputs from each step are used as inputs for subsequent steps.

1. **[Pheno_processing_100000_tr.R](Pheno_processing_100000_tr.R)**  
   *Processes Sentinel-2 NDVI time series and detects green-up day (GUD).*  
   **Outputs:**  
   - hundred_thousand_NDVI_time_series.csv  
   - `GUP.csv`

2. **[Rain_Data.R](Rain_Data.R)**  
   *Processes CHIRPS rainfall data and computes start of the rainy season (SRS).*  
   **Inputs:**  
   - `GUP.csv`  
   **Outputs:**  
   - `Pheno_rain_annual_100000_trees.csv`

3. **[Seasonality_screen_pheno_metrics.R](Seasonality_screen_pheno_metrics.R)**  
   *Screens NDVI time series for bimodality and computes phenological metrics.*  
   **Inputs:**  
   - hundred_thousand_NDVI_time_series.csv  
   - `Pheno_rain_annual_100000_trees.csv`  
   **Outputs:**  
   - `Full_pheno_params_cleaned.csv`

4. **[Data_downsampling.R](Data_downsampling.R)**  
   *Integrates environmental variables and downsamples the dataset for spatial modeling.*  
   **Inputs:**  
   - `Pheno_rain_annual_100000_trees.csv`  
   - `Full_pheno_params_cleaned.csv`  
   - `soil_pca_0_60_data.csv`  
   - `DEMData_100000_trees.csv`  
   **Outputs:**  
   - `Env_data_annual.csv`  
   - `Final_8500_pts.csv`  
   - `subset_dNDVI_pts.csv`  

5. **[Spatial_modelling.R](Spatial_modelling.R)**  
   *Fits spatial autoregressive models for green-up timing and pre-rain productivity.*  
   **Inputs:**  
   - `Final_8500_pts.csv`  
   - `subset_dNDVI_pts.csv`  
   **Outputs:**  
   - Model objects saved as `.rds` files (see `Model Outputs/`)






