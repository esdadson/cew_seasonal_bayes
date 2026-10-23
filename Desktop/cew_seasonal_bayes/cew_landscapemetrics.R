library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(landscapemetrics)
library(rnaturalearth)
library(parallel)
library(future)
library(future.apply)

# Load CEW data
data_full <- read.csv("/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/doi_10_5061_dryad_m0cfxpp5x__v20221206/Lawton_etal_data.csv") %>%
  filter(longitude >= -105)

# States vector
states_sf <- ne_states(country = 'united states of america', returnclass = 'sf') %>%
  select(postal, geometry)  # keep only what's needed

# Reclassification of cropland data layer
rcl <- matrix(c(
  # Corn
  1,1, 12,1, 13,1,
  # Cotton
  2,2,
  # Soybean
  5,5, 25,5, 26,5,
  # Peanuts
  10,10, 
  # Tobacco
  11,11,
  # Other crops
  3,44, 4,44, 6,44, 14,44, 21,44, 22,44, 23,44, 24,44, 25,44,
  27,44, 28,44, 29,44, 30,44, 31,44, 32,44, 33,44, 35,44, 36,44,
  39,44, 41,44, 42,44, 43,44, 44,44, 45,44, 46,44, 47,44, 48,44,
  49,44, 50,44, 53,44, 54,44, 55,44, 57,44, 59,44, 60,44, 61,44,
  67,44, 68,44, 69,44, 70,44, 71,44, 72,44, 74,44, 76,44, 77,44,
  204,44, 205,44, 206,44, 207,44, 208,44, 209,44, 211,44,
  # Non crop
  37,63, 58,63, 63,63, 64,63, 65,63, 82,63, 83,63, 87,63, 92,63,
  111,63, 121,63, 122,63, 123,63, 124,63, 131,63, 141,63, 142,63,
  143,63, 152,63, 176,63, 190,63, 195,63
), ncol = 2, byrow = TRUE)

# Reclassify everything else to non crop
unmapped <- setdiff(1:255, rcl[,1])
rcl_complete <- rbind(rcl, cbind(unmapped, 63L))

# Set class names table
class_names <- data.frame(
  value = c(1, 2, 5, 10, 11, 44, 63),
  class = c("Corn", "Cotton", "Soybeans", "Peanuts", "Tobacco", "Other_Crops", "Non_Crop")
)

# ============================================================================
# 4. FUNCTION TO PROCESS ONE YEAR
# ============================================================================
process_year <- function(yr, data_full, states_sf, rcl_complete, class_names) {
  
  message("Processing year: ", yr)
  
  # --- 4a. Subset points for this year ---
  locs <- data_full %>%
    filter(year == yr) %>%
    select(location2, longitude, latitude) %>%
    distinct(longitude, latitude, .keep_all = TRUE)
  
  if (nrow(locs) == 0) {
    message("  No points for year ", yr, " — skipping.")
    return(NULL)
  }
  
  # --- 4b. Identify only the states that contain points this year ---
  points_sf_4326 <- st_as_sf(locs, coords = c("longitude", "latitude"), crs = 4326)
  
  # fix: filter to only points that fall within a US state (drops points outside the US)
  points_with_state <- st_join(points_sf_4326, states_sf, join = st_within)
  points_in_us <- points_with_state %>% filter(!is.na(postal))
  
  if (nrow(points_in_us) == 0) {
    message("No points fell within any state for year ", yr, " — skipping.")
    return(NULL)
  }
  
  active_states <- unique(points_in_us$postal)
  message("  States with points (", yr, "): ", paste(sort(active_states), collapse = ", "))
  
  states_mask_4326 <- states_sf %>%
    filter(postal %in% active_states) %>%
    st_union() %>%
    st_as_sf()
  
  # --- 4c. Load raster, reproject mask, crop to active states only ---
  cdl_path <- paste0("/Volumes/cmjone25/Data/Raster/USA/landcover/", yr, "_30m_cdls/", yr, "_30m_cdls.tif")
  r        <- rast(cdl_path)
  
  states_mask_native <- vect(st_transform(states_mask_4326, crs(r)))
  r <- crop(r, states_mask_native, mask = TRUE)
  r <- classify(r, rcl_complete)
  
  # --- 4d. Reproject the US-filtered points to raster CRS ---
  # fix: use only the in-US points, and snap to exact same CRS object from raster
  points_in_us <- st_transform(points_in_us, crs(r))
  points_mat   <- st_coordinates(points_in_us)
  
  # --- 4e. Extract landscape metrics ---
  metrics_raw <- tryCatch(
    sample_lsm(
      landscape = r,
      y         = points_mat,
      plot_id   = points_in_us$location2,   # use unique plot_id
      size      = 1000,
      shape     = "circle",
      what      = c("lsm_c_pland", "lsm_c_lpi", "lsm_c_pd", "lsm_c_ed"),
      verbose   = FALSE,
      progress  = FALSE
    ),
    error = function(e) {
      message("  Error in year ", yr, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(metrics_raw)) return(NULL)
  
  # --- 4f. Reshape, then join plot_id back to original location name ---
  metrics_wide <- metrics_raw %>%
    merge(class_names, by.x = "class", by.y = "value") %>%
    select(location2, class.y, metric, value) %>%
    pivot_wider(
      id_cols     = location2,
      names_from  = c(class.y, metric),
      values_from = value,
      values_fill = 0
    ) %>%
    # # fix: recover original location name by joining on the unique plot_id
    # left_join(st_drop_geometry(points_in_us) %>% select(location2, plot_id), by = "plot_id") %>%
    # select(-plot_id) %>%
    mutate(year = yr)
  
  data_full %>%
    filter(year == yr) %>%
    left_join(metrics_wide, by = c("location2", "year"))
}

# ============================================================================
# 5. PARALLEL LOOP OVER YEARS
# ============================================================================
years   <- 2008:2021

results_list <- lapply(
  years,
  process_year,
  data_full    = data_full,
  states_sf    = states_sf,
  rcl_complete = rcl_complete,
  class_names  = class_names
)

# ============================================================================
# 6. COMBINE
# ============================================================================
d_new <- bind_rows(Filter(Negate(is.null), results_list))


# ============================================================================
# CLIMATE COVARIATE EXTRACTION
# ============================================================================
states <- vect(states_sf)
d_new$tmax <- NA
d_new$tmin <- NA
d_new$gdd <- NA
# Create spatial points from your data
d_spatial <- vect(d_new, geom = c('longitude', 'latitude'), crs = "EPSG:4326")
for (i in seq(2008, 2021)) {
  # Load raster data
  tmax_rast <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmax/daymet_v4_daily_na_tmax_", i, ".nc"))
  tmin_rast <- rast(paste0("/Volumes/cmjone25/Data/Original/Daymet/tmin/daymet_v4_daily_na_tmin_", i, ".nc"))
  agdds_rast <- rast(paste0("/Volumes/cmjone25/Data/Raster/USA/pops_casestudies/corn_earworm/gdd/cum_gdd_", i, ".tiff"))
  states_sf <- st_transform(states_sf, crs(tmax_rast))
  b = st_bbox(states_sf)
  #tmax_rast <- crop(tmax_rast, b)
  #tmin_rast <- crop(tmin_rast, b)
  
  unique_doys <- unique(d_new[d_new$year == i, ]$doy)
  unique_doys <- unique_doys[!is.na(unique_doys)]
  # Extract TMAX at location and date
  tmax_values <- rep(NA_real_, nrow(d_new[d_new$year == i, ]))
  for (doy in unique_doys) {
    # Skip if DOY is outside raster range
    if(doy > nlyr(tmax_rast) || doy < 1) next
    # Find all rows with this DOY
    rows_for_doy <- which(d_new[d_new$year == i, ]$doy == doy)
    # Get corresponding spatial points
    points_for_doy <- d_spatial[d_spatial$year == i,][rows_for_doy,]
    # Extract for all points at once
    extracted <- terra::extract(tmax_rast[[doy]], points_for_doy)[, 2]
    # Assign to result vector
    tmax_values[rows_for_doy] <- extracted
  }
  d_new[d_new$year == i, ]$tmax <- tmax_values
  # Extract TMIN at location and date
  tmin_values <- rep(NA_real_, nrow(d_new[d_new$year == i, ]))
  # Loop over unique DOYs instead of all rows
  for(doy in unique_doys) {
    # Skip if DOY is outside raster range
    if(doy > nlyr(tmax_rast) || doy < 1) next
    # Find all rows with this DOY
    rows_for_doy <- which(d_new[d_new$year ==i,]$doy == doy)
    # Get corresponding spatial points
    points_for_doy <- d_spatial[d_spatial$year == i,][rows_for_doy,]
    # Extract for all points at once
    extracted <- terra::extract(tmin_rast[[doy]], points_for_doy)[, 2]
    # Assign to result vector
    tmin_values[rows_for_doy] <- extracted
  }
  d_new[d_new$year ==i,]$tmin <- tmin_values
  
  # Extract AGDD at location and date
  agdd_values <- rep(NA_real_, nrow(d_new[d_new$year == i, ]))
  # Loop over unique DOYs instead of all rows
  for(doy in unique_doys) {
    # Skip if DOY is outside raster range
    if(doy > nlyr(agdds_rast) || doy < 1) next
    # Find all rows with this DOY
    rows_for_doy <- which(d_new[d_new$year == i,]$doy == doy)
    # Get corresponding spatial points
    points_for_doy <- d_spatial[d_spatial$year == i,][rows_for_doy,]
    # Extract for all points at once
    extracted <- terra::extract(agdds_rast[[doy]], points_for_doy)[, 2]
    # Assign to result vector
    agdd_values[rows_for_doy] <- extracted
  }
  d_new[d_new$year == i,]$agdd <- agdd_values
  rm(tmax_rast)
  rm(tmin_rast)
  rm(agdds_rast)
}
t_base <- 10
d_new <- d_new %>%
  mutate(
    gdd = pmax(((tmax + tmin) / 2) - t_base, 0)  # vectorized, no ifelse needed
  )
write.csv(d_new, file = "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/cew_landscape_metrics.csv", row.names = F)
d_new <- bind_rows(Filter(Negate(is.null), results_climate))


d_new <- read.csv("/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/cew_landscape_metrics.csv")
d_new$doy <- as.numeric(format(as.Date(paste(d_new$year, d_new$woy * 7, sep = "-"), 
                                       format = "%Y-%j"), "%j"))