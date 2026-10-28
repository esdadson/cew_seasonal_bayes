# =============================================================================
# step1/step1_data.R
# Modelling data for step 1, built from the covariate extraction
# (cew_covariates_hpc.R -> covariates/): ISO-week trap counts with weekly
# Daymet GDD, CDL crop shares and MODIS EVI for every trapped site-year
# 2008-2020, plus a prediction grid - weeks 1-52 of every site-year that was
# monitored - which step 1 predicts and step 2 uses.
#
# A cleaned copy of sections 0-4 of cew_stage1_data.R (that file is left as
# it was): every year is in the fit (no held-out 2020), and winter severity,
# tmin/precipitation and scaled lon/lat, which steps 1-2 don't use, are gone.
# Run from the project directory. Creates sites, d_all (observed rows +
# prediction grid), d_obs (observed rows), yrs, my_vars and anom_vars.
# =============================================================================

library(dplyr); library(tidyr); library(lubridate); library(zoo)

yrs    <- 2008:2020          # CDL coverage + dense Southern trapping
t_base <- 10                 # GDD base temperature (deg C)

# -----------------------------------------------------------------------------
# Helpers (from cew_source_sink_functions.R)
# -----------------------------------------------------------------------------
# Data root (contains Original/): $CEW_DATA if set, else the first that exists
# of Hazel research storage, the Mac mount, the Windows drive
default_data_root <- function() {
  env <- Sys.getenv("CEW_DATA")
  if (nzchar(env)) return(env)
  cand <- c("/rs1/researchers/c/cmjone25/Data", "/Volumes/cmjone25/Data", "Z:/Data")
  hit <- cand[dir.exists(cand)]
  if (!length(hit)) stop("Data root not found; set CEW_DATA. Tried: ", paste(cand, collapse = ", "))
  hit[1]
}

# ONE site table: location2 by sorted coordinates (the ids of covariates/),
# plus CONUS Albers km for the GP - degrees of longitude shrink with
# latitude, so lon/lat is not isotropic
make_site_table <- function(counts_raw) {
  sites <- counts_raw %>%
    group_by(longitude, latitude) %>%
    summarise(location = first(location), zone_30y = first(zone_30y), .groups = "drop") %>%
    arrange(longitude, latitude) %>%
    mutate(location2 = row_number())
  xy <- albers_km(sites$longitude, sites$latitude)
  if (anyNA(xy)) stop(sum(!complete.cases(xy)), " of ", nrow(xy), " sites could not be projected to EPSG:5070")
  sites$x_km <- xy[, 1]; sites$y_km <- xy[, 2]
  sites
}

# CONUS Albers equal-area conic (EPSG:5070: GRS80, standard parallels 29.5 and
# 45.5, origin 23N 96W), km. Snyder (1987) eqs. 14-3 to 14-15; agrees with
# PROJ to ~1 m. Closed form because PROJ on Hazel's conda env returned NA.
albers_km <- function(lon, lat) {
  a <- 6378137; f <- 1 / 298.257222101; e2 <- 2 * f - f^2; e <- sqrt(e2); rad <- pi / 180
  q <- function(p) { s <- sin(p); (1 - e2) * (s / (1 - e2 * s^2) - log((1 - e * s) / (1 + e * s)) / (2 * e)) }
  m <- function(p) cos(p) / sqrt(1 - e2 * sin(p)^2)
  p0 <- 23 * rad; p1 <- 29.5 * rad; p2 <- 45.5 * rad; l0 <- -96 * rad
  n   <- (m(p1)^2 - m(p2)^2) / (q(p2) - q(p1)); C <- m(p1)^2 + n * q(p1)
  rho <- function(p) a * sqrt(C - n * q(p)) / n
  th  <- n * (lon * rad - l0)
  cbind(rho(lat * rad) * sin(th), rho(p0) - rho(lat * rad) * cos(th)) / 1000
}

# Standardise with means/sds from reference rows. map: new column = source column
make_scaler <- function(ref, map) {
  m <- vapply(map, function(v) mean(ref[[v]], na.rm = TRUE), numeric(1))
  s <- vapply(map, function(v) sd(ref[[v]], na.rm = TRUE), numeric(1))
  function(df) {
    for (k in names(map)) df[[k]] <- (df[[map[[k]]]] - m[[k]]) / s[[k]]
    df
  }
}

# Fourier basis on ISO week (period = 52 weeks)
add_fourier <- function(df, H = 2, week_col = "woy", period = 52) {
  w <- df[[week_col]]
  for (h in seq_len(H)) {
    df[[paste0("s", h)]] <- sin(2 * pi * h * w / period)
    df[[paste0("c", h)]] <- cos(2 * pi * h * w / period)
  }
  df
}

root    <- file.path(default_data_root(), "Original")
cov_dir <- Sys.getenv("CEW_OUT", "covariates")   # output of cew_covariates_hpc.R

# =============================================================================
# 0. Site table -----------------------------------------------------------------
#    Every join (counts, Daymet, CDL, EVI) uses this location2.
# =============================================================================
counts_raw <- read.csv(file.path(root, "pest-occurrence/corn_earworm/doi_10_5061_dryad_m0cfxpp5x__v20221206/Lawton_etal_data.csv")) %>%
  filter(longitude >= -105)

sites <- make_site_table(counts_raw)
# The covariates are keyed by the same table; stop if they have drifted
sites_cov <- read.csv(file.path(cov_dir, "cew_sites.csv"))
stopifnot(identical(sites_cov$location2, sites$location2),
          isTRUE(all.equal(sites_cov$longitude, sites$longitude)),
          isTRUE(all.equal(sites_cov$latitude, sites$latitude)))

# =============================================================================
# 1. Trap counts on ISO weeks ---------------------------------------------------
#    ISO year/week from the date, so counts and Daymet weeks align. Site-years
#    with >= 5 sampled weeks, at sites with >= 3 such years.
# =============================================================================
counts <- counts_raw %>%
  mutate(date = as.Date(substr(date, 1, 10)),
         year = isoyear(date), woy = isoweek(date),
         CEW  = suppressWarnings(as.numeric(CEW_sum))) %>%
  filter(!is.na(CEW), year %in% yrs, woy <= 52) %>%
  left_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude")) %>%
  group_by(location2, year, woy) %>%
  summarise(CEW = mean(CEW), trap_type = first(trap_type), .groups = "drop")

keep_ly <- counts %>% count(location2, year, name = "n_wk") %>% filter(n_wk >= 5)
counts  <- counts %>% semi_join(keep_ly, by = c("location2", "year")) %>%
  group_by(location2) %>% filter(n_distinct(year) >= 3) %>% ungroup()

# =============================================================================
# 2. Weekly Daymet GDD ------------------------------------------------------------
#    covariates/daymet_<var>_<year>.rds: all days of each extracted year
#    (the previous year too, for accumulated GDD from Jan 1).
#    gdd_wk = weekly GDD sum (step 1); the anomalies from the site x week
#    climatology are step 2's local weather controls.
# =============================================================================
read_daymet <- function(yr) {
  f <- file.path(cov_dir, sprintf("daymet_%s_%d.rds", c("tmax", "tmin"), yr))
  if (!all(file.exists(f))) stop("Missing Daymet extraction: ", paste(f[!file.exists(f)], collapse = ", "))
  Reduce(function(a, b) left_join(a, b, by = c("location2", "date")), lapply(f, readRDS))
}

weekly_clim <- lapply((min(yrs) - 1):max(yrs), read_daymet) %>% bind_rows() %>%
  arrange(location2, date) %>%
  mutate(gdd_day = pmax((tmax + tmin) / 2 - t_base, 0), cal_year = year(date)) %>%
  group_by(location2, cal_year) %>% mutate(agdd = cumsum(gdd_day)) %>% ungroup() %>%
  mutate(year = isoyear(date), woy = isoweek(date)) %>%
  filter(year %in% yrs, woy <= 52) %>%
  group_by(location2, year, woy) %>%
  summarise(gdd_wk = sum(gdd_day), agdd_wk = last(agdd), n_days = n(), .groups = "drop") %>%
  # Drop weeks cut off where a site's extracted years end (sums would be short);
  # Daymet's missing Dec 31 in leap years shortens that week for every site alike
  group_by(year, woy) %>% filter(n_days == max(n_days)) %>% ungroup() %>% select(-n_days) %>%
  group_by(location2, woy) %>%
  mutate(across(c(gdd_wk, agdd_wk), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_anom")) %>%
  ungroup()

# =============================================================================
# 3. CDL crop shares and EVI ------------------------------------------------------
#    covariates/cdl_<year>.rds: <Crop>_pland = % of a 1 km radius.
#    EVI (MCD43): gaps linearly interpolated, then the site x week mean
#    (EVI_filled in cew_fourier_covariates.R); evi_anom = departure from it.
# =============================================================================
landscape <- bind_rows(lapply(yrs, function(yr) readRDS(file.path(cov_dir, sprintf("cdl_%d.rds", yr))))) %>%
  # a class absent around every site in one year is a missing column in that
  # year's file; across years it is 0, not NA
  mutate(across(-c(location2, year), ~ coalesce(.x, 0)))

# EVI points carry the site ids of cew_fourier_covariates.R (groups of the
# full-timeseries file sorted by lon/lat) - re-key them by coordinates
evi_key <- read.csv(file.path(root, "pest-occurrence/corn_earworm/cew_landscape_metrics_full_timeseries.csv")) %>%
  distinct(longitude, latitude) %>% arrange(longitude, latitude) %>% mutate(id = row_number()) %>%
  inner_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude"))

evi <- read.csv(file.path(root, "pest-occurrence/corn_earworm/MCD43_EVI_points.csv")) %>%
  inner_join(evi_key %>% select(id, location2), by = "id") %>%
  mutate(date = as.Date(date), year = isoyear(date), woy = isoweek(date)) %>%
  filter(woy <= 52) %>%
  group_by(location2, year, woy) %>% summarise(EVI = mean(EVI, na.rm = TRUE), .groups = "drop") %>%
  arrange(location2, year, woy) %>% group_by(location2) %>%
  mutate(EVI = na.approx(EVI, na.rm = FALSE)) %>%
  group_by(location2, woy) %>%
  mutate(evi_clim = mean(EVI, na.rm = TRUE), evi_wk = coalesce(EVI, evi_clim),
         evi_anom = evi_wk - evi_clim) %>% ungroup()

# =============================================================================
# 4. Modelling data ---------------------------------------------------------------
#    Observed site-weeks plus the prediction grid: the unsampled weeks 1-52 of
#    every monitored site-year (Daymet/CDL exist for trapped site-years only).
#    Only observed rows enter the fit. Rows missing a covariate (e.g. sites
#    without EVI) are dropped.
# =============================================================================
# Model covariates (names as in cew_fourier_covariates.R) = source columns.
# gdd there is GDD on the trap date; here the weekly sum.
my_vars <- c(EVI_filled = "evi_wk", gdd = "gdd_wk",
             corn = "Corn_pland", cotton = "Cotton_pland", soybeans = "Soybeans_pland",
             tobacco = "Tobacco_pland", peanuts = "Peanuts_pland", other = "Other_Crops_pland")
anom_vars <- c("gdd_wk_anom", "agdd_wk_anom", "evi_anom")   # step-2 local weather controls

grid_rows <- counts %>%
  group_by(location2, year) %>%
  reframe(woy = 1:52, trap_type = first(trap_type)) %>%
  anti_join(counts, by = c("location2", "year", "woy"))

d_all <- bind_rows(counts %>% mutate(observed = TRUE),
                   grid_rows %>% mutate(observed = FALSE, CEW = NA_real_)) %>%
  left_join(sites, by = "location2") %>%
  left_join(weekly_clim, by = c("location2", "year", "woy")) %>%
  left_join(landscape, by = c("location2", "year")) %>%
  left_join(evi %>% select(location2, year, woy, evi_wk, evi_anom), by = c("location2", "year", "woy")) %>%
  mutate(log_CEW   = log(CEW + 1),
         zone_fac  = factor(zone_30y),
         trap_type = factor(trap_type),
         year_f    = factor(year),
         loc_year  = paste(location2, year, sep = "_")) %>%   # character: site-year key
  add_fourier(H = 2) %>%
  filter(if_all(all_of(c(unname(my_vars), anom_vars)), ~ !is.na(.x))) %>%
  arrange(location2, year, woy)
print(count(d_all, zone_30y, observed))

# Standardise on the observed rows (anomalies: scale only)
scale_my <- make_scaler(d_all %>% filter(observed), my_vars)
d_all <- d_all %>% scale_my() %>%
  mutate(across(all_of(anom_vars), ~ .x / sd(.x[observed])))

# =============================================================================
# 5. Seasonal basis -----------------------------------------------------------
#    The 2-harmonic Fourier basis (s1, c1, s2, c2) centred and orthogonalised
#    over the observed rows, with unit sd: f = B_c R^-1 sqrt(n - 1), where
#    B_c = Q R. Over the trapping season (weeks ~20-38) the raw columns barely
#    vary and sit far from 0 (cos(2 pi 29/52) = -0.94), so each coefficient is
#    nearly collinear with the zone intercept, which NUTS samples poorly. The
#    seasonal curves - and peaks - are the same functions of week either way.
# =============================================================================
four  <- c("s1", "c1", "s2", "c2")
B_obs <- as.matrix(d_all[d_all$observed, four])
B_mu  <- colMeans(B_obs)
B_R   <- qr.R(qr(sweep(B_obs, 2, B_mu))) / sqrt(nrow(B_obs) - 1)
d_all[paste0("f", 1:4)] <- sweep(as.matrix(d_all[, four]), 2, B_mu) %*% solve(B_R)

d_obs <- d_all %>% filter(observed)
