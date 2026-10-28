# =============================================================================
# cew_stage1_data.R
# Sections 0-4 of cew_source_sink_pipeline.R (site table, counts, Daymet, CDL,
# EVI, modelling data) and the stage-1 model (section 5's formula, priors and
# brm call). Sourced by the pipeline and by cew_stage1_hpc.R, so the fit made
# on the HPC and the pipeline use identical data.
# =============================================================================

library(dplyr); library(tidyr); library(lubridate); library(zoo)
library(brms); library(posterior); library(ggplot2)
options(brms.backend = "cmdstanr")
source("cew_source_sink_functions.R")

yrs       <- 2008:2020       # CDL coverage + dense Southern trapping
yr_test   <- 2020            # held out: fit and select on the other years
yrs_train <- setdiff(yrs, yr_test)
t_base    <- 10              # GDD base temperature (deg C) - your choice
fill_sources <- TRUE         # fill unsampled Southern weeks with stage-1 predictions
root    <- file.path(default_data_root(), "Original")
cov_dir <- Sys.getenv("CEW_OUT", "covariates")   # output of cew_covariates_hpc.R

# =============================================================================
# 0. ONE site table -----------------------------------------------------------
#    Every join (counts, Daymet, CDL, EVI) uses this location2.
# =============================================================================
counts_raw <- read.csv(file.path(root, "pest-occurrence/corn_earworm/doi_10_5061_dryad_m0cfxpp5x__v20221206/Lawton_etal_data.csv")) %>%
  filter(longitude >= -105)

sites <- make_site_table(counts_raw)
# The HPC covariates are keyed by the same table; stop if they have drifted
sites_hpc <- read.csv(file.path(cov_dir, "cew_sites.csv"))
stopifnot(identical(sites_hpc$location2, sites$location2),
          isTRUE(all.equal(sites_hpc$longitude, sites$longitude)),
          isTRUE(all.equal(sites_hpc$latitude, sites$latitude)))
write.csv(sites, "cew_sites.csv", row.names = FALSE)

# =============================================================================
# 1. Trap counts on ISO weeks ---------------------------------------------------
#    Recompute ISO year/week from the date so counts and Daymet weeks align.
# =============================================================================
counts <- counts_raw %>%
  mutate(date = as.Date(substr(date, 1, 10)),
         year = isoyear(date), woy = isoweek(date),
         CEW  = suppressWarnings(as.numeric(CEW_sum))) %>%
  filter(!is.na(CEW), year %in% yrs, woy <= 52) %>%
  left_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude")) %>%
  group_by(location2, year, woy) %>%
  summarise(CEW = mean(CEW), trap_type = first(trap_type), .groups = "drop")

# Sampling filter. NOTE: the earlier ">= 50 moths per location-year" rule selects
# on the outcome and removes low years - exactly the variation stage 2 needs.
keep_ly <- counts %>% count(location2, year, name = "n_wk") %>% filter(n_wk >= 5)
counts  <- counts %>% semi_join(keep_ly, by = c("location2", "year")) %>%
  group_by(location2) %>% filter(n_distinct(year) >= 3) %>% ungroup()

# =============================================================================
# 2. Daymet climate covariates ------------------------------------------------
#    All days (not just trap dates) of each extracted year. With CEW_SITES=active
#    (default) a site is extracted only around the years it was trapped (its
#    trapped years +/- 1), so its climatology is the mean over those years; with
#    CEW_SITES=all every site has every year. covariates/daymet_<var>_<year>.rds.
# =============================================================================
read_daymet <- function(yr) {                               # previous year for winter
  f <- file.path(cov_dir, sprintf("daymet_%s_%d.rds", c("tmax", "tmin", "prcp"), yr))
  if (!all(file.exists(f))) stop("Missing Daymet extraction: ", paste(f[!file.exists(f)], collapse = ", "))
  Reduce(function(a, b) left_join(a, b, by = c("location2", "date")), lapply(f, readRDS))
}

daily <- lapply((min(yrs) - 1):max(yrs), read_daymet) %>% bind_rows() %>%
  arrange(location2, date) %>%
  mutate(gdd_day = pmax((tmax + tmin) / 2 - t_base, 0),
         cal_year = year(date)) %>%
  group_by(location2, cal_year) %>% mutate(agdd = cumsum(gdd_day)) %>% ungroup()

weekly_clim <- daily %>%
  mutate(year = isoyear(date), woy = isoweek(date)) %>%
  filter(year %in% yrs, woy <= 52) %>%
  group_by(location2, year, woy) %>%
  summarise(gdd_wk  = sum(gdd_day), tmin_wk = mean(tmin), prcp_wk = sum(prcp),
            agdd_wk = last(agdd),            # accumulated GDD at end of week
            n_days  = n(), .groups = "drop") %>%
  # Drop weeks cut off where a site's extracted years end (sums would be short);
  # Daymet's missing Dec 31 in leap years shortens that week for every site alike
  group_by(year, woy) %>% filter(n_days == max(n_days)) %>% ungroup() %>% select(-n_days) %>%
  group_by(location2, year) %>% arrange(woy, .by_group = TRUE) %>%
  mutate(gdd_wk_lag1 = lag(gdd_wk)) %>% ungroup()

# Climatology (site x week mean over yrs) and anomalies. Stage 1 uses the raw
# weekly GDD (as in cew_fourier_covariates.R); its climatology defines the
# expectation in section 6. The anomalies are the local weather controls of
# stage 2. (Covariates only - no counts - so using 2020 weather here is no leak.)
weekly_clim <- weekly_clim %>%
  group_by(location2, woy) %>%
  mutate(across(c(gdd_wk, gdd_wk_lag1, agdd_wk, tmin_wk, prcp_wk),
                ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_anom")) %>%
  ungroup() %>%
  mutate(gdd_wk_clim = gdd_wk - gdd_wk_anom)

# Winter severity preceding each season (Dec of t-1 to Feb of t), a proxy for
# overwintering success; matters most for Southern/Transitional sources.
winter <- daily %>%
  mutate(season = ifelse(month(date) == 12, cal_year + 1, cal_year)) %>%
  filter(month(date) %in% c(12, 1, 2), season %in% yrs) %>%
  group_by(location2, year = season) %>%
  summarise(winter_tmin_min = min(tmin), winter_cold_days = sum(tmin < -10),
            n_days = n(), .groups = "drop") %>%
  filter(n_days >= 89) %>% select(-n_days) %>%    # complete Dec-Feb only (89-90 Daymet days)
  group_by(location2) %>%
  mutate(across(c(winter_tmin_min, winter_cold_days), ~ .x - mean(.x), .names = "{.col}_anom")) %>%
  ungroup()

# =============================================================================
# 3. CDL landscape covariates and EVI -------------------------------------------
#    Extracted on the HPC (covariates/cdl_<year>.rds) with the corrected
#    reclassification in cew_covariates_hpc.R. Columns per crop:
#      <Crop>_pland (1 km radius), <Crop>_ed (1 km), <Crop>_pland16k (16.1 km).
#    Each year covers only that year's trapped sites (CEW_SITES=active).
# =============================================================================
landscape <- bind_rows(lapply(yrs, function(yr) readRDS(file.path(cov_dir, sprintf("cdl_%d.rds", yr))))) %>%
  # a class absent around every site in one year is a missing column in that
  # year's file; across years it is 0, not NA
  mutate(across(-c(location2, year), ~ coalesce(.x, 0)))

# EVI points were extracted with the site ids of cew_fourier_covariates.R
# (location2 there = groups of the full-timeseries file sorted by lon/lat),
# NOT this table's location2 - re-key them by coordinates.
evi_key <- read.csv(file.path(root, "pest-occurrence/corn_earworm/cew_landscape_metrics_full_timeseries.csv")) %>%
  distinct(longitude, latitude) %>% arrange(longitude, latitude) %>% mutate(id = row_number()) %>%
  inner_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude"))

# Gaps: linear interpolation, then the site x week mean (as EVI_filled in
# cew_fourier_covariates.R); evi_anom = departure from that mean
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
# 4. Assemble modelling data --------------------------------------------------
#    Observed site-weeks plus a prediction grid: the unsampled weeks 1-52 of
#    every trapped site-year (only trapped site-years have Daymet/CDL with
#    CEW_SITES=active). in_season marks weeks between the site-year's first
#    and last sampled week; the Southern source field fills only those.
#    The grid rows never enter the fit (d_train = observed rows).
# =============================================================================
crops <- c("Corn_pland", "Cotton_pland", "Soybeans_pland", "Peanuts_pland",
           "Tobacco_pland", "Other_Crops_pland")             # Non_Crop = reference
anom_vars <- c("gdd_wk_anom", "agdd_wk_anom", "evi_anom")   # stage-2 local weather controls

# Stage-1 covariates under their names in cew_fourier_covariates.R = source
# column here. gdd there is GDD on the trap date; here the weekly sum.
my_vars <- c(EVI_filled = "evi_wk", gdd = "gdd_wk",
             corn = "Corn_pland", cotton = "Cotton_pland", soybeans = "Soybeans_pland",
             tobacco = "Tobacco_pland", peanuts = "Peanuts_pland", other = "Other_Crops_pland",
             lon_scaled = "longitude", lat_scaled = "latitude")

grid_rows <- counts %>%
  group_by(location2, year) %>%
  reframe(in_season = 1:52 >= min(woy) & 1:52 <= max(woy), woy = 1:52,
          trap_type = first(trap_type)) %>%
  anti_join(counts, by = c("location2", "year", "woy"))

d_all <- bind_rows(counts %>% mutate(observed = TRUE, in_season = TRUE),
                   grid_rows %>% mutate(observed = FALSE, CEW = NA_real_)) %>%
  left_join(sites, by = "location2") %>%
  left_join(weekly_clim, by = c("location2", "year", "woy")) %>%
  left_join(winter, by = c("location2", "year")) %>%
  left_join(landscape, by = c("location2", "year")) %>%
  left_join(evi %>% select(location2, year, woy, evi_wk, evi_clim, evi_anom),
            by = c("location2", "year", "woy")) %>%
  mutate(log_CEW   = log(CEW + 1),
         train     = year %in% yrs_train,
         zone_fac  = factor(zone_30y),
         trap_type = factor(trap_type),
         year_f    = factor(year)) %>%
  add_fourier(H = 2) %>%
  filter(if_all(all_of(c(unname(my_vars), anom_vars)), ~ !is.na(.x)))   # e.g. sites without EVI
print(count(d_all, observed, zone_30y, train))

# Standardise with the training years' means/sds (anomalies: scale only)
scale_my <- make_scaler(d_all %>% filter(observed, train), my_vars)
d_all <- d_all %>% scale_my() %>%
  mutate(across(all_of(anom_vars), ~ .x / sd(.x[observed & train])))

d_train <- d_all %>% filter(observed, train)

# =============================================================================
# 5. Stage 1: seasonal model of cew_fourier_covariates.R (fit_gp2.1) ------------
#
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)
#  mu_itw = alpha + x_itw' beta (EVI, GDD, crops)
#           + sum_h [ u^s_zh sin(2 pi h w/52) + u^c_zh cos(2 pi h w/52) ]   (h = 1, 2)
#           + f(lon_i, lat_i) + a_t + e_trap
#  u_zh ~ N(0, tau_h^2) per zone, f ~ GP (exponential kernel),
#  a_t ~ N(0, s_a^2), e_trap ~ N(0, s_e^2)
#  Fit to the training years only.
# =============================================================================
f_my <- bf(log_CEW ~
             EVI_filled + gdd +
             corn + cotton + soybeans + tobacco + peanuts + other +
             (0 + s1 + c1 + s2 + c2 || zone_fac) +
             gp(lon_scaled, lat_scaled, cov = "exponential") +
             (1 | year_f) +
             (1 | trap_type))

pri_my <- c(prior(normal(0, 1), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(student_t(3, 0, 2.5), class = "sd"),
            prior(exponential(1), class = "sdgp"),
            prior(inv_gamma(5, 5), class = "lscale", coef = "gplon_scaledlat_scaled"))

# Sampler defaults = the long run (20000 iterations; the first 3100 took ~13 h
# on a laptop, hence the HPC: hpc/submit_stage1.sh). The pipeline only reads
# the saved fit (file_refit = "never"); delete the .rds to refit.
stage1_file <- "fit1_gp21_train"
fit_stage1 <- function(data, chains = 4, iter = 20000, warmup = 19000, thin = 4,
                       file = stage1_file, file_refit = "never", ...) {
  brm(f_my, data = data, family = gaussian(), prior = pri_my,
      chains = chains, cores = chains, iter = iter, warmup = warmup, thin = thin,
      seed = 42, control = list(adapt_delta = 0.95),
      file = file, file_refit = file_refit, ...)
}
