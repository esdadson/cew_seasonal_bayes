# =============================================================================
# cew_source_sink_pipeline.R
# H. zea: seasonal curves by overwintering zone (stage 1) and
#          Southern -> Transitional/Northern source-sink dynamics (stage 2)
#
# Sections
#   0. Setup and ONE site table (location2 defined once)
#   1. Trap counts on ISO weeks
#   2. Daymet climate covariates (daily -> ISO week, climatology, anomalies)
#   3. CDL landscape covariates (reclassification fixes, multiple radii)
#   4. Assemble modelling data
#   5. Stage 1: Gaussian Fourier model of log(CEW + 1)  -> peak timing
#   6. Climatological expectation and anomalies
#   7. Stage 2: source pressure, kernel/lag screening, Bayesian sink model
#   8. Falsification tests
#   9. Outputs for "when" and "where"
#
# Sections 2-3 read the raster extractions made on the HPC by
# cew_covariates_hpc.R (see hpc/submit_covariates.sh); copy its output
# directory here as covariates/. Sections 5-7 (brms) could not be run here
# (no Stan). The helper functions in cew_source_sink_functions.R were tested
# on cew_landscape_vars.csv.
# =============================================================================

library(dplyr); library(tidyr); library(lubridate); library(zoo)
library(sf)
library(brms); library(posterior); library(ggplot2)
options(brms.backend = "cmdstanr")
source("cew_source_sink_functions.R")

yrs    <- 2008:2020          # CDL coverage + dense Southern trapping
t_base <- 10                 # GDD base temperature (deg C) - your choice
root   <- file.path(default_data_root(), "Original")
cov_dir <- Sys.getenv("CEW_OUT", "covariates")   # output of cew_covariates_hpc.R

# =============================================================================
# 0. ONE site table -----------------------------------------------------------
#    Every join (counts, Daymet, CDL, EVI) uses this location2. If your EVI
#    extraction used an older location2 definition, re-key it to this table
#    by coordinates (or re-extract) before joining.
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
#    ALL days for ALL sites and years (not just trap dates), so each site has
#    a full climatology even in years it wasn't trapped. Extracted on the HPC:
#    covariates/daymet_<var>_<year>.rds, one file per variable and year.
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
            .groups = "drop") %>%
  group_by(location2, year) %>% arrange(woy, .by_group = TRUE) %>%
  mutate(gdd_wk_lag1 = lag(gdd_wk)) %>% ungroup()

# Climatology (site x week mean over yrs) and anomalies.
# Climatology is collinear with the Fourier + site terms, so only anomalies
# enter the model: the Fourier curve then IS the climatological seasonal curve.
weekly_clim <- weekly_clim %>%
  group_by(location2, woy) %>%
  mutate(across(c(gdd_wk, gdd_wk_lag1, agdd_wk, tmin_wk, prcp_wk),
                ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_anom")) %>%
  ungroup()

# Winter severity preceding each season (Dec of t-1 to Feb of t), a proxy for
# overwintering success; matters most for Southern/Transitional sources.
winter <- daily %>%
  mutate(season = ifelse(month(date) == 12, cal_year + 1, cal_year)) %>%
  filter(month(date) %in% c(12, 1, 2), season %in% yrs) %>%
  group_by(location2, year = season) %>%
  summarise(winter_tmin_min = min(tmin), winter_cold_days = sum(tmin < -10), .groups = "drop") %>%
  group_by(location2) %>%
  mutate(across(c(winter_tmin_min, winter_cold_days), ~ .x - mean(.x), .names = "{.col}_anom")) %>%
  ungroup()

# =============================================================================
# 3. CDL landscape covariates -------------------------------------------------
#    Extracted on the HPC (covariates/cdl_<year>.rds) with the corrected
#    reclassification in cew_covariates_hpc.R. Columns per crop:
#      <Crop>_pland (1 km radius), <Crop>_ed (1 km), <Crop>_pland16k (16.1 km).
#    Studies summarised in your 2024 review used 1.61 and 16.1 km; choose the
#    radius per crop by LOO in stage 1 (or keep both if not collinear).
# =============================================================================
landscape <- bind_rows(lapply(yrs, function(yr) readRDS(file.path(cov_dir, sprintf("cdl_%d.rds", yr)))))

# EVI: same anomaly treatment (site x week climatology removed)
evi <- read.csv(file.path(root, "pest-occurrence/corn_earworm/MCD43_EVI_points.csv")) %>%
  rename(location2 = id) %>%                             # must match the site table!
  mutate(date = as.Date(date), year = isoyear(date), woy = isoweek(date)) %>%
  group_by(location2, year, woy) %>% summarise(EVI = mean(EVI, na.rm = TRUE), .groups = "drop") %>%
  arrange(location2, year, woy) %>% group_by(location2) %>%
  mutate(EVI = na.approx(EVI, na.rm = FALSE)) %>%
  group_by(location2, woy) %>% mutate(evi_anom = EVI - mean(EVI, na.rm = TRUE)) %>% ungroup()

# =============================================================================
# 4. Assemble modelling data --------------------------------------------------
# =============================================================================
crops <- c("Corn_pland", "Cotton_pland", "Soybeans_pland", "Peanuts_pland",
           "Tobacco_pland", "Other_Crops_pland")             # Non_Crop = reference
anom_vars <- c("gdd_wk_anom", "agdd_wk_anom", "evi_anom")   # weekly anomalies in the model

d_fit <- counts %>%
  left_join(sites, by = "location2") %>%
  left_join(weekly_clim, by = c("location2", "year", "woy")) %>%
  left_join(winter, by = c("location2", "year")) %>%
  left_join(landscape, by = c("location2", "year")) %>%
  left_join(evi %>% select(location2, year, woy, evi_anom), by = c("location2", "year", "woy")) %>%
  mutate(log_CEW   = log(CEW + 1),
         zone_fac  = factor(zone_30y),
         trap_type = factor(trap_type),
         year_f    = factor(year),
         loc_year  = factor(paste(location2, year, sep = "_")),
         location2 = factor(location2),
         evi_anom  = coalesce(evi_anom, 0)) %>%              # missing EVI -> climatology
  add_fourier(H = 3) %>%
  mutate(across(all_of(crops), ~ as.numeric(scale(.x)), .names = "{.col}_z"),
         across(all_of(anom_vars), ~ .x / sd(.x, na.rm = TRUE)),   # anomalies: scale only
         gx = x_km / 1000, gy = y_km / 1000) %>%
  filter(if_all(all_of(anom_vars), ~ !is.na(.x)))

# =============================================================================
# 5. Stage 1: seasonal curves ---------------------------------------------------
#
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)
#  mu_itw = alpha_z + sum_h [ (beta^s_zh + u^s_ih) sin(2 pi h w/52)
#                           + (beta^c_zh + u^c_ih) cos(2 pi h w/52) ]
#           + x_it' gamma (landscape) + c_itw' delta (weather anomalies)
#           + eta_trap + f(s_i) + a_t + b_it
#  u_ih ~ N(0, tau_h^2), f ~ GP(0, k(.)), a_t ~ N(0, s_a^2), b_it ~ N(0, s_b^2)
# =============================================================================
crop_z <- paste0(crops, "_z")

f_H2 <- bf(as.formula(paste(
  "log_CEW ~ 0 + zone_fac + zone_fac:(s1 + c1 + s2 + c2) + trap_type +",
  paste(c(crop_z, anom_vars), collapse = " + "),
  "+ (0 + s1 + c1 + s2 + c2 || location2)",
  "+ gp(gx, gy, k = 25, c = 5/4, cov = 'matern32')",        # HSGP; see ?gp for supported cov
  "+ (1 | year_f) + (1 | loc_year)")))

f_H3 <- bf(as.formula(paste(
  "log_CEW ~ 0 + zone_fac + zone_fac:(s1 + c1 + s2 + c2 + s3 + c3) + trap_type +",
  paste(c(crop_z, anom_vars), collapse = " + "),
  "+ (0 + s1 + c1 + s2 + c2 + s3 + c3 || location2)",
  "+ gp(gx, gy, k = 25, c = 5/4, cov = 'matern32')",
  "+ (1 | year_f) + (1 | loc_year)")))

pri <- c(prior(normal(0, 2), class = "b"),
         prior(student_t(3, 0, 1), class = "sd"),
         prior(exponential(1), class = "sdgp"),
         prior(student_t(3, 0, 1), class = "sigma"))

fit_args <- list(data = d_fit, family = gaussian(), prior = pri, chains = 4, cores = 4,
                 iter = 2000, warmup = 1000, seed = 42, control = list(adapt_delta = 0.95))
fit1_H2 <- do.call(brm, c(list(formula = f_H2, file = "fit1_H2"), fit_args))
fit1_H3 <- do.call(brm, c(list(formula = f_H3, file = "fit1_H3"), fit_args))
loo_compare(loo(fit1_H2), loo(fit1_H3))   # Southern bimodality may need H = 3
fit1 <- fit1_H2                            # <- set to the preferred model
hz <- c("s1", "c1", "s2", "c2")            # add "s3","c3" if H = 3

# ---- Peak timing -------------------------------------------------------------
# Weather enters only as anomalies, so the Fourier curve alone is the
# climatological seasonal curve: zone coefficients + site deviations.
fx <- as_draws_df(fit1)
re <- ranef(fit1, summary = FALSE)$location2[, , hz]      # draws x sites x coefs
site_ids <- dimnames(re)[[2]]
site_z   <- sites$zone_30y[match(as.integer(site_ids), sites$location2)]

stopifnot(all(paste0("b_zone_fac", rep(0:2, each = length(hz)), ":", hz) %in% names(fx)))
zone_arr <- array(NA, c(dim(re)[1], 3, length(hz)), dimnames = list(NULL, c("0", "1", "2"), hz))
for (z in 0:2) zone_arr[, as.character(z), ] <- as.matrix(fx[, paste0("b_zone_fac", z, ":", hz)])
site_arr <- re
for (z in 0:2) {
  idx <- which(site_z == z)
  site_arr[, idx, ] <- re[, idx, , drop = FALSE] +
    aperm(array(zone_arr[, as.character(z), ], c(dim(re)[1], length(hz), length(idx))), c(1, 3, 2))
}

# Search only the sampled season (5th-95th percentile of sampled weeks);
# outside it the curve is extrapolation.
win <- range(quantile(d_fit$woy, c(0.05, 0.95)))
pk_zone <- peak_from_draws(zone_arr, window = win[1]:win[2])
pk_site <- peak_from_draws(site_arr, window = win[1]:win[2])

lag_N <- pk_zone$peak[, "0"] - pk_zone$peak[, "2"]
lag_T <- pk_zone$peak[, "1"] - pk_zone$peak[, "2"]
rbind(N_minus_S = quantile(lag_N, c(.025, .5, .975)), T_minus_S = quantile(lag_T, c(.025, .5, .975)))
c(P_N_later = mean(lag_N > 0), P_T_later = mean(lag_T > 0))

peak_map <- data.frame(location2 = as.integer(site_ids),
                       peak_med = apply(pk_site$peak, 2, median),
                       peak_lo  = apply(pk_site$peak, 2, quantile, 0.05),
                       peak_hi  = apply(pk_site$peak, 2, quantile, 0.95),
                       at_edge  = pk_site$at_edge) %>%
  left_join(sites, by = "location2")
# Check at_edge and bimodal peak posteriors (esp. Southern sites) before mapping.

# =============================================================================
# 6. Climatological expectation and anomalies ---------------------------------
#    E[y] with weather anomalies = 0 and no year / location-year effects.
#    (re_formula controls group-level terms only; the GP is always included.)
# =============================================================================
nd <- d_fit %>% mutate(across(all_of(anom_vars), ~ 0))
mu_clim <- colMeans(posterior_linpred(
  fit1, newdata = nd, ndraws = 200,
  re_formula = as.formula(paste("~ (0 +", paste(hz, collapse = " + "), "|| location2)"))))
d_fit$anom <- d_fit$log_CEW - mu_clim

# Quick first look (Lawton-style): yearly zone mean anomalies
d_fit %>% group_by(zone_30y, year) %>% summarise(a = mean(anom), .groups = "drop") %>%
  pivot_wider(names_from = zone_30y, values_from = a) %>%
  summarise(r_S_T = cor(`2`, `1`), r_S_N = cor(`2`, `0`))

# =============================================================================
# 7. Stage 2: source pressure -------------------------------------------------
# =============================================================================
dd   <- d_fit %>% mutate(location2 = as.integer(as.character(location2)))
tw   <- make_tw_index(dd$year)
src  <- sites %>% filter(location2 %in% dd$location2, zone_30y == 2)
snk  <- sites %>% filter(location2 %in% dd$location2, zone_30y != 2)
A    <- anomaly_matrix(dd %>% filter(zone_30y == 2), src$location2, tw)   # sources x year-week
geom <- pair_geometry(snk, src)
obs  <- dd %>% filter(zone_30y != 2)

base_terms <- anom_vars                     # local weather controls in screening

# ---- 7a. Distance decay (where) ----------------------------------------------
lambdas <- c(100, 200, 400, 800, 1600, Inf)
P_dist  <- setNames(lapply(lambdas, function(l) source_pressure(kernel_weights(geom, l), A, tw$woy)),
                    paste0("lambda_", lambdas))
scr_N <- grid_screen(obs %>% filter(zone_30y == 0), P_dist, tw, base_terms = base_terms)
scr_T <- grid_screen(obs %>% filter(zone_30y == 1), P_dist, tw, base_terms = base_terms)
scr_N; scr_T
lam_best <- 400                              # <- set from the screens

# ---- 7b. Direction (where) ---------------------------------------------------
dir_grid <- expand.grid(theta0 = c(0, 30, 60, 90), kappa = c(0, 1, 3))
P_dir <- setNames(lapply(seq_len(nrow(dir_grid)), function(g)
  source_pressure(kernel_weights(geom, lam_best, dir_grid$kappa[g], dir_grid$theta0[g]), A, tw$woy)),
  paste0("theta", dir_grid$theta0, "_k", dir_grid$kappa))
grid_screen(obs, P_dir, tw, base_terms = c("zone_fac", base_terms))

# ---- 7c. Travelling wave: lag grows with distance (when + speed) -------------
speeds <- c(100, 200, 400, 800)              # km per week
P_wave <- setNames(lapply(speeds, function(v)
  source_pressure(kernel_weights(geom, lam_best), A, tw$woy,
                  lag_mat = round(geom$D / v), max_lag = 8)), paste0("v_", speeds))
grid_screen(obs, P_wave, tw, lags = 0:4, df = 3, base_terms = c("zone_fac", base_terms))

# ---- 7d. Bayesian sink model at the chosen kernel (when) ---------------------
lags  <- 0:8
Blag  <- lag_basis(lags, df = 4)
P_best <- P_dist[[paste0("lambda_", lam_best)]]
obs_b <- attach_lags(obs, P_best, tw, lags = c(-3:-1, lags)) %>% add_Q(Blag) %>%
  filter(if_all(all_of(colnames(Blag)), ~ !is.na(.x))) %>%
  mutate(zone_fac = droplevels(zone_fac), location2 = factor(location2), year_f = factor(year))

# Own-site anomaly ~1 generation earlier: separates local recruitment from immigration
own <- anomaly_matrix(obs, snk$location2, tw)
own4 <- shift_within_year(own, 4, tw$woy)
obs_b$own_lag4 <- own4[cbind(match(as.integer(as.character(obs_b$location2)), snk$location2), obs_b$tw)]
obs_b <- obs_b %>% mutate(own_lag4_miss = as.integer(is.na(own_lag4)), own_lag4 = coalesce(own_lag4, 0))

Qs <- colnames(Blag)
f2 <- bf(as.formula(paste(
  "anom ~ 0 + zone_fac + zone_fac:(", paste(Qs, collapse = " + "), ") +",
  paste(anom_vars, collapse = " + "), "+ own_lag4 + own_lag4_miss",
  "+ (1 | year_f) + (1 | location2)")))
# (1 | year_f) => lag effects identified from within-season and between-site
# variation only (conservative). Refit without it to include year-level
# Southern booms; report both.
fit2 <- brm(f2, data = obs_b, family = gaussian(), chains = 4, cores = 4,
            prior = c(prior(normal(0, 0.5), class = "b"), prior(student_t(3, 0, 1), class = "sd")),
            seed = 42, file = "fit2_sink")

# Lag profile beta_k = Blag %*% theta, per draw and sink zone
fx2 <- as_draws_df(fit2)
lag_profile <- lapply(c("0", "1"), function(z) {
  th <- as.matrix(fx2[, paste0("b_zone_fac", z, ":", Qs)])
  bk <- th %*% t(Blag)                                     # draws x lags
  data.frame(zone = z, lag = lags, med = apply(bk, 2, median),
             lo = apply(bk, 2, quantile, .025), hi = apply(bk, 2, quantile, .975),
             p_pos = colMeans(bk > 0))
}) %>% bind_rows()

# When in the season: let the combined pressure effect vary by week
theta_hat <- sapply(c("0", "1"), function(z) colMeans(as.matrix(fx2[, paste0("b_zone_fac", z, ":", Qs)])))
obs_b$Pstar <- rowSums(as.matrix(obs_b[, Qs]) * t(theta_hat[, as.character(obs_b$zone_fac)]))
# (two-step: Pstar uses fit2's posterior means, so treat fit2_season / fit2_site as descriptive)
fit2_season <- brm(bf(anom ~ zone_fac + gdd_wk_anom + agdd_wk_anom + evi_anom +
                        s(woy, by = Pstar, k = 6) + (1 | year_f) + (1 | location2)),
                   data = obs_b, family = gaussian(), chains = 4, cores = 4, seed = 42,
                   file = "fit2_season")

# Where on the sink side: site-specific sensitivity to source pressure
fit2_site <- brm(bf(anom ~ zone_fac + gdd_wk_anom + agdd_wk_anom + evi_anom + Pstar +
                      (1 + Pstar || location2) + (1 | year_f)),
                 data = obs_b, family = gaussian(), chains = 4, cores = 4, seed = 42,
                 file = "fit2_site")

# =============================================================================
# 8. Falsification tests ------------------------------------------------------
# =============================================================================
# (i) Leads: future Southern anomalies should not predict sinks
loyo_gain(obs_b, c("P_lead1", "P_lead2", "P_lead3"), c("zone_fac", anom_vars))

# (ii) Flat kernel (lambda = Inf) should do worse than the chosen kernel: see 7a

# (iii) Year-permutation null for the CV gain
obs_raw <- obs                                              # before lag/Q columns
real_gain <- grid_screen(obs_raw, list(real = P_best), tw, base_terms = c("zone_fac", anom_vars))$gain
null_gain <- replicate(200, {
  Pp <- source_pressure(kernel_weights(geom, lam_best), permute_source_years(A, tw), tw$woy)
  grid_screen(obs_raw, list(perm = Pp), tw, base_terms = c("zone_fac", anom_vars))$gain
})
mean(null_gain >= real_gain)                                 # permutation p-value

# (iv) Reverse direction: Transitional/Northern anomalies -> Southern sinks
A_rev    <- anomaly_matrix(obs, snk$location2, tw)
geom_rev <- pair_geometry(src, snk)
P_rev    <- source_pressure(kernel_weights(geom_rev, lam_best), A_rev, tw$woy)
grid_screen(dd %>% filter(zone_30y == 2), list(reverse = P_rev), tw, base_terms = anom_vars)

# =============================================================================
# 9. Outputs ------------------------------------------------------------------
# =============================================================================
# When: lag profile by sink zone
ggplot(lag_profile, aes(lag, med, ymin = lo, ymax = hi, colour = zone)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_pointrange(position = position_dodge(0.3)) +
  labs(x = "Lag (weeks) behind Southern anomaly", y = "Effect on sink log-count anomaly") +
  theme_minimal()
conditional_effects(fit2_season, effects = "woy:Pstar")     # when in the season

# Where (sink side): site multipliers
site_sens <- coef(fit2_site)$location2[, , "Pstar"] %>% as.data.frame() %>%
  mutate(location2 = as.integer(rownames(.))) %>% left_join(sites, by = "location2")

# Where (source side): average share of sink pressure contributed by each source
Wb <- kernel_weights(geom, lam_best); Wn <- Wb / rowSums(Wb)
src_share <- data.frame(location2 = src$location2, share = colMeans(Wn),
                        anom_sd = apply(A, 1, sd, na.rm = TRUE)) %>%
  mutate(influence = share * anom_sd) %>% left_join(sites, by = "location2")
