# =============================================================================
# cew_peak_granger.R
# H. zea weekly log(CEW + 1), 2008-2020
#   Step 1: peak timing - the week of peak expected log(CEW + 1) for every
#           site-year, from every posterior draw of the weekly spatial model
#           (cew_fourier_covariates.R, restructured in cew_peak_model.R),
#           summarised by zone and year.
#   Step 2: do Southern populations (zone_30y == 2) lead Transitional (1) and
#           Northern (0) ones? The step-1 filled weekly populations (observed
#           weeks + model fills of unsampled weeks) drive regional clustering
#           and Granger causality: Northern/Transitional populations modelled
#           from lagged Southern ones at lags of 1-8 weeks, scored in-sample
#           (F test) and out of sample (leave-one-year-out, held-out years
#           filled by the step-1 fit without them), with a year-permutation null.
#   Replaces the Bayesian source-sink pressure model of
#   cew_source_sink_pipeline.R for step 2.
#
# Sections
#   0. Data and step-1 model (cew_peak_model.R)
#   1. Step-1 fit and peak draws (fit on the HPC: cew_peak_hpc.R)
#   2. Peak week per site-year, zone-year and site
#   3. Step 1 out of sample: leave-one-year-out (cew_peak_hpc.R loyo)
#   4. Filled weekly populations and regional clustering
#   5. Regional weekly series
#   6. Granger causality between zones (6b: step-1 uncertainty in the fills)
#   7. Granger causality between regions
#   8. Save outputs for CEW_JIPM_Figures.R
#
# Sections 4-8 need the step-1 fills (fill_<model>.rds, loyo/<model>/) to
# use filled populations; without them they fall back to observed counts
# only, so they can be tried while the step-1 fits are on the HPC. With the
# peak draws present, section 7 adds the step-1 peak delay of each region pair.
# =============================================================================

source("cew_peak_model.R")       # cew_stage1_data.R + model, peak and LOYO functions
source("cew_granger_functions.R")
zl <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")
as_int <- function(f) as.integer(as.character(f))
q95 <- function(x) c(med = median(x), lo = quantile(x, 0.025, names = FALSE),
                     hi = quantile(x, 0.975, names = FALSE))

# =============================================================================
# 1. Step-1 fit and peak draws ---------------------------------------------------
#    ./hpc/submit_peak.sh fit (hpc/README.md step 11) fits each candidate
#    model (cew_peak_model.R) on all years and writes peak_draws_<model>.rds
#    (peak week of every site-year in every draw) and loo_peak_<model>.rds,
#    which are all sections 2-3 need: copy them here. Given only a fit
#    (fit_peak_<model>.rds), the peaks are computed from it. To fit locally
#    instead (hours): fit_peak(d_obs, file = peak_files()$fit, threads = 2)
#    peak_use = the candidate whose peaks are summarised; choose it from the
#    comparison in section 3 (Sys.setenv(PEAK_MODEL = ...) before sourcing).
# =============================================================================
peak_use <- peak_model
pf <- peak_files(peak_use)
pk <- NULL
if (file.exists(pf$draws)) {
  pk <- readRDS(pf$draws)
} else if (file.exists(paste0(pf$fit, ".rds"))) {
  fit_pk <- fit_peak(d_obs, file = pf$fit)                   # reads the saved fit
  pk <- peak_draws(fit_pk, peak_grid(d_all))
  saveRDS(pk, pf$draws)
} else {
  message("No ", pf$draws, " or ", pf$fit, ".rds yet: skipping sections 1-2")
}

# =============================================================================
# 2. Peak week per site-year, zone-year and site --------------------------------
#    Every summary is computed draw by draw, then summarised, so its interval
#    carries the uncertainty of every site-year peak in it.
#    Zone-year = mean over the zone's trapped site-years that year: the set of
#    sites changes between years, so part of a year-to-year change can be
#    which sites were trapped.
# =============================================================================
if (!is.null(pk)) {
  P  <- pk$draws
  sy <- pk$site_years %>% mutate(location2 = as_int(location2), loc_year = as.character(loc_year))

  peak_sy <- summarise_peaks(list(draws = P, site_years = sy)) %>%
    left_join(sites %>% select(location2, location, longitude, latitude, x_km, y_km), by = "location2")
  # CI widths, and the share of site-years whose peak most likely lies within
  # their own sampled weeks (the rest are placed by pooling)
  print(peak_sy %>% group_by(zone = zl[as.character(zone_30y)]) %>%
          summarise(n = n(), median_ci_width = median(ci_width), share_in_sampled = mean(in_sampled > 0.5)))

  # Zone-year mean peak, draws x zone-years
  ZY <- vapply(split(seq_len(ncol(P)), paste(sy$zone_30y, sy$year)),
               function(j) rowMeans(P[, j, drop = FALSE]), numeric(nrow(P)))
  peak_zone_year <- data.frame(key = colnames(ZY), t(apply(ZY, 2, q95))) %>%
    tidyr::separate(key, c("zone_30y", "year"), convert = TRUE) %>%
    left_join(count(sy, zone_30y, year, name = "n_sites"), by = c("zone_30y", "year")) %>%
    mutate(zone = zl[as.character(zone_30y)])

  # Delay of the Transitional and Northern zones behind the Southern zone:
  # per year, and the mean over years (P_later = posterior P(delay > 0))
  zone_delay <- bind_rows(lapply(c(1, 0), function(z) {
    ys <- intersect(sy$year[sy$zone_30y == z], sy$year[sy$zone_30y == 2])
    Dz <- sapply(ys, function(y) ZY[, paste(z, y)] - ZY[, paste(2, y)])
    colnames(Dz) <- ys
    Dz <- cbind(Dz, all_years = rowMeans(Dz))
    data.frame(zone = zl[as.character(z)], year = colnames(Dz), t(apply(Dz, 2, q95)),
               P_later = colMeans(Dz > 0), row.names = NULL)
  }))
  print(zone_delay %>% filter(year == "all_years"))

  # Site peak, averaged over the site's trapped years (per draw): for maps
  peak_site <- data.frame(location2 = sy$location2, j = seq_len(ncol(P))) %>%
    group_by(location2) %>% summarise(j = list(j), n_years = n()) %>%
    mutate(s = lapply(j, function(jj) q95(rowMeans(P[, jj, drop = FALSE])))) %>%
    tidyr::unnest_wider(s) %>% select(-j) %>%
    left_join(sites, by = "location2")
}

# =============================================================================
# 3. Step 1 out of sample: leave-one-year-out -----------------------------------
#    ./hpc/submit_peak.sh loyo refits the model 13 times, each without one
#    year, and predicts that year's weekly counts and peaks (its year and
#    zone-year effects drawn from their fitted distributions).
#    Copy loyo/ here. Scored: weekly log(CEW + 1) (RMSE, 90% coverage, ELPD)
#    and the peak week - against the full fit's estimate and against the
#    observed peak week of well-sampled site-years (the maximum not at the
#    first/last sampled week).
# =============================================================================
# ---- 3a. Candidate models: site vs zone_year -----------------------------------
#    PSIS-LOO (one week left out): do zone-year shifts improve the fit to the
#    observed seasons? LOYO (one year left out): do they help predict a new
#    year, or overfit? A difference within ~2 se_diff is not a preference;
#    prefer the simpler model then.
loo_have <- Filter(file.exists, sapply(names(peak_terms), function(m) peak_files(m)$loo, simplify = FALSE))
if (length(loo_have) >= 2) print(loo::loo_compare(lapply(loo_have, readRDS)))

loyo_lpd <- lapply(names(peak_terms), function(m) {
  f <- list.files(peak_files(m)$loyo_dir, "^loyo_peak_\\d{4}\\.rds$", full.names = TRUE)
  if (length(f)) bind_rows(lapply(f, function(x) readRDS(x)$pred)) %>%
    transmute(location2, year, woy, model = m, lpd)
})
loyo_lpd <- bind_rows(loyo_lpd)
if (n_distinct(loyo_lpd$model) >= 2) {
  w <- loyo_lpd %>% tidyr::pivot_wider(names_from = model, values_from = lpd) %>% na.omit()
  print(loyo_compare(as.matrix(w[, unique(loyo_lpd$model)]), w$year))   # same rows for all models
}

# ---- 3b. Out-of-sample skill of the model in use (peak_use) --------------------
loyo_files <- list.files(pf$loyo_dir, "^loyo_peak_\\d{4}\\.rds$", full.names = TRUE)
if (length(loyo_files)) {
  lo <- lapply(loyo_files, readRDS)
  names(lo) <- sub("^loyo_peak_(\\d{4})\\.rds$", "\\1", basename(loyo_files))
  loyo_pred <- bind_rows(lapply(lo, `[[`, "pred"))
  loyo_scores <- loyo_pred %>%
    group_by(zone = zl[as.character(zone_30y)]) %>%
    group_modify(~ data.frame(n = nrow(.x), rmse = sqrt(mean((.x$log_CEW - .x$pred_mean)^2)),
                              bias = mean(.x$pred_mean - .x$log_CEW),
                              cover90 = mean(.x$log_CEW >= .x$pred_lo & .x$log_CEW <= .x$pred_hi),
                              elpd_per_obs = mean(.x$lpd)))
  print(loyo_scores)
  print(bind_rows(lapply(lo, `[[`, "scores"), .id = "year") %>% filter(group == "all"))   # by held-out year

  # Observed peak week of each site-year (week of the largest count)
  obs_peak <- d_obs %>% group_by(loc_year = as.character(loc_year)) %>%
    summarise(obs_peak = woy[which.max(log_CEW)], interior = obs_peak > min(woy) & obs_peak < max(woy))

  loyo_peaks <- bind_rows(lapply(lo, `[[`, "peaks")) %>%
    transmute(loc_year = as.character(loc_year), zone_30y, year,
              loyo_med = peak_med, loyo_lo = peak_lo, loyo_hi = peak_hi) %>%
    left_join(obs_peak, by = "loc_year")
  if (!is.null(pk)) loyo_peaks <- loyo_peaks %>% left_join(peak_sy %>% select(loc_year, fit_med = peak_med), by = "loc_year")

  # Peak week: out-of-sample vs full fit, and both vs the observed peak
  print(loyo_peaks %>% group_by(zone = zl[as.character(zone_30y)]) %>%
          summarise(n = n(),
                    mae_vs_fit = if (!is.null(pk)) mean(abs(loyo_med - fit_med), na.rm = TRUE) else NA,
                    mae_vs_obs = mean(abs(loyo_med - obs_peak)[interior]),
                    cover95_obs = mean((obs_peak >= loyo_lo & obs_peak <= loyo_hi)[interior]),
                    mae_fit_vs_obs = if (!is.null(pk)) mean(abs(fit_med - obs_peak)[interior], na.rm = TRUE) else NA))

  # Zone-year mean peak predicted for the held-out year vs the full fit
  loyo_zone_year <- bind_rows(lapply(lo, function(r) {
    s <- r$peaks
    Z <- vapply(split(seq_len(ncol(r$peak_draws)), paste(s$zone_30y, s$year)),
                function(j) rowMeans(r$peak_draws[, j, drop = FALSE]), numeric(nrow(r$peak_draws)))
    data.frame(key = colnames(Z), t(apply(Z, 2, q95)))
  })) %>% tidyr::separate(key, c("zone_30y", "year"), convert = TRUE)
  if (!is.null(pk))
    loyo_zone_year <- loyo_zone_year %>%
      left_join(peak_zone_year %>% select(zone_30y, year, fit_med = med), by = c("zone_30y", "year"))
}

# =============================================================================
# 4. Filled weekly populations and regional clustering ---------------------------
#    Step-1 fills (fill_<model>.rds, from cew_peak_hpc.R fit): E[log(CEW + 1)]
#    for every week of every trapped site-year, sampled or not.
#    fill_mode
#      "gaps"  : observed log(CEW + 1) where the trap was sampled, the step-1
#                expectation in unsampled weeks (default)
#      "model" : the step-1 expectation everywhere. The model has no site-year
#                terms, so its week-to-week changes come only from EVI/GDD and
#                the shared year effect: Granger would then test the model's
#                covariates rather than the moths. For comparison only.
#    Held-out years (LOYO, section 6): their unsampled weeks are filled by the
#    step-1 fit without that year (loyo/<model>/), so a held-out year never
#    informs its own predictors. Without fills, observed counts only are used.
#    Clustering: single traps are sampled a few weeks a year, so sites are
#    grouped into regions within their zone. Ward clustering on km coordinates
#    and each site's mean filled curve over its years (season weeks), mixed
#    by space_weight; k per zone = the largest with mean silhouette >= 0.5
#    and every region >= min_sites sites (cluster_regions()). Check the
#    silhouette table; fix k with k_by_zone (e.g. c(`2` = 4)) if wanted.
# =============================================================================
fill_mode    <- "gaps"
space_weight <- 0.5
wk_season <- min(peak_win$wk_lo):max(peak_win$wk_hi)            # zone windows
wk_range  <- (min(peak_win$wk_lo) - 8):max(peak_win$wk_hi)      # + 8 lag weeks
ctrl <- d_all %>% select(location2, year, woy, all_of(anom_vars))
# The step-1 model is Gaussian on log(CEW + 1), so its expectation can dip
# below 0 where the curve extrapolates (season edges); log(CEW + 1) cannot,
# so fills are floored at 0 (the raw values stay in fill_<model>.rds)
fill_at <- function(observed, log_CEW, fill) {
  fill <- pmax(fill, 0)
  if (fill_mode == "gaps") ifelse(observed, log_CEW, fill) else fill
}
with_fill <- function(f) f %>%
  mutate(value = fill_at(observed, log_CEW, fill_mean),
         obs_value = ifelse(observed, log_CEW, NA_real_), location2 = as_int(location2)) %>%
  left_join(ctrl %>% mutate(location2 = as_int(location2)), by = c("location2", "year", "woy"))

fl <- if (file.exists(pf$fill)) readRDS(pf$fill) else NULL
if (!is.null(fl)) {
  d_g <- with_fill(fl$fill)
} else {
  message("No ", pf$fill, " (cew_peak_hpc.R fit): step 2 uses observed counts only")
  d_g <- d_obs %>% mutate(location2 = as_int(location2), value = log_CEW, obs_value = log_CEW)
}
loyo_fill_files <- file.path(pf$loyo_dir, sprintf("loyo_peak_%d.rds", yrs))
d_g_oos <- NULL
if (!is.null(fl)) {
  if (all(file.exists(loyo_fill_files))) {
    d_g_oos <- with_fill(bind_rows(lapply(loyo_fill_files, function(f) readRDS(f)$fill)))
  } else message("LOYO fills incomplete: held-out years use the full-fit fills (their ",
                 "unsampled weeks then carry information from the held-out year)")
}

# Each site's mean filled curve over its years (season weeks) for clustering
curves <- NULL
if (!is.null(fl)) {
  cw <- d_g %>% filter(woy %in% wk_season) %>% group_by(location2, woy) %>%
    summarise(v = mean(value), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = woy, values_from = v)
  curves <- as.matrix(cw[, -1]); rownames(curves) <- cw$location2
  miss <- which(is.na(curves), arr.ind = TRUE)                   # rare grid gaps
  curves[miss] <- colMeans(curves, na.rm = TRUE)[miss[, "col"]]
}
sites_g <- sites %>% filter(location2 %in% as_int(d_obs$location2))
cl <- cluster_regions(sites_g, curves = curves, space_weight = space_weight,
                      k_range = 2:6, min_sites = 10, sil_min = 0.5)
cl$silhouette
regions <- cl$sites
reg_info <- region_info(regions)

# =============================================================================
# 5. Regional weekly series -----------------------------------------------------
#    Y      = mean of `value` (observed + fills) over the region's trapped
#             site-years: the lagged predictors, complete within the season
#    Y_resp = mean over the sites actually trapped that week: the response,
#             so every Granger regression is scored on observed moths
#    Controls = the region's mean local weather/EVI anomalies (anom_vars).
#    Seasonality is removed inside each Granger fit (section 6), not here,
#    so a held-out year never informs its own climatology.
# =============================================================================
add_region <- function(d) d %>% inner_join(regions %>% select(location2, region), by = "location2")
series <- function(d, zone_level = FALSE) {
  if (zone_level) d <- d %>% mutate(region = zone_code[as.character(zone_30y)])
  region_series(d, value = "value", resp = "obs_value", controls = anom_vars,
                years = yrs, weeks = wk_range)
}
rs_reg  <- series(add_region(d_g))
rs_zone <- series(d_g, zone_level = TRUE)
rs_reg_oos  <- if (!is.null(d_g_oos)) series(add_region(d_g_oos))
rs_zone_oos <- if (!is.null(d_g_oos)) series(d_g_oos, zone_level = TRUE)

# Response weeks: in the sink zone's season window, >= min_obs trapped sites
min_obs <- 1
zone_of <- function(r) as.integer(names(zone_code)[zone_code == substr(r, 1, 1)])
resp_ok <- function(rs, r) {
  w  <- peak_win %>% filter(zone_30y == zone_of(r))
  ok <- !is.na(rs$Y_resp[r, , ]) & rs$n_obs[r, , ] >= min_obs
  ok[, setdiff(1:52, w$wk_lo:w$wk_hi)] <- FALSE
  ok
}

# Coverage: share of in-window weeks with a trapped site, sites per such week
cov_tab <- function(rs) {
  bind_rows(lapply(rs$regions, function(r) {
    w  <- peak_win %>% filter(zone_30y == zone_of(r))
    nw <- rs$n_obs[r, , w$wk_lo:w$wk_hi]; nw[is.na(nw)] <- 0
    data.frame(region = r, weeks_with_data = mean(nw > 0), sites_per_week = mean(nw[nw > 0]))
  }))
}
cov_tab(rs_zone)
# Regions trapped in fewer than min_cover of their in-window weeks give too
# few responses to test; they stay on the maps but get no Granger pairs
min_cover <- 0.5
reg_info <- reg_info %>% left_join(cov_tab(rs_reg), by = "region") %>%
  mutate(tested = weeks_with_data >= min_cover)
reg_info

# =============================================================================
# 6. Granger causality between zones ----------------------------------------------
#    For sink zone j and source zone i, with a_ = log(CEW + 1) minus its
#    seasonal climatology (Fourier, H = 3, fit on the training years):
#      R: a_j(w) = c + sum_l phi_l a_j(w - l) + controls_j(w) + e
#      U: R + sum_k gamma_k a_i(w - k)
#    Southern lags "Granger-cause" the sink if U predicts better than R.
#    Lags stay within a season (no carry-over across winter).
#    Settings
#      lag1-lag8   : one Southern lag k (own lags 1-2): when does the South
#                    best predict the sink? (the arrival delay)
#      order1-6    : the standard Granger test of order p (own and source lags 1..p)
#      lead1-3     : FUTURE Southern anomalies - falsification: should not help
#    Scores
#      p_F       : in-sample F test of the gammas (nominal; serially
#                  correlated errors make it optimistic)
#      loyo_gain : 1 - SSE_U / SSE_R over held-out years (> 0: out-of-sample gain)
#      p_perm    : year-permutation p-value of loyo_gain; p_perm_max for the
#                  best lag across lag1-lag8 (accounts for picking that lag)
#    Reverse pairs (Transitional/Northern -> Southern) are a falsification:
#    a northward flight should not show up as strongly in reverse.
# =============================================================================
set.seed(42)
n_perm <- 199
own_p  <- 1:2
set_lag  <- setNames(lapply(1:8, function(k) list(own = own_p, src = k)), paste0("lag", 1:8))
set_ord  <- setNames(lapply(1:6, function(p) list(own = 1:p, src = 1:p)), paste0("order", 1:6))
set_lead <- setNames(lapply(1:3, function(k) list(own = own_p, src = -k)), paste0("lead", 1:3))
set_all  <- c(set_lag, set_ord, set_lead)

# One sink/source pair: lags from the filled series (held-out years from the
# LOYO fills when rs_oos is given), responses from trapped sites only
run_pair <- function(rs, source, sink, settings, rs_oos = NULL, perm_settings = set_lag,
                     n_perm = 0, keep_pred = FALSE) {
  a <- list(Ys = rs$Y[sink, , ], Yk = rs$Y[source, , ], X = lapply(rs$X, function(x) x[sink, , ]),
            Ys_resp = rs$Y_resp[sink, , ], resp_ok = resp_ok(rs, sink),
            Ys_oos = if (!is.null(rs_oos)) rs_oos$Y[sink, , ], Yk_oos = if (!is.null(rs_oos)) rs_oos$Y[source, , ])
  g  <- granger_pair(a$Ys, a$Yk, a$X, settings, Ys_resp = a$Ys_resp, resp_ok = a$resp_ok,
                     Ys_oos = a$Ys_oos, Yk_oos = a$Yk_oos, keep_pred = keep_pred)
  tab <- g$table %>% mutate(source = source, sink = sink, .before = 1)
  if (n_perm > 0) {
    null <- granger_perm(a$Ys, a$Yk, a$X, perm_settings, n_perm = n_perm, Yk_oos = a$Yk_oos,
                         Ys_resp = a$Ys_resp, resp_ok = a$resp_ok, Ys_oos = a$Ys_oos)
    pp <- perm_p(tab$loyo_gain[match(names(perm_settings), tab$setting)], null,
                 max_cols = names(perm_settings) %in% names(set_lag))
    tab$p_perm <- unname(pp$per_setting[match(tab$setting, names(perm_settings))])
    tab$p_perm_max <- pp$max
  }
  pred <- if (keep_pred) g$pred %>% mutate(source = source, sink = sink, year = rs$years[year_i]) else NULL
  list(table = tab, pred = pred)
}

zone_pairs <- data.frame(source = c("S", "S", "T", "N"), sink = c("T", "N", "S", "S"),
                         direction = c("forward", "forward", "reverse", "reverse"))
gz <- lapply(seq_len(nrow(zone_pairs)), function(i)
  run_pair(rs_zone, zone_pairs$source[i], zone_pairs$sink[i], set_all, rs_oos = rs_zone_oos,
           perm_settings = c(set_lag, set_ord), n_perm = n_perm, keep_pred = TRUE))
granger_zone <- bind_rows(lapply(gz, `[[`, "table")) %>%
  left_join(zone_pairs, by = c("source", "sink"))
granger_pred_zone <- bind_rows(lapply(gz, `[[`, "pred"))

# Lag profile and order tests, forward pairs
granger_zone %>% filter(direction == "forward") %>%
  select(source, sink, setting, n, p_F, loyo_gain, years_better, years_tested, p_perm) %>%
  print(n = Inf)
# Falsification: leads and reverse direction should show no gain
granger_zone %>% filter(grepl("^lead", setting) | direction == "reverse") %>%
  select(source, sink, setting, p_F, loyo_gain, years_better, p_perm) %>% print(n = Inf)

# Best Southern lag per sink zone
best_zone <- granger_zone %>% filter(direction == "forward", setting %in% names(set_lag)) %>%
  group_by(sink) %>% slice_max(loyo_gain, n = 1) %>% ungroup()
best_zone

# ---- 6b. Step-1 uncertainty in the fills ---------------------------------------
#    The tests above use the posterior mean fills. Rerun the forward zone lag
#    profiles with n_fill posterior draws of the fills (fill_<model>.rds keeps
#    them): the spread of loyo_gain is the step-1 uncertainty carried into
#    step 2. Observed weeks are the same in every draw.
n_fill <- 50
granger_fill_unc <- NULL
if (!is.null(fl) && ncol(fl$draws) > 0) {
  m_use <- seq_len(min(n_fill, ncol(fl$draws)))
  unc <- lapply(m_use, function(m) {
    d_m <- d_g %>% mutate(value = fill_at(observed, log_CEW, fl$draws[, m]))
    rs_m <- series(d_m, zone_level = TRUE)
    bind_rows(lapply(c("T", "N"), function(snk)
      run_pair(rs_m, "S", snk, set_lag, rs_oos = rs_zone_oos)$table)) %>% mutate(draw = m)
  })
  granger_fill_unc <- bind_rows(unc) %>% group_by(source, sink, setting) %>%
    summarise(gain_med = median(loyo_gain), gain_lo = quantile(loyo_gain, 0.05),
              gain_hi = quantile(loyo_gain, 0.95), P_gain_pos = mean(loyo_gain > 0), .groups = "drop") %>%
    left_join(granger_zone %>% select(source, sink, setting, gain_at_mean = loyo_gain),
              by = c("source", "sink", "setting"))
  print(granger_fill_unc, n = Inf)
}

# =============================================================================
# 7. Granger causality between regions --------------------------------------------
#    Every tested Southern region -> every tested Transitional/Northern
#    region (coverage >= min_cover, section 5), lags 1-8.
#    A link = the best lag's LOYO gain > 0 with p_perm_max (that best-lag
#    statistic against its permutation null) < 0.05 after Benjamini-Hochberg
#    across pairs. Geometry: centroid distance (km) and bearing; distance /
#    best lag = apparent northward speed. With step-1 peaks: the posterior
#    delay of the sink's mean peak behind the source's, over common years.
# =============================================================================
tested <- reg_info$region[reg_info$tested]
reg_pairs <- expand.grid(source = tested[startsWith(tested, "S")],
                         sink = tested[!startsWith(tested, "S")],
                         stringsAsFactors = FALSE)
t0 <- Sys.time()
gr <- lapply(seq_len(nrow(reg_pairs)), function(i)
  run_pair(rs_reg, reg_pairs$source[i], reg_pairs$sink[i], set_lag, rs_oos = rs_reg_oos, n_perm = n_perm)$table)
message(sprintf("Region pairs: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
granger_reg <- bind_rows(gr)

ri <- reg_info %>% select(region, longitude, latitude)
granger_pairs <- granger_reg %>%
  group_by(source, sink) %>% slice_max(loyo_gain, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(source, sink, best_lag = as.integer(sub("lag", "", setting)), n, best_gain = loyo_gain,
            years_better, years_tested, p_F, p_perm_max) %>%
  mutate(p_adj = p.adjust(p_perm_max, "BH"),
         link = best_gain > 0 & p_adj < 0.05) %>%
  left_join(ri %>% rename_with(~ paste0(.x, "_src"), -region), by = c(source = "region")) %>%
  left_join(ri %>% rename_with(~ paste0(.x, "_snk"), -region), by = c(sink = "region")) %>%
  mutate(dist_km = geosphere::distHaversine(cbind(longitude_src, latitude_src),
                                            cbind(longitude_snk, latitude_snk)) / 1000,
         bearing = geosphere::bearing(cbind(longitude_src, latitude_src),
                                      cbind(longitude_snk, latitude_snk)),
         speed_km_wk = dist_km / best_lag,
         sink_zone = ifelse(startsWith(sink, "T"), "Transitional", "Northern"))

# Step-1 peak delay of each pair (sink - source, weeks)
if (!is.null(pk)) {
  RP <- region_year_peaks(P, sy %>% left_join(regions %>% select(location2, region), by = "location2"))
  granger_pairs <- granger_pairs %>%
    bind_cols(bind_rows(lapply(seq_len(nrow(.)), function(i)
      peak_delay(RP, granger_pairs$source[i], granger_pairs$sink[i]))))
}
granger_pairs %>% arrange(sink, desc(best_gain)) %>% print(n = Inf)

# Does the predictive lag grow with distance (a travelling front)? Over the
# links, or over all pairs with a positive gain when there are too few links
lag_dist <- granger_pairs %>% filter(if (sum(link) >= 5) link else best_gain > 0)
if (nrow(lag_dist) >= 3) summary(lm(best_lag ~ dist_km, data = lag_dist))

# =============================================================================
# 8. Save outputs for CEW_JIPM_Figures.R --------------------------------------------
# =============================================================================
write.csv(regions %>% select(location2, location, zone_30y, longitude, latitude, x_km, y_km, region),
          "cew_regions.csv", row.names = FALSE)
write.csv(granger_zone, "granger_zone.csv", row.names = FALSE)
write.csv(granger_reg, "granger_region_lags.csv", row.names = FALSE)
write.csv(granger_pairs, "granger_pairs.csv", row.names = FALSE)
# Each region's mean seasonal curve (observed + fills, averaged over years)
region_curves <- add_region(d_g) %>% filter(woy %in% wk_season) %>%
  group_by(region, zone_30y, woy) %>% summarise(value = mean(value), .groups = "drop")
out_step2 <- list(regions = regions, reg_info = reg_info, silhouette = cl$silhouette,
                  region_curves = region_curves,
                  granger_zone = granger_zone, granger_pred_zone = granger_pred_zone,
                  granger_reg = granger_reg, granger_pairs = granger_pairs,
                  granger_fill_unc = granger_fill_unc,
                  settings = list(n_perm = n_perm, own_lags = own_p, years = yrs,
                                  filled = !is.null(fl), fill_mode = fill_mode,
                                  loyo_fills = !is.null(d_g_oos), space_weight = space_weight,
                                  peak_model = peak_use))
if (!is.null(granger_fill_unc)) write.csv(granger_fill_unc, "granger_fill_uncertainty.csv", row.names = FALSE)
if (!is.null(pk)) {
  write.csv(peak_sy, "peak_site_year.csv", row.names = FALSE)
  write.csv(peak_zone_year, "peak_zone_year.csv", row.names = FALSE)
  write.csv(zone_delay, "peak_zone_delay.csv", row.names = FALSE)
  write.csv(peak_site, "peak_site.csv", row.names = FALSE)
}
if (length(loyo_files)) {
  write.csv(loyo_scores, "loyo_peak_scores.csv", row.names = FALSE)
  write.csv(loyo_peaks, "loyo_peak_site_year.csv", row.names = FALSE)
  write.csv(loyo_zone_year, "loyo_peak_zone_year.csv", row.names = FALSE)
}
saveRDS(out_step2, "granger_results.rds")
