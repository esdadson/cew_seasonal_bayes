# =============================================================================
# step2/step2_granger.R
# Step 2: do Southern populations (zone_30y == 2) lead Transitional (1) and
# Northern (0) ones?
#   1. The step-1 weekly populations (observed weeks + the selected model's
#      predictions for unsampled weeks) are clustered spatially into regions
#      within each zone.
#   2. Granger causality: Transitional/Northern regional series modelled from
#      their own history plus lagged Southern series, at lags of 1-8 weeks -
#      in-sample F test, out-of-sample gain over leave-one-year-out folds,
#      and a year-permutation null. The most predictive lag is the apparent
#      arrival delay; it is compared with the step-1 peak delay.
# Run from the project directory, after step 1 (step1/step1_hpc.R compare)
# has written results/step1/. Writes results/step2/ for CEW_JIPM_Figures.R.
#
# Sections
#   0. Step-1 outputs
#   1. Weekly populations and regional clustering
#   2. Regional weekly series
#   3. Granger causality between zones (3b: step-1 uncertainty)
#   4. Granger causality between regions, lag times and peak delays
#   5. Save outputs
# =============================================================================

library(dplyr); library(tidyr)
source("step2/step2_functions.R")
step1_dir <- "results/step1"
step2_dir <- "results/step2"
dir.create(step2_dir, recursive = TRUE, showWarnings = FALSE)
zl <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")

# =============================================================================
# 0. Step-1 outputs ----------------------------------------------------------------
#    Copy results/step1/ back from Hazel (hpc/README.md step 11):
#    model_data.rds (modelling data), predictions.rds (E[log(CEW + 1)] for
#    weeks 1-52 of every monitored site-year: mean + draws) and
#    peak_draws.rds (peak week of every site-year in every draw).
# =============================================================================
need <- file.path(step1_dir, c("model_data.rds", "predictions.rds", "peak_draws.rds"))
if (!all(file.exists(need)))
  stop("Missing step-1 outputs: ", paste(need[!file.exists(need)], collapse = ", "),
       ". Run ./hpc/submit_step1.sh compare and copy results/step1/ here.")
md <- readRDS(need[1]); pr <- readRDS(need[2]); pk <- readRDS(need[3])
d_all <- md$d_all; sites <- md$sites; anom_vars <- md$anom_vars; yrs <- md$yrs; peak_win <- md$peak_win
d_obs <- d_all %>% filter(observed)
message("Step-1 model: ", pr$model)

# =============================================================================
# 1. Weekly populations and regional clustering ------------------------------------
#    fill_mode
#      "gaps"  : observed log(CEW + 1) where the trap was sampled, the step-1
#                prediction in unsampled weeks (default)
#      "model" : the step-1 prediction everywhere. The model has no site-year
#                terms, so its week-to-week changes within a site-year come
#                only from EVI/GDD: Granger would then test the model's
#                covariates rather than the moths. For comparison only.
#    The step-1 fit used every year, so a held-out year's unsampled weeks are
#    filled by a model that saw that year's counts (through its year effect);
#    the Granger responses are observed counts only (section 2).
#    Clustering: single traps are sampled a few weeks a year, so sites are
#    grouped into regions within their zone. Ward clustering on km coordinates
#    and each site's mean weekly curve over its years (season weeks), mixed by
#    space_weight; k per zone = the largest with mean silhouette >= 0.5 and
#    every region >= min_sites sites (cluster_regions()). Check the
#    silhouette table; fix k with k_by_zone (e.g. c(`2` = 4)) if wanted.
# =============================================================================
fill_mode    <- "gaps"
space_weight <- 0.5
wk_season <- min(peak_win$wk_lo):max(peak_win$wk_hi)            # zone windows
wk_range  <- (min(peak_win$wk_lo) - 8):max(peak_win$wk_hi)      # + 8 lag weeks
ctrl <- d_all %>% select(location2, year, woy, all_of(anom_vars))
# The step-1 model is Gaussian on log(CEW + 1), so its expectation can dip
# below 0 where the curve extrapolates (season edges); log(CEW + 1) cannot,
# so predictions are floored at 0 (the raw values stay in predictions.rds)
fill_at <- function(observed, log_CEW, fill) {
  fill <- pmax(fill, 0)
  if (fill_mode == "gaps") ifelse(observed, log_CEW, fill) else fill
}
d_g <- pr$pred %>%
  mutate(value = fill_at(observed, log_CEW, pred_mean),
         obs_value = ifelse(observed, log_CEW, NA_real_)) %>%
  left_join(ctrl, by = c("location2", "year", "woy"))

# Each site's mean weekly curve over its years (season weeks)
cw <- d_g %>% filter(woy %in% wk_season) %>% group_by(location2, woy) %>%
  summarise(v = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = woy, values_from = v)
curves <- as.matrix(cw[, -1]); rownames(curves) <- cw$location2
miss <- which(is.na(curves), arr.ind = TRUE)                     # rare grid gaps
curves[miss] <- colMeans(curves, na.rm = TRUE)[miss[, "col"]]

sites_g <- sites %>% filter(location2 %in% d_obs$location2)
cl <- cluster_regions(sites_g, curves = curves, space_weight = space_weight,
                      k_range = 2:6, min_sites = 10, sil_min = 0.5)
cl$silhouette
regions <- cl$sites
reg_info <- region_info(regions)

# =============================================================================
# 2. Regional weekly series ----------------------------------------------------------
#    Y      = mean of `value` (observed + predictions) over the region's
#             monitored site-years: the lagged predictors, complete in season
#    Y_resp = mean over the sites actually trapped that week: the response,
#             so every Granger regression is scored on observed moths
#    Controls = the region's mean local weather/EVI anomalies (anom_vars).
#    Seasonality is removed inside each Granger fit (section 3), not here,
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
# 3. Granger causality between zones -------------------------------------------------
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

# One sink/source pair: lags from the filled series, responses from trapped sites only
run_pair <- function(rs, source, sink, settings, perm_settings = set_lag, n_perm = 0, keep_pred = FALSE) {
  a <- list(Ys = rs$Y[sink, , ], Yk = rs$Y[source, , ], X = lapply(rs$X, function(x) x[sink, , ]),
            Ys_resp = rs$Y_resp[sink, , ], resp_ok = resp_ok(rs, sink))
  g  <- granger_pair(a$Ys, a$Yk, a$X, settings, Ys_resp = a$Ys_resp, resp_ok = a$resp_ok,
                     keep_pred = keep_pred)
  tab <- g$table %>% mutate(source = source, sink = sink, .before = 1)
  if (n_perm > 0) {
    null <- granger_perm(a$Ys, a$Yk, a$X, perm_settings, n_perm = n_perm,
                         Ys_resp = a$Ys_resp, resp_ok = a$resp_ok)
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
  run_pair(rs_zone, zone_pairs$source[i], zone_pairs$sink[i], set_all,
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

# ---- 3b. Step-1 uncertainty ----------------------------------------------------
#    The tests above use the posterior mean predictions. Rerun the forward
#    zone lag profiles with n_fill posterior draws of them (predictions.rds
#    keeps STEP1_PRED_DRAWS): the spread of loyo_gain is the step-1
#    uncertainty carried into step 2. Observed weeks are the same in every draw.
n_fill <- 50
granger_fill_unc <- NULL
if (ncol(pr$draws) > 0) {
  m_use <- seq_len(min(n_fill, ncol(pr$draws)))
  unc <- lapply(m_use, function(m) {
    d_m <- d_g %>% mutate(value = fill_at(observed, log_CEW, pr$draws[, m]))
    rs_m <- series(d_m, zone_level = TRUE)
    bind_rows(lapply(c("T", "N"), function(snk) run_pair(rs_m, "S", snk, set_lag)$table)) %>%
      mutate(draw = m)
  })
  granger_fill_unc <- bind_rows(unc) %>% group_by(source, sink, setting) %>%
    summarise(gain_med = median(loyo_gain), gain_lo = quantile(loyo_gain, 0.05),
              gain_hi = quantile(loyo_gain, 0.95), P_gain_pos = mean(loyo_gain > 0), .groups = "drop") %>%
    left_join(granger_zone %>% select(source, sink, setting, gain_at_mean = loyo_gain),
              by = c("source", "sink", "setting"))
  print(granger_fill_unc, n = Inf)
}

# =============================================================================
# 4. Granger causality between regions, lag times and peak delays -------------------
#    Every tested Southern region -> every tested Transitional/Northern
#    region (coverage >= min_cover, section 2), lags 1-8.
#    A link = the best lag's LOYO gain > 0 with p_perm_max (that best-lag
#    statistic against its permutation null) < 0.05 after Benjamini-Hochberg
#    across pairs. Geometry: centroid distance (km) and bearing; distance /
#    best lag = apparent northward speed. Step-1 peak delay: the posterior
#    delay of the sink's mean peak week behind the source's, over common years.
# =============================================================================
tested <- reg_info$region[reg_info$tested]
reg_pairs <- expand.grid(source = tested[startsWith(tested, "S")],
                         sink = tested[!startsWith(tested, "S")],
                         stringsAsFactors = FALSE)
if (!nrow(reg_pairs))
  stop("No Southern -> Transitional/Northern pair of tested regions (coverage >= min_cover = ",
       min_cover, "): see reg_info; lower min_cover or change the clustering (section 1)")
t0 <- Sys.time()
gr <- lapply(seq_len(nrow(reg_pairs)), function(i)
  run_pair(rs_reg, reg_pairs$source[i], reg_pairs$sink[i], set_lag, n_perm = n_perm)$table)
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
sy <- pk$site_years %>% left_join(regions %>% select(location2, region), by = "location2")
RP <- region_year_peaks(pk$draws[, !is.na(sy$region), drop = FALSE], sy %>% filter(!is.na(region)))
granger_pairs <- granger_pairs %>%
  bind_cols(bind_rows(lapply(seq_len(nrow(.)), function(i)
    peak_delay(RP, granger_pairs$source[i], granger_pairs$sink[i]))))
granger_pairs %>% arrange(sink, desc(best_gain)) %>% print(n = Inf)

# Does the predictive lag grow with distance (a travelling front)? Over the
# links, or over all pairs with a positive gain when there are too few links
lag_dist <- granger_pairs %>% filter(if (sum(link) >= 5) link else best_gain > 0)
if (nrow(lag_dist) >= 3) print(summary(lm(best_lag ~ dist_km, data = lag_dist)))

# =============================================================================
# 5. Save outputs (results/step2/, read by CEW_JIPM_Figures.R) -----------------------
# =============================================================================
write.csv(regions %>% select(location2, location, zone_30y, longitude, latitude, x_km, y_km, region),
          file.path(step2_dir, "cew_regions.csv"), row.names = FALSE)
write.csv(granger_zone, file.path(step2_dir, "granger_zone.csv"), row.names = FALSE)
write.csv(granger_reg, file.path(step2_dir, "granger_region_lags.csv"), row.names = FALSE)
write.csv(granger_pairs, file.path(step2_dir, "granger_pairs.csv"), row.names = FALSE)
if (!is.null(granger_fill_unc))
  write.csv(granger_fill_unc, file.path(step2_dir, "granger_fill_uncertainty.csv"), row.names = FALSE)
# Each region's mean seasonal curve (observed + predictions, averaged over years)
region_curves <- add_region(d_g) %>% filter(woy %in% wk_season) %>%
  group_by(region, zone_30y, woy) %>% summarise(value = mean(value), .groups = "drop")
saveRDS(list(regions = regions, reg_info = reg_info, silhouette = cl$silhouette,
             region_curves = region_curves,
             granger_zone = granger_zone, granger_pred_zone = granger_pred_zone,
             granger_reg = granger_reg, granger_pairs = granger_pairs,
             granger_fill_unc = granger_fill_unc,
             settings = list(n_perm = n_perm, own_lags = own_p, years = yrs, fill_mode = fill_mode,
                             space_weight = space_weight, step1_model = pr$model)),
        file.path(step2_dir, "granger_results.rds"))
