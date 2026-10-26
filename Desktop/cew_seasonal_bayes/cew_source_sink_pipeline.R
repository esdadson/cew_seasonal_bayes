# =============================================================================
# cew_source_sink_pipeline.R
# H. zea: seasonal curves by overwintering zone (stage 1) and
#          Southern -> Transitional/Northern source-sink dynamics (stage 2)
#
# Sections
#   0. Setup and ONE site table (location2 defined once)
#   1. Trap counts on ISO weeks
#   2. Daymet climate covariates (daily -> ISO week, climatology, anomalies)
#   3. CDL landscape covariates and EVI
#   4. Assemble modelling data (+ unsampled Southern weeks to fill)
#   5. Stage 1: the seasonal model of cew_fourier_covariates.R (fit_gp2.1)
#   6. Climatological expectation, anomalies, filled Southern weeks
#   7. Stage 2: source pressure, kernel/lag screening, Bayesian sink model
#   8. Falsification tests
#   9. Outputs for "when" and "where"
#  10. Out-of-sample evaluation: 2020
#
# Stage 1 is the model in cew_fourier_covariates.R (log(CEW + 1) ~ EVI + gdd +
# crops + zone Fourier curves + exponential GP + year + trap type). Here it is
# refit on the training years only (yrs_train), and it is used to
#   - define each site-week's climatological expectation (anomaly = observed -
#     expectation) and
#   - fill unsampled weeks at Southern traps, so the source pressure on
#     Northern/Transitional sinks uses a complete weekly Southern field.
# 2020 (yr_test) is held out of every fit and every screening/selection step,
# then predicted in section 10.
#
# Sections 2-3 read the raster extractions made on the HPC by
# cew_covariates_hpc.R (see hpc/submit_covariates.sh); copy its output
# directory here as covariates/.
# =============================================================================

library(dplyr); library(tidyr); library(lubridate); library(zoo)
library(brms); library(posterior); library(ggplot2)

# Sections 0-4 (data) and the stage-1 model spec (f_my, pri_my, fit_stage1)
source("cew_stage1_data.R")

# =============================================================================
# 5. Stage 1: seasonal model of cew_fourier_covariates.R (fit_gp2.1) ------------
#    Model and priors: cew_stage1_data.R. Fit on the HPC (hpc/submit_stage1.sh,
#    see hpc/README.md) and copy fit1_gp21_train.rds here; this only reads it.
# =============================================================================
fit1 <- fit_stage1(d_train)
hz <- c("s1", "c1", "s2", "c2")

# ---- Peak timing -------------------------------------------------------------
# The zone curves are the zone_fac Fourier effects (there is no population-level
# Fourier term). EVI and GDD are raw covariates, so their seasonal shape is
# NOT in these curves: the peaks are of the curve net of EVI/GDD seasonality.
zone_arr <- ranef(fit1, summary = FALSE)$zone_fac[, , hz]  # draws x zones x coefs

# Search only the sampled season (5th-95th percentile of sampled weeks);
# outside it the curve is extrapolation.
win <- range(quantile(d_train$woy, c(0.05, 0.95)))
pk_zone <- peak_from_draws(zone_arr, window = win[1]:win[2])

lag_N <- pk_zone$peak[, "0"] - pk_zone$peak[, "2"]
lag_T <- pk_zone$peak[, "1"] - pk_zone$peak[, "2"]
rbind(N_minus_S = quantile(lag_N, c(.025, .5, .975)), T_minus_S = quantile(lag_T, c(.025, .5, .975)))
c(P_N_later = mean(lag_N > 0), P_T_later = mean(lag_T > 0))
pk_zone$at_edge                                            # check before reporting

# =============================================================================
# 6. Climatological expectation, anomalies, filled Southern weeks ---------------
#    Expectation: EVI and GDD at their site x week climatology, no year effect
#    (the GP is always included). anomaly = observed - expectation.
#    Filled (unsampled) Southern weeks get the anomaly the model predicts from
#    that week's own EVI/GDD. The year effect is left out of the fill: it is
#    shared by all zones and partly estimated from sink traps, so including it
#    would pass sink information into the source field.
#    2020 rows are predicted from the training fit (no new-level terms).
# =============================================================================
re_clim <- ~ (0 + s1 + c1 + s2 + c2 || zone_fac) + (1 | trap_type)   # all but year
set.seed(42)
ids <- sample(ndraws(fit1), 200)
lin <- function(nd) colMeans(posterior_linpred(fit1, newdata = nd, re_formula = re_clim,
                                               draw_ids = ids, allow_new_levels = TRUE))
mu_clim <- lin(d_all %>% mutate(evi_wk = evi_clim, gdd_wk = gdd_wk_clim) %>% scale_my())
mu_wx   <- lin(d_all)
d_all <- d_all %>%
  mutate(mu_clim = mu_clim,
         anom = ifelse(observed, log_CEW - mu_clim, mu_wx - mu_clim))
d_fit <- d_all %>% filter(observed)

# Full prediction grid: weeks 1-52 of every trapped site-year, all effects
# (2020's year effect is a new level, drawn from its fitted distribution).
# pred_* = posterior predictive mean and 90% interval of log(CEW + 1).
pp_grid <- posterior_predict(fit1, newdata = d_all, draw_ids = ids,
                             allow_new_levels = TRUE, sample_new_levels = "gaussian")
pred_grid <- d_all %>%
  transmute(location2, location, zone_30y, longitude, latitude, year, woy, observed, in_season,
            log_CEW, mu_clim, anom,
            pred_mean = colMeans(pp_grid),
            pred_lo   = apply(pp_grid, 2, quantile, 0.05),
            pred_hi   = apply(pp_grid, 2, quantile, 0.95),
            pred_CEW  = exp(pred_mean) - 1) %>%
  arrange(location2, year, woy)
write.csv(pred_grid, "stage1_prediction_grid.csv", row.names = FALSE)
rm(pp_grid)

# Quick first look (Lawton-style): yearly zone mean anomalies, training years
d_fit %>% filter(train) %>% group_by(zone_30y, year) %>% summarise(a = mean(anom), .groups = "drop") %>%
  pivot_wider(names_from = zone_30y, values_from = a) %>%
  summarise(r_S_T = cor(`2`, `1`), r_S_N = cor(`2`, `0`))

# =============================================================================
# 7. Stage 2: source pressure -------------------------------------------------
#    Every screen and fit here uses the training years only; 2020 columns of
#    the source field (observed 2020 Southern anomalies + fills) are built now
#    but used only in section 10.
# =============================================================================
dd   <- d_fit
tw   <- make_tw_index(dd$year)
src  <- sites %>% filter(location2 %in% dd$location2, zone_30y == 2)
snk  <- sites %>% filter(location2 %in% dd$location2, zone_30y != 2)
A_obs <- anomaly_matrix(dd %>% filter(zone_30y == 2), src$location2, tw)   # sources x year-week
A     <- A_obs
if (fill_sources) {
  # in-season gaps only; drop in_season to fill all 52 weeks of the grid
  A_fill <- anomaly_matrix(d_all %>% filter(!observed), src$location2, tw)
  A[is.na(A)] <- A_fill[is.na(A)]
}
c(share_sampled = mean(!is.na(A_obs)), share_after_fill = mean(!is.na(A)))
geom <- pair_geometry(snk, src)
obs    <- dd %>% filter(zone_30y != 2)
obs_tr <- obs %>% filter(train)

base_terms <- anom_vars                     # local weather controls in screening

# ---- 7a. Distance decay (where) ----------------------------------------------
lambdas <- c(100, 200, 400, 800, 1600, Inf)
P_dist  <- setNames(lapply(lambdas, function(l) source_pressure(kernel_weights(geom, l), A, tw$woy)),
                    paste0("lambda_", lambdas))
scr_N <- grid_screen(obs_tr %>% filter(zone_30y == 0), P_dist, tw, base_terms = base_terms)
scr_T <- grid_screen(obs_tr %>% filter(zone_30y == 1), P_dist, tw, base_terms = base_terms)
scr_N; scr_T
lam_best <- 400                              # <- set from the screens

# ---- 7b. Direction (where) ---------------------------------------------------
dir_grid <- expand.grid(theta0 = c(0, 30, 60, 90), kappa = c(0, 1, 3))
P_dir <- setNames(lapply(seq_len(nrow(dir_grid)), function(g)
  source_pressure(kernel_weights(geom, lam_best, dir_grid$kappa[g], dir_grid$theta0[g]), A, tw$woy)),
  paste0("theta", dir_grid$theta0, "_k", dir_grid$kappa))
grid_screen(obs_tr, P_dir, tw, base_terms = c("zone_fac", base_terms))

# ---- 7c. Travelling wave: lag grows with distance (when + speed) -------------
speeds <- c(100, 200, 400, 800)              # km per week
P_wave <- setNames(lapply(speeds, function(v)
  source_pressure(kernel_weights(geom, lam_best), A, tw$woy,
                  lag_mat = round(geom$D / v), max_lag = 8)), paste0("v_", speeds))
grid_screen(obs_tr, P_wave, tw, lags = 0:4, df = 3, base_terms = c("zone_fac", base_terms))

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
obs_b_tr <- obs_b %>% filter(train)

Qs <- colnames(Blag)
f2 <- bf(as.formula(paste(
  "anom ~ 0 + zone_fac + zone_fac:(", paste(Qs, collapse = " + "), ") +",
  paste(anom_vars, collapse = " + "), "+ own_lag4 + own_lag4_miss",
  "+ (1 | year_f) + (1 | location2)")))
# Same model without Southern pressure: the reference for section 10
f2_base <- bf(as.formula(paste(
  "anom ~ 0 + zone_fac +", paste(anom_vars, collapse = " + "), "+ own_lag4 + own_lag4_miss",
  "+ (1 | year_f) + (1 | location2)")))
# (1 | year_f) => lag effects identified from within-season and between-site
# variation only (conservative). Refit without it to include year-level
# Southern booms; report both.
pri2 <- c(prior(normal(0, 0.5), class = "b"), prior(student_t(3, 0, 1), class = "sd"))
fit2 <- brm(f2, data = obs_b_tr, family = gaussian(), chains = 4, cores = 4,
            prior = pri2, seed = 42, file = "fit2_sink_train", file_refit = "on_change")
fit2_base <- brm(f2_base, data = obs_b_tr, family = gaussian(), chains = 4, cores = 4,
                 prior = pri2, seed = 42, file = "fit2_base_train", file_refit = "on_change")

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
obs_b_tr$Pstar <- rowSums(as.matrix(obs_b_tr[, Qs]) * t(theta_hat[, as.character(obs_b_tr$zone_fac)]))
# (two-step: Pstar uses fit2's posterior means, so treat fit2_season / fit2_site as descriptive)
fit2_season <- brm(bf(anom ~ zone_fac + gdd_wk_anom + agdd_wk_anom + evi_anom +
                        s(woy, by = Pstar, k = 6) + (1 | year_f) + (1 | location2)),
                   data = obs_b_tr, family = gaussian(), chains = 4, cores = 4, seed = 42,
                   file = "fit2_season", file_refit = "on_change")

# Where on the sink side: site-specific sensitivity to source pressure
fit2_site <- brm(bf(anom ~ zone_fac + gdd_wk_anom + agdd_wk_anom + evi_anom + Pstar +
                      (1 + Pstar || location2) + (1 | year_f)),
                 data = obs_b_tr, family = gaussian(), chains = 4, cores = 4, seed = 42,
                 file = "fit2_site", file_refit = "on_change")

# =============================================================================
# 8. Falsification tests (training years) ---------------------------------------
# =============================================================================
# (i) Leads: future Southern anomalies should not predict sinks
loyo_gain(obs_b_tr, c("P_lead1", "P_lead2", "P_lead3"), c("zone_fac", anom_vars))

# (ii) Flat kernel (lambda = Inf) should do worse than the chosen kernel: see 7a

# (iii) Year-permutation null for the CV gain (training years permuted among themselves)
tw_tr <- tw %>% filter(year %in% yrs_train)
real_gain <- grid_screen(obs_tr, list(real = P_best), tw, base_terms = c("zone_fac", anom_vars))$gain
null_gain <- replicate(200, {
  Pp <- source_pressure(kernel_weights(geom, lam_best), permute_source_years(A, tw_tr), tw$woy)
  grid_screen(obs_tr, list(perm = Pp), tw, base_terms = c("zone_fac", anom_vars))$gain
})
mean(null_gain >= real_gain)                                 # permutation p-value

# (iv) Reverse direction: Transitional/Northern anomalies -> Southern sinks
A_rev    <- anomaly_matrix(obs, snk$location2, tw)
geom_rev <- pair_geometry(src, snk)
P_rev    <- source_pressure(kernel_weights(geom_rev, lam_best), A_rev, tw$woy)
grid_screen(dd %>% filter(zone_30y == 2, train), list(reverse = P_rev), tw, base_terms = anom_vars)

# (v) Model-filled vs observed-only Southern field (same rows)
grid_screen(obs_tr, list(source_field = P_best,
                         observed_only = source_pressure(kernel_weights(geom, lam_best), A_obs, tw$woy)),
            tw, base_terms = c("zone_fac", anom_vars))

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
                        anom_sd = apply(A[, tw$year %in% yrs_train], 1, sd, na.rm = TRUE)) %>%
  mutate(influence = share * anom_sd) %>% left_join(sites, by = "location2")

# =============================================================================
# 10. Out-of-sample evaluation: 2020 ---------------------------------------------
#   No 2020 count entered any fit, screen or kernel choice above. 2020
#   Southern traps enter only as observed source anomalies (relative to the
#   training-fit expectation); 2020 Northern/Transitional counts are only
#   predicted. The 2020 year effect is a new level, drawn from its fitted
#   distribution, and new sites likewise (sample_new_levels = "gaussian").
#   Note: P_lag0 is the same-week Southern anomaly and own_lag4 is the sink's
#   own count 4 weeks earlier, so this scores a nowcast with recent data, not a
#   season-ahead forecast.
# =============================================================================
zone_lab <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")
new_lv   <- list(allow_new_levels = TRUE, sample_new_levels = "gaussian")

# ---- 10a. Stage 1 alone: 2020 weekly counts, all zones --------------------------
d_test <- d_fit %>% filter(year == yr_test)
pp1  <- do.call(posterior_predict, c(list(fit1, newdata = d_test), new_lv))
lpd1 <- lpd_pointwise(do.call(log_lik, c(list(fit1, newdata = d_test), new_lv)))
oos_stage1 <- score_by(d_test$log_CEW, pp1, lpd1, zone_lab[as.character(d_test$zone_30y)])
oos_stage1

# ---- 10b. Sinks: does Southern source pressure improve 2020 predictions? -----
#   stage1   : seasonal model alone
#   local    : + stage 2 without Southern pressure (local weather, own lag)
#   southern : + Southern source pressure (fit2)
# Scored on the anomaly scale; log_CEW = mu_clim + anomaly is a fixed shift, so
# the scores on the count scale are identical. Same rows for all three.
test_b <- obs_b %>% filter(year == yr_test)
i1 <- match(paste(test_b$location2, test_b$woy), paste(d_test$location2, d_test$woy))
pred_b <- list(
  stage1   = sweep(pp1[, i1, drop = FALSE], 2, test_b$mu_clim),
  local    = do.call(posterior_predict, c(list(fit2_base, newdata = test_b), new_lv)),
  southern = do.call(posterior_predict, c(list(fit2, newdata = test_b), new_lv)))
lpd_b <- list(
  stage1   = lpd1[i1],
  local    = lpd_pointwise(do.call(log_lik, c(list(fit2_base, newdata = test_b), new_lv))),
  southern = lpd_pointwise(do.call(log_lik, c(list(fit2, newdata = test_b), new_lv))))
grp_b <- zone_lab[as.character(test_b$zone_30y)]
oos_sink <- bind_rows(lapply(names(pred_b), function(m)
  score_by(test_b$anom, pred_b[[m]], lpd_b[[m]], grp_b) %>% mutate(model = m, .before = 1))) %>%
  arrange(group, model)
oos_sink

# ELPD gain from Southern pressure. 2020 observations share the year and site
# effects, so this se (independent obs) is optimistic.
dl <- lpd_b$southern - lpd_b$local
bind_rows(lapply(c(split(seq_along(dl), grp_b), list(all = seq_along(dl))), function(i)
  data.frame(n = length(i), elpd_diff = sum(dl[i]), se = sqrt(length(i)) * sd(dl[i]))), .id = "group")

# 2020 weekly zone means: observed vs predicted (log(CEW + 1))
oos_wk <- bind_cols(test_b %>% select(zone_30y, woy, observed = log_CEW),
                    lapply(pred_b, function(p) test_b$mu_clim + colMeans(p))) %>%
  pivot_longer(-c(zone_30y, woy), names_to = "series") %>%
  group_by(zone = zone_lab[as.character(zone_30y)], woy, series) %>%
  summarise(value = mean(value), .groups = "drop")
ggplot(oos_wk, aes(woy, value, colour = series, linetype = series == "observed")) +
  geom_line() + facet_wrap(~ zone) + guides(linetype = "none") +
  labs(x = "ISO week (2020)", y = "Mean log(CEW + 1)", colour = NULL) +
  theme_minimal()
