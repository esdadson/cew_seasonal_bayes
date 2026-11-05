# =============================================================================
# step1/step1_models.R
# Step 1: weekly spatial Bayesian regression of log(CEW + 1) - four nested
# candidate models fit to 2008-2019 and compared with loo_compare, R^2 and
# MAE, and validated out of sample on the held-out year 2020. The selected
# one is predicted over weeks 1-52 of every monitored site-year, and the peak
# week of every site-year is taken from each posterior draw.
# Sourced by step1/step1_hpc.R (the Slurm driver: fits, loo, metrics,
# predictions). Run from the project directory.
#
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)
#  fixed        : mu = a_z + trap + x_itw' beta + u_t
#                 (zone and trap type fixed; x = weekly EVI and GDD, crop
#                 shares; u_t ~ N(0, sd_year), so a held-out year can be predicted)
#  spatial      : fixed + f(x_i, y_i), a spatial error term: approximate
#                 (Hilbert-space) GP, exponential kernel, on CONUS Albers km
#  fourier      : spatial + sum_k F_k(w) (b_k + b_zk)
#                 F_1..4 = the 2-harmonic Fourier basis of week (orthogonalised,
#                 step1_data.R section 5); b_k = seasonal curve shared by all
#                 sites, b_zk ~ N(0, sd_k) = zone deviations (|| = uncorrelated)
#  fourier_site : fourier + sum_k F_k(w) v_ik, v_ik ~ N(0, sd'_k) = site curve
#                 offsets (site-specific). Every 2020 site was trapped in
#                 2008-2019, so its offsets are estimated.
#  No site x year terms: thousands of coefficients from ~15 weeks each would
#  overparameterise the model. Each added block must earn its place in
#  loo_compare (the simplest model within 2 SE of the best is selected).
# =============================================================================

source("step1/step1_data.R")
library(brms); library(posterior)
options(brms.backend = "cmdstanr")

step1_dir <- "results/step1"          # every step-1 output
dir.create(step1_dir, recursive = TRUE, showWarnings = FALSE)

# Out-of-sample validation year: the models are fit to the other years and
# predict this one (STEP1_HOLDOUT=none fits every year, with no validation)
ho <- Sys.getenv("STEP1_HOLDOUT", "2020")
holdout <- if (ho %in% c("", "none")) NA_integer_ else as.integer(ho)
d_fit  <- if (is.na(holdout)) d_obs else d_obs %>% filter(year != holdout)   # training rows
d_test <- if (is.na(holdout)) d_obs[0, ] else d_obs %>% filter(year == holdout)

# =============================================================================
# Candidate models ----------------------------------------------------------------
# =============================================================================
gp_k <- as.integer(Sys.getenv("STEP1_GP_K", 15))    # basis functions per dimension (k^2 in all)
# Exponential (Matern 1/2) kernel, as in fit_gp2.1. Its spectral density
# decays slowly, so short length-scales need many basis functions: raise
# STEP1_GP_K if lscale piles up near the bottom of its prior.
gp_term <- paste0("gp(x_km, y_km, k = ", gp_k, ", c = 5/4, cov = 'exponential')")

year_term <- "(1 | year_f)"
terms_fixed   <- c("0 + zone_fac + trap_type", year_term,
                   "EVI_filled + gdd + corn + cotton + soybeans + tobacco + peanuts + other")
terms_spatial <- c(terms_fixed, gp_term)
terms_fourier <- c(terms_spatial, "f1 + f2 + f3 + f4", "(0 + f1 + f2 + f3 + f4 || zone_fac)")
terms_fourier_site <- c(terms_fourier, "(0 + f1 + f2 + f3 + f4 || location2)")
step1_terms <- list(fixed = terms_fixed, spatial = terms_spatial, fourier = terms_fourier,
                    fourier_site = terms_fourier_site)       # simplest first

step1_formula <- function(model)
  bf(as.formula(paste("log_CEW ~", paste(step1_terms[[model]], collapse = " + "))))
has_gp <- function(model) gp_term %in% step1_terms[[model]]
has_term <- function(model, term) term %in% step1_terms[[model]]

# brms drops rows with any NA in a formula variable; fail here, naming the
# column, rather than at compile time with "all observations removed"
na_vars <- unique(unlist(lapply(names(step1_terms), function(m) all.vars(step1_formula(m)$formula))))
na_vars <- na_vars[colSums(is.na(d_obs[intersect(na_vars, names(d_obs))])) > 0]
if (length(na_vars)) stop("NA in model variables of the observed rows: ", paste(na_vars, collapse = ", "))

# brms scales the GP inputs so the largest distance between sites (~3200 km)
# is 1: a length-scale of 0.1 is ~320 km. inv_gamma(4, 0.5) puts 90% of its
# mass between ~0.06 and ~0.36 (~190-1150 km). The lscale prior needs its
# coef: a class-only lscale prior is silently replaced by brms' default.
# Seasonal curve: the f's have unit sd, so normal(0, 2) allows a seasonal
# swing of several log units. Zone deviations: three zones barely inform a
# standard deviation, so its prior keeps them to fractions of the shared
# curve (half-t, scale 0.5) - weakly informative, against overfitting.
# brms rejects priors on classes a model doesn't have, so each model gets
# only its own.
step1_priors <- function(model) {
  p <- c(prior(normal(0, 1), class = "b"),
         prior(normal(0, 3), class = "b", coef = "zone_fac0"),
         prior(normal(0, 3), class = "b", coef = "zone_fac1"),
         prior(normal(0, 3), class = "b", coef = "zone_fac2"),
         prior(student_t(3, 0, 1), class = "sigma"),
         prior(student_t(3, 0, 1), class = "sd"))                     # year (and sites)
  if (has_gp(model)) p <- c(p, prior(exponential(1), class = "sdgp"),
                            prior(inv_gamma(4, 0.5), class = "lscale", coef = "gpx_kmy_km"))
  if (has_term(model, "f1 + f2 + f3 + f4"))
    p <- c(p, prior(normal(0, 2), class = "b", coef = "f1"), prior(normal(0, 2), class = "b", coef = "f2"),
           prior(normal(0, 2), class = "b", coef = "f3"), prior(normal(0, 2), class = "b", coef = "f4"),
           prior(student_t(3, 0, 0.5), class = "sd", group = "zone_fac"))
  p
}

# Output files
step1_files <- function(model)
  list(fit = file.path(step1_dir, paste0("fit_", model)),          # brms adds .rds
       loo = file.path(step1_dir, paste0("loo_", model, ".rds")),
       metrics = file.path(step1_dir, paste0("metrics_", model, ".rds")),
       csv = file.path(step1_dir, "stan_csv", model))

# =============================================================================
# Sampler -------------------------------------------------------------------------
#   15000 iterations per chain, the first 7500 warmup, every 30th draw kept:
#   4 chains x 250 = 1000 posterior draws. Thinning keeps the saved fit, loo
#   and predictions small; the long chains are for convergence. Override with
#   STEP1_ITER / STEP1_WARMUP / STEP1_THIN (step1/step1_hpc.R).
#   threads > 1 splits each chain's likelihood over that many cores.
#   max_treedepth 12 allows up to 4095 leapfrog steps per iteration (10: 1023):
#   fewer iterations hit the cap, but each can take longer.
# =============================================================================
step1_control <- list(adapt_delta = 0.95, max_treedepth = 12)
fit_step1 <- function(model, data = d_fit, chains = 4, iter = 15000, warmup = 7500, thin = 30,
                      threads = 1, file = step1_files(model)$fit, file_refit = "never", ...) {
  brm(step1_formula(model), data = droplevels(data), family = gaussian(), prior = step1_priors(model),
      chains = chains, cores = chains, iter = iter, warmup = warmup, thin = thin,
      threads = if (threads > 1) threading(threads) else NULL,
      seed = 42, init = 0, control = step1_control,
      file = file, file_refit = file_refit, ...)
}

# =============================================================================
# Predicting a year the fit never saw ------------------------------------------------
#   Its year effect is unknown: per posterior draw, u_new ~ N(0, sd_year),
#   shared by every site that year. (brms' sample_new_levels would redraw it
#   in every posterior_epred() call, i.e. per chunk of sites.) Rows of the
#   training years use their estimated year effects.
# =============================================================================
new_year_draws <- function(fit, seed = 42) {
  s <- as_draws_df(fit, variable = "sd_year_f__Intercept")$sd_year_f__Intercept
  set.seed(seed); s * rnorm(length(s))
}
# Group-level terms of a model without the year term (re_formula for new years)
re_no_year <- function(model) {
  re <- setdiff(grep("|", step1_terms[[model]], fixed = TRUE, value = TRUE), year_term)
  if (length(re)) as.formula(paste("~", paste(re, collapse = " + "))) else NA
}
# E[log(CEW + 1)] (draws x rows) for rows of seen and unseen years
epred_rows <- function(fit, model, nd, u_new = new_year_draws(fit)) {
  seen <- nd$year_f %in% levels(fit$data$year_f)
  out  <- matrix(NA_real_, ndraws(fit), nrow(nd))
  # droplevels: brms reads unused factor levels as new
  if (any(seen)) out[, seen] <- posterior_epred(fit, newdata = droplevels(nd[seen, ]))
  if (any(!seen))
    out[, !seen] <- posterior_epred(fit, newdata = droplevels(nd[!seen, ]), re_formula = re_no_year(model),
                                    allow_new_levels = TRUE) + u_new     # u_new[s] added to draw s
  out
}

# =============================================================================
# Fit and out-of-sample metrics ----------------------------------------------------
#   Training years (in-sample and leave-one-out):
#     bayes_R2 = Bayesian R^2 (Gelman et al. 2019), var(mu) / (var(mu) + sigma^2)
#     loo_R2   = the same with leave-one-out predictions (overfitting-robust)
#     MAE, RMSE of the posterior mean prediction; mae_loo = with LOO means
#   Held-out year (every 2020 trap-week):
#     elpd_test = sum log mean_s p(y | draw s): the out-of-sample counterpart of
#                 elpd_loo (larger = better)
#     MAE, RMSE, R2 = 1 - SSE / SST around the held-out mean
#     cov90 = share of observations inside the 90% posterior predictive interval
#   lo: loo(fit, save_psis = TRUE). Returns list(metrics = one row,
#   test_rows = held-out predictions, test_zone_week = zone weekly means).
# =============================================================================
step1_metrics <- function(fit, model, lo, d_test) {
  y   <- fit$data$log_CEW
  ep  <- posterior_epred(fit)
  mu  <- colMeans(ep)
  br2 <- bayes_R2(fit, summary = TRUE)
  lr2 <- loo_R2(fit, summary = TRUE)
  mu_loo <- loo_predict(fit, type = "mean", psis_object = lo$psis_object)
  est <- lo$estimates
  m <- data.frame(model = model, n_train = length(y),
                  elpd_loo = est["elpd_loo", 1], elpd_loo_se = est["elpd_loo", 2],
                  p_loo = est["p_loo", 1], looic = est["looic", 1],
                  pareto_k_gt_0.7 = sum(lo$diagnostics$pareto_k > 0.7),
                  bayes_R2 = br2[1, "Estimate"], bayes_R2_lo = br2[1, "Q2.5"], bayes_R2_hi = br2[1, "Q97.5"],
                  loo_R2 = lr2[1, "Estimate"], loo_R2_lo = lr2[1, "Q2.5"], loo_R2_hi = lr2[1, "Q97.5"],
                  mae_train = mean(abs(y - mu)), rmse_train = sqrt(mean((y - mu)^2)),
                  mae_loo = mean(abs(y - mu_loo)))
  rm(ep)
  if (!nrow(d_test)) return(list(metrics = m, test_rows = NULL, test_zone_week = NULL))

  # Held-out year
  et  <- epred_rows(fit, model, d_test)
  sig <- as_draws_df(fit, variable = "sigma")$sigma
  yt  <- d_test$log_CEW
  ll  <- dnorm(matrix(yt, nrow(et), ncol(et), byrow = TRUE), et, sig, log = TRUE)
  lpd <- matrixStats::colLogSumExps(ll) - log(nrow(ll))
  yrep <- et + sig * matrix(rnorm(length(et)), nrow(et))          # posterior predictive
  pi90 <- matrixStats::colQuantiles(yrep, probs = c(0.05, 0.95))
  mt <- colMeans(et)
  m <- cbind(m, data.frame(n_test = length(yt), elpd_test = sum(lpd), elpd_test_se = sd(lpd) * sqrt(length(lpd)),
                           mae_test = mean(abs(yt - mt)), rmse_test = sqrt(mean((yt - mt)^2)),
                           R2_test = 1 - sum((yt - mt)^2) / sum((yt - mean(yt))^2),
                           cov90_test = mean(yt >= pi90[, 1] & yt <= pi90[, 2])))
  test_rows <- d_test %>% select(location2, loc_year, zone_30y, year, woy, log_CEW) %>%
    mutate(model = model, pred_mean = mt, pi_lo = pi90[, 1], pi_hi = pi90[, 2])
  # Zone weekly mean over the sites trapped that week, per draw
  g  <- paste(d_test$zone_30y, d_test$woy)
  Z  <- sapply(split(seq_along(g), g), function(j) rowMeans(et[, j, drop = FALSE]))
  zw <- data.frame(key = colnames(Z), t(apply(Z, 2, q95))) %>%
    separate(key, c("zone_30y", "woy"), convert = TRUE) %>%
    left_join(d_test %>% group_by(zone_30y, woy) %>% summarise(obs_mean = mean(log_CEW), n_sites = n(), .groups = "drop"),
              by = c("zone_30y", "woy")) %>%
    mutate(model = model, year = holdout)
  list(metrics = m, test_rows = test_rows, test_zone_week = zw)
}

# =============================================================================
# Predictions over weeks 1-52 and peak weeks ---------------------------------------
#   E[log(CEW + 1)] for every row of d_all: weeks 1-52 of every monitored
#   site-year, sampled or not, with each week's own EVI and GDD - including
#   the held-out year, predicted out of sample.
#   Peak week of a site-year = the week of its largest E[log(CEW + 1)], per
#   draw, searched within the zone's 5th-95th percentile of sampled weeks
#   (outside it the curve is extrapolation). in_sampled = the peak falls
#   within the site-year's own first-last sampled week; low values mean the
#   peak is placed by pooling, not by that site-year's counts.
# =============================================================================
peak_win <- d_obs %>% group_by(zone_30y) %>%
  summarise(wk_lo = floor(quantile(woy, 0.05)), wk_hi = ceiling(quantile(woy, 0.95)))
zone_lab <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")

# Returns
#   pred       : d's keys + pred_mean, pred_sd, pred_lo / pred_med / pred_hi (95%)
#   draws      : rows x n_keep posterior draws (step 2: step-1 uncertainty)
#   peaks      : list(draws = draws x site-years peak weeks, site_years)
#   zone_curve : draws x zone-year-weeks, mean E[log(CEW + 1)] over the zone's site-years
# Chunked over sites (each site-year's 52 weeks in one call); memory ~ draws x chunk rows.
# The held-out year is predicted out of sample (epred_rows: its year effect is
# drawn, shared by all its sites).
predict_step1 <- function(fit, model, d = d_all, n_keep = 100, sites_per_chunk = 25) {
  u_new <- new_year_draws(fit)
  d <- d %>% arrange(location2, year, woy) %>%
    group_by(loc_year) %>%
    mutate(samp_lo = min(woy[observed]), samp_hi = max(woy[observed])) %>% ungroup() %>%
    left_join(peak_win, by = "zone_30y")
  nd   <- ndraws(fit)
  keep <- sort(sample(nd, min(n_keep, nd)))
  sy   <- d %>% distinct(loc_year, location2, year, zone_30y, samp_lo, samp_hi)
  zw   <- paste(d$zone_30y, d$year, d$woy)                  # zone-year-week of each row
  zkey <- sort(unique(zw))
  z_sum <- matrix(0, nd, length(zkey), dimnames = list(NULL, zkey))
  z_n   <- table(factor(zw, zkey))
  P <- matrix(NA_integer_, nd, nrow(sy), dimnames = list(NULL, sy$loc_year))
  s <- matrix(NA_real_, nrow(d), 5)
  D <- matrix(NA_real_, nrow(d), length(keep))
  locs <- unique(d$location2)
  for (b in split(locs, ceiling(seq_along(locs) / sites_per_chunk))) {
    i  <- which(d$location2 %in% b)
    ep <- epred_rows(fit, model, d[i, ], u_new)
    s[i, ] <- cbind(colMeans(ep), matrixStats::colSds(ep),
                    matrixStats::colQuantiles(ep, probs = c(0.025, 0.5, 0.975)))
    D[i, ] <- t(ep[keep, , drop = FALSE])
    zi <- factor(zw[i], zkey)
    z_sum <- z_sum + ep %*% Matrix::sparseMatrix(seq_along(i), as.integer(zi), x = 1,
                                                 dims = c(length(i), length(zkey)))
    di <- d[i, ]
    for (k in unique(di$loc_year)) {
      j <- which(di$loc_year == k & di$woy >= di$wk_lo & di$woy <= di$wk_hi)
      if (!length(j)) next
      P[, k] <- as.integer(di$woy[j][max.col(ep[, j, drop = FALSE], ties.method = "first")])
    }
  }
  colnames(s) <- c("pred_mean", "pred_sd", "pred_lo", "pred_med", "pred_hi")
  list(pred = bind_cols(d %>% select(location2, loc_year, zone_30y, year, woy, observed, log_CEW), as.data.frame(s)),
       draws = D,
       peaks = list(draws = P, site_years = sy),
       zone_curve = as.matrix(z_sum) / rep(as.numeric(z_n), each = nd))
}

q95 <- function(x) c(med = median(x), lo = quantile(x, 0.025, names = FALSE),
                     hi = quantile(x, 0.975, names = FALSE))

# Peak summaries for maps and tables. Every summary is computed draw by draw,
# then summarised, so its interval carries the uncertainty of every site-year
# peak in it. Zone-year = mean over the zone's monitored site-years that year
# (the set of sites changes between years, so part of a year-to-year change
# can be which sites were trapped).
peak_summaries <- function(pk) {
  P <- pk$draws; sy <- pk$site_years
  site_year <- sy %>%
    mutate(peak_med = apply(P, 2, median), peak_mean = colMeans(P),
           peak_lo = apply(P, 2, quantile, 0.025), peak_hi = apply(P, 2, quantile, 0.975),
           peak_week_sd = apply(P, 2, sd), ci_width = peak_hi - peak_lo,
           in_sampled = colMeans(sweep(P, 2, samp_lo, ">=") & sweep(P, 2, samp_hi, "<="))) %>%
    left_join(sites %>% select(location2, location, longitude, latitude, x_km, y_km), by = "location2")

  ZY <- vapply(split(seq_len(ncol(P)), paste(sy$zone_30y, sy$year)),
               function(j) rowMeans(P[, j, drop = FALSE]), numeric(nrow(P)))
  zone_year <- data.frame(key = colnames(ZY), t(apply(ZY, 2, q95))) %>%
    separate(key, c("zone_30y", "year"), convert = TRUE) %>%
    left_join(count(sy, zone_30y, year, name = "n_sites"), by = c("zone_30y", "year")) %>%
    mutate(zone = zone_lab[as.character(zone_30y)])

  # Delay of the Transitional and Northern zones behind the Southern zone, per
  # year and averaged over years (P_later = posterior P(delay > 0))
  zone_delay <- bind_rows(lapply(c(1, 0), function(z) {
    ys <- intersect(sy$year[sy$zone_30y == z], sy$year[sy$zone_30y == 2])
    Dz <- sapply(ys, function(y) ZY[, paste(z, y)] - ZY[, paste(2, y)])
    colnames(Dz) <- ys
    Dz <- cbind(Dz, all_years = rowMeans(Dz))
    data.frame(zone = zone_lab[as.character(z)], year = colnames(Dz), t(apply(Dz, 2, q95)),
               P_later = colMeans(Dz > 0), row.names = NULL)
  }))

  # Site peak, averaged over the site's monitored years (per draw)
  site <- data.frame(location2 = sy$location2, j = seq_len(ncol(P))) %>%
    group_by(location2) %>% summarise(j = list(j), n_years = n()) %>%
    mutate(s = lapply(j, function(jj) q95(rowMeans(P[, jj, drop = FALSE])))) %>%
    unnest_wider(s) %>% select(-j) %>%
    left_join(sites, by = "location2")

  list(site_year = site_year, zone_year = zone_year, zone_delay = zone_delay, site = site)
}

# Zone mean predicted curve per year and week (95% CrI), from predict_step1()$zone_curve
zone_curve_summary <- function(Z) {
  data.frame(key = colnames(Z), t(apply(Z, 2, q95))) %>%
    separate(key, c("zone_30y", "year", "woy"), convert = TRUE) %>%
    mutate(zone = zone_lab[as.character(zone_30y)])
}
