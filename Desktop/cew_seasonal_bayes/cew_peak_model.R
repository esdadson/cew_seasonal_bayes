# =============================================================================
# cew_peak_model.R
# Step 1 of cew_peak_granger.R: the weekly spatial model of log(CEW + 1) from
# cew_fourier_covariates.R, restructured so that
#   (a) the peak week of every site-year is estimable from the posterior, and
#   (b) a fit takes hours, not the ~164 h per chain of the exact-GP model.
# Sourced by cew_peak_granger.R and by cew_peak_hpc.R (the Slurm driver), so
# the HPC fits and the analysis use identical data.
#
# Data: sections 0-4 of cew_stage1_data.R (site table, ISO-week counts,
# Daymet, CDL, EVI, prediction grid), with every year 2008-2020 in the fit -
# a peak is wanted for every site-year. Out-of-sample skill comes from
# leave-one-year-out refits (cew_peak_hpc.R loyo), not from a held-out 2020.
#
# Changes from fit_gp2.1, and why
#   - Exact GP -> Hilbert-space approximate GP (gp(..., k = )) on projected
#     km. The exact GP factorises a 422 x 422 covariance matrix at every
#     gradient (up to 1023 per iteration); the approximation's cost is linear
#     in sites x basis functions. Degrees of longitude shrink with latitude,
#     so lon/lat is not isotropic; x_km/y_km (CONUS Albers) are.
#   - Zone curves as population-level effects. A variance across 3 zones is
#     not estimable (the old (s1 + ... || zone_fac) term).
#   - Zone shares seasonality over space; location adds site-specific
#     curve offsets (|| = uncorrelated). No location x year terms (thousands
#     of coefficients from ~15 weeks each would overparameterise the model):
#     a site-year's peak differs from the site's other years through its own
#     weekly EVI/GDD and, in the zone_year candidate, a zone-year curve shift.
#     The candidates are compared with loo_compare (see Model below).
#   - Trap type as a fixed effect (a variance from 3 groups is unidentified).
#   - The Fourier basis is centred and orthogonalised over the observed
#     weeks (f1-f4 = a fixed linear transform of s1, c1, s2, c2). Over the
#     trapping season (weeks ~20-38) the raw columns barely vary and sit far
#     from 0 (cos(2 pi 29/52) = -0.94), so each coefficient is nearly
#     collinear with the intercept at its level: NUTS then needs the maximum
#     1023 tiny steps per iteration (both the exact-GP model and the raw-basis
#     version of this one did). The curves - and so the peaks - are the same
#     functions of week either way; only the parameterisation changes.
# =============================================================================

source("cew_stage1_data.R")

# Standardise on every observed year (the stage-1 data used 2008-2019 only).
# Covariates only, so the LOYO refits can share this scaling without leaking
# the held-out year's counts.
# Site id stays integer and loc_year (site-year key for the peaks) character:
# brms reads the unused levels of a factor in newdata as new levels.
scale_my <- make_scaler(d_all %>% filter(observed), my_vars)
d_all <- d_all %>% scale_my() %>%
  mutate(loc_year = paste(location2, year, sep = "_"))

# Centred, orthogonal Fourier basis with unit sd over the observed rows:
# B_c = Q R (QR of the centred basis), f = B_c R^-1 sqrt(n - 1). Computed
# from weeks only (no counts), so the LOYO refits can share it.
four  <- c("s1", "c1", "s2", "c2")
B_obs <- as.matrix(d_all[d_all$observed, four])
B_mu  <- colMeans(B_obs)
B_R   <- qr.R(qr(sweep(B_obs, 2, B_mu))) / sqrt(nrow(B_obs) - 1)
add_orth_fourier <- function(d) {
  d[paste0("f", 1:4)] <- sweep(as.matrix(d[, four]), 2, B_mu) %*% solve(B_R)
  d
}
d_all <- add_orth_fourier(d_all)
d_obs <- d_all %>% filter(observed)

# =============================================================================
# Model ----------------------------------------------------------------------
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)
#  mu_itw = a_z + x_itw' beta + f(x_i, y_i) + a_t
#           + sum_k F_k(w) (b_zk + v_ik [+ u_ztk])
#  F_1..4(w) = f1-f4: the 2-harmonic Fourier basis (sin/cos(2 pi h w / 52),
#              h = 1, 2), centred and orthogonalised (above)
#  b_z: zone curve; v_i: site offset; u_zt: zone-year shift (zone_year only)
#  f ~ approximate GP (exponential kernel) on km coordinates
#
# Candidates (PEAK_MODEL), each fit with cew_peak_hpc.R and compared in
# cew_peak_granger.R section 3:
#   site      : zone curves + site offsets + GP + year level
#   zone_year : site + zone-year curve shifts (39 zone-years x 4 coefficients),
#               letting a zone's peak move between years beyond EVI/GDD.
#   loo_compare on PSIS-LOO (leave one week out) asks whether the zone-year
#   shifts improve the fit to the observed seasons; the same comparison on
#   the leave-one-year-out folds asks whether they overfit - they cannot
#   inform a year that was not in the fit.
# =============================================================================
gp_k <- as.integer(Sys.getenv("PEAK_GP_K", 15))    # basis functions per dimension (k^2 in all)
# Exponential (Matern 1/2) kernel, as in fit_gp2.1. Its spectral density
# decays slowly, so short length-scales need many basis functions: raise
# PEAK_GP_K if lscale_1 piles up near the bottom of its prior.
gp_term <- paste0("gp(x_km, y_km, k = ", gp_k, ", c = 5/4, cov = 'exponential')")

peak_terms <- list(
  site      = "(0 + f1 + f2 + f3 + f4 || location2)",
  zone_year = "(0 + f1 + f2 + f3 + f4 || location2) + (0 + f1 + f2 + f3 + f4 || zone_fac:year_f)")
make_f_peak <- function(model) bf(as.formula(paste(
  "log_CEW ~ 0 + zone_fac + zone_fac:(f1 + f2 + f3 + f4) +",
  "EVI_filled + gdd + corn + cotton + soybeans + tobacco + peanuts + other + trap_type +",
  gp_term, "+ (1 | year_f) +", peak_terms[[model]])))

peak_model <- Sys.getenv("PEAK_MODEL", "site")
stopifnot(peak_model %in% names(peak_terms))
f_peak <- make_f_peak(peak_model)
# brms drops rows with any NA in the formula's variables; fail here, naming
# the column, rather than at compile time with "all observations removed"
na_vars <- all.vars(f_peak$formula)
na_vars <- na_vars[colSums(is.na(d_obs[intersect(na_vars, names(d_obs))])) > 0]
if (length(na_vars)) stop("NA in model variables of the observed rows: ", paste(na_vars, collapse = ", "))
# Output files of the chosen candidate
peak_files <- function(model = peak_model)
  list(fit = paste0("fit_peak_", model), draws = paste0("peak_draws_", model, ".rds"),
       loo = paste0("loo_peak_", model, ".rds"), fill = paste0("fill_", model, ".rds"),
       loyo_dir = file.path("loyo", model))

# brms scales the GP inputs so the largest distance between sites (~3200 km)
# is 1: a length-scale of 0.1 is ~320 km. inv_gamma(4, 0.5) puts 90% of its
# mass between ~0.06 and ~0.36 (~190-1150 km); much shorter than 0.06 the
# k = 15 basis can't resolve the field anyway (raise PEAK_GP_K for that).
# The lscale prior needs its coef: brms gives the GP its own default prior,
# and a class-only lscale prior is silently dropped (as happened to the
# inv_gamma(5, 5) of fit_gp2.1 - its Stan code shows brms' default).
pri_peak <- c(prior(normal(0, 1), class = "b"),
              prior(normal(0, 3), class = "b", coef = "zone_fac0"),
              prior(normal(0, 3), class = "b", coef = "zone_fac1"),
              prior(normal(0, 3), class = "b", coef = "zone_fac2"),
              prior(student_t(3, 0, 1), class = "sd"),
              prior(student_t(3, 0, 1), class = "sigma"),
              prior(exponential(1), class = "sdgp"),
              prior(inv_gamma(4, 0.5), class = "lscale", coef = "gpx_kmy_km"))

# Sampler: every post-warmup draw is kept (the peaks use them all). Set the
# run length with PEAK_ITER / PEAK_WARMUP (cew_peak_hpc.R). threads > 1 splits
# each chain's likelihood over that many cores (reduce_sum).
peak_control <- list(adapt_delta = 0.95, max_treedepth = 12)
fit_peak <- function(data, file, chains = 4, iter = 2000, warmup = 1000, threads = 1,
                     file_refit = "never", ...) {
  brm(f_peak, data = droplevels(data), family = gaussian(), prior = pri_peak,   # LOYO: drop the held-out year_f level
      chains = chains, cores = chains, iter = iter, warmup = warmup,
      threads = if (threads > 1) threading(threads) else NULL,
      seed = 42, init = 0, control = peak_control,
      file = file, file_refit = file_refit, ...)
}

# =============================================================================
# Peak week per site-year, from every posterior draw --------------------------
#   E[log(CEW + 1)] is evaluated at each week of a site-year (its own weekly
#   EVI and GDD included, so this is the peak of the expected count, not of
#   the Fourier curve alone) and the argmax is taken per draw.
#   Search window: the zone's 5th-95th percentile of sampled weeks (outside
#   it the curve is extrapolation). in_sampled = the peak falls within the
#   site-year's own first-last sampled week; low values mean the peak is
#   placed by pooling, not by that site-year's counts.
# =============================================================================
peak_win <- d_obs %>% group_by(zone_30y) %>%
  summarise(wk_lo = floor(quantile(woy, 0.05)), wk_hi = ceiling(quantile(woy, 0.95)))

# Site-years and their weeks to evaluate: every week of the zone window
# (prediction-grid rows supply the covariates of unsampled weeks)
peak_grid <- function(d) {
  d %>% group_by(loc_year) %>%
    mutate(samp_lo = min(woy[observed]), samp_hi = max(woy[observed])) %>%
    ungroup() %>%
    inner_join(peak_win, by = "zone_30y") %>%
    filter(woy >= wk_lo, woy <= wk_hi) %>%
    arrange(loc_year, woy)
}

# Returns list(draws = draws x site-years integer matrix of peak weeks,
#              site_years = one row per site-year, in the columns' order).
# Chunked over site-years so memory stays ~ draws x chunk rows.
peak_draws <- function(fit, grid, chunk = 150, ...) {
  sy <- grid %>% distinct(loc_year, location2, year, zone_30y, samp_lo, samp_hi)
  keys <- sy$loc_year
  out <- matrix(NA_integer_, ndraws(fit), length(keys), dimnames = list(NULL, keys))
  for (b in split(keys, ceiling(seq_along(keys) / chunk))) {
    g  <- droplevels(grid[grid$loc_year %in% b, ])   # unused factor levels read as new
    lp <- posterior_linpred(fit, newdata = g, ...)
    for (k in b) {
      j <- which(g$loc_year == k)
      out[, k] <- as.integer(g$woy[j][max.col(lp[, j, drop = FALSE], ties.method = "first")])
    }
  }
  list(draws = out, site_years = sy)
}

# Posterior summaries per site-year
summarise_peaks <- function(pk) {
  P <- pk$draws
  pk$site_years %>%
    mutate(peak_med = apply(P, 2, median),
           peak_mean = colMeans(P),
           peak_lo = apply(P, 2, quantile, 0.025),
           peak_hi = apply(P, 2, quantile, 0.975),
           peak_week_sd = apply(P, 2, sd),
           ci_width = peak_hi - peak_lo,
           in_sampled = colMeans(sweep(P, 2, samp_lo, ">=") & sweep(P, 2, samp_hi, "<=")))
}

# =============================================================================
# Filled weekly populations (input to step 2) -----------------------------------
#   E[log(CEW + 1)] for every row of d (the prediction grid: weeks 1-52 of
#   every trapped site-year, sampled or not): posterior mean and sd, plus
#   n_keep draws (rows x draws) so step 2 can carry step-1 uncertainty.
#   Chunked over rows; with new levels (LOYO) use chunk = Inf so a held-out
#   year's zone-year shift is one draw for all its rows.
# =============================================================================
fill_draws <- function(fit, d, chunk = 5000, n_keep = 100, ...) {
  keep <- sort(sample(ndraws(fit), min(n_keep, ndraws(fit))))
  idx  <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / chunk))
  res  <- lapply(idx, function(i) {
    lp <- posterior_linpred(fit, newdata = droplevels(d[i, ]), ...)
    list(mean = colMeans(lp), sd = apply(lp, 2, sd), draws = t(lp[keep, , drop = FALSE]))
  })
  list(fill = d %>% select(location2, loc_year, zone_30y, year, woy, observed, log_CEW) %>%
         mutate(fill_mean = unlist(lapply(res, `[[`, "mean"), use.names = FALSE),
                fill_sd = unlist(lapply(res, `[[`, "sd"), use.names = FALSE)),
       draws = do.call(rbind, lapply(res, `[[`, "draws")))
}

# =============================================================================
# Leave-one-year-out evaluation (used by cew_peak_hpc.R loyo) -----------------
#   Fit on every other year; predict the held-out year's weekly counts. Its
#   year and zone-year effects (and sites trapped only that year)
#   are new levels, drawn from their fitted distributions, so this is a
#   season-ahead prediction from covariates, zone/site curves and the GP.
#   Saves scores, per-observation predictions, predicted peak draws and the
#   held-out year's filled weekly series (fill_all = all its grid rows; step
#   2 uses these for that year's predictors) - not the fit, which is large.
# =============================================================================
# loo_compare-style table from leave-one-year-out pointwise log predictive
# densities (lpd: observations x models, same rows; year: held-out year of
# each row). The SE of a difference is blocked by year (weeks of one season
# are not independent): sqrt(K) x sd of the K yearly sums, as kfold() does.
loyo_compare <- function(lpd, year) {
  el <- colSums(lpd); best <- which.max(el)
  yd <- rowsum(lpd - lpd[, best], year)
  out <- data.frame(model = colnames(lpd), elpd_loyo = el, elpd_diff = el - el[best],
                    se_diff = sqrt(nrow(yd)) * apply(yd, 2, sd),
                    years_better = colSums(yd > 0), row.names = NULL)
  out[order(-out$elpd_loyo), ]
}

loyo_eval <- function(fit, d_test, grid_test, fill_all) {
  d_test <- droplevels(d_test)
  nl <- list(allow_new_levels = TRUE, sample_new_levels = "gaussian")
  pp  <- do.call(posterior_predict, c(list(fit, newdata = d_test), nl))
  lpd <- lpd_pointwise(do.call(log_lik, c(list(fit, newdata = d_test), nl)))
  zl  <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")
  # one chunk: new levels are drawn per call, and the held-out zone-year
  # shift must be the same draw for every site-year of that zone
  pk  <- do.call(peak_draws, c(list(fit, grid_test, chunk = Inf), nl))
  fl  <- do.call(fill_draws, c(list(fit, fill_all, chunk = Inf, n_keep = 0), nl))
  list(scores = score_by(d_test$log_CEW, pp, lpd, zl[as.character(d_test$zone_30y)]),
       pred = d_test %>% select(location2, loc_year, zone_30y, year, woy, log_CEW) %>%
         mutate(pred_mean = colMeans(pp), pred_lo = apply(pp, 2, quantile, 0.05),
                pred_hi = apply(pp, 2, quantile, 0.95), lpd = lpd),
       peaks = summarise_peaks(pk), peak_draws = pk$draws, fill = fl$fill)
}
