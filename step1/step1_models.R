# =============================================================================
# step1/step1_models.R
# Step 1: weekly spatial Bayesian regression of log(CEW + 1), 2008-2020 -
# three nested candidate models compared with loo_compare; the selected one
# is predicted over weeks 1-52 of every monitored site-year, and the peak
# week of every site-year is taken from each posterior draw.
# Sourced by step1/step1_hpc.R (the Slurm driver: fits, loo, predictions).
# Run from the project directory.
#
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)
#  fixed   : mu = a_z + a_t + trap + x_itw' beta
#            (zone, year and trap type as fixed effects; x = weekly EVI and
#             GDD, crop shares)
#  spatial : fixed + f(x_i, y_i), a spatial error term: approximate
#            (Hilbert-space) GP, exponential kernel, on CONUS Albers km
#  fourier : spatial + sum_k F_k(w) (b_zk + v_ik)
#            F_1..4 = the 2-harmonic Fourier basis of week (orthogonalised,
#            step1_data.R section 5); b_z = zone seasonal curve (shared over
#            space), v_i = site curve offset (site-specific; || = uncorrelated)
#  simple  : (not nested in the others; from cew_fourier_covariates.R)
#            mu = a + x_itw' beta + sum_k F_k(w) b_zk + f(x_i, y_i) + a_t + trap
#            zone curves as uncorrelated zone-level slopes, year and trap type
#            as random intercepts, no zone intercepts or site curves
#  No site x year terms: thousands of coefficients from ~15 weeks each would
#  overparameterise the model.
# =============================================================================

source("step1/step1_data.R")
library(brms); library(posterior)
options(brms.backend = "cmdstanr")

step1_dir <- "results/step1"          # every step-1 output
dir.create(step1_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Candidate models ----------------------------------------------------------------
# =============================================================================
gp_k <- as.integer(Sys.getenv("STEP1_GP_K", 15))    # basis functions per dimension (k^2 in all)
# Exponential (Matern 1/2) kernel, as in fit_gp2.1. Its spectral density
# decays slowly, so short length-scales need many basis functions: raise
# STEP1_GP_K if lscale piles up near the bottom of its prior.
gp_term <- paste0("gp(x_km, y_km, k = ", gp_k, ", c = 5/4, cov = 'exponential')")

terms_fixed   <- c("0 + zone_fac + year_f + trap_type",
                   "EVI_filled + gdd + corn + cotton + soybeans + tobacco + peanuts + other")
terms_spatial <- c(terms_fixed, gp_term)
terms_fourier <- c(terms_spatial, "zone_fac:(f1 + f2 + f3 + f4)", "(0 + f1 + f2 + f3 + f4 || location2)")
# fit_gp's formula with f1-f4 for s1, c1, s2, c2 (the same curves; section 5
# of step1_data.R) and the approximate GP for its exact one
terms_simple  <- c("EVI_filled + gdd + corn + cotton + soybeans + tobacco + peanuts + other",
                   "(0 + f1 + f2 + f3 + f4 || zone_fac)", gp_term, "(1 | year_f) + (1 | trap_type)")
step1_terms <- list(fixed = terms_fixed, spatial = terms_spatial, fourier = terms_fourier,
                    simple = terms_simple)

step1_formula <- function(model)
  bf(as.formula(paste("log_CEW ~", paste(step1_terms[[model]], collapse = " + "))))
has_gp <- function(model) gp_term %in% step1_terms[[model]]
has_sd <- function(model) any(grepl("|", step1_terms[[model]], fixed = TRUE))
has_zone <- function(model) any(startsWith(step1_terms[[model]], "0 + zone_fac"))

# brms drops rows with any NA in a formula variable; fail here, naming the
# column, rather than at compile time with "all observations removed"
na_vars <- unique(unlist(lapply(names(step1_terms), function(m) all.vars(step1_formula(m)$formula))))
na_vars <- na_vars[colSums(is.na(d_obs[intersect(na_vars, names(d_obs))])) > 0]
if (length(na_vars)) stop("NA in model variables of the observed rows: ", paste(na_vars, collapse = ", "))

# brms scales the GP inputs so the largest distance between sites (~3200 km)
# is 1: a length-scale of 0.1 is ~320 km. inv_gamma(4, 0.5) puts 90% of its
# mass between ~0.06 and ~0.36 (~190-1150 km). The lscale prior needs its
# coef: a class-only lscale prior is silently replaced by brms' default.
# brms rejects priors on classes a model doesn't have, so each model gets
# only its own.
step1_priors <- function(model) {
  p <- c(prior(normal(0, 1), class = "b"), prior(student_t(3, 0, 1), class = "sigma"))
  if (has_zone(model)) p <- c(p, prior(normal(0, 3), class = "b", coef = "zone_fac0"),
                              prior(normal(0, 3), class = "b", coef = "zone_fac1"),
                              prior(normal(0, 3), class = "b", coef = "zone_fac2"))
  if (has_gp(model)) p <- c(p, prior(exponential(1), class = "sdgp"),
                            prior(inv_gamma(4, 0.5), class = "lscale", coef = "gpx_kmy_km"))
  if (has_sd(model)) p <- c(p, prior(student_t(3, 0, 1), class = "sd"))
  p
}

# Output files
step1_files <- function(model)
  list(fit = file.path(step1_dir, paste0("fit_", model)),          # brms adds .rds
       loo = file.path(step1_dir, paste0("loo_", model, ".rds")),
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
fit_step1 <- function(model, data = d_obs, chains = 4, iter = 15000, warmup = 7500, thin = 30,
                      threads = 1, file = step1_files(model)$fit, file_refit = "never", ...) {
  brm(step1_formula(model), data = droplevels(data), family = gaussian(), prior = step1_priors(model),
      chains = chains, cores = chains, iter = iter, warmup = warmup, thin = thin,
      threads = if (threads > 1) threading(threads) else NULL,
      seed = 42, init = 0, control = step1_control,
      file = file, file_refit = file_refit, ...)
}

# =============================================================================
# Predictions over weeks 1-52 and peak weeks ---------------------------------------
#   E[log(CEW + 1)] for every row of d_all: weeks 1-52 of every monitored
#   site-year, sampled or not, with each week's own EVI and GDD.
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
predict_step1 <- function(fit, d = d_all, n_keep = 100, sites_per_chunk = 25) {
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
    ep <- posterior_epred(fit, newdata = droplevels(d[i, ]))   # brms reads unused levels as new
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
