# =============================================================================
# step1_stan/zone_harmonic.R
# Stan fit of the zone harmonic model (zone_harmonic.stan; model written out
# in zone_harmonic_common.R):
#   zone-specific harmonic seasonal curves (fixed effects)
#   + year intercept
#   + site effect = short-range Matern field + site nugget, on raw site locations
#   + (ZH_SVC = 1) site-level phase / amplitude shifts of the first harmonic
#   every site-level effect sums to zero within each zone
# No temporal knots, no splines, no merging of sites. The brms version of the
# same model is zone_harmonic_brms.R.
#
# Modes:
#   Rscript step1_stan/zone_harmonic.R dryrun    # one short chain -> time per gradient,
#                                                # projected run time
#   Rscript step1_stan/zone_harmonic.R fit       # full fit -> ZH_OUT/stan/draws.rds
#   Rscript step1_stan/zone_harmonic.R summary   # tables and figures (with brms if fitted)
#
# Environment (plus the ZH_* settings in zone_harmonic_common.R):
#   STAN_CHAINS, STAN_ITER, STAN_WARMUP, STAN_THIN   4, 2000, 1000, 2 (iter includes warmup)
#   STAN_ADAPT_DELTA, STAN_MAX_TREEDEPTH             0.9, 10
#   ZH_CENTER_SITE   1 = centered site effect w (default: ~90 rows per site)
#   ZH_OUT           output folder                                  results/zone_harmonic
#   DRY_WARMUP, DRY_SAMPLING                         dry-run chain length (150 + 50)
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "summary"))

source("step1_stan/zone_harmonic_common.R")

out_dir <- file.path(env("ZH_OUT", "results/zone_harmonic"), "stan")
csv_dir <- file.path(out_dir, "stan_csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
smp <- list(chains = env_num("STAN_CHAINS", 4), iter = env_num("STAN_ITER", 2000),
            warmup = env_num("STAN_WARMUP", 1000), thin = env_num("STAN_THIN", 2),
            adapt_delta = env_num("STAN_ADAPT_DELTA", 0.9), max_td = env_num("STAN_MAX_TREEDEPTH", 10))

# =============================================================================
# 1. Stan data -------------------------------------------------------------------
#    Priors: beta ~ N(0, 10^2) zone intercepts, N(0, 2.5^2) zone harmonics,
#    N(0, 1) trap / covariates (standardised); random-effect sds ~
#    exponential(3) (P(sd > 1) = 0.05); range lognormal (98% in
#    [ZH_RANGE_LO, ZH_RANGE_HI] km); residual sd ~ half-Student-t(3, 0, 2.5)
# =============================================================================
stan_data <- list(
  N = length(fit_idx), P = ncol(X), T = max(d$year_id), S = nrow(site_xy), Z = nlevels(d$zone),
  y = d$log_CEW[fit_idx], X = X[fit_idx, ], week = d$woy[fit_idx], year = d$year_id[fit_idx],
  site = d$site_id[fit_idx],
  xy = lapply(seq_len(nrow(site_xy)), function(i) c(site_xy$x_km[i], site_xy$y_km[i])),
  zone_site = as.integer(site_xy$zone),
  svc = svc, centered_site = as.integer(env_num("ZH_CENTER_SITE", 1)), jitter = 1e-6,
  beta_sd = beta_sd, range_mu = range_mu, range_sd = range_sd,
  lambda_sd = -log(0.05), sigma_scale = 2.5,
  N_all = nrow(d), X_all = X, week_all = d$woy, year_all = d$year_id, site_all = d$site_id,
  win_lo = win[1], win_hi = win[2],
  N_test = length(test_idx), test_idx = test_idx, y_test = d$log_CEW[test_idx])

inits <- function() list(
  beta = rep(0, ncol(X)) + ifelse(is_zone, 1.5, 0), u_z = rep(0, stan_data$T), sigma_year = 0.3,
  w_par = rep(0, stan_data$S), sigma_field = 0.3, sigma_site = 0.3,
  a_z = rep(0, svc * stan_data$S), b_z = rep(0, svc * stan_data$S), sigma_svc = as.array(rep(0.2, svc)),
  range_km = exp(range_mu) * runif(1, 0.8, 1.25), sigma = 1.2)

hyper <- c("range_km", "sigma_field", "sigma_site", if (svc) "sigma_svc[1]", "sigma_year", "sigma")

# =============================================================================
# 2. Dry run / fit -------------------------------------------------------------------
# =============================================================================
if (mode %in% c("dryrun", "fit")) {
  library(cmdstanr)
  t0  <- Sys.time()
  mod <- cmdstan_model("step1_stan/zone_harmonic.stan")
  t_compile <- mins(t0)
}

if (mode == "dryrun") {
  w <- env_num("DRY_WARMUP", 150); s <- env_num("DRY_SAMPLING", 50)
  fit <- mod$sample(data = stan_data, init = inits, chains = 1, iter_warmup = w, iter_sampling = s,
                    seed = 42, adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                    save_warmup = TRUE, refresh = 10, output_dir = csv_dir)
  tm <- fit$time()$chains
  dg <- posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))
  lf <- dg$n_leapfrog__; in_w <- seq_along(lf) <= w
  sec_grad <- tm$total / sum(lf)
  steps <- mean(lf[!in_w])
  h_obs <- smp$iter * steps * sec_grad / 3600
  message("\n--- Dry run (1 chain, ", w, " warmup + ", s, " sampling iterations) ---")
  message(sprintf("Compile %.1f min. Warmup %.1f s/iteration, sampling %.1f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps, 2^smp$max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * mean(dg$treedepth__[!in_w] >= smp$max_td), sum(dg$divergent__[!in_w])))
  print(fit$summary(hyper), digits = 3)
  message(sprintf("Projected full run (%d iterations): %.1f h per chain", smp$iter, h_obs))
}

if (mode == "fit") {
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d)", smp$chains, smp$iter, smp$warmup, smp$thin))
  t0  <- Sys.time()
  fit <- mod$sample(data = stan_data, init = inits, chains = smp$chains, parallel_chains = smp$chains,
                    iter_warmup = smp$warmup, iter_sampling = smp$iter - smp$warmup, thin = smp$thin,
                    adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                    seed = 42, refresh = 50, output_dir = csv_dir)
  message(sprintf("Sampling: %.2f h", mins(t0) / 60))
  saveRDS(list(draws = fit$draws(c(hyper, "beta", "u", "w", "alpha", "beta_s", "peak", "eta_test"),
                                 format = "draws_matrix"),
               diag = fit$diagnostic_summary(), time = fit$time(), smp = smp,
               setup = list(H = H, svc = svc, range = c(range_lo, range_hi), win = win, xcols = colnames(X),
                            n_sub = n_sub)),
          file.path(out_dir, "draws.rds"))
  print(fit$summary(hyper), digits = 3)
  message("Saved ", file.path(out_dir, "draws.rds"))
}

# =============================================================================
# 3. Summary -----------------------------------------------------------------------
#    holdout.csv, zone_stats.csv (site-year peaks), zone_curve_peaks.csv
#    (peak of each zone's harmonic curve), hyper.csv, beta.csv,
#    peak_hist.png, zone_curves.png, site maps. Includes the brms fit when
#    ZH_OUT/brms/draws.rds exists.
# =============================================================================
if (mode == "summary") {
  res <- readRDS(file.path(out_dir, "draws.rds"))
  dr  <- matrix(res$draws, nrow(res$draws), dimnames = dimnames(res$draws))
  col <- function(p) dr[, grep(paste0("^", p, "\\["), colnames(dr)), drop = FALSE]

  hyp <- posterior::summarise_draws(posterior::subset_draws(res$draws, variable = hyper))
  print(hyp, digits = 3)
  write.csv(hyp, file.path(out_dir, "hyper.csv"), row.names = FALSE)
  print(res$diag)

  Bz <- col("beta"); colnames(Bz) <- colnames(X)
  bt <- posterior::summarise_draws(posterior::subset_draws(res$draws, variable = "beta"))
  bt$variable <- colnames(X)
  write.csv(bt, file.path(out_dir, "beta.csv"), row.names = FALSE)

  fits <- list(Stan = list(PK = col("peak"), ET = col("eta_test"), sigma = dr[, "sigma"], Bz = Bz))
  brms_file <- file.path(env("ZH_OUT", "results/zone_harmonic"), "brms", "draws.rds")
  if (file.exists(brms_file)) {
    rb <- readRDS(brms_file)
    if (identical(rb$sy_loc_year, sy$loc_year)) fits$brms <- rb[c("PK", "ET", "sigma", "Bz")]
    else message("brms draws were fitted to different data (subset?) - not compared")
  }

  ho <- bind_rows(lapply(names(fits), function(m) holdout_metrics(m, fits[[m]]$ET, fits[[m]]$sigma)))
  zs <- bind_rows(lapply(names(fits), function(m) zone_peak_stats(m, fits[[m]]$PK)))
  zc <- lapply(fits, function(f) zone_curves(f$Bz))
  zcp <- bind_rows(lapply(names(zc), function(m) data.frame(
    model = m, peak_Southern = qs(zc[[m]]$Southern$peak), peak_Transitional = qs(zc[[m]]$Transitional$peak),
    peak_Northern = qs(zc[[m]]$Northern$peak),
    delay_Transitional = qs(zc[[m]]$Transitional$peak - zc[[m]]$Southern$peak),
    delay_Northern = qs(zc[[m]]$Northern$peak - zc[[m]]$Southern$peak))))
  message("\nHeld-out ", yr_test, ":"); print(ho, digits = 4)
  message("\nSite-year peak weeks (linear predictor incl. covariates) and delays behind Southern:"); print(zs)
  message("\nZone harmonic curve peaks (covariates at 0):"); print(zcp)
  write.csv(ho, file.path(out_dir, "holdout.csv"), row.names = FALSE)
  write.csv(zs, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)
  write.csv(zcp, file.path(out_dir, "zone_curve_peaks.csv"), row.names = FALSE)

  plot_peak_hist(lapply(fits, `[[`, "PK"), file.path(out_dir, "peak_hist.png"))
  plot_zone_curves(zc$Stan, file.path(out_dir, "zone_curves.png"),
                   sprintf("Stan: zone seasonal curves, %d harmonic(s) (95%% CrI)", H))
  plot_site_map(colMeans(col("w")), file.path(out_dir, "map_site_effect.png"),
                "Site effect w (posterior mean; centred within zone)")
  if (svc) {
    A <- col("alpha"); Bs <- col("beta_s")
    # first-harmonic shift as a peak-timing change in weeks: the phase shift
    # of b1 cos + c1 sin when (alpha, beta) is added, using the zone's b1, c1
    z <- as.character(site_xy$zone)
    b1 <- Bz[, paste0(z, ":c1"), drop = FALSE]; c1 <- Bz[, paste0(z, ":s1"), drop = FALSE]
    shift <- (atan2(c1 + Bs, b1 + A) - atan2(c1, b1)) * 52 / (2 * pi)
    shift <- (shift + 26) %% 52 - 26
    plot_site_map(colMeans(shift), file.path(out_dir, "map_peak_shift.png"),
                  "First-harmonic peak shift (weeks; + = later than the zone curve)")
    message(sprintf("Site peak shifts from the svc fields: sd across sites %.2f weeks, range %.1f to %.1f",
                    sd(colMeans(shift)), min(colMeans(shift)), max(colMeans(shift))))
  }
  message("Summary written to ", out_dir)
}
