# =============================================================================
# step1_stan/spacetime_simple.R
# Stan fit of the zone-free space-time model from standard parts
# (spacetime_simple.stan; written out in spacetime_simple_common.R):
#
#   log(CEW+1) = intercept + shared harmonic seasonal curve + trap + covariates
#              + (1 | year)                                year random intercept
#              + (1 | site), spatial + independent          site random intercept
#              + (0 + cos1 + sin1 | site), spatial          site random slopes
#              + AR(1) residual within site-year
#
# The brms version is spacetime_simple_brms.R.
#
# Modes:
#   Rscript step1_stan/spacetime_simple.R dryrun    # one short chain: time per gradient,
#                                                   # projected run time
#   Rscript step1_stan/spacetime_simple.R fit       # full fit -> ST_OUT/stan/draws.rds
#   Rscript step1_stan/spacetime_simple.R summary   # tables and figures (with brms if fitted)
#
# Environment (plus the ST_* settings in spacetime_simple_common.R):
#   STAN_CHAINS, STAN_ITER, STAN_WARMUP, STAN_THIN   4, 2000, 1000, 2 (iter includes warmup)
#   STAN_ADAPT_DELTA, STAN_MAX_TREEDEPTH             0.9, 10
#   ST_CENTER_SITE   1 = centered independent site intercept v (0 = non-centered, default)
#   ST_OUT           output folder                                results/spacetime_simple
#   DRY_WARMUP, DRY_SAMPLING                         dry-run chain length (150 + 50)
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "summary"))

source("step1_stan/spacetime_simple_common.R")

out_root <- env("ST_OUT", "results/spacetime_simple")
out_dir  <- file.path(out_root, "stan")
csv_dir  <- file.path(out_dir, "stan_csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
smp <- list(chains = env_num("STAN_CHAINS", 4), iter = env_num("STAN_ITER", 2000),
            warmup = env_num("STAN_WARMUP", 1000), thin = env_num("STAN_THIN", 2),
            adapt_delta = env_num("STAN_ADAPT_DELTA", 0.9), max_td = env_num("STAN_MAX_TREEDEPTH", 10))

# =============================================================================
# 1. Stan data -------------------------------------------------------------------
#    Priors: intercept N(0, 10^2), harmonics N(0, 2.5^2), trap / covariates
#    N(0, 1); random-effect sds exponential(3) (P(sd > 1) = 0.05); range PC
#    prior (P(range < ST_RANGE0) = 0.1); phi uniform(-1, 1); residual sd
#    half-Student-t(3, 0, 2.5)
# =============================================================================
stan_data <- list(
  N = length(fit_idx), P = ncol(X), T = max(d$year_id), S = nrow(site_xy),
  y = d$log_CEW[fit_idx], X = X[fit_idx, ], week = d$woy[fit_idx], year = d$year_id[fit_idx],
  site = d$site_id[fit_idx],
  xy = lapply(seq_len(nrow(site_xy)), function(i) c(site_xy$x_km[i], site_xy$y_km[i])),
  slopes = slopes, ar = use_ar, centered_site = as.integer(env_num("ST_CENTER_SITE", 0)), jitter = 1e-6,
  N_first = length(idx_first), idx_first = idx_first,
  N_next = length(idx_next), idx_next = idx_next, idx_prev = idx_prev, gap = gap,
  beta_sd = beta_sd, lambda_sd = -log(0.05), lambda_range = lambda_range, sigma_scale = 2.5,
  N_all = nrow(d), X_all = X, week_all = d$woy, year_all = d$year_id, site_all = d$site_id,
  win_lo = win[1], win_hi = win[2],
  N_test = length(test_idx), test_idx = test_idx, y_test = d$log_CEW[test_idx])

inits <- function() list(
  beta = c(1.5, rep(0, ncol(X) - 1)), u_z = rep(0, stan_data$T), sigma_year = 0.3,
  z_g = rep(0, stan_data$S), sigma_g = 0.5, v_par = rep(0, stan_data$S), sigma_v = 0.4,
  z_a = rep(0, slopes * stan_data$S), z_b = rep(0, slopes * stan_data$S),
  sigma_a = as.array(rep(0.3, slopes)), sigma_b = as.array(rep(0.3, slopes)),
  range_km = runif(1, 100, 300), phi = as.array(rep(0.5, use_ar)), sigma = 1.2)

hyper <- c("range_km", "sigma_g", "sigma_v", if (slopes) c("sigma_a[1]", "sigma_b[1]"), "sigma_year",
           if (use_ar) "phi[1]", "sigma", "sigma_innov")

# =============================================================================
# 2. Dry run / fit -------------------------------------------------------------------
# =============================================================================
if (mode %in% c("dryrun", "fit")) {
  library(cmdstanr)
  t0  <- Sys.time()
  mod <- cmdstan_model("step1_stan/spacetime_simple.stan")
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
  message("\n--- Dry run (1 chain, ", w, " warmup + ", s, " sampling iterations) ---")
  message(sprintf("Compile %.1f min. Warmup %.1f s/iteration, sampling %.1f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps, 2^smp$max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * mean(dg$treedepth__[!in_w] >= smp$max_td), sum(dg$divergent__[!in_w])))
  print(fit$summary(hyper), digits = 3)
  message(sprintf("Projected full run (%d iterations): %.1f h per chain at this step count; %.1f h if every iteration hits max depth",
                  smp$iter, smp$iter * steps * sec_grad / 3600, smp$iter * (2^smp$max_td - 1) * sec_grad / 3600))
}

if (mode == "fit") {
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d)", smp$chains, smp$iter, smp$warmup, smp$thin))
  t0  <- Sys.time()
  fit <- mod$sample(data = stan_data, init = inits, chains = smp$chains, parallel_chains = smp$chains,
                    iter_warmup = smp$warmup, iter_sampling = smp$iter - smp$warmup, thin = smp$thin,
                    adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                    seed = 42, refresh = 50, output_dir = csv_dir)
  message(sprintf("Sampling: %.2f h", mins(t0) / 60))
  saveRDS(list(draws = fit$draws(c(hyper, "beta", "u", "g", "v", "a", "b", "peak", "eta_test"),
                                 format = "draws_matrix"),
               diag = fit$diagnostic_summary(), time = fit$time(), smp = smp,
               setup = list(H = H, slopes = slopes, ar = use_ar, range0 = range0, win = win,
                            xcols = colnames(X), n_sub = n_sub)),
          file.path(out_dir, "draws.rds"))
  print(fit$summary(hyper), digits = 3)
  message("Saved ", file.path(out_dir, "draws.rds"))
}

# =============================================================================
# 3. Summary -----------------------------------------------------------------------
#    hyper.csv, beta.csv, holdout.csv, zone_stats.csv, peak_hist.png,
#    season_curve.png, map_site_intercept.png, map_site_peak.png. Includes
#    the brms fit when ST_OUT/brms/draws.rds exists.
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
  brms_file <- file.path(out_root, "brms", "draws.rds")
  if (file.exists(brms_file)) {
    rb <- readRDS(brms_file)
    if (identical(rb$sy_loc_year, sy$loc_year)) fits$brms <- rb[c("PK", "ET", "sigma", "Bz")]
    else message("brms draws were fitted to different data (subset?) - not compared")
  }

  ho <- bind_rows(lapply(names(fits), function(m) holdout_metrics(m, fits[[m]]$ET, fits[[m]]$sigma)))
  zs <- bind_rows(lapply(names(fits), function(m) zone_peak_stats(m, fits[[m]]$PK)))
  message("\nHeld-out ", yr_test, " (marginal predictive, no AR conditioning):"); print(ho, digits = 4)
  message("\nPeak weeks by zone (zones assigned after fitting) and delays behind Southern:"); print(zs)
  write.csv(ho, file.path(out_dir, "holdout.csv"), row.names = FALSE)
  write.csv(zs, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)

  plot_peak_hist(lapply(fits, `[[`, "PK"), file.path(out_dir, "peak_hist.png"))
  plot_curves(lapply(fits, function(f) shared_curve(f$Bz)), file.path(out_dir, "season_curve.png"))
  plot_site_map(colMeans(col("g")) + colMeans(col("v")), file.path(out_dir, "map_site_intercept.png"),
                "Site random intercept g + v (posterior mean)")
  sp <- colMeans(site_peaks(fits$Stan$PK))
  plot_site_map(sp, file.path(out_dir, "map_site_peak.png"),
                "Mean peak week of each site (posterior mean over its site-years)", mid = median(sp))
  message("Summary written to ", out_dir)
}
