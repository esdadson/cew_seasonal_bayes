# =============================================================================
# step1_stan/step1_stan.R
# Fits and compares three nested Stan models of weekly log(CEW + 1)
# (definitions in stan_common.R and in each .stan file):
#
#   fixed           fixed.stan           intercept, trap type, covariates
#   fourier         fourier.stan         + shared harmonic seasonal curve (<= 2 harmonics)
#                                        + year random intercept
#   spatiotemporal  spatiotemporal.stan  + site intercept (spatial Matern + independent)
#                                        + spatial site slopes on cos1 / sin1
#                                        + AR(1) residual within site-year
#
# 2020 is held out and predicted out of sample.
#
# Modes:
#   Rscript step1_stan/step1_stan.R dryrun    # one short chain per model: seconds per
#                                             # gradient, steps per iteration, projected time
#                                             # of the full run -> results/step1_stan/timing.csv
#   Rscript step1_stan/step1_stan.R fit       # full fits -> results/step1_stan/<model>/
#   Rscript step1_stan/step1_stan.R summary   # comparison of the fitted models
#
# Environment (plus STAN_HARM, STAN_RANGE0, STAN_SUBSET_SITES, STEP1_HOLDOUT,
# STEP1_DATA_CACHE in stan_common.R):
#   STAN_MODEL                 fixed, fourier, spatiotemporal or all (comma list ok)   all
#   STAN_CHAINS                4
#   STAN_ITER                  iterations per chain, warmup included                  15000
#   STAN_WARMUP                7500
#   STAN_THIN                  30  -> 4 x (15000 - 7500) / 30 = 1000 posterior draws
#   STAN_ADAPT_DELTA, STAN_MAX_TREEDEPTH                                              0.9, 10
#   STAN_OUT                   output folder                                          results/step1_stan
#   DRY_WARMUP, DRY_SAMPLING   dry-run chain length                                   150, 50
# Thinning only shrinks the saved output; every iteration still has to run.
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "summary"))

source("step1_stan/stan_common.R")

models <- strsplit(env("STAN_MODEL", "all"), ",")[[1]]
if (identical(models, "all")) models <- MODELS
stopifnot(all(models %in% MODELS))
out_root <- env("STAN_OUT", "results/step1_stan")
smp <- list(chains = env_num("STAN_CHAINS", 4), iter = env_num("STAN_ITER", 15000),
            warmup = env_num("STAN_WARMUP", 7500), thin = env_num("STAN_THIN", 30),
            adapt_delta = env_num("STAN_ADAPT_DELTA", 0.9), max_td = env_num("STAN_MAX_TREEDEPTH", 10))
n_draws <- smp$chains * floor((smp$iter - smp$warmup) / smp$thin)
message(sprintf("Sampler: %d chains x %d iterations (%d warmup), thin %d -> %d posterior draws",
                smp$chains, smp$iter, smp$warmup, smp$thin, n_draws))

model_dir <- function(m) { p <- file.path(out_root, m); dir.create(file.path(p, "stan_csv"), recursive = TRUE, showWarnings = FALSE); p }
compile   <- function(m) cmdstanr::cmdstan_model(file.path("step1_stan", paste0(m, ".stan")))
saved_vars <- function(m) c(hyper_for(m), "beta", if (m != "fixed") "u",
                            if (m == "spatiotemporal") c("g", "v", "a", "b"), "peak", "eta_test")

# =============================================================================
# 1. Dry run: time per gradient and projected run time ----------------------------------
# =============================================================================
if (mode == "dryrun") {
  w <- env_num("DRY_WARMUP", 150); s <- env_num("DRY_SAMPLING", 50)
  timing <- list()
  for (m in models) {
    mod <- compile(m)
    fit <- mod$sample(data = stan_data_for(m), init = inits_for(m), chains = 1,
                      iter_warmup = w, iter_sampling = s, seed = 42,
                      adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                      save_warmup = TRUE, refresh = 25, output_dir = file.path(model_dir(m), "stan_csv"))
    tm <- fit$time()$chains
    dg <- posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))
    lf <- dg$n_leapfrog__; in_w <- seq_along(lf) <= w
    sec_grad <- tm$total / sum(lf)            # includes generated quantities: an upper bound
    steps    <- mean(lf[!in_w])
    timing[[m]] <- data.frame(
      model = m, sec_per_gradient = sec_grad, leapfrog_steps = steps,
      pct_max_treedepth = 100 * mean(dg$treedepth__[!in_w] >= smp$max_td),
      divergent = sum(dg$divergent__[!in_w]),
      warmup_s_per_iter = tm$warmup / w, sampling_s_per_iter = tm$sampling / s,
      iter = smp$iter,
      hours_per_chain = smp$iter * steps * sec_grad / 3600,
      hours_if_max_depth = smp$iter * (2^smp$max_td - 1) * sec_grad / 3600)
    message(sprintf("\n--- %s: dry run (1 chain, %d warmup + %d sampling) ---", m, w, s))
    message(sprintf("Seconds per gradient %.4f; leapfrog steps %.0f; max tree depth %.0f%%; divergent %d",
                    sec_grad, steps, timing[[m]]$pct_max_treedepth, timing[[m]]$divergent))
    print(fit$summary(hyper_for(m)), digits = 3)
    message(sprintf("Projected %d iterations: %.1f h per chain (chains run in parallel); %.1f h if every iteration hits max depth",
                    smp$iter, timing[[m]]$hours_per_chain, timing[[m]]$hours_if_max_depth))
  }
  timing <- bind_rows(timing)
  print(timing, digits = 3)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  write.csv(timing, file.path(out_root, "timing.csv"), row.names = FALSE)
  message(sprintf("All models in one job, %d chains in parallel: about %.1f h (worst case %.1f h). Request ~1.5x that.",
                  smp$chains, sum(timing$hours_per_chain), sum(timing$hours_if_max_depth)))
}

# =============================================================================
# 2. Fit ---------------------------------------------------------------------------------
#    Saves draws.rds (hyperparameters, beta, effects, peak weeks, held-out
#    means), loo.rds (PSIS-LOO) and diagnostics. The Stan CSVs stay in
#    <model>/stan_csv/.
# =============================================================================
if (mode == "fit") {
  for (m in models) {
    mod <- compile(m)
    message(sprintf("\n=== %s ===", m))
    t0  <- Sys.time()
    fit <- mod$sample(data = stan_data_for(m), init = inits_for(m),
                      chains = smp$chains, parallel_chains = smp$chains,
                      iter_warmup = smp$warmup, iter_sampling = smp$iter - smp$warmup, thin = smp$thin,
                      adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                      seed = 42, refresh = 500, output_dir = file.path(model_dir(m), "stan_csv"))
    hours <- mins(t0) / 60
    dg  <- fit$diagnostic_summary()
    sm  <- fit$summary(c(hyper_for(m), "beta"))
    diag <- data.frame(model = m, hours = hours, draws = n_draws,
                       divergent = sum(dg$num_divergent), max_treedepth_hits = sum(dg$num_max_treedepth),
                       min_ebfmi = min(dg$ebfmi), max_rhat = max(sm$rhat, na.rm = TRUE),
                       min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE), min_ess_tail = min(sm$ess_tail, na.rm = TRUE))
    print(sm, n = 40, digits = 3); print(diag, digits = 3)
    saveRDS(list(draws = fit$draws(saved_vars(m), format = "draws_matrix"), diag = diag, smp = smp,
                 setup = list(H = H, range0 = range0, win = win, xcols = colnames(X_list[[m]]), n_sub = n_sub,
                              holdout = yr_test, sy_loc_year = sy$loc_year)),
            file.path(model_dir(m), "draws.rds"))
    saveRDS(fit$loo(cores = smp$chains), file.path(model_dir(m), "loo.rds"))
    message(sprintf("%s: %.2f h. Saved %s", m, hours, model_dir(m)))
  }
}

# =============================================================================
# 3. Summary: compare whichever models have been fitted --------------------------------------
#    comparison.csv  PSIS-LOO (elpd_loo, p_loo, Pareto k > 0.7, elpd_diff vs
#                    best), held-out 2020 (elpd, R2, RMSE, MAE, 90% coverage),
#                    sampler diagnostics
#    zone_stats.csv  peak weeks by zone, delays behind Southern, eta^2
#    <model>/hyper.csv, <model>/beta.csv
#    peak_hist.png, season_curve.png, map_site_intercept.png, map_site_peak.png
# =============================================================================
if (mode == "summary") {
  done <- models[file.exists(file.path(out_root, models, "draws.rds"))]
  if (!length(done)) stop("No fitted models in ", out_root)
  res <- lapply(setNames(done, done), function(m) readRDS(file.path(out_root, m, "draws.rds")))
  for (m in done) if (!identical(res[[m]]$setup$sy_loc_year, sy$loc_year))
    stop(m, " was fitted to different data (subset or holdout?) than the current settings")
  dm  <- lapply(res, function(r) matrix(r$draws, nrow(r$draws), dimnames = dimnames(r$draws)))
  col <- function(m, p) dm[[m]][, grep(paste0("^", p, "\\["), colnames(dm[[m]])), drop = FALSE]

  for (m in done) {
    hyp <- posterior::summarise_draws(posterior::subset_draws(res[[m]]$draws, variable = hyper_for(m)))
    bt  <- posterior::summarise_draws(posterior::subset_draws(res[[m]]$draws, variable = "beta"))
    bt$variable <- res[[m]]$setup$xcols
    message("\n", m); print(hyp, digits = 3); print(bt, digits = 3)
    write.csv(hyp, file.path(out_root, m, "hyper.csv"), row.names = FALSE)
    write.csv(bt, file.path(out_root, m, "beta.csv"), row.names = FALSE)
  }

  loos <- lapply(setNames(done, done), function(m) readRDS(file.path(out_root, m, "loo.rds")))
  lc   <- loo::loo_compare(loos)
  cmp <- bind_rows(lapply(done, function(m) {
    est <- loos[[m]]$estimates
    cbind(data.frame(model = m,
                     elpd_loo = est["elpd_loo", "Estimate"], se_elpd_loo = est["elpd_loo", "SE"],
                     p_loo = est["p_loo", "Estimate"],
                     n_pareto_k_gt_0.7 = sum(loo::pareto_k_values(loos[[m]]) > 0.7),
                     elpd_diff = lc[m, "elpd_diff"], se_diff = lc[m, "se_diff"]),
          holdout_metrics(col(m, "eta_test"), dm[[m]][, "sigma"]),
          res[[m]]$diag[, c("hours", "divergent", "max_treedepth_hits", "max_rhat", "min_ess_bulk")])
  }))
  zs <- bind_rows(lapply(done, function(m) cbind(model = m, zone_peak_stats(col(m, "peak")))))
  message("\nModel comparison (elpd_diff: against the best LOO model; held-out ", yr_test,
          " scored with the marginal predictive):"); print(cmp, digits = 4)
  message("\nPeak weeks by zone (zones assigned after fitting), delays behind Southern:"); print(zs)
  write.csv(cmp, file.path(out_root, "comparison.csv"), row.names = FALSE)
  write.csv(zs, file.path(out_root, "zone_stats.csv"), row.names = FALSE)

  plot_peak_hist(lapply(setNames(done, done), function(m) col(m, "peak")), file.path(out_root, "peak_hist.png"))
  seas <- intersect(done, c("fourier", "spatiotemporal"))
  if (length(seas)) plot_curves(lapply(setNames(seas, seas), function(m) {
    Bz <- col(m, "beta"); colnames(Bz) <- res[[m]]$setup$xcols; shared_curve(Bz) }),
    file.path(out_root, "season_curve.png"))
  if ("spatiotemporal" %in% done) {
    m <- "spatiotemporal"
    plot_site_map(colMeans(col(m, "g")) + colMeans(col(m, "v")), file.path(out_root, "map_site_intercept.png"),
                  "Spatiotemporal: site intercept g + v (posterior mean)")
    sp <- colMeans(site_peaks(col(m, "peak")))
    plot_site_map(sp, file.path(out_root, "map_site_peak.png"),
                  "Spatiotemporal: mean peak week of each site (posterior mean over its site-years)", mid = median(sp))
  }
  message("Summary written to ", out_root)
}
