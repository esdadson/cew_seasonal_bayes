# =============================================================================
# step1/step1_hpc.R
# Step 1 (step1/step1_models.R) as Slurm jobs on Hazel. Step 2 and the
# figures only read what these jobs write to results/step1/.
#
# Always run inside a Slurm job from the project directory (the script
# refuses otherwise on Linux):
#   STEP1_MODEL=fourier sbatch hpc/step1_dryrun.sbatch   # dry run -> suggested resources
#   ./hpc/submit_step1.sh all       # fit the three models, then compare + predict
# Modes:
#   Rscript step1/step1_hpc.R dryrun    # compile + one short chain of STEP1_MODEL
#   Rscript step1/step1_hpc.R fit       # fit STEP1_MODEL -> fit_<model>.rds, loo_<model>.rds
#   Rscript step1/step1_hpc.R compare   # loo_compare of the fitted models; predict
#                                       # the selected one over weeks 1-52 + peaks
#
# Environment:
#   STEP1_MODEL     fixed | spatial | fourier | simple (dryrun, fit; default fourier)
#   STEP1_CHAINS, STEP1_ITER, STEP1_WARMUP, STEP1_THIN   (defaults: fit_step1(),
#                   4 x 15000 iterations, 7500 warmup, thin 30 -> 1000 draws)
#   STEP1_THREADS   threads per chain (the job needs STEP1_CHAINS x STEP1_THREADS cores)
#   STEP1_GP_K      GP basis functions per dimension (step1_models.R)
#   STEP1_FROM_CSV  1 = fit mode rebuilds fit_<model>.rds from the run's saved
#                   CmdStan CSVs (results/step1/stan_csv/<model>/), no sampling
#   STEP1_SELECT    model to predict in compare mode (default: the best ELPD)
#   STEP1_PRED_DRAWS  posterior draws of the 52-week predictions kept for step 2 (100)
#   DRY_WARMUP, DRY_SAMPLING   dry-run chain length (default 150 + 50)
#   CEW_DATA, CEW_OUT          data root and covariates/ (step1_data.R)
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with `sbatch hpc/step1_dryrun.sbatch` or ",
       "`./hpc/submit_step1.sh` - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "compare"))

t_start <- Sys.time()
source("step1/step1_models.R")
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("Data: %d observed rows (%s), %d sites, %d site-years, %d prediction rows; %.1f min",
                nrow(d_obs), paste(range(d_obs$year), collapse = "-"), n_distinct(d_obs$location2),
                n_distinct(d_obs$loc_year), nrow(d_all), mins(t_start)))
message("CmdStan ", cmdstanr::cmdstan_version(), " at ", cmdstanr::cmdstan_path())

env_int <- function(name, default) as.integer(Sys.getenv(name, default))
def  <- formals(fit_step1)
full <- list(chains = env_int("STEP1_CHAINS", def$chains), iter = env_int("STEP1_ITER", def$iter),
             warmup = env_int("STEP1_WARMUP", def$warmup), thin = env_int("STEP1_THIN", def$thin),
             threads = env_int("STEP1_THREADS", def$threads))
max_td <- step1_control$max_treedepth
model  <- Sys.getenv("STEP1_MODEL", "fourier")
if (mode != "compare") stopifnot(model %in% names(step1_terms))

peak_mem_gb <- function() {          # Linux only: peak resident memory of this R process
  s <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
  kb <- as.numeric(gsub("\\D", "", grep("^VmHWM", s, value = TRUE)))
  if (length(kb)) kb / 1024^2 else NA_real_
}

if (mode == "dryrun") {
  # Same Stan program and data brm() builds, sampled directly so compile,
  # warmup and sampling are timed separately. One chain with the full run's
  # threads per chain: chains run in parallel, so this is one chain's speed.
  t0  <- Sys.time()
  thr <- if (full$threads > 1) threading(full$threads) else NULL
  dd  <- droplevels(d_obs)
  code  <- make_stancode(step1_formula(model), data = dd, family = gaussian(), prior = step1_priors(model), threads = thr)
  sdata <- make_standata(step1_formula(model), data = dd, family = gaussian(), prior = step1_priors(model), threads = thr)
  mod   <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code),
                                   cpp_options = if (full$threads > 1) list(stan_threads = TRUE) else list())
  t_compile <- mins(t0)

  w <- env_int("DRY_WARMUP", 150); s <- env_int("DRY_SAMPLING", 50)
  fit <- mod$sample(data = sdata, chains = 1, iter_warmup = w, iter_sampling = s,
                    threads_per_chain = if (full$threads > 1) full$threads else NULL,
                    seed = 42, init = 0, adapt_delta = step1_control$adapt_delta, max_treedepth = max_td,
                    save_warmup = TRUE, refresh = 10)
  tm <- fit$time()$chains
  dg <- posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))
  lf <- dg$n_leapfrog__; in_w <- seq_along(lf) <= w
  sec_grad <- tm$total / sum(lf)

  # Projections for one full-run chain, in hours (thinning saves draws, not time)
  steps_obs <- mean(lf[!in_w])
  h_obs <- full$iter * steps_obs * sec_grad / 3600
  h_max <- full$iter * (2^max_td - 1) * sec_grad / 3600
  h_post <- 1                                 # + data, compile, loo
  # padded: adaptation in the dry run is short (the full spatial run took
  # 1.9 x h_obs); if most iterations still hit max tree depth, size for the
  # worst case
  at_max <- mean(dg$treedepth__[!in_w] >= max_td)
  h_req <- (if (at_max > 0.5) 1.15 * h_max else 2.5 * h_obs) + h_post
  qos   <- if (h_req <= 4 * 24) "normal" else if (h_req <= 10 * 24) "long" else NA
  fmt_t <- function(h) { m <- ceiling(min(h, 240) * 60)
    sprintf("%d-%02d:%02d:00", m %/% 1440, (m %% 1440) %/% 60, m %% 60) }
  # Memory: the saved draws of every parameter (incl. brms' transformed r_*),
  # read into R ~3 times over while CmdStan's CSVs are parsed, plus loo's
  # draws x observations log-likelihood
  n_draws <- full$chains * (full$iter - full$warmup) %/% full$thin
  n_var <- ncol(posterior::as_draws_df(fit$draws()))
  gb_draws <- n_var * n_draws * 8 / 1024^3
  gb_loo <- 2 * n_draws * nrow(d_obs) * 8 / 1024^3
  mem <- peak_mem_gb()

  message("\n--- Dry run of '", model, "' (1 chain, ", w, " warmup + ", s, " sampling iterations, ",
          full$threads, " thread(s)) ---")
  message(sprintf("Compile: %.1f min. Warmup %.2f s/iteration, sampling %.2f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps_obs, 2^max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * at_max, sum(dg$divergent__[!in_w])))
  message(sprintf("Parameters saved per draw: %d; %d draws -> %.2f GB of draws, %.1f GB for loo",
                  n_var, n_draws, gb_draws, gb_loo))
  message("Peak memory of R in the dry run: ", round(mem, 2), " GB")

  message(sprintf("\n--- Full run: %d chains x %d iterations (%d warmup, thin %d), %d thread(s) per chain ---",
                  full$chains, full$iter, full$warmup, full$thin, full$threads))
  message(sprintf("Per chain at the dry run's tree depth: %.1f h; if every iteration hits max depth: %.0f h",
                  h_obs, h_max))
  message("\n--- Suggested settings for hpc/submit_step1.sh (", model, ") ---")
  if (is.na(qos)) message("Exceeds Hazel's 10-day limit (qos=long): lower STEP1_ITER or raise STEP1_THREADS.")
  message("STEP1_TIME=", fmt_t(h_req), "   STEP1_QOS=", ifelse(is.na(qos), "long", qos),
          "   STEP1_MEM=", sprintf("%dG", max(16, ceiling(3 * gb_draws + gb_loo + 2 * ifelse(is.na(mem), 4, mem)))),
          "   # cores = STEP1_CHAINS x STEP1_THREADS")
  message("Check seff <jobid> before the full run.")

} else if (mode == "fit") {
  pf <- step1_files(model)
  message(sprintf("Model '%s': fitting %d chains x %d iterations (%d warmup, thin %d), %d thread(s) per chain, on %s cores",
                  model, full$chains, full$iter, full$warmup, full$thin, full$threads,
                  Sys.getenv("SLURM_CPUS_PER_TASK", "?")))
  print(step1_formula(model))
  # CmdStan's CSVs are kept in results/step1/stan_csv/<model>/: if R fails
  # after sampling, STEP1_FROM_CSV=1 rebuilds the fit from them instead of refitting
  dir.create(pf$csv, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  if (Sys.getenv("STEP1_FROM_CSV") == "1") {
    # newest run only: one run's chains share the id after the last "-"
    csvs <- list.files(pf$csv, pattern = "[.]csv$", full.names = TRUE)
    if (!length(csvs)) stop("No CmdStan CSVs in ", pf$csv)
    run  <- sub(".*-", "", tools::file_path_sans_ext(basename(csvs)))
    csvs <- csvs[run == run[which.max(file.mtime(csvs))]]
    message("Rebuilding the fit from ", length(csvs), " CSVs:\n  ", paste(basename(csvs), collapse = "\n  "))
    if (length(csvs) != full$chains) stop("Expected ", full$chains, " chains, found ", length(csvs))
    fit <- fit_step1(model, chains = full$chains, iter = full$iter, warmup = full$warmup, thin = full$thin,
                     threads = full$threads, file = NULL, empty = TRUE)
    fit$fit <- brms:::read_csv_as_stanfit(csvs, exclude = brms:::exclude_pars(fit))
    fit <- rename_pars(fit)
    n_exp <- full$chains * (full$iter - full$warmup) %/% full$thin
    if (ndraws(fit) != n_exp) stop("The CSVs hold ", ndraws(fit), " draws, expected ", n_exp, " (a chain did not finish)")
    saveRDS(fit, paste0(pf$fit, ".rds"))
  } else {
    fit <- fit_step1(model, chains = full$chains, iter = full$iter, warmup = full$warmup, thin = full$thin,
                     threads = full$threads, file_refit = "always", output_dir = pf$csv, refresh = 500)
  }
  message(sprintf("Sampling + saving: %.1f h. Saved %s.rds", mins(t0) / 60, pf$fit))
  print(summary(fit))
  np <- nuts_params(fit)
  message(sprintf("Divergent transitions: %d; at max tree depth: %.1f%%; max Rhat: %.3f; min bulk ESS: %.0f; min tail ESS: %.0f",
                  sum(np$Value[np$Parameter == "divergent__"]),
                  100 * mean(np$Value[np$Parameter == "treedepth__"] >= max_td),
                  max(rhat(fit), na.rm = TRUE),
                  min(summarise_draws(as_draws_df(fit), "ess_bulk")$ess_bulk, na.rm = TRUE),
                  min(summarise_draws(as_draws_df(fit), "ess_tail")$ess_tail, na.rm = TRUE)))

  # PSIS-LOO (leave one week out) for loo_compare in compare mode
  t0 <- Sys.time()
  lo <- loo(fit)
  saveRDS(lo, pf$loo)
  print(lo)
  message(sprintf("PSIS-LOO: %.1f min. Saved %s", mins(t0), pf$loo))
  message("Peak memory of R: ", round(peak_mem_gb(), 2), " GB")

} else if (mode == "compare") {
  # ---- Model comparison -------------------------------------------------------
  # elpd_diff = expected log predictive density relative to the best model;
  # a difference within ~2 se_diff is not a clear preference - prefer the
  # simpler model then (STEP1_SELECT overrides the choice).
  loo_have <- Filter(file.exists, sapply(names(step1_terms), function(m) step1_files(m)$loo))
  if (!length(loo_have)) stop("No loo_<model>.rds in ", step1_dir, ": run the fit jobs first")
  los <- lapply(loo_have, readRDS)
  for (m in names(los)) { message("\n", m, ": Pareto k diagnostics"); print(loo::pareto_k_table(los[[m]])) }
  if (length(los) >= 2) {
    lc <- loo::loo_compare(los)
    print(lc, simplify = FALSE)
    lc_tab <- data.frame(model = rownames(lc), as.data.frame(unclass(lc)), row.names = NULL)
  } else {
    lc_tab <- data.frame(model = names(los), t(los[[1]]$estimates[, "Estimate"]), row.names = NULL)
  }
  write.csv(lc_tab, file.path(step1_dir, "loo_compare.csv"), row.names = FALSE)
  sel <- Sys.getenv("STEP1_SELECT", lc_tab$model[1])
  stopifnot(sel %in% names(los))
  message(sprintf("\nSelected model: '%s'%s", sel,
                  if (nzchar(Sys.getenv("STEP1_SELECT"))) " (STEP1_SELECT)" else " (best ELPD)"))

  # ---- Predictions over weeks 1-52 and peaks -----------------------------------
  t0 <- Sys.time(); set.seed(42)
  fit <- readRDS(paste0(step1_files(sel)$fit, ".rds"))
  pr  <- predict_step1(fit, d_all, n_keep = env_int("STEP1_PRED_DRAWS", 100))
  saveRDS(list(model = sel, pred = pr$pred, draws = pr$draws), file.path(step1_dir, "predictions.rds"))
  saveRDS(c(pr$peaks, model = sel), file.path(step1_dir, "peak_draws.rds"))
  message(sprintf("Predictions (%d rows x %d draws; %d kept) and peaks (%d site-years): %.1f min",
                  nrow(pr$pred), ndraws(fit), ncol(pr$draws), ncol(pr$peaks$draws), mins(t0)))

  ps <- peak_summaries(pr$peaks)
  print(ps$site_year %>% group_by(zone = zone_lab[as.character(zone_30y)]) %>%
          summarise(n = n(), median_ci_width = median(ci_width), share_in_sampled = mean(in_sampled > 0.5)))
  print(ps$zone_delay %>% filter(year == "all_years"))
  write.csv(ps$site_year, file.path(step1_dir, "peak_site_year.csv"), row.names = FALSE)
  write.csv(ps$zone_year, file.path(step1_dir, "peak_zone_year.csv"), row.names = FALSE)
  write.csv(ps$zone_delay, file.path(step1_dir, "peak_zone_delay.csv"), row.names = FALSE)
  write.csv(ps$site, file.path(step1_dir, "peak_site.csv"), row.names = FALSE)
  write.csv(zone_curve_summary(pr$zone_curve), file.path(step1_dir, "pred_zone_week.csv"), row.names = FALSE)

  # Modelling data for step 2 (so it doesn't rebuild it from the covariates)
  saveRDS(list(d_all = d_all, sites = sites, anom_vars = anom_vars, yrs = yrs, peak_win = peak_win),
          file.path(step1_dir, "model_data.rds"))
  message("Saved to ", step1_dir, ": loo_compare.csv, predictions.rds, peak_draws.rds, peak_*.csv, ",
          "pred_zone_week.csv, model_data.rds")
  message("Peak memory of R: ", round(peak_mem_gb(), 2), " GB")
}
