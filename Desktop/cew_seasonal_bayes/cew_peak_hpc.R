# =============================================================================
# cew_peak_hpc.R
# The step-1 peak model (cew_peak_model.R) as Slurm jobs on Hazel. The
# analysis (cew_peak_granger.R) only reads what these jobs write.
#
# Always run inside a Slurm job (the script refuses otherwise on Linux):
#   sbatch hpc/peak_dryrun.sbatch           # dry run -> suggested resources
#   ./hpc/submit_peak.sh fit                # full fit
#   ./hpc/submit_peak.sh loyo               # leave-one-year-out refits (array)
# Modes (for the candidate model PEAK_MODEL = site | zone_year, cew_peak_model.R):
#   Rscript cew_peak_hpc.R dryrun   # compile + one short chain on all years
#   Rscript cew_peak_hpc.R fit      # all years -> fit_peak_<model>.rds, loo_peak_<model>.rds,
#                                   #   peak_draws_<model>.rds, fill_<model>.rds
#   Rscript cew_peak_hpc.R loyo     # year yrs[SLURM_ARRAY_TASK_ID + 1] held out
#                                   # -> loyo/<model>/loyo_peak_<year>.rds
#
# Environment:
#   PEAK_MODEL      candidate model (default site)
#   PEAK_CHAINS, PEAK_ITER, PEAK_WARMUP   full-run settings (defaults: fit_peak())
#   PEAK_THREADS    threads per chain (within-chain parallelism; the job needs
#                   PEAK_CHAINS x PEAK_THREADS cores)
#   PEAK_GP_K       GP basis functions per dimension (cew_peak_model.R)
#   PEAK_SAVE_LOYO_FIT=1   also keep each LOYO fit (large; not needed)
#   PEAK_FILL_DRAWS posterior draws of the filled weekly series kept for step 2 (100)
#   DRY_WARMUP, DRY_SAMPLING   dry-run chain length (default 150 + 50)
#   CEW_DATA, CEW_OUT          as for cew_covariates_hpc.R
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with `sbatch hpc/peak_dryrun.sbatch` or ",
       "`./hpc/submit_peak.sh` - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "loyo"))

t_start <- Sys.time()
source("cew_peak_model.R")
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("Data: %d observed rows (%s), %d sites, %d site-years; %.1f min",
                nrow(d_obs), paste(range(d_obs$year), collapse = "-"),
                n_distinct(d_obs$location2), n_distinct(d_obs$loc_year), mins(t_start)))
message("CmdStan ", cmdstanr::cmdstan_version(), " at ", cmdstanr::cmdstan_path())

env_int <- function(name, default) as.integer(Sys.getenv(name, default))
def  <- formals(fit_peak)
full <- list(chains = env_int("PEAK_CHAINS", def$chains), iter = env_int("PEAK_ITER", def$iter),
             warmup = env_int("PEAK_WARMUP", def$warmup), threads = env_int("PEAK_THREADS", def$threads))
max_td <- peak_control$max_treedepth           # cew_peak_model.R

peak_mem_gb <- function() {          # Linux only: peak resident memory of this R process
  s <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
  kb <- as.numeric(gsub("\\D", "", grep("^VmHWM", s, value = TRUE)))
  if (length(kb)) kb / 1024^2 else NA_real_
}
diagnostics <- function(fit) {
  np <- nuts_params(fit)
  message(sprintf("Divergent transitions: %d; at max tree depth: %.1f%%; max Rhat: %.3f; min bulk ESS: %.0f",
                  sum(np$Value[np$Parameter == "divergent__"]),
                  100 * mean(np$Value[np$Parameter == "treedepth__"] >= max_td),
                  max(rhat(fit), na.rm = TRUE),
                  min(posterior::summarise_draws(as_draws_df(fit), "ess_bulk")$ess_bulk, na.rm = TRUE)))
}
fit_args <- function(data, file, refit = "always")
  list(data = data, file = file, chains = full$chains, iter = full$iter, warmup = full$warmup,
       threads = full$threads, file_refit = refit, output_dir = "stan_csv", refresh = 100)

if (mode == "dryrun") {
  # Same Stan program and data brm() builds, sampled directly so compile,
  # warmup and sampling are timed separately. One chain with the full run's
  # threads per chain: chains run in parallel, so this is one chain's speed.
  t0  <- Sys.time()
  thr <- if (full$threads > 1) threading(full$threads) else NULL
  code  <- make_stancode(f_peak, data = d_obs, family = gaussian(), prior = pri_peak, threads = thr)
  sdata <- make_standata(f_peak, data = d_obs, family = gaussian(), prior = pri_peak, threads = thr)
  mod   <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code),
                                   cpp_options = if (full$threads > 1) list(stan_threads = TRUE) else list())
  t_compile <- mins(t0)

  w <- env_int("DRY_WARMUP", 150); s <- env_int("DRY_SAMPLING", 50)
  fit <- mod$sample(data = sdata, chains = 1, iter_warmup = w, iter_sampling = s,
                    threads_per_chain = if (full$threads > 1) full$threads else NULL,
                    seed = 42, init = 0, adapt_delta = peak_control$adapt_delta, max_treedepth = max_td,
                    save_warmup = TRUE, refresh = 10)
  tm <- fit$time()$chains
  dg <- posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))
  lf <- dg$n_leapfrog__; in_w <- seq_along(lf) <= w
  sec_grad <- tm$total / sum(lf)

  # Projections for one full-run chain, in hours (the LOYO folds are the same size)
  steps_obs <- mean(lf[!in_w])
  h_obs <- full$iter * steps_obs * sec_grad / 3600
  h_max <- full$iter * (2^max_td - 1) * sec_grad / 3600
  h_post <- 1                                 # + data, compile, peak draws / LOYO predictions
  # padded: adaptation in the dry run is short; if most iterations still hit
  # max tree depth, size for the worst case
  at_max <- mean(dg$treedepth__[!in_w] >= max_td)
  h_req <- (if (at_max > 0.5) 1.15 * h_max else 1.5 * h_obs) + h_post
  qos   <- if (h_req <= 4 * 24) "normal" else if (h_req <= 10 * 24) "long" else NA
  fmt_t <- function(h) { m <- ceiling(min(h, 240) * 60)
    sprintf("%d-%02d:%02d:00", m %/% 1440, (m %% 1440) %/% 60, m %% 60) }
  # Memory: all draws of every parameter (incl. brms' transformed r_*), read
  # into R ~3 times over while CmdStan's CSVs are parsed
  n_var <- ncol(posterior::as_draws_df(fit$draws()))
  gb_draws <- n_var * full$chains * (full$iter - full$warmup) * 8 / 1024^3
  mem <- peak_mem_gb()

  message("\n--- Dry run (1 chain, ", w, " warmup + ", s, " sampling iterations, ",
          full$threads, " thread(s)) ---")
  message(sprintf("Compile: %.1f min. Warmup %.1f s/iteration, sampling %.1f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps_obs, 2^max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * mean(dg$treedepth__[!in_w] >= max_td), sum(dg$divergent__[!in_w])))
  message(sprintf("Parameters saved per draw: %d -> %.1f GB of draws for the full run", n_var, gb_draws))
  message("Peak memory of R in the dry run: ", round(mem, 2), " GB")

  message(sprintf("\n--- Full run: %d chains x %d iterations (%d warmup), %d thread(s) per chain ---",
                  full$chains, full$iter, full$warmup, full$threads))
  message(sprintf("Per chain at the dry run's tree depth: %.1f h; if every iteration hits max depth: %.0f h",
                  h_obs, h_max))
  message("\n--- Suggested settings for hpc/submit_peak.sh ---")
  if (is.na(qos)) message("Exceeds Hazel's 10-day limit (qos=long): lower PEAK_ITER or raise PEAK_THREADS.")
  message("PEAK_TIME=", fmt_t(h_req), "   PEAK_QOS=", ifelse(is.na(qos), "long", qos),
          "   PEAK_MEM=", sprintf("%dG", max(16, ceiling(3 * gb_draws + 2 * ifelse(is.na(mem), 4, mem)))),
          "   # cores = PEAK_CHAINS x PEAK_THREADS")
  message("Check seff <jobid> before the full run.")

} else if (mode == "fit") {
  pf <- peak_files()
  message(sprintf("Model '%s': fitting %d chains x %d iterations (%d warmup), %d thread(s) per chain, on %s cores",
                  peak_model, full$chains, full$iter, full$warmup, full$threads,
                  Sys.getenv("SLURM_CPUS_PER_TASK", "?")))
  # CmdStan's CSVs are kept in stan_csv/<model>/: if R fails after sampling,
  # read them with cmdstanr::as_cmdstan_fit() instead of refitting
  csv_dir <- file.path("stan_csv", peak_model); dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  a <- fit_args(d_obs, pf$fit); a$output_dir <- csv_dir
  fit <- do.call(fit_peak, a)
  message(sprintf("Sampling + saving: %.1f h. Saved %s.rds", mins(t0) / 60, pf$fit))
  print(summary(fit))
  diagnostics(fit)

  # PSIS-LOO for loo_compare (cew_peak_granger.R section 3); the pointwise
  # log-likelihood is draws x observations (~1.2 GB for 4000 x 38568)
  t0 <- Sys.time()
  lo <- loo(fit)
  saveRDS(lo, pf$loo)
  print(lo)
  message(sprintf("PSIS-LOO: %.1f min. Saved %s", mins(t0), pf$loo))

  t0 <- Sys.time()
  pk <- peak_draws(fit, peak_grid(d_all))
  saveRDS(pk, pf$draws)
  message(sprintf("Peak draws (%d draws x %d site-years): %.1f min. Saved %s",
                  nrow(pk$draws), ncol(pk$draws), mins(t0), pf$draws))

  # Filled weekly populations for step 2 (cew_peak_granger.R sections 4-6)
  t0 <- Sys.time(); set.seed(42)
  fl <- fill_draws(fit, d_all, n_keep = env_int("PEAK_FILL_DRAWS", 100))
  saveRDS(fl, pf$fill)
  message(sprintf("Filled weekly series (%d rows, %d draws kept): %.1f min. Saved %s",
                  nrow(fl$fill), ncol(fl$draws), mins(t0), pf$fill))
  message("Peak memory of R: ", round(peak_mem_gb(), 2), " GB")

} else if (mode == "loyo") {
  task <- env_int("SLURM_ARRAY_TASK_ID", NA)
  if (is.na(task)) stop("loyo mode runs as a Slurm array: ./hpc/submit_peak.sh loyo")
  y_out <- yrs[task + 1]
  ldir  <- peak_files()$loyo_dir
  message(sprintf("Model '%s', leave-one-year-out: holding out %d (%d of %d observed rows)",
                  peak_model, y_out, sum(d_obs$year == y_out), nrow(d_obs)))
  dir.create(ldir, recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path("stan_csv", peak_model, paste0("loyo_", y_out))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  keep_fit <- Sys.getenv("PEAK_SAVE_LOYO_FIT") == "1"

  t0 <- Sys.time()
  a <- fit_args(d_obs %>% filter(year != y_out),
                if (keep_fit) file.path(ldir, paste0("fit_peak_loyo_", y_out)) else NULL)
  a$output_dir <- csv_dir
  fit <- do.call(fit_peak, a)
  message(sprintf("Sampling: %.1f h", mins(t0) / 60))
  diagnostics(fit)

  t0 <- Sys.time()
  res <- loyo_eval(fit, d_obs %>% filter(year == y_out), peak_grid(d_all %>% filter(year == y_out)),
                   d_all %>% filter(year == y_out))
  res$year <- y_out
  saveRDS(res, file.path(ldir, sprintf("loyo_peak_%d.rds", y_out)))
  print(res$scores)
  message(sprintf("Predictions: %.1f min. Saved %s/loyo_peak_%d.rds", mins(t0), ldir, y_out))
  if (!keep_fit) unlink(csv_dir, recursive = TRUE)    # the CSVs are the size of the fit
  message("Peak memory of R: ", round(peak_mem_gb(), 2), " GB")
}
