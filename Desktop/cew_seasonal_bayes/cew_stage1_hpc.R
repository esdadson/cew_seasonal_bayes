# =============================================================================
# cew_stage1_hpc.R
# The stage-1 seasonal model (cew_stage1_data.R) as a Slurm job on Hazel.
# The pipeline then only reads the saved fit (copy fit1_gp21_train.rds back).
#
# Always run inside a Slurm job (the script refuses otherwise on Linux):
#   sbatch hpc/stage1_dryrun.sbatch                  # dry run
#   ./hpc/submit_stage1.sh                           # full fit
# Modes:
#   Rscript cew_stage1_hpc.R dryrun   # compile + one short chain on the full
#                                     # training data -> suggested resources
#   Rscript cew_stage1_hpc.R fit      # full fit -> fit1_gp21_train.rds
#
# Environment:
#   STAGE1_CHAINS, STAGE1_ITER, STAGE1_WARMUP, STAGE1_THIN
#             full-run settings (defaults: fit_stage1() in cew_stage1_data.R);
#             the dry run sizes resources for these
#   DRY_WARMUP, DRY_SAMPLING   dry-run chain length (default 150 + 50)
#   CEW_DATA, CEW_OUT          as for cew_covariates_hpc.R
#
# Chains run in parallel, one core each, so wall time ~ one chain's time.
# Time per iteration = leapfrog steps x seconds per gradient. The dry run
# measures seconds per gradient; a NUTS iteration takes at most
# 2^max_treedepth - 1 = 1023 steps, which bounds the full run's time.
# =============================================================================

# Hazel AUP: no code on login nodes. On Linux, refuse to run outside a Slurm job.
if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with `sbatch hpc/stage1_dryrun.sbatch` or ",
       "`./hpc/submit_stage1.sh` - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit"))

t_start <- Sys.time()
source("cew_stage1_data.R")
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("Data: %d training rows (%s), %d sites; %.1f min",
                nrow(d_train), paste(range(d_train$year), collapse = "-"),
                nrow(distinct(d_train, longitude, latitude)), mins(t_start)))
message("CmdStan ", cmdstanr::cmdstan_version(), " at ", cmdstanr::cmdstan_path())

env_int <- function(name, default) as.integer(Sys.getenv(name, default))
def  <- formals(fit_stage1)
full <- list(chains = env_int("STAGE1_CHAINS", def$chains), iter = env_int("STAGE1_ITER", def$iter),
             warmup = env_int("STAGE1_WARMUP", def$warmup), thin = env_int("STAGE1_THIN", def$thin))
max_td <- 10                                  # NUTS default; fit_stage1() doesn't change it

peak_mem_gb <- function() {          # Linux only: peak resident memory of this R process
  s <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
  kb <- as.numeric(gsub("\\D", "", grep("^VmHWM", s, value = TRUE)))
  if (length(kb)) kb / 1024^2 else NA_real_
}

if (mode == "dryrun") {
  # Same Stan program and data brm() would build, sampled directly with cmdstanr
  # so compile time, warmup and sampling are timed separately
  t0 <- Sys.time()
  code  <- make_stancode(f_my, data = d_train, family = gaussian(), prior = pri_my)
  sdata <- make_standata(f_my, data = d_train, family = gaussian(), prior = pri_my)
  mod   <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  t_compile <- mins(t0)

  w <- env_int("DRY_WARMUP", 150); s <- env_int("DRY_SAMPLING", 50)
  fit <- mod$sample(data = sdata, chains = 1, iter_warmup = w, iter_sampling = s,
                    seed = 42, adapt_delta = 0.95, max_treedepth = max_td,
                    save_warmup = TRUE, refresh = 10)
  tm <- fit$time()$chains
  dg <- posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))
  lf <- dg$n_leapfrog__; in_w <- seq_along(lf) <= w
  sec_grad <- tm$total / sum(lf)              # seconds per leapfrog step (one gradient)

  # Projections for one full-run chain, in hours
  steps_obs <- mean(lf[!in_w])                # once adapted (short adaptation: rough)
  h_obs <- full$iter * steps_obs * sec_grad / 3600
  h_max <- full$iter * (2^max_td - 1) * sec_grad / 3600
  h_req <- 1.15 * h_max + 1                   # + data, compile, reading the draws
  qos   <- if (h_req <= 4 * 24) "normal" else if (h_req <= 10 * 24) "long" else NA
  fmt_t <- function(h) { m <- ceiling(min(h, 240) * 60)
    sprintf("%d-%02d:%02d:00", m %/% 1440, (m %% 1440) %/% 60, m %% 60) }
  mem <- peak_mem_gb()

  message("\n--- Dry run (1 chain, ", w, " warmup + ", s, " sampling iterations) ---")
  message(sprintf("Compile: %.1f min. Warmup %.1f s/iteration, sampling %.1f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps_obs, 2^max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * mean(dg$treedepth__[!in_w] >= max_td), sum(dg$divergent__[!in_w])))
  message("Peak memory of R (data + draws): ", round(mem, 2), " GB")

  message(sprintf("\n--- Full run: %d chains x %d iterations (%d warmup, thin %d) ---",
                  full$chains, full$iter, full$warmup, full$thin))
  message(sprintf("Per chain at this run's tree depth: %.0f h; if every iteration hits max depth: %.0f h",
                  h_obs, h_max))
  message("\n--- Suggested settings for hpc/submit_stage1.sh ---")
  if (is.na(qos)) {
    it_fit <- floor((10 * 24 - 1) / 1.15 * 3600 / ((2^max_td - 1) * sec_grad))
    message("Worst case exceeds Hazel's 10-day limit (qos=long). Lower STAGE1_ITER to <= ",
            it_fit, " (warmup included), or use the tree-depth estimate above at your own risk.")
  }
  message("STAGE1_TIME=", fmt_t(h_req), "   STAGE1_QOS=", ifelse(is.na(qos), "long", qos),
          "   STAGE1_MEM=", if (is.na(mem)) "16G" else sprintf("%dG", max(8, ceiling(2 * mem) + full$chains)),
          "   # --cpus-per-task = STAGE1_CHAINS")
  message("Check seff <jobid> (Memory Utilized) before the full run.")

} else if (mode == "fit") {
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d) on %s cores",
                  full$chains, full$iter, full$warmup, full$thin, Sys.getenv("SLURM_CPUS_PER_TASK", "?")))
  # CmdStan's CSV files are kept in stan_csv/: if R fails after sampling, the
  # draws are not lost (read them with cmdstanr::as_cmdstan_fit())
  dir.create("stan_csv", showWarnings = FALSE)
  t0 <- Sys.time()
  fit1 <- fit_stage1(d_train, chains = full$chains, iter = full$iter, warmup = full$warmup,
                     thin = full$thin, file_refit = "always", output_dir = "stan_csv", refresh = 100)
  message(sprintf("Sampling + saving: %.1f h. Saved %s.rds", mins(t0) / 60, stage1_file))

  print(summary(fit1))
  np <- nuts_params(fit1)
  message(sprintf("Divergent transitions: %d; at max tree depth: %.0f%%; max Rhat: %.3f",
                  sum(np$Value[np$Parameter == "divergent__"]),
                  100 * mean(np$Value[np$Parameter == "treedepth__"] >= max_td),
                  max(rhat(fit1), na.rm = TRUE)))
  message("Peak memory of R: ", round(peak_mem_gb(), 2), " GB")
}
