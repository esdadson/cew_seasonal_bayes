# =============================================================================
# step1_stan/spacetime_simple_brms.R
# brms version of the zone-free space-time model from standard parts
# (spacetime_simple.stan; written out in spacetime_simple_common.R):
#
#   log_CEW ~ c1 + s1 + c2 + s2 + c3 + s3 + trap_type + covariates
#           + gp(x_km, y_km, cov = "matern32")            spatial site intercept
#           + gp(x_km, y_km, by = c1) + gp(x_km, y_km, by = s1)   spatial site slopes
#           + (1 | site_f) + (1 | year_f)                 independent site, year
#           + ar(time = woy, gr = loc_year, p = 1)        AR(1) within site-year
#
# Same data, priors and scoring as the Stan version. Differences brms cannot
# express:
#   - each gp() has its own range; Stan shares one range across the three
#   - ar() treats consecutive observed weeks as one week apart (92% are);
#     Stan uses phi^gap for the 8% with gaps
#   - brms sigma is the AR innovation sd (Stan: sigma_innov); the first
#     observation of a site-year is not given the stationary variance
# gp() gives one latent value per distinct coordinate pair (gr = TRUE), i.e.
# per raw site; scale = FALSE keeps distances in km.
#
# Modes:
#   Rscript step1_stan/spacetime_simple_brms.R code      # print the generated Stan code
#   Rscript step1_stan/spacetime_simple_brms.R fit       # fit -> ST_OUT/brms/draws.rds
#   Rscript step1_stan/spacetime_simple_brms.R summary   # brms tables; spacetime_simple.R
#                                                        # summary compares with Stan
# Environment (plus the ST_* settings in spacetime_simple_common.R):
#   BRMS_CHAINS, BRMS_ITER, BRMS_WARMUP, BRMS_THIN   4, 2000, 1000, 2
#   BRMS_ADAPT_DELTA, BRMS_MAX_TREEDEPTH             0.9, 10
#   BRMS_NDRAWS     posterior draws used for peaks / held-out scores    400
#   ST_OUT          output folder                                results/spacetime_simple
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "code"
stopifnot(mode %in% c("code", "fit", "summary"))

source("step1_stan/spacetime_simple_common.R")
suppressPackageStartupMessages(library(brms))

out_dir <- file.path(env("ST_OUT", "results/spacetime_simple"), "brms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
smp <- list(chains = env_num("BRMS_CHAINS", 4), iter = env_num("BRMS_ITER", 2000),
            warmup = env_num("BRMS_WARMUP", 1000), thin = env_num("BRMS_THIN", 2),
            adapt_delta = env_num("BRMS_ADAPT_DELTA", 0.9), max_td = env_num("BRMS_MAX_TREEDEPTH", 10))

# =============================================================================
# 1. Formula and priors (as spacetime_simple.R) ----------------------------------------
# =============================================================================
gp_term <- function(by = NULL)
  paste0("gp(x_km, y_km, cov = 'matern32', scale = FALSE", if (!is.null(by)) paste0(", by = ", by), ")")
rhs <- c(harm_terms, "trap_type", covs, gp_term(), if (slopes) c(gp_term("c1"), gp_term("s1")),
         "(1 | site_f)", "(1 | year_f)", if (use_ar) "ar(time = woy, gr = loc_year, p = 1)")
f_brms <- bf(as.formula(paste("log_CEW ~", paste(rhs, collapse = " + "))))
print(f_brms)

d_fit <- d[fit_idx, ]
pr <- get_prior(f_brms, data = d_fit)
b_coefs <- pr$coef[pr$class == "b" & nzchar(pr$coef)]
# PC prior on range: range ~ inv_gamma(1, lambda_range); length-scale = range / 2
# ~ inv_gamma(1, lambda_range / 2)
priors <- c(
  set_prior("normal(0, 10)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  do.call(c, lapply(b_coefs[grepl("^[cs][0-9]+$", b_coefs)],
                    function(k) set_prior("normal(0, 2.5)", class = "b", coef = k))),
  do.call(c, lapply(unique(pr$coef[pr$class == "lscale" & nzchar(pr$coef)]), function(k)
    set_prior(sprintf("inv_gamma(1, %.4f)", lambda_range / 2), class = "lscale", coef = k))),
  set_prior("exponential(3)", class = "sdgp"),
  set_prior("exponential(3)", class = "sd"),
  set_prior("student_t(3, 0, 2.5)", class = "sigma"),
  if (use_ar) set_prior("uniform(-1, 1)", class = "ar", lb = -1, ub = 1))   # brms default: flat, unbounded

if (mode == "code") {
  print(priors)
  cat(make_stancode(f_brms, data = d_fit, prior = priors))
}

# =============================================================================
# 2. Fit ------------------------------------------------------------------------------
# =============================================================================
if (mode == "fit") {
  n_gp <- 1 + 2 * slopes
  inits <- function() c(
    setNames(lapply(1:n_gp, function(i) matrix(runif(1, 50, 150), 1, 1)), paste0("lscale_", 1:n_gp)),
    setNames(lapply(1:n_gp, function(i) array(0.3, 1)), paste0("sdgp_", 1:n_gp)),
    list(sigma = 0.8), if (use_ar) list(ar = array(0.5, 1)))
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d)", smp$chains, smp$iter, smp$warmup, smp$thin))
  t0 <- Sys.time()
  fit <- brm(f_brms, data = d_fit, prior = priors, backend = "cmdstanr",
             chains = smp$chains, cores = smp$chains, iter = smp$iter, warmup = smp$warmup, thin = smp$thin,
             init = inits, control = list(adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td),
             seed = 42, refresh = 50, file = file.path(out_dir, "fit_brms"), file_refit = "always")
  t_fit <- mins(t0)
  message(sprintf("Sampling: %.2f h", t_fit / 60))

  # Mean (no AR residual) for every row of every site-year; held-out year and
  # held-out-only sites drawn from their priors / kriged from neighbours
  nd  <- min(env_num("BRMS_NDRAWS", 400), ndraws(fit))
  ids <- sort(sample(ndraws(fit), nd))
  # ar() groups in new data must be site-years brms has seen, even with
  # incl_autocor = FALSE (the AR term is then unused). Held-out site-years
  # therefore go in a second call, each under a different training site-year
  # label so week numbers stay unique within a group.
  nd_all <- d %>% mutate(log_CEW = ifelse(is.na(log_CEW), 0, log_CEW))   # response unused
  seen <- nd_all$loc_year %in% unique(d_fit$loc_year)
  if (use_ar && any(!seen)) {
    new_sy <- unique(nd_all$loc_year[!seen])
    stand_in <- unique(d_fit$loc_year)[seq_along(new_sy)]
    nd_all$loc_year[!seen] <- stand_in[match(nd_all$loc_year[!seen], new_sy)]
  }
  lp <- function(rows) posterior_linpred(fit, newdata = nd_all[rows, ], allow_new_levels = TRUE,
                                         sample_new_levels = "gaussian", incl_autocor = FALSE, draw_ids = ids)
  t0 <- Sys.time()
  ETA <- matrix(NA_real_, nd, nrow(d))
  ETA[, seen] <- lp(which(seen))
  if (any(!seen)) ETA[, !seen] <- lp(which(!seen))
  message(sprintf("Linear predictor for %d rows x %d draws: %.1f min", nrow(d), nd, mins(t0)))
  PK <- t(apply(ETA, 1, function(e) win[1] - 1 + apply(matrix(e, 52)[win[1]:win[2], , drop = FALSE], 2, which.max)))
  dr <- as_draws_matrix(fit)[ids, ]
  bcols <- grep("^b_", colnames(dr), value = TRUE)
  Bz <- matrix(dr[, bcols], nd, dimnames = list(NULL, sub("^trap_type", "", sub("^b_", "", bcols))))
  Bz <- Bz[, colnames(X)]
  # marginal residual sd = innovation sd / sqrt(1 - phi^2), as Stan's sigma
  sig <- as.vector(dr[, "sigma"])
  if (use_ar) sig <- sig / sqrt(1 - as.vector(dr[, "ar[1]"])^2)
  saveRDS(list(PK = PK, ET = ETA[, test_idx, drop = FALSE], sigma = sig, Bz = Bz,
               sy_loc_year = sy$loc_year, time_min = t_fit,
               setup = list(H = H, slopes = slopes, ar = use_ar, range0 = range0, win = win, n_sub = n_sub)),
          file.path(out_dir, "draws.rds"))
  print(summary(fit), digits = 3)
  message("Saved ", file.path(out_dir, "draws.rds"))
}

# =============================================================================
# 3. Summary --------------------------------------------------------------------------
# =============================================================================
if (mode == "summary") {
  fit <- readRDS(file.path(out_dir, "fit_brms.rds"))
  rb  <- readRDS(file.path(out_dir, "draws.rds"))
  print(summary(fit), digits = 3)
  ls <- posterior::summarise_draws(posterior::subset_draws(as_draws_df(fit), variable = "^lscale_", regex = TRUE))
  message("Ranges (km) = 2 x lscale:")
  print(transform(ls[, c("variable", "median", "q5", "q95")], median = 2 * median, q5 = 2 * q5, q95 = 2 * q95))
  ho <- holdout_metrics("brms", rb$ET, rb$sigma)
  zs <- zone_peak_stats("brms", rb$PK)
  print(ho, digits = 4); print(zs)
  write.csv(ho, file.path(out_dir, "holdout.csv"), row.names = FALSE)
  write.csv(zs, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)
  message("Summary written to ", out_dir)
}
