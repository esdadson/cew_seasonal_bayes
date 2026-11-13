# =============================================================================
# step1_stan/zone_harmonic_brms.R
# brms version of the zone harmonic model (zone_harmonic.stan; written out in
# zone_harmonic_common.R):
#
#   log_CEW ~ 0 + zone + zone:(c1 + s1 + c2 + s2)        zone seasonal curves
#           + trap_type + covariates
#           + gp(x_km, y_km, cov = "matern32")           short-range site field
#           + gp(x_km, y_km, by = c1) + gp(x_km, y_km, by = s1)   (ZH_SVC = 1)
#           + (1 | year_f) + (1 | site_f)                year, site nugget
#
# Same data, priors and scoring as the Stan version, with three differences
# brms cannot express:
#   - each gp() has its own range (lscale) and sd; Stan shares one range
#     across fields and one sd between the cos / sin fields
#   - no sum-to-zero within zone: the fields can drift towards zone means;
#     the short-range prior is what discourages it
#   - field and nugget are separate terms (same model, other parameterisation)
# gp() gives one latent value per distinct coordinate pair (gr = TRUE), i.e.
# per raw site; scale = FALSE keeps distances in km so the range prior is in km.
#
# Modes:
#   Rscript step1_stan/zone_harmonic_brms.R code      # print the generated Stan code
#   Rscript step1_stan/zone_harmonic_brms.R fit       # fit -> ZH_OUT/brms/draws.rds
#   Rscript step1_stan/zone_harmonic_brms.R summary   # brms tables; zone_harmonic.R
#                                                     # summary compares with Stan
# Environment (plus the ZH_* settings in zone_harmonic_common.R):
#   BRMS_CHAINS, BRMS_ITER, BRMS_WARMUP, BRMS_THIN   4, 2000, 1000, 2
#   BRMS_ADAPT_DELTA, BRMS_MAX_TREEDEPTH             0.9, 10
#   BRMS_NDRAWS     posterior draws used for peaks / held-out scores    400
#   ZH_OUT          output folder                                  results/zone_harmonic
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "code"
stopifnot(mode %in% c("code", "fit", "summary"))

source("step1_stan/zone_harmonic_common.R")
suppressPackageStartupMessages(library(brms))

out_dir <- file.path(env("ZH_OUT", "results/zone_harmonic"), "brms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
smp <- list(chains = env_num("BRMS_CHAINS", 4), iter = env_num("BRMS_ITER", 2000),
            warmup = env_num("BRMS_WARMUP", 1000), thin = env_num("BRMS_THIN", 2),
            adapt_delta = env_num("BRMS_ADAPT_DELTA", 0.9), max_td = env_num("BRMS_MAX_TREEDEPTH", 10))

# =============================================================================
# 1. Formula and priors (as zone_harmonic.R) ------------------------------------------
# =============================================================================
gp_term <- function(by = NULL)
  paste0("gp(x_km, y_km, cov = 'matern32', scale = FALSE", if (!is.null(by)) paste0(", by = ", by), ")")
rhs <- c("0 + zone", paste0("zone:(", harm_terms, ")"), "trap_type", covs, gp_term(),
         if (svc) c(gp_term("c1"), gp_term("s1")), "(1 | year_f)", "(1 | site_f)")
f_brms <- bf(as.formula(paste("log_CEW ~", paste(rhs, collapse = " + "))))
print(f_brms)

d_fit <- d[fit_idx, ]
pr <- get_prior(f_brms, data = d_fit)
b_coefs <- pr$coef[pr$class == "b" & nzchar(pr$coef)]
lscale_mu <- range_mu - log(2)                       # range = 2 x length-scale
priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  do.call(c, lapply(b_coefs[b_coefs %in% paste0("zone", levels(d$zone))],
                    function(k) set_prior("normal(0, 10)", class = "b", coef = k))),
  do.call(c, lapply(b_coefs[grepl(":", b_coefs)],
                    function(k) set_prior("normal(0, 2.5)", class = "b", coef = k))),
  # per gp() term: brms otherwise keeps its own data-based inverse-gamma default
  do.call(c, lapply(unique(pr$coef[pr$class == "lscale" & nzchar(pr$coef)]), function(k)
    set_prior(sprintf("lognormal(%.4f, %.4f)", lscale_mu, range_sd), class = "lscale", coef = k))),
  set_prior("exponential(3)", class = "sdgp"),
  set_prior("exponential(3)", class = "sd"),
  set_prior("student_t(3, 0, 2.5)", class = "sigma"))

if (mode == "code") {
  print(priors)
  cat(make_stancode(f_brms, data = d_fit, prior = priors))
}

# =============================================================================
# 2. Fit ------------------------------------------------------------------------------
# =============================================================================
if (mode == "fit") {
  n_gp <- 1 + 2 * svc
  inits <- function() c(
    setNames(lapply(1:n_gp, function(i) matrix(exp(lscale_mu) * runif(1, 0.8, 1.25), 1, 1)),
             paste0("lscale_", 1:n_gp)),
    setNames(lapply(1:n_gp, function(i) array(0.3, 1)), paste0("sdgp_", 1:n_gp)),
    list(sigma = 1.2))
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d)", smp$chains, smp$iter, smp$warmup, smp$thin))
  t0 <- Sys.time()
  fit <- brm(f_brms, data = d_fit, prior = priors, backend = "cmdstanr",
             chains = smp$chains, cores = smp$chains, iter = smp$iter, warmup = smp$warmup, thin = smp$thin,
             init = inits, control = list(adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td),
             seed = 42, refresh = 50, file = file.path(out_dir, "fit_brms"), file_refit = "always")
  t_fit <- mins(t0)
  message(sprintf("Sampling: %.2f h", t_fit / 60))

  # Linear predictor for every row of every site-year (held-out year and any
  # held-out-only sites drawn from their priors / kriged from neighbours)
  nd <- min(env_num("BRMS_NDRAWS", 400), ndraws(fit))
  t0 <- Sys.time()
  ids <- sort(sample(ndraws(fit), nd))
  ETA <- posterior_linpred(fit, newdata = d, allow_new_levels = TRUE, sample_new_levels = "gaussian",
                           draw_ids = ids)
  message(sprintf("Linear predictor for %d rows x %d draws: %.1f min", nrow(d), nd, mins(t0)))
  PK <- t(apply(ETA, 1, function(e) win[1] - 1 + apply(matrix(e, 52)[win[1]:win[2], , drop = FALSE], 2, which.max)))
  dr <- as_draws_matrix(fit)[ids, ]
  bcols <- grep("^b_", colnames(dr), value = TRUE)
  Bz <- matrix(dr[, bcols], nd, dimnames = list(NULL, gsub("^b_|zone|trap_type", "", bcols)))
  Bz <- Bz[, colnames(X)]                            # same column order as the Stan fit
  saveRDS(list(PK = PK, ET = ETA[, test_idx, drop = FALSE], sigma = as.vector(dr[, "sigma"]), Bz = Bz,
               sy_loc_year = sy$loc_year, time_min = t_fit,
               setup = list(H = H, svc = svc, range = c(range_lo, range_hi), win = win, n_sub = n_sub)),
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
  # lscale is a length-scale: range (correlation ~0.14) = 2 x lscale
  ls <- posterior::summarise_draws(posterior::subset_draws(as_draws_df(fit), variable = "^lscale_", regex = TRUE))
  message("Ranges (km) = 2 x lscale:"); print(transform(ls[, c("variable", "median", "q5", "q95")],
                                                    median = 2 * median, q5 = 2 * q5, q95 = 2 * q95))
  ho <- holdout_metrics("brms", rb$ET, rb$sigma)
  zs <- zone_peak_stats("brms", rb$PK)
  zc <- zone_curves(rb$Bz)
  print(ho, digits = 4); print(zs)
  write.csv(ho, file.path(out_dir, "holdout.csv"), row.names = FALSE)
  write.csv(zs, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)
  plot_zone_curves(zc, file.path(out_dir, "zone_curves.png"),
                   sprintf("brms: zone seasonal curves, %d harmonic(s) (95%% CrI)", H))
  message("Summary written to ", out_dir)
}
