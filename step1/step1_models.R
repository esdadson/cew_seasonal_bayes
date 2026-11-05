# =============================================================================
# step1/step1_models.R
# Step 1 candidate models of weekly log(CEW + 1), their priors, the sampler,
# and the comparison metrics. Built up one block at a time; each block must
# earn its place (loo_compare + the held-out year) before the next is added.
# Sourced by step1/step1_compare.R. Run from the project directory.
#
#  y_itw = log(CEW_itw + 1) ~ Normal(mu_itw, sigma)        i site, t year, w ISO week
#  fixed    : mu = a_zone + b_trap + x_itw' beta + u_t
#             x = weekly EVI, weekly GDD, six crop shares (standardised);
#             u_t ~ N(0, sd_year), so a held-out year can be predicted
#  seasonal : fixed + sum_k f_k(w) (g_k + g_zk)
#             f_1..f_4 = 2-harmonic Fourier basis of week, orthogonalised
#             (step1_data.R section 5: f1, f2 span the first harmonic, f3, f4
#             add the second); g_k = seasonal curve shared by all sites,
#             g_zk ~ N(0, sd_k) = zone deviations (|| = uncorrelated)
#  Later (not yet fit): seasonal sigma (sigma ~ f1 + f2), site and spatial
#  site intercepts. A model is list(mu = <terms>, sigma = <terms or NULL>).
# =============================================================================

source("step1/step1_data.R")
library(brms); library(posterior)
options(brms.backend = "cmdstanr")

step1_dir <- Sys.getenv("STEP1_OUT", "results/step1")
dir.create(step1_dir, recursive = TRUE, showWarnings = FALSE)

# Held-out validation year (STEP1_HOLDOUT=none: fit every year)
ho      <- Sys.getenv("STEP1_HOLDOUT", "2020")
holdout <- if (ho %in% c("", "none")) NA_integer_ else as.integer(ho)
d_fit   <- if (is.na(holdout)) d_obs else filter(d_obs, year != holdout)
d_test  <- if (is.na(holdout)) d_obs[0, ] else filter(d_obs, year == holdout)

# =============================================================================
# Candidate models -------------------------------------------------------------
# =============================================================================
year_term <- "(1 | year_f)"
covs      <- "EVI_filled + gdd + corn + cotton + soybeans + tobacco + peanuts + other"
four_mu   <- c("f1 + f2 + f3 + f4", "zone_fac:(f1 + f2 + f3 + f4)")

terms_fixed <- c("0 + zone_fac + trap_type", year_term, covs)
step1_specs <- list(
  fixed    = list(mu = terms_fixed),
  seasonal = list(mu = c(terms_fixed, four_mu))
)

has_term <- function(model, term) term %in% step1_specs[[model]]$mu
step1_formula <- function(model) {
  s <- step1_specs[[model]]
  f <- as.formula(paste("log_CEW ~", paste(s$mu, collapse = " + ")))
  if (is.null(s$sigma)) bf(f) else bf(f, as.formula(paste("sigma ~", paste(s$sigma, collapse = " + "))))
}

# brms silently drops rows with NA in a formula variable; fail early instead
na_vars <- unique(unlist(lapply(names(step1_specs), function(m) all.vars(step1_formula(m)$formula))))
na_vars <- na_vars[colSums(is.na(d_obs[intersect(na_vars, names(d_obs))])) > 0]
if (length(na_vars)) stop("NA in model variables of the observed rows: ", paste(na_vars, collapse = ", "))

# Priors (weakly informative; covariates and f's have unit sd)
#   b ~ N(0, 1): one sd of a covariate moves log(CEW + 1) by < ~2
#   zone intercepts ~ N(0, 3); sigma, sd_year ~ half-t(3, 0, 1)
#   seasonal curve g_k ~ N(0, 2): allows a seasonal swing of several log units
#   zone deviations sd_k ~ half-t(3, 0, 0.5): three zones barely inform an sd,
#   so deviations are kept to fractions of the shared curve
step1_priors <- function(model) {
  s <- step1_specs[[model]]
  p <- c(prior(normal(0, 1), class = "b"),
         prior(normal(0, 3), class = "b", coef = "zone_fac0"),
         prior(normal(0, 3), class = "b", coef = "zone_fac1"),
         prior(normal(0, 3), class = "b", coef = "zone_fac2"),
         prior(student_t(3, 0, 1), class = "sd", group = "year_f"))
  p <- c(p, if (is.null(s$sigma)) prior(student_t(3, 0, 1), class = "sigma") else
            c(prior(normal(0, 1), class = "Intercept", dpar = "sigma"),
              prior(normal(0, 0.5), class = "b", dpar = "sigma")))
  if (has_term(model, four_mu[1]))
    p <- c(p, prior(normal(0, 2), class = "b", coef = "f1"), prior(normal(0, 2), class = "b", coef = "f2"),
           prior(normal(0, 2), class = "b", coef = "f3"), prior(normal(0, 2), class = "b", coef = "f4")
  p
}

# =============================================================================
# Sampler -------------------------------------------------------------------------
#   Defaults: 4 chains x 2000 iterations (1000 warmup), thin 4 = 1000 draws.
#   Enough for these models (no site terms); override with STEP1_ITER etc.
# =============================================================================
env_int <- function(k, d) as.integer(Sys.getenv(k, d))
step1_sampler <- list(chains = env_int("STEP1_CHAINS", 4), iter = env_int("STEP1_ITER", 2000),
                      warmup = env_int("STEP1_WARMUP", 1000), thin = env_int("STEP1_THIN", 4),
                      threads = env_int("STEP1_THREADS", 1))
step1_control <- list(adapt_delta = 0.95, max_treedepth = 12)

fit_file <- function(model) file.path(step1_dir, paste0("fit_", model))   # brms adds .rds

fit_step1 <- function(model, data = d_fit, smp = step1_sampler) {
  brm(step1_formula(model), data = droplevels(data), family = gaussian(),
      prior = step1_priors(model), chains = smp$chains, cores = smp$chains,
      iter = smp$iter, warmup = smp$warmup, thin = smp$thin,
      threads = if (smp$threads > 1) threading(smp$threads) else NULL,
      seed = 42, init = 0, control = step1_control,
      refresh = max(smp$iter %/% 10, 1),
      file = fit_file(model), file_refit = "on_change")
}

# Sampling time and convergence
sampler_stats <- function(fit, model, smp = step1_sampler) {
  np <- nuts_params(fit)
  v  <- function(p) np$Value[np$Parameter == p]
  et <- rstan::get_elapsed_time(fit$fit)                       # chains x (warmup, sample), s
  sm <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  steps <- mean(v("n_leapfrog__"))
  data.frame(model = model,
             warmup_min = mean(et[, 1]) / 60, sampling_min = mean(et[, 2]) / 60,
             steps_per_iter = steps,
             pct_max_treedepth = 100 * mean(v("treedepth__") >= step1_control$max_treedepth),
             divergent = sum(v("divergent__")),
             ms_per_grad = 1000 * mean(et[, 2]) / ((smp$iter - smp$warmup) * steps),
             max_rhat = max(sm$rhat, na.rm = TRUE),
             min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
             min_ess_tail = min(sm$ess_tail, na.rm = TRUE),
             worst = sm$variable[which.min(sm$ess_bulk)])
}

# =============================================================================
# Predicting a year the fit never saw ------------------------------------------
#   Its year effect is unknown: per draw, u_new ~ N(0, sd_year), shared by
#   every site that year. Training years use their estimated effects.
# =============================================================================
new_year_draws <- function(fit, seed = 42) {
  s <- as_draws_df(fit, variable = "sd_year_f__Intercept")$sd_year_f__Intercept
  set.seed(seed); s * rnorm(length(s))
}
re_no_year <- function(model) {
  re <- setdiff(grep("|", step1_specs[[model]]$mu, fixed = TRUE, value = TRUE), year_term)
  if (length(re)) as.formula(paste("~", paste(re, collapse = " + "))) else NA
}
# E[log(CEW + 1)], draws x rows
epred_rows <- function(fit, model, nd, u_new = new_year_draws(fit)) {
  seen <- nd$year_f %in% levels(fit$data$year_f)
  out  <- matrix(NA_real_, ndraws(fit), nrow(nd))
  if (any(seen)) out[, seen] <- posterior_epred(fit, newdata = droplevels(nd[seen, ]))
  if (any(!seen))
    out[, !seen] <- posterior_epred(fit, newdata = droplevels(nd[!seen, ]), re_formula = re_no_year(model),
                                    allow_new_levels = TRUE) + u_new
  out
}
# sigma, draws x rows (constant, or modelled from fixed terms only)
sigma_rows <- function(fit, model, nd) {
  if (is.null(step1_specs[[model]]$sigma))
    return(matrix(as_draws_df(fit, variable = "sigma")$sigma, ndraws(fit), nrow(nd)))
  posterior_linpred(fit, newdata = droplevels(nd), dpar = "sigma", transform = TRUE, re_formula = NA)
}

# =============================================================================
# Metrics ------------------------------------------------------------------------
#   Training years: elpd_loo (PSIS-LOO), bayes_R2, loo_R2, MAE / RMSE of the
#   posterior mean, MAE of the LOO mean. Ordinary LOO leaves out one trap-week;
#   weeks of a site-year are autocorrelated, so it is optimistic - the
#   held-out year is the stricter test.
#   Held-out year: elpd_test (sum log mean_s p(y | s)), MAE, RMSE, R2 around
#   the held-out mean, cov90 = share inside the 90% predictive interval,
#   overall and by season block (cov90 per block shows whether a constant
#   sigma is too wide off-season and too narrow at the peak).
# =============================================================================
season_block <- function(woy) cut(woy, c(0, 22, 27, 32, 37, 53),
                                  labels = c("<=22", "23-27", "28-32", "33-37", ">=38"))

step1_metrics <- function(fit, model, lo) {
  y   <- fit$data$log_CEW
  mu  <- colMeans(posterior_epred(fit))
  br2 <- bayes_R2(fit, summary = TRUE)
  lr2 <- loo_R2(fit, summary = TRUE)
  mu_loo <- loo_predict(fit, type = "mean", psis_object = lo$psis_object)
  est <- lo$estimates
  m <- data.frame(model = model, n_train = length(y),
                  elpd_loo = est["elpd_loo", 1], elpd_loo_se = est["elpd_loo", 2],
                  p_loo = est["p_loo", 1], pareto_k_gt_0.7 = sum(lo$diagnostics$pareto_k > 0.7),
                  bayes_R2 = br2[1, "Estimate"], loo_R2 = lr2[1, "Estimate"],
                  mae_train = mean(abs(y - mu)), rmse_train = sqrt(mean((y - mu)^2)),
                  mae_loo = mean(abs(y - mu_loo)))
  if (!nrow(d_test)) return(list(metrics = m, test_rows = NULL))

  et  <- epred_rows(fit, model, d_test)
  sg  <- sigma_rows(fit, model, d_test)
  yt  <- d_test$log_CEW
  Y   <- matrix(yt, nrow(et), ncol(et), byrow = TRUE)
  lpd <- matrixStats::colLogSumExps(dnorm(Y, et, sg, log = TRUE)) - log(nrow(et))
  set.seed(42)
  pi90 <- matrixStats::colQuantiles(et + sg * matrix(rnorm(length(et)), nrow(et)), probs = c(0.05, 0.95))
  mt  <- colMeans(et)
  m <- cbind(m, data.frame(n_test = length(yt), elpd_test = sum(lpd),
                           elpd_test_se = sd(lpd) * sqrt(length(lpd)),
                           mae_test = mean(abs(yt - mt)), rmse_test = sqrt(mean((yt - mt)^2)),
                           R2_test = 1 - sum((yt - mt)^2) / sum((yt - mean(yt))^2),
                           cov90_test = mean(yt >= pi90[, 1] & yt <= pi90[, 2])))
  test_rows <- d_test %>%
    select(location2, loc_year, zone_30y, year, woy, log_CEW) %>%
    mutate(model = model, pred_mean = mt, pi_lo = pi90[, 1], pi_hi = pi90[, 2], lpd = lpd,
           block = season_block(woy))
  list(metrics = m, test_rows = test_rows)
}
