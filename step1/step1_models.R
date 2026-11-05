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
#             add the second); g_k = the Northern-limits curve, g_zk (fixed,
#             z = Transitional, Southern) = each zone's difference from it:
#             one curve per zone, 12 coefficients
#  spatial  : seasonal + v_i, a spatially correlated site effect (exact GP
#             over the monitored sites with a fixed range):
#               v ~ N(0, sd_site^2 K),  K = rho C + (1 - rho) I,
#               C_ij = exp(-d_ij / phi),  d = CONUS Albers km between sites
#             rho (spatial share) and phi (range) are estimated once from the
#             correlogram of site-mean residuals (section "Spatial site
#             effect") and held fixed, so the model adds ONE parameter
#             (sd_site), samples like an independent site intercept, and an
#             unmonitored location can be kriged from its lon/lat:
#               v0 | v ~ N(rho c0' K^-1 v, sd_site^2 (1 - rho^2 c0' K^-1 c0))
#  spatial_gp : seasonal + gp(longitude, latitude, gr = TRUE, iso = TRUE,
#             scale = TRUE), brms' exact GP (no k): computed over the unique
#             site coordinates (gr), one length-scale for both axes (iso),
#             coordinates scaled so the largest distance is 1 (scale);
#             squared-exponential kernel (brms default). Estimates sd_gp and
#             the length-scale (brms' default priors) plus one latent value
#             per site. Cost per gradient grows with sites^3.
#             Note: in degrees, 1 deg of longitude (~80-95 km here) is shorter
#             than 1 deg of latitude (111 km), so iso = TRUE is only
#             approximately isotropic; x_km, y_km would be exact.
#  Later (not yet fit): seasonal sigma (sigma ~ f1 + f2).
#  A model is list(mu = <terms>, sigma = <terms or NULL>).
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

site_term <- "(1 | gr(location2, cov = K))"
gp_term   <- "gp(longitude, latitude, gr = TRUE, iso = TRUE, scale = TRUE)"

terms_fixed <- c("0 + zone_fac + trap_type", year_term, covs)
step1_specs <- list(
  fixed    = list(mu = terms_fixed),
  seasonal = list(mu = c(terms_fixed, four_mu)),
  spatial  = list(mu = c(terms_fixed, four_mu, site_term)),
  spatial_gp = list(mu = c(terms_fixed, four_mu, gp_term))
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

# =============================================================================
# Spatial site effect: range and spatial share ---------------------------------
#   Trend = the seasonal model's terms with years as fixed effects, by OLS on
#   the fitting rows. Site-mean residuals, standardised; their products
#   averaged by distance bin give the correlogram; rho exp(-h / phi) is fit by
#   weighted least squares (weights = site pairs per bin). STEP1_RHO and
#   STEP1_PHI (km) override the estimates (e.g. a range sensitivity run).
#   Earlier estimate (2008-2019): rho 0.76, phi 46 km.
# =============================================================================
env_num  <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) as.numeric(v) else d }
dist_km  <- function(a, b) sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
corr_exp <- function(D, phi) exp(-D / phi)

correlogram <- function(r, d, breaks = c(0, 5, 10, 15, 20, 30, 40, 50, 65, 80, 100, 130, 160, 200, 300)) {
  sm <- d %>% mutate(r = r) %>% group_by(location2) %>%
    summarise(r = mean(r), x = first(x_km), y = first(y_km), .groups = "drop")
  D  <- dist_km(as.matrix(sm[, c("x", "y")]), as.matrix(sm[, c("x", "y")]))
  z  <- (sm$r - mean(sm$r)) / sd(sm$r)
  ut <- upper.tri(D)
  tibble(dist = D[ut], prod = outer(z, z)[ut]) %>%
    filter(dist > 0, dist <= max(breaks)) %>%
    mutate(bin = cut(dist, breaks)) %>% group_by(bin) %>%
    summarise(h = mean(dist), pairs = n(), cor = mean(prod), .groups = "drop")
}

site_xy <- d_fit %>% distinct(location2, x_km, y_km) %>% arrange(location2)
stopifnot(!anyDuplicated(site_xy$location2))
trend <- lm(log_CEW ~ 0 + zone_fac + trap_type + year_f + EVI_filled + gdd + corn + cotton + soybeans +
              tobacco + peanuts + other + zone_fac:(f1 + f2 + f3 + f4), data = d_fit)
cg_trend <- correlogram(resid(trend), d_fit)
cfit <- tryCatch(nls(cor ~ rho * exp(-h / phi), data = cg_trend, weights = pairs,
                     start = list(rho = 0.8, phi = 40)),
                 error = function(e) { message("Correlogram fit failed (", conditionMessage(e),
                                               "); using rho 0.76, phi 46 km"); NULL })
sp_est  <- if (is.null(cfit)) c(rho = 0.76, phi = 46) else coef(cfit)
sp_rho  <- env_num("STEP1_RHO", sp_est[["rho"]])
sp_phi  <- env_num("STEP1_PHI", sp_est[["phi"]])
stopifnot(sp_rho > 0, sp_rho < 1, sp_phi > 0)
message(sprintf("Spatial site effect: %d sites; correlogram rho = %.2f, phi = %.0f km; using rho = %.2f, phi = %.0f km",
                nrow(site_xy), sp_est[["rho"]], sp_est[["phi"]], sp_rho, sp_phi))

# K over the fitting sites, rows/columns named by location2 (brms matches the
# names to the levels of location2)
K_sites <- local({
  XY <- as.matrix(site_xy[, c("x_km", "y_km")])
  K  <- sp_rho * corr_exp(dist_km(XY, XY), sp_phi) + (1 - sp_rho) * diag(nrow(XY))
  dimnames(K) <- list(site_xy$location2, site_xy$location2)
  K
})

# Held-out sites never fit would get an independent N(0, sd_site) draw from
# brms (allow_new_levels ignores K); flag them
new_test_sites <- setdiff(unique(d_test$location2), site_xy$location2)
if (length(new_test_sites))
  message(length(new_test_sites), " held-out site(s) not in the fitting years: their spatial ",
          "effect is drawn independently, not kriged")

# Priors (weakly informative; covariates and f's have unit sd)
#   b ~ N(0, 1): one sd of a covariate moves log(CEW + 1) by < ~2
#   zone intercepts ~ N(0, 3); sigma, sd_year ~ half-t(3, 0, 1)
#   seasonal curve g_k ~ N(0, 2): allows a seasonal swing of several log units.
#   zone_fac:(f1..f4) is coded against the reference zone (0, Northern limits),
#   because f1..f4 are in the model too: f_k = the Northern curve, and
#   zone_fac1:f_k, zone_fac2:f_k = how the Transitional / Southern curves
#   differ from it (the default b ~ N(0, 1), shrinking toward a shared shape).
#   Fixed effects, not random: three zones can't inform a standard deviation,
#   and the random version (|| zone_fac) sampled very slowly.
#   sd_site ~ half-t(3, 0, 1), like sd_year.
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
           prior(normal(0, 2), class = "b", coef = "f3"), prior(normal(0, 2), class = "b", coef = "f4"))
  if (has_term(model, site_term))
    p <- c(p, prior(student_t(3, 0, 1), class = "sd", group = "location2"))
  p
}

# =============================================================================
# Sampler -------------------------------------------------------------------------
#   Defaults: 4 chains x 2000 iterations (1000 warmup), thin 4 = 1000 draws;
#   override with STEP1_ITER etc. Check min_ess_bulk in comparison.csv: the
#   spatial model (one effect per site) may need longer chains.
# =============================================================================
env_int <- function(k, d) as.integer(Sys.getenv(k, d))
step1_sampler <- list(chains = env_int("STEP1_CHAINS", 4), iter = env_int("STEP1_ITER", 2000),
                      warmup = env_int("STEP1_WARMUP", 1000), thin = env_int("STEP1_THIN", 4),
                      threads = env_int("STEP1_THREADS", 1))
# dense_e: a full mass matrix, adapted in warmup - the seasonal terms, GDD and
# EVI are correlated, which a diagonal metric handles poorly
step1_control <- list(adapt_delta = 0.95, max_treedepth = 12, metric = "dense_e")

fit_file <- function(model) file.path(step1_dir, paste0("fit_", model))   # brms adds .rds

fit_step1 <- function(model, data = d_fit, smp = step1_sampler) {
  brm(step1_formula(model), data = droplevels(data), family = gaussian(),
      data2 = if (has_term(model, site_term)) list(K = K_sites) else NULL,
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
