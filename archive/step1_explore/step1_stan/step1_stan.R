# =============================================================================
# step1_stan/step1_stan.R
# Stan mirror of the R-INLA zone-free space-time model
# (step1_inla/inla_spacetime_free.R, model "spacetime"), with the spatial
# range and sd ESTIMATED and the AR1 correlation across knots kept:
#
#   y = X beta + s(w) + u_t + v_i + g(c_i, w) + e      (spacetime_free.stan)
#
# What matches INLA: covariates, cyclic RW2 seasonal curve (scaled), year and
# site random intercepts, Matern field per temporal knot with AR1 across
# knots, the same quadratic B-spline temporal basis, the same PC priors, the
# same held-out year, peak-week window and peak definition.
# What differs:
#   - the field is an exact Matern GP over field locations = sites merged
#     within STAN_CUTOFF km (INLA: SPDE on a triangular mesh with the same
#     kind of cutoff). Sites closer than a few metres made the site-level
#     covariance singular; 15 km gives 239 locations, well conditioned.
#   - full MCMC over every hyperparameter (INLA's spacetime fit used "eb",
#     hyperparameters fixed at their mode), so intervals can be wider.
#   - the seasonal curve's sum-to-zero constraint is soft (sd 0.001 per week)
#   - residual sd: half-Student-t(3, 0, 2.5) (INLA: its vague default)
#
# Modes (like cew_stage1_hpc.R):
#   Rscript step1_stan/step1_stan.R dryrun    # compile + one short chain -> time per
#                                             # gradient, projected run time, resources
#   Rscript step1_stan/step1_stan.R fit       # full fit -> STAN_OUT/draws.rds
#   Rscript step1_stan/step1_stan.R summary   # tables and figures from draws.rds,
#                                             # compared with INLA when its results exist
# On Hazel through Slurm only (hpc/step1_stan.sbatch); on Windows / Mac any mode.
#
# Environment:
#   STAN_CHAINS, STAN_ITER, STAN_WARMUP, STAN_THIN   4, 2000, 1000, 2 (iter includes warmup)
#   STAN_ADAPT_DELTA, STAN_MAX_TREEDEPTH             0.9, 10
#   STAN_NU             1 (Matern nu = 1, as INLA) or 1.5 (built-in, faster)
#   STAN_CENTER_SITE    1 = centered site effects (default: many rows per site)
#   STAN_CENTER_FIELD   0 = non-centered field (default)
#   STAN_CUTOFF         km; sites closer share a field location        15
#   STAN_KNOT_WEEKS     weeks between temporal knots                   4
#   STAN_RANGE0         P(range < STAN_RANGE0 km) = 0.1                20
#   STAN_SUBSET_SITES   fit a random subset of this many sites (quick tests)
#   STAN_OUT            output folder                                  results/step1_stan
#   STEP1_HOLDOUT       held-out year                                  2020
#   STEP1_DATA_CACHE    an .rds of list(d_all, sites) (else step1/step1_data.R)
#   DRY_WARMUP, DRY_SAMPLING                         dry-run chain length (150 + 50)
#   INLA_RESULTS        INLA results to compare with   results/inla_spacetime_free
# =============================================================================

if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with sbatch hpc/step1_stan.sbatch - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"
stopifnot(mode %in% c("dryrun", "fit", "summary"))

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
env     <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
env_num <- function(k, d) as.numeric(env(k, d))
mins    <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

out_dir <- env("STAN_OUT", "results/step1_stan")
csv_dir <- file.path(out_dir, "stan_csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
yr_test <- as.integer(env("STEP1_HOLDOUT", "2020"))
smp <- list(chains = env_num("STAN_CHAINS", 4), iter = env_num("STAN_ITER", 2000),
            warmup = env_num("STAN_WARMUP", 1000), thin = env_num("STAN_THIN", 2),
            adapt_delta = env_num("STAN_ADAPT_DELTA", 0.9), max_td = env_num("STAN_MAX_TREEDEPTH", 10))
zone_lab <- c(`0` = "Northern limits", `1` = "Transitional zone", `2` = "Southern range")

# =============================================================================
# 1. Data (as step1_inla/inla_spacetime_free.R) -----------------------------------
# =============================================================================
cache <- Sys.getenv("STEP1_DATA_CACHE")
if (nzchar(cache) && file.exists(cache)) {
  d_all <- readRDS(cache)$d_all
} else {
  source("step1/step1_data.R")
  if (nzchar(cache)) saveRDS(list(d_all = d_all, sites = sites), cache)
}

n_sub <- env_num("STAN_SUBSET_SITES", 0)
if (n_sub > 0) {                                  # quick tests only
  set.seed(42)
  keep <- sample(unique(d_all$location2), n_sub)
  d_all <- d_all %>% filter(location2 %in% keep)
  message("SUBSET: ", n_sub, " sites (test only)")
}

covs <- c("EVI_filled", "gdd", "corn", "cotton", "soybeans", "tobacco", "peanuts", "other")
d <- d_all %>%
  arrange(location2, year, woy) %>%
  mutate(fit_row  = observed & year != yr_test,
         test_row = observed & year == yr_test,
         year_id  = match(year, sort(unique(year))),
         site_id  = match(location2, sort(unique(location2))))
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))
fit_idx  <- which(d$fit_row)
test_idx <- which(d$test_row)

trap_X <- model.matrix(~ trap_type, d)[, -1, drop = FALSE]
colnames(trap_X) <- sub("^trap_type", "trap_", colnames(trap_X))
X <- cbind(intercept = 1, trap_X, as.matrix(d[, covs]))

site_xy <- d %>% distinct(site_id, location2, longitude, latitude, x_km, y_km, zone_30y) %>% arrange(site_id)
sy <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "site_id", "year", "zone_30y")]

# Field locations: sites merged within the cutoff (complete linkage, so every
# cluster is at most cutoff km across), at the cluster's mean position
cutoff  <- env_num("STAN_CUTOFF", 15)
cluster <- cutree(hclust(dist(site_xy[, c("x_km", "y_km")]), "complete"), h = cutoff)
xy_loc  <- as.matrix(aggregate(site_xy[, c("x_km", "y_km")], list(cluster), mean)[, -1])

# Temporal basis: quadratic B-splines with knots every STAN_KNOT_WEEKS weeks,
# identical to INLA's fm_mesh_1d(knots, degree = 2, boundary = "free")
knot_wk <- env_num("STAN_KNOT_WEEKS", 4)
knots   <- unique(c(seq(1, 52, by = knot_wk), 52))
h1 <- diff(knots)[1]; h2 <- tail(diff(knots), 1)
B_wk <- splines::splineDesign(c(knots[1] - 2 * h1, knots[1] - h1, knots, 52 + h2, 52 + 2 * h2),
                              1:52, ord = 3, outer.ok = TRUE)

peak_win <- quantile(d$woy[fit_idx], c(0.05, 0.95), names = FALSE)
win <- c(floor(peak_win[1]), ceiling(peak_win[2]))

message(sprintf("Training rows %d; held-out %d rows %d; sites %d; field locations %d (cutoff %g km); knots %d; peak window %d-%d",
                length(fit_idx), yr_test, length(test_idx), nrow(site_xy), nrow(xy_loc), cutoff,
                ncol(B_wk), win[1], win[2]))

# =============================================================================
# 2. Priors (INLA's PC priors) ---------------------------------------------------
#    pc.prec(u, a) on a precision  ==  sd ~ exponential(-log(a) / u)
#    pcmatern range: P(range < r0) = a  ==  1 / range ~ exponential(-log(a) r0)
#    pc.cor1(u, a): P(rho > u) = a, exponential prior on sqrt(1 - rho) over
#    (-1, 1); theta solved here as INLA does
#    scale.model: the RW2 precision is multiplied by the geometric mean of
#    the marginal variances of the unscaled structure (generalised inverse)
# =============================================================================
D2 <- matrix(0, 52, 52)                                  # cyclic second differences
for (w in 1:52) D2[w, c((w - 2) %% 52 + 1, w, w %% 52 + 1)] <- c(1, -2, 1)
R_rw2 <- crossprod(D2)
eg <- eigen(R_rw2, symmetric = TRUE); pos <- eg$values > 1e-8
rw2_scale <- exp(mean(log(diag(eg$vectors[, pos] %*% diag(1 / eg$values[pos]) %*% t(eg$vectors[, pos])))))

theta_rho <- uniroot(function(th) (1 - exp(-th * sqrt(1 - 0.5))) / (1 - exp(-th * sqrt(2))) - 0.9,
                     c(1e-6, 100))$root

stan_data <- list(
  N = length(fit_idx), P = ncol(X), T = max(d$year_id), S = nrow(site_xy), M = nrow(xy_loc), K = ncol(B_wk),
  y = d$log_CEW[fit_idx], X = X[fit_idx, ], week = d$woy[fit_idx], year = d$year_id[fit_idx],
  site = d$site_id[fit_idx], cluster = as.integer(cluster), xy = lapply(seq_len(nrow(xy_loc)), function(i) xy_loc[i, ]),
  B = B_wk,
  nu_code = if (env_num("STAN_NU", 1) == 1) 1L else 2L,
  centered_site = as.integer(env_num("STAN_CENTER_SITE", 1)),
  centered_field = as.integer(env_num("STAN_CENTER_FIELD", 0)),
  jitter = 1e-6,
  rw2_scale = rw2_scale,
  beta_sd = ifelse(colnames(X) == "intercept", 10, 1),
  lambda_season = -log(0.05), lambda_year = -log(0.01), lambda_site = -log(0.01),
  lambda_field_sd = -log(0.01), lambda_range = -log(0.1) * env_num("STAN_RANGE0", 20),
  theta_rho = theta_rho, sigma_scale = 2.5,
  N_all = nrow(d), X_all = X, week_all = d$woy, year_all = d$year_id, site_all = d$site_id,
  win_lo = win[1], win_hi = win[2],
  N_test = length(test_idx), test_idx = test_idx, y_test = d$log_CEW[test_idx])

# Start every chain from moderate values (a GP sampler started at extreme
# ranges can spend its warmup getting back)
inits <- function() list(
  beta = c(1, rep(0, ncol(X) - 1)), s = rep(0, 52), sigma_season = 1,
  u_z = rep(0, stan_data$T), sigma_year = 0.3,
  v_par = rep(0, stan_data$S), sigma_site = 0.5,
  G_par = matrix(0, stan_data$M, stan_data$K),
  range_km = runif(1, 150, 300), sigma_field = 0.5, rho = 0.5, sigma = 1.2)

hyper <- c("range_km", "sigma_field", "rho", "sigma_season", "sigma_year", "sigma_site", "sigma")

# =============================================================================
# 3. Dry run / fit -------------------------------------------------------------------
# =============================================================================
if (mode %in% c("dryrun", "fit")) {
  library(cmdstanr)
  message("CmdStan ", cmdstan_version(), " | nu = ", c("1", "3/2")[stan_data$nu_code],
          " | centered site ", stan_data$centered_site, ", field ", stan_data$centered_field)
  t0  <- Sys.time()
  mod <- cmdstan_model("step1_stan/spacetime_free.stan")
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
  h_max <- smp$iter * (2^smp$max_td - 1) * sec_grad / 3600
  h_req <- 1.15 * min(h_max, 3 * h_obs) + 1
  fmt_t <- function(h) { m <- ceiling(min(h, 240) * 60)
    sprintf("%d-%02d:%02d:00", m %/% 1440, (m %% 1440) %/% 60, m %% 60) }

  message("\n--- Dry run (1 chain, ", w, " warmup + ", s, " sampling iterations) ---")
  message(sprintf("Compile %.1f min. Warmup %.1f s/iteration, sampling %.1f s/iteration",
                  t_compile, tm$warmup / w, tm$sampling / s))
  message(sprintf("Seconds per gradient: %.4f. Leapfrog steps per sampling iteration: %.0f (max %d)",
                  sec_grad, steps, 2^smp$max_td - 1))
  message(sprintf("Sampling iterations at max tree depth: %.0f%%; divergent: %d",
                  100 * mean(dg$treedepth__[!in_w] >= smp$max_td), sum(dg$divergent__[!in_w])))
  print(fit$summary(hyper), digits = 3)
  message(sprintf("\n--- Full run: %d chains x %d iterations (%d warmup) ---", smp$chains, smp$iter, smp$warmup))
  message(sprintf("Per chain at this tree depth: %.1f h; if every iteration hits max depth: %.1f h", h_obs, h_max))
  message("Suggested: sbatch --time=", fmt_t(h_req), if (h_req > 96) " --qos=long", " --cpus-per-task=", smp$chains,
          " --mem=16G hpc/step1_stan.sbatch   (dry-run estimates are rough: short adaptation)")
}

if (mode == "fit") {
  message(sprintf("Fitting %d chains x %d iterations (%d warmup, thin %d)", smp$chains, smp$iter, smp$warmup, smp$thin))
  t0  <- Sys.time()
  fit <- mod$sample(data = stan_data, init = inits, chains = smp$chains, parallel_chains = smp$chains,
                    iter_warmup = smp$warmup, iter_sampling = smp$iter - smp$warmup, thin = smp$thin,
                    adapt_delta = smp$adapt_delta, max_treedepth = smp$max_td,
                    seed = 42, refresh = 50, output_dir = csv_dir)
  message(sprintf("Sampling: %.1f h", mins(t0) / 60))
  # keep what the summary needs (the CSVs in stan_csv/ hold everything)
  saveRDS(list(draws = fit$draws(c(hyper, "beta", "s", "u", "v", "G", "peak", "eta_test", "log_lik_test"),
                                 format = "draws_matrix"),
               diag = fit$diagnostic_summary(), time = fit$time(), smp = smp,
               setup = list(stan_data = stan_data[c("nu_code", "centered_site", "centered_field", "M", "K")],
                            cutoff = cutoff, knots = knots, win = win, xcols = colnames(X))),
          file.path(out_dir, "draws.rds"))
  print(fit$summary(hyper), digits = 3)
  message("Saved ", file.path(out_dir, "draws.rds"))
}

# =============================================================================
# 4. Summary (local or HPC) ------------------------------------------------------------
# =============================================================================
if (mode == "summary") {
  res <- readRDS(file.path(out_dir, "draws.rds"))
  dr  <- matrix(res$draws, nrow(res$draws), dimnames = dimnames(res$draws))   # plain matrix: draws x variables
  n_dr <- nrow(dr)
  col  <- function(p) dr[, grep(paste0("^", p, "\\["), colnames(dr)), drop = FALSE]

  # Hyperparameters with convergence diagnostics
  hyp <- posterior::summarise_draws(posterior::subset_draws(res$draws, variable = hyper))
  print(hyp, digits = 3)
  write.csv(hyp, file.path(out_dir, "hyper.csv"), row.names = FALSE)

  # Held-out year
  y_test <- d$log_CEW[test_idx]
  LL <- col("log_lik_test"); ET <- col("eta_test")
  mu <- colMeans(ET)
  sig <- dr[, "sigma"]
  yrep <- ET + matrix(rnorm(length(ET)), n_dr) * sig
  q <- apply(yrep, 2, quantile, c(0.05, 0.95))
  ho <- data.frame(model = "stan_spacetime",
                   elpd_test = sum(matrixStats::colLogSumExps(LL) - log(n_dr)),
                   R2_test = 1 - sum((y_test - mu)^2) / sum((y_test - mean(y_test))^2),
                   mae_test = mean(abs(y_test - mu)),
                   cov90_test = mean(y_test >= q[1, ] & y_test <= q[2, ]))

  # Peak weeks: draws x site-years, zones assigned after fitting
  PK <- col("peak")
  stopifnot(ncol(PK) == nrow(sy))
  zmean <- sapply(c(`2` = 2, `1` = 1, `0` = 0), function(z) rowMeans(PK[, sy$zone_30y == z, drop = FALSE]))
  Ws <- Matrix::sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)
  Ws <- Ws %*% Matrix::Diagonal(x = 1 / Matrix::colSums(Ws))
  SP <- as.matrix(PK %*% Ws)
  eta2 <- apply(SP, 1, function(x) {
    zm <- tapply(x, site_xy$zone_30y, mean)[as.character(site_xy$zone_30y)]
    sum((zm - mean(x))^2) / sum((x - mean(x))^2) })
  qs <- function(v) sprintf("%.2f (%.2f, %.2f)", median(v), quantile(v, 0.025), quantile(v, 0.975))
  zone_stats <- data.frame(model = "stan_spacetime",
                           peak_Southern = qs(zmean[, "2"]), peak_Transitional = qs(zmean[, "1"]),
                           peak_Northern = qs(zmean[, "0"]),
                           delay_Transitional = qs(zmean[, "1"] - zmean[, "2"]),
                           delay_Northern = qs(zmean[, "0"] - zmean[, "2"]), eta2 = qs(eta2))

  # INLA results for comparison, if inla_spacetime_free.R has been run with the same data
  inla_dir <- env("INLA_RESULTS", "results/inla_spacetime_free")
  pk_long <- data.frame(model = "Stan", zone_30y = rep(sy$zone_30y, each = n_dr), week = as.vector(PK))
  if (file.exists(file.path(inla_dir, "comparison.csv")) && n_sub == 0) {
    inla_cmp <- read.csv(file.path(inla_dir, "comparison.csv")) %>% filter(model == "spacetime")
    ho <- bind_rows(ho, inla_cmp %>% transmute(model = "inla_spacetime", elpd_test, R2_test, mae_test, cov90_test))
    zone_stats <- bind_rows(zone_stats, read.csv(file.path(inla_dir, "zone_stats.csv")) %>%
                              filter(model == "spacetime") %>% mutate(model = "inla_spacetime"))
    ip <- readRDS(file.path(inla_dir, "posterior_draws.rds"))
    if (identical(ip$sy$loc_year, sy$loc_year))
      pk_long <- bind_rows(pk_long, data.frame(model = "INLA", zone_30y = rep(sy$zone_30y, each = nrow(ip$post$spacetime$PK)),
                                               week = as.vector(ip$post$spacetime$PK)))
    fi <- file.path(inla_dir, "fit_spacetime.rds")
    if (file.exists(fi)) {
      h <- readRDS(fi)$summary.hyperpar
      message("\nINLA spacetime hyperparameters (precisions as sd):")
      prec <- grepl("^Precision", rownames(h))
      h[prec, c("0.5quant", "0.025quant", "0.975quant")] <- 1 / sqrt(h[prec, c("0.5quant", "0.975quant", "0.025quant")])
      print(round(h[, c("0.5quant", "0.025quant", "0.975quant")], 3))
    }
  }
  print(ho, digits = 4); print(zone_stats)
  write.csv(ho, file.path(out_dir, "holdout.csv"), row.names = FALSE)
  write.csv(zone_stats, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)

  pk_long <- pk_long %>% mutate(zone = factor(zone_lab[as.character(zone_30y)], zone_lab))
  p <- ggplot(pk_long, aes(week, fill = zone, colour = zone)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                   position = "identity", alpha = 0.35, linewidth = 0.3) +
    facet_wrap(~ model, ncol = 1) +
    scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2") +
    labs(x = "Peak week", y = "Density", fill = NULL, colour = NULL,
         subtitle = "Posterior peak-week draws, all site-years, by zone (zones assigned after fitting)") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "peak_hist.png"), p, width = 7, height = 2 + 2.5 * n_distinct(pk_long$model), dpi = 300)

  # Seasonal curve
  S_dr <- col("s")
  sc <- data.frame(week = 1:52, mean = colMeans(S_dr), lo = apply(S_dr, 2, quantile, 0.025),
                   hi = apply(S_dr, 2, quantile, 0.975))
  p <- ggplot(sc, aes(week, mean)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25) + geom_line() +
    labs(x = "Week of year", y = "s(w)", subtitle = "Stan: shared seasonal curve (95% CrI)") + theme_bw(base_size = 9)
  ggsave(file.path(out_dir, "season_curve.png"), p, width = 6, height = 3.5, dpi = 300)
  message("Summary written to ", out_dir)
}
