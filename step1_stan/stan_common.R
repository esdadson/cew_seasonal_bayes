# =============================================================================
# step1_stan/stan_common.R
# Shared by step1_stan.R: settings, data, design matrices, Stan data for the
# three models, and the scoring / plotting helpers.
#
# Three nested Gaussian models of y = log(CEW + 1), site i, year t, week w:
#
#   fixed          y = b0 + trap + x' beta + e
#   fourier        fixed + sum_h [c_h cos(2 pi h w/52) + s_h sin(2 pi h w/52)]
#                        + u_t                                   (year intercept)
#   spatiotemporal fourier + g_i + v_i                           (site intercept:
#                                                                 spatial + independent)
#                        + a_i cos(2 pi w/52) + b_i sin(2 pi w/52)  (spatial site slopes:
#                                                                 each site's peak shifts)
#                        + e AR(1) over weeks within each site-year
#
#   x   weekly EVI and GDD, six crop shares (standardised); trap = BL / BL_PH / PH
#   Harmonics h = 1..STAN_HARM (at most 2). g, a, b ~ MVN(0, sigma^2 C), C a
#   Matern 3/2 correlation between the raw site locations with one range.
#
# Zones are not in any model; they are used only to summarise peak weeks.
#
# Environment:
#   STAN_HARM         harmonics in the shared seasonal curve (1 or 2)     2
#   STAN_RANGE0       PC prior on the range: P(range < STAN_RANGE0 km) = 0.1   20
#   STAN_SUBSET_SITES fit a random subset of this many sites (tests only)
#   STEP1_HOLDOUT     held-out year                                       2020
#   STEP1_DATA_CACHE  an .rds of list(d_all, sites) (else step1/step1_data.R)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
env     <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
env_num <- function(k, d) as.numeric(env(k, d))
mins    <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

MODELS  <- c("fixed", "fourier", "spatiotemporal")
yr_test <- as.integer(env("STEP1_HOLDOUT", "2020"))
H       <- as.integer(env_num("STAN_HARM", 2))
range0  <- env_num("STAN_RANGE0", 20)
stopifnot(H %in% 1:2)
zone_lab <- c(`2` = "Southern", `1` = "Transitional", `0` = "Northern")   # zone_30y codes

# =============================================================================
# Data ---------------------------------------------------------------------------
#    d_all (step1/step1_data.R): observed weeks plus every other week 1-52 of
#    each monitored site-year, so full seasonal curves and peak weeks can be
#    predicted. Only observed rows outside the held-out year enter the fit.
# =============================================================================
cache <- Sys.getenv("STEP1_DATA_CACHE")
if (nzchar(cache) && file.exists(cache)) {
  d_all <- readRDS(cache)$d_all
} else {
  source("step1/step1_data.R")
  if (nzchar(cache)) saveRDS(list(d_all = d_all, sites = sites), cache)
}

n_sub <- env_num("STAN_SUBSET_SITES", 0)
if (n_sub > 0) {
  set.seed(42)
  keep  <- sample(unique(d_all$location2), n_sub)
  d_all <- d_all %>% filter(location2 %in% keep)
  message("SUBSET: ", n_sub, " sites (test only)")
}

covs <- c("EVI_filled", "gdd", "corn", "cotton", "soybeans", "tobacco", "peanuts", "other")
d <- d_all %>%
  arrange(location2, year, woy) %>%
  mutate(fit_row  = observed & year != yr_test,
         test_row = observed & year == yr_test,
         year_id  = match(year, sort(unique(year))),
         site_id  = match(location2, sort(unique(location2))),
         zone     = factor(zone_lab[as.character(zone_30y)], zone_lab))
for (h in 1:H) {
  d[[paste0("c", h)]] <- cos(2 * pi * h * d$woy / 52)
  d[[paste0("s", h)]] <- sin(2 * pi * h * d$woy / 52)
}
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))
fit_idx  <- which(d$fit_row)
test_idx <- which(d$test_row)

site_xy <- d %>% distinct(site_id, location2, x_km, y_km, zone_30y, zone) %>% arrange(site_id)
stopifnot(nrow(site_xy) == max(d$site_id))
sy <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "site_id", "year", "zone_30y", "zone")]

# Peak weeks are searched inside the weeks that hold 90% of the observations
pw  <- quantile(d$woy[fit_idx], c(0.05, 0.95), names = FALSE)
win <- c(floor(pw[1]), ceiling(pw[2]))

# =============================================================================
# Design matrices ------------------------------------------------------------------
#    fixed: Intercept, trap type (BL_PH, PH against BL), covariates
#    fourier / spatiotemporal: the same plus c1, s1 (, c2, s2)
# =============================================================================
harm <- c(rbind(paste0("c", 1:H), paste0("s", 1:H)))
design <- function(terms) {
  X <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), d)
  colnames(X) <- sub("^trap_type", "", sub("^\\(Intercept\\)$", "Intercept", colnames(X)))
  X
}
X_list <- list(fixed = design(c("trap_type", covs)))
X_list$fourier <- X_list$spatiotemporal <- design(c(harm, "trap_type", covs))
beta_prior_sd <- function(X)       # intercept N(0, 10^2), harmonics N(0, 2.5^2), others N(0, 1)
  ifelse(colnames(X) == "Intercept", 10, ifelse(grepl("^[cs][0-9]$", colnames(X)), 2.5, 1))

# AR(1) structure of the training rows: each row's previous observed row in
# the same site-year, and the gap between them in weeks
fit_d     <- d[fit_idx, c("loc_year", "woy")]
first     <- c(TRUE, fit_d$loc_year[-1] != fit_d$loc_year[-nrow(fit_d)])
idx_first <- which(first)
idx_next  <- which(!first)
idx_prev  <- idx_next - 1
gap       <- fit_d$woy[idx_next] - fit_d$woy[idx_prev]
stopifnot(all(gap >= 1))

message(sprintf("Training rows %d (%d site-years); held-out %d rows %d; sites %d; peak window weeks %d-%d; %d harmonic(s)",
                length(fit_idx), length(idx_first), yr_test, length(test_idx), nrow(site_xy), win[1], win[2], H))

# =============================================================================
# Stan data -------------------------------------------------------------------------
#    Priors: random-effect sds exponential(3) (P(sd > 1) = 0.05); Matern range
#    PC prior, range ~ inv_gamma(1, -log(0.1) * range0) (P(range < range0) =
#    0.1); phi uniform(-1, 1); residual sd half-Student-t(3, 0, 2.5)
# =============================================================================
stan_data_for <- function(model) {
  X <- X_list[[model]]
  sd_ <- list(
    N = length(fit_idx), P = ncol(X), y = d$log_CEW[fit_idx], X = X[fit_idx, , drop = FALSE],
    beta_sd = beta_prior_sd(X), sigma_scale = 2.5,
    N_all = nrow(d), X_all = X, win_lo = win[1], win_hi = win[2],
    N_test = length(test_idx), test_idx = test_idx)
  if (model != "fixed") sd_ <- c(sd_, list(
    T = max(d$year_id), year = d$year_id[fit_idx], year_all = d$year_id, lambda_sd = -log(0.05)))
  if (model == "spatiotemporal") sd_ <- c(sd_, list(
    S = nrow(site_xy), site = d$site_id[fit_idx], site_all = d$site_id,
    week = d$woy[fit_idx], week_all = d$woy,
    xy = lapply(seq_len(nrow(site_xy)), function(i) c(site_xy$x_km[i], site_xy$y_km[i])),
    jitter = 1e-6, lambda_range = -log(0.1) * range0,
    N_first = length(idx_first), idx_first = idx_first,
    N_next = length(idx_next), idx_next = idx_next, idx_prev = idx_prev, gap = gap))
  sd_
}

inits_for <- function(model) {
  P <- ncol(X_list[[model]]); T <- max(d$year_id); S <- nrow(site_xy)
  base <- list(beta = c(1.5, rep(0, P - 1)), sigma = 1.2)
  function() {
    i <- base
    if (model != "fixed") i <- c(i, list(u_z = rep(0, T), sigma_year = 0.3))
    if (model == "spatiotemporal") i <- c(i, list(
      z_g = rep(0, S), z_a = rep(0, S), z_b = rep(0, S), v_z = rep(0, S),
      sigma_g = 0.5, sigma_a = 0.3, sigma_b = 0.3, sigma_v = 0.4,
      range_km = runif(1, 100, 300), phi = 0.5))
    i
  }
}

hyper_for <- function(model) switch(model,
  fixed          = "sigma",
  fourier        = c("sigma_year", "sigma"),
  spatiotemporal = c("range_km", "sigma_g", "sigma_a", "sigma_b", "sigma_v", "sigma_year",
                     "phi", "sigma", "sigma_innov"))

# =============================================================================
# Scoring and plotting helpers (draws matrices: rows = posterior draws) -----------------
# =============================================================================
qs <- function(v) sprintf("%.2f (%.2f, %.2f)", median(v), quantile(v, 0.025), quantile(v, 0.975))

# Held-out year scored with the marginal predictive N(eta, sigma): for the
# spatiotemporal model sigma is the marginal residual sd and earlier 2020
# observations are not conditioned on, so all three models are scored alike
holdout_metrics <- function(ET, sigma) {
  y  <- d$log_CEW[test_idx]
  LL <- dnorm(matrix(y, nrow(ET), ncol(ET), byrow = TRUE), ET, sigma, log = TRUE)
  mu <- colMeans(ET)
  yrep <- ET + matrix(rnorm(length(ET)), nrow(ET)) * sigma
  q <- apply(yrep, 2, quantile, c(0.05, 0.95))
  data.frame(elpd_test = sum(matrixStats::colLogSumExps(LL) - log(nrow(ET))),
             R2_test   = 1 - sum((y - mu)^2) / sum((y - mean(y))^2),
             rmse_test = sqrt(mean((y - mu)^2)),
             mae_test  = mean(abs(y - mu)),
             cov90_test = mean(y >= q[1, ] & y <= q[2, ]))
}

# Mean peak week of each site (draws x sites), from the site-year peaks
site_peaks <- function(PK) {
  Ws <- Matrix::sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)
  as.matrix(PK %*% (Ws %*% Matrix::Diagonal(x = 1 / Matrix::colSums(Ws))))
}

# Zone mean peak weeks, delays behind Southern, and eta^2 (share of the
# between-site variance in peak week explained by zone)
zone_peak_stats <- function(PK) {
  zm <- sapply(levels(sy$zone), function(z) rowMeans(PK[, sy$zone == z, drop = FALSE]))
  SP <- site_peaks(PK)
  eta2 <- apply(SP, 1, function(x) {
    m <- tapply(x, site_xy$zone, mean)[as.character(site_xy$zone)]
    v <- sum((x - mean(x))^2)
    if (v > 0) sum((m - mean(x))^2) / v else NA })
  data.frame(peak_Southern = qs(zm[, "Southern"]), peak_Transitional = qs(zm[, "Transitional"]),
             peak_Northern = qs(zm[, "Northern"]),
             delay_Transitional = qs(zm[, "Transitional"] - zm[, "Southern"]),
             delay_Northern = qs(zm[, "Northern"] - zm[, "Southern"]),
             eta2 = if (all(is.na(eta2))) NA else qs(na.omit(eta2)))
}

# Shared seasonal curve from the fixed effects (Bz: draws x columns of X)
shared_curve <- function(Bz) {
  wk <- 1:52
  cv <- Bz[, "Intercept"] %o% rep(1, 52)
  for (h in 1:H) cv <- cv + Bz[, paste0("c", h)] %o% cos(2 * pi * h * wk / 52) +
                            Bz[, paste0("s", h)] %o% sin(2 * pi * h * wk / 52)
  data.frame(week = wk, mean = colMeans(cv), lo = apply(cv, 2, quantile, 0.025),
             hi = apply(cv, 2, quantile, 0.975))
}

plot_peak_hist <- function(pk_list, file) {
  pk_long <- bind_rows(lapply(names(pk_list), function(m)
    data.frame(model = m, zone = rep(sy$zone, each = nrow(pk_list[[m]])), week = as.vector(pk_list[[m]]))))
  pk_long$model <- factor(pk_long$model, names(pk_list))
  p <- ggplot(pk_long, aes(week, fill = zone, colour = zone)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                   position = "identity", alpha = 0.35, linewidth = 0.3) +
    facet_wrap(~ model, ncol = 1) +
    scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2") +
    labs(x = "Peak week", y = "Density", fill = NULL, colour = NULL,
         subtitle = "Posterior peak-week draws, all site-years, by zone (zones assigned after fitting)") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file, p, width = 7, height = 2 + 2.5 * length(pk_list), dpi = 300)
}

plot_curves <- function(cv_list, file) {
  cd <- bind_rows(lapply(names(cv_list), function(m) cbind(model = m, cv_list[[m]])))
  p <- ggplot(cd, aes(week, mean, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) + geom_line() +
    annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf, alpha = 0.06) +
    labs(x = "Week of year", y = "Shared curve (log scale, covariates at 0)", colour = NULL, fill = NULL,
         subtitle = sprintf("Shared seasonal curve, %d harmonic(s) (95%% CrI); shaded = peak window", H)) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file, p, width = 6.5, height = 4, dpi = 300)
}

plot_site_map <- function(v, file, title, mid = 0) {
  p <- ggplot(cbind(site_xy, v = v), aes(x_km, y_km, colour = v)) +
    geom_point(size = 1.4) + coord_equal() +
    scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = mid) +
    labs(x = "x (km)", y = "y (km)", colour = NULL, subtitle = title) + theme_bw(base_size = 9)
  ggsave(file, p, width = 7, height = 4.5, dpi = 300)
}
