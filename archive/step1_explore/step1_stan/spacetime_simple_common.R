# =============================================================================
# step1_stan/spacetime_simple_common.R
# Shared by spacetime_simple.R (Stan) and spacetime_simple_brms.R (brms):
# settings, data, design matrix, priors and summary helpers, so both fits use
# the same data and are scored the same way.
#
# The model: the zone-free space-time model rebuilt from standard parts
#
#   y_n = x_n' theta                              intercept, shared seasonal curve
#                                                 (harmonics), trap type, covariates
#       + u_t                                     year random intercept
#       + g_i + v_i                               site random intercept:
#                                                 spatial (Matern) + independent
#       + a_i cos(2 pi w/52) + b_i sin(2 pi w/52) site random slopes on the first
#                                                 harmonic, spatially correlated
#                                                 (the spatial-seasonal deviation)
#       + e_n,  e AR(1) over weeks within each site-year
#
# Zones are not in the model; they are used only to summarise peak weeks.
#
# Environment (both scripts):
#   ST_HARM          harmonics in the shared seasonal curve                3
#   ST_SLOPES        1 = site random slopes on the first harmonic          1
#   ST_AR            1 = AR(1) residual within site-year                   1
#   ST_RANGE0        PC prior on the range: P(range < ST_RANGE0 km) = 0.1   20
#   ST_SUBSET_SITES  fit a random subset of this many sites (tests)
#   STEP1_HOLDOUT    held-out year                                          2020
#   STEP1_DATA_CACHE an .rds of list(d_all, sites) (else step1/step1_data.R)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
env     <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
env_num <- function(k, d) as.numeric(env(k, d))
mins    <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

yr_test <- as.integer(env("STEP1_HOLDOUT", "2020"))
H       <- as.integer(env_num("ST_HARM", 3))
slopes  <- as.integer(env_num("ST_SLOPES", 1))
use_ar  <- as.integer(env_num("ST_AR", 1))
range0  <- env_num("ST_RANGE0", 20)
zone_lab <- c(`2` = "Southern", `1` = "Transitional", `0` = "Northern")   # zone_30y codes

# =============================================================================
# Data --------------------------------------------------------------------------
# =============================================================================
cache <- Sys.getenv("STEP1_DATA_CACHE")
if (nzchar(cache) && file.exists(cache)) {
  d_all <- readRDS(cache)$d_all
} else {
  source("step1/step1_data.R")
  if (nzchar(cache)) saveRDS(list(d_all = d_all, sites = sites), cache)
}

n_sub <- env_num("ST_SUBSET_SITES", 0)
if (n_sub > 0) {                                  # quick tests only
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
         zone     = factor(zone_lab[as.character(zone_30y)], zone_lab),
         year_f   = factor(year), site_f = factor(location2))
for (h in 1:H) {                                  # harmonics of week of year
  d[[paste0("c", h)]] <- cos(2 * pi * h * d$woy / 52)
  d[[paste0("s", h)]] <- sin(2 * pi * h * d$woy / 52)
}
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))
fit_idx  <- which(d$fit_row)
test_idx <- which(d$test_row)

site_xy <- d %>% distinct(site_id, location2, x_km, y_km, zone_30y, zone) %>% arrange(site_id)
stopifnot(nrow(site_xy) == max(d$site_id))
sy <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "site_id", "year", "zone_30y", "zone")]

peak_win <- quantile(d$woy[fit_idx], c(0.05, 0.95), names = FALSE)
win <- c(floor(peak_win[1]), ceiling(peak_win[2]))

# Fixed effects: intercept, shared seasonal curve (H harmonics), trap type, covariates
harm_terms <- paste(c(rbind(paste0("c", 1:H), paste0("s", 1:H))), collapse = " + ")
f_fixed <- as.formula(paste("~", harm_terms, "+ trap_type +", paste(covs, collapse = " + ")))
X <- model.matrix(f_fixed, d)
colnames(X) <- sub("^trap_type", "", sub("^\\(Intercept\\)$", "Intercept", colnames(X)))
is_harm <- grepl("^[cs][0-9]+$", colnames(X))
beta_sd <- ifelse(colnames(X) == "Intercept", 10, ifelse(is_harm, 2.5, 1))

# AR(1) structure of the training rows: each row's previous observed row in
# the same site-year and the gap between them in weeks
fit_d   <- d[fit_idx, c("loc_year", "woy")]
first   <- c(TRUE, fit_d$loc_year[-1] != fit_d$loc_year[-nrow(fit_d)])
idx_first <- which(first)
idx_next  <- which(!first)
idx_prev  <- idx_next - 1
gap       <- fit_d$woy[idx_next] - fit_d$woy[idx_prev]
stopifnot(all(gap >= 1))

# PC prior on the Matern range in 2-d: 1 / range ~ exponential(lambda_range),
# P(range < range0) = 0.1  <=>  range ~ inv_gamma(1, lambda_range)
lambda_range <- -log(0.1) * range0

message(sprintf("Training rows %d (%d site-years, %.0f%% of steps 1 week); held-out %d rows %d; sites %d; peak window %d-%d",
                length(fit_idx), length(idx_first), 100 * mean(gap == 1), yr_test, length(test_idx),
                nrow(site_xy), win[1], win[2]))
message(sprintf("Shared curve %d harmonic(s); site slopes: %s; AR(1): %s; P(range < %g km) = 0.1",
                H, if (slopes) "yes" else "no", if (use_ar) "yes" else "no", range0))

# =============================================================================
# Summary helpers (draws matrices: rows = posterior draws) ---------------------------
# =============================================================================
qs <- function(v) sprintf("%.2f (%.2f, %.2f)", median(v), quantile(v, 0.025), quantile(v, 0.975))

# Held-out year scored with the marginal predictive N(eta, sigma_marginal):
# no conditioning on earlier 2020 observations, so it is comparable with
# models without AR(1)
holdout_metrics <- function(model, ET, sigma) {
  y  <- d$log_CEW[test_idx]
  LL <- dnorm(matrix(y, nrow(ET), ncol(ET), byrow = TRUE), ET, sigma, log = TRUE)
  mu <- colMeans(ET)
  yrep <- ET + matrix(rnorm(length(ET)), nrow(ET)) * sigma
  q <- apply(yrep, 2, quantile, c(0.05, 0.95))
  data.frame(model = model,
             elpd_test = sum(matrixStats::colLogSumExps(LL) - log(nrow(ET))),
             R2_test   = 1 - sum((y - mu)^2) / sum((y - mean(y))^2),
             mae_test  = mean(abs(y - mu)),
             cov90_test = mean(y >= q[1, ] & y <= q[2, ]))
}

# Peak weeks (PK: draws x site-years) -> zone means, delays behind Southern,
# eta^2 (share of between-site peak variance due to zone). Zones assigned
# after fitting.
zone_peak_stats <- function(model, PK) {
  zm <- sapply(levels(sy$zone), function(z) rowMeans(PK[, sy$zone == z, drop = FALSE]))
  Ws <- Matrix::sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)
  Ws <- Ws %*% Matrix::Diagonal(x = 1 / Matrix::colSums(Ws))
  SP <- as.matrix(PK %*% Ws)
  eta2 <- apply(SP, 1, function(x) {
    m <- tapply(x, site_xy$zone, mean)[as.character(site_xy$zone)]
    sum((m - mean(x))^2) / sum((x - mean(x))^2) })
  data.frame(model = model,
             peak_Southern = qs(zm[, "Southern"]), peak_Transitional = qs(zm[, "Transitional"]),
             peak_Northern = qs(zm[, "Northern"]),
             delay_Transitional = qs(zm[, "Transitional"] - zm[, "Southern"]),
             delay_Northern = qs(zm[, "Northern"] - zm[, "Southern"]), eta2 = qs(eta2))
}

# Mean peak week of each site (draws x sites), from the site-year peaks
site_peaks <- function(PK) {
  Ws <- Matrix::sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)
  Ws <- Ws %*% Matrix::Diagonal(x = 1 / Matrix::colSums(Ws))
  as.matrix(PK %*% Ws)
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
         subtitle = sprintf("Shared seasonal curve, %d harmonics (95%% CrI); shaded = peak window", H)) +
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
