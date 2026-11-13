# =============================================================================
# step1_stan/zone_harmonic_common.R
# Shared by zone_harmonic.R (Stan) and zone_harmonic_brms.R (brms):
# settings, data, design matrices, priors and the summary helpers, so the two
# fits use identical data and are scored identically.
#
# The model (both versions):
#
#   y = X beta                                    trap type + covariates
#     + a_z + sum_h [b_zh cos(2 pi h w/52)         zone-specific seasonal curve
#                  + c_zh sin(2 pi h w/52)]         (fixed effects, h = 1..H)
#     + u_t                                        year intercept
#     + w_i                                        site effect: short-range Matern
#                                                  field + site nugget, raw sites
#     + alpha_i cos(2 pi w/52) + beta_i sin(2 pi w/52)   (ZH_SVC = 1) site-level
#                                                  phase / amplitude shift of the
#                                                  first harmonic: short-range fields
#     + e
#   y = log(CEW + 1), w = week of year, i = site, t = year, z = zone of site i
#
# Environment (both scripts):
#   ZH_HARM          harmonics in the zone curves (H)                    2
#   ZH_SVC           1 = add the site-level first-harmonic fields         1
#   ZH_RANGE_LO/HI   98% prior interval of the Matern range (km)         10, 150
#                    (lognormal; range = distance where correlation ~0.14)
#   ZH_SUBSET_SITES  fit a random subset of this many sites (tests)
#   STEP1_HOLDOUT    held-out year                                       2020
#   STEP1_DATA_CACHE an .rds of list(d_all, sites) (else step1/step1_data.R)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
env     <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
env_num <- function(k, d) as.numeric(env(k, d))
mins    <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

yr_test  <- as.integer(env("STEP1_HOLDOUT", "2020"))
H        <- as.integer(env_num("ZH_HARM", 2))
svc      <- as.integer(env_num("ZH_SVC", 1))
range_lo <- env_num("ZH_RANGE_LO", 10)
range_hi <- env_num("ZH_RANGE_HI", 150)
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

n_sub <- env_num("ZH_SUBSET_SITES", 0)
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
for (h in 1:max(H, 1)) {                          # harmonics of week of year
  d[[paste0("c", h)]] <- cos(2 * pi * h * d$woy / 52)
  d[[paste0("s", h)]] <- sin(2 * pi * h * d$woy / 52)
}
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))
fit_idx  <- which(d$fit_row)
test_idx <- which(d$test_row)

# One row per site (all sites, including any seen only in the held-out year:
# their effects are predicted from neighbours) and per site-year
site_xy <- d %>% distinct(site_id, location2, x_km, y_km, zone_30y, zone) %>% arrange(site_id)
stopifnot(nrow(site_xy) == max(d$site_id))
sy <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "site_id", "year", "zone_30y", "zone")]

# Peak-week search window: 5th-95th percentile of observed training weeks
peak_win <- quantile(d$woy[fit_idx], c(0.05, 0.95), names = FALSE)
win <- c(floor(peak_win[1]), ceiling(peak_win[2]))

# Fixed effects: zone intercepts, zone x harmonics, trap type, covariates.
# No global intercept (each zone has its own).
harm_terms <- paste(c(rbind(paste0("c", 1:H), paste0("s", 1:H))), collapse = " + ")
f_fixed <- as.formula(paste("~ 0 + zone + zone:(", harm_terms, ") + trap_type +", paste(covs, collapse = " + ")))
X <- model.matrix(f_fixed, d)
colnames(X) <- gsub("zone|trap_type", "", colnames(X))          # "Southern", "Southern:c1", "BL_PH", ...
is_zone <- colnames(X) %in% levels(d$zone)
is_harm <- grepl(":", colnames(X))
beta_sd <- ifelse(is_zone, 10, ifelse(is_harm, 2.5, 1))

# Matern range prior: lognormal with 98% of its mass in [range_lo, range_hi]
range_mu <- log(sqrt(range_lo * range_hi))
range_sd <- log(range_hi / range_lo) / (2 * qnorm(0.99))

message(sprintf("Training rows %d; held-out %d rows %d; sites %d; site-years %d; peak window %d-%d",
                length(fit_idx), yr_test, length(test_idx), nrow(site_xy), nrow(sy), win[1], win[2]))
message(sprintf("Zone curves with %d harmonic(s); site-level harmonic fields: %s; range prior 98%% in [%g, %g] km",
                H, if (svc) "yes" else "no", range_lo, range_hi))

# =============================================================================
# Summary helpers (draws matrices: rows = posterior draws) ---------------------------
# =============================================================================
qs <- function(v) sprintf("%.2f (%.2f, %.2f)", median(v), quantile(v, 0.025), quantile(v, 0.975))

# Held-out year: elpd from pointwise log densities, R2 / MAE of the posterior
# mean, coverage of 90% predictive intervals
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

# Peak weeks per site-year (PK: draws x site-years) -> zone means, delays
# behind Southern, and eta^2 = share of between-site peak variance due to zone
zone_peak_stats <- function(model, PK) {
  zm <- sapply(levels(sy$zone), function(z) rowMeans(PK[, sy$zone == z, drop = FALSE]))
  Ws <- Matrix::sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)
  Ws <- Ws %*% Matrix::Diagonal(x = 1 / Matrix::colSums(Ws))
  SP <- as.matrix(PK %*% Ws)                                   # draws x sites
  eta2 <- apply(SP, 1, function(x) {
    m <- tapply(x, site_xy$zone, mean)[as.character(site_xy$zone)]
    sum((m - mean(x))^2) / sum((x - mean(x))^2) })
  data.frame(model = model,
             peak_Southern = qs(zm[, "Southern"]), peak_Transitional = qs(zm[, "Transitional"]),
             peak_Northern = qs(zm[, "Northern"]),
             delay_Transitional = qs(zm[, "Transitional"] - zm[, "Southern"]),
             delay_Northern = qs(zm[, "Northern"] - zm[, "Southern"]), eta2 = qs(eta2))
}

# Zone harmonic curves from the fixed effects (Bz: draws x columns of X).
# Curve = zone intercept + harmonics; its peak week is the zone's seasonal
# peak with covariates (gdd, EVI, crops) held at zero.
zone_curves <- function(Bz) {
  wk <- 1:52
  out <- lapply(levels(d$zone), function(z) {
    cv <- Bz[, z] %o% rep(1, 52)
    for (h in 1:H) {
      cv <- cv + Bz[, paste0(z, ":c", h)] %o% cos(2 * pi * h * wk / 52) +
                 Bz[, paste0(z, ":s", h)] %o% sin(2 * pi * h * wk / 52)
    }
    pk <- win[1] - 1 + apply(cv[, win[1]:win[2], drop = FALSE], 1, which.max)
    list(curve = data.frame(zone = z, week = wk, mean = colMeans(cv),
                            lo = apply(cv, 2, quantile, 0.025), hi = apply(cv, 2, quantile, 0.975)),
         peak = pk)
  })
  names(out) <- levels(d$zone)
  out
}

# Peak-week histogram by zone, one panel per model
plot_peak_hist <- function(pk_list, file) {
  pk_long <- bind_rows(lapply(names(pk_list), function(m)
    data.frame(model = m, zone = rep(sy$zone, each = nrow(pk_list[[m]])), week = as.vector(pk_list[[m]]))))
  p <- ggplot(pk_long, aes(week, fill = zone, colour = zone)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                   position = "identity", alpha = 0.35, linewidth = 0.3) +
    facet_wrap(~ model, ncol = 1) +
    scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2") +
    labs(x = "Peak week", y = "Density", fill = NULL, colour = NULL,
         subtitle = "Posterior peak-week draws, all site-years, by zone") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file, p, width = 7, height = 2 + 2.5 * length(pk_list), dpi = 300)
}

plot_zone_curves <- function(zc, file, title) {
  cd <- bind_rows(lapply(zc, `[[`, "curve")) %>% mutate(zone = factor(zone, levels(d$zone)))
  p <- ggplot(cd, aes(week, mean, colour = zone, fill = zone)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) + geom_line() +
    annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf, alpha = 0.05) +
    scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2") +
    labs(x = "Week of year", y = "Zone curve (log scale, covariates at 0)", colour = NULL, fill = NULL,
         subtitle = title) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file, p, width = 6.5, height = 4, dpi = 300)
}

# Site effects on a map (vector over site_xy rows)
plot_site_map <- function(v, file, title) {
  p <- ggplot(cbind(site_xy, v = v), aes(x_km, y_km, colour = v)) +
    geom_point(size = 1.4) + coord_equal() +
    scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b") +
    labs(x = "x (km)", y = "y (km)", colour = NULL, subtitle = title) + theme_bw(base_size = 9)
  ggsave(file, p, width = 7, height = 4.5, dpi = 300)
}
