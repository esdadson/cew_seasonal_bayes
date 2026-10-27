# =============================================================================
# cew_source_sink_functions.R
# Helper functions for the H. zea seasonal-curve (stage 1) and source-sink
# (stage 2) workflow. Source this file from cew_source_sink_pipeline.R.
# Only base R + dplyr/tidyr + geosphere + splines are needed here.
# =============================================================================

library(dplyr)
library(tidyr)
library(geosphere)
library(splines)

# -----------------------------------------------------------------------------
# Data root (contains Original/ and Raster/): $CEW_DATA if set, else the first
# that exists of Hazel research storage, the Mac mount, the Windows drive.
# -----------------------------------------------------------------------------
default_data_root <- function() {
  env <- Sys.getenv("CEW_DATA")
  if (nzchar(env)) return(env)
  cand <- c("/rs1/researchers/c/cmjone25/Data", "/Volumes/cmjone25/Data", "Z:/Data")
  hit <- cand[dir.exists(cand)]
  if (!length(hit)) stop("Data root not found; set CEW_DATA. Tried: ", paste(cand, collapse = ", "))
  hit[1]
}

# -----------------------------------------------------------------------------
# ONE site table: location2 by sorted coordinates, plus CONUS Albers km.
# Shared by the pipeline and cew_covariates_hpc.R so every join uses the same ids.
# -----------------------------------------------------------------------------
make_site_table <- function(counts_raw) {
  sites <- counts_raw %>%
    group_by(longitude, latitude) %>%
    summarise(location = first(location), zone_30y = first(zone_30y), .groups = "drop") %>%
    arrange(longitude, latitude) %>%
    mutate(location2 = row_number())
  # Projected coordinates (km) for the GP - degrees of longitude shrink with
  # latitude, so scaled lon/lat is not isotropic
  # (closed-form, not terra::project: on Hazel's conda env PROJ returned NA
  # for every site, which left the GP with no usable rows)
  xy <- albers_km(sites$longitude, sites$latitude)
  if (anyNA(xy)) stop(sum(!complete.cases(xy)), " of ", nrow(xy), " sites could not be projected to EPSG:5070")
  sites$x_km <- xy[, 1]; sites$y_km <- xy[, 2]
  sites
}

# CONUS Albers equal-area conic (EPSG:5070: GRS80, standard parallels 29.5 and
# 45.5, origin 23N 96W), km. Snyder (1987) eqs. 14-3 to 14-15; agrees with
# PROJ to ~1 m (lon/lat taken as NAD83).
albers_km <- function(lon, lat) {
  a <- 6378137; f <- 1 / 298.257222101; e2 <- 2 * f - f^2; e <- sqrt(e2); rad <- pi / 180
  q <- function(p) { s <- sin(p); (1 - e2) * (s / (1 - e2 * s^2) - log((1 - e * s) / (1 + e * s)) / (2 * e)) }
  m <- function(p) cos(p) / sqrt(1 - e2 * sin(p)^2)
  p0 <- 23 * rad; p1 <- 29.5 * rad; p2 <- 45.5 * rad; l0 <- -96 * rad
  n   <- (m(p1)^2 - m(p2)^2) / (q(p2) - q(p1)); C <- m(p1)^2 + n * q(p1)
  rho <- function(p) a * sqrt(C - n * q(p)) / n
  th  <- n * (lon * rad - l0)
  cbind(rho(lat * rad) * sin(th), rho(p0) - rho(lat * rad) * cos(th)) / 1000
}

# -----------------------------------------------------------------------------
# Standardise with means/sds from reference rows (the training years), so
# held-out and prediction rows are put on the scale the model was fit on.
#   map: named character vector, new column = source column
# Returns a function df -> df that (re)computes the new columns.
# -----------------------------------------------------------------------------
make_scaler <- function(ref, map) {
  m <- vapply(map, function(v) mean(ref[[v]], na.rm = TRUE), numeric(1))
  s <- vapply(map, function(v) sd(ref[[v]], na.rm = TRUE), numeric(1))
  function(df) {
    for (k in names(map)) df[[k]] <- (df[[map[[k]]]] - m[[k]]) / s[[k]]
    df
  }
}

# -----------------------------------------------------------------------------
# Fourier basis on ISO week (period = 52 weeks)
# -----------------------------------------------------------------------------
add_fourier <- function(df, H = 3, week_col = "woy", period = 52) {
  w <- df[[week_col]]
  for (h in seq_len(H)) {
    df[[paste0("s", h)]] <- sin(2 * pi * h * w / period)
    df[[paste0("c", h)]] <- cos(2 * pi * h * w / period)
  }
  df
}

fourier_matrix <- function(weeks, H = 2, period = 52) {
  B <- do.call(cbind, lapply(seq_len(H), function(h)
    cbind(sin(2 * pi * h * weeks / period), cos(2 * pi * h * weeks / period))))
  colnames(B) <- as.vector(rbind(paste0("s", 1:H), paste0("c", 1:H)))
  B
}

# -----------------------------------------------------------------------------
# Peak week from Fourier coefficient draws
#   coef_draws: draws x sites x 2H array (columns ordered s1,c1,s2,c2,...)
#   window:     weeks to search (restrict to the sampled season!)
# Returns list(peak = draws x sites matrix, at_edge = share of draws at window edge)
# -----------------------------------------------------------------------------
peak_from_draws <- function(coef_draws, window = 17:39, step = 0.1) {
  H <- dim(coef_draws)[3] / 2
  wk <- seq(min(window), max(window), by = step)
  B <- fourier_matrix(wk, H)
  peak <- t(vapply(seq_len(dim(coef_draws)[1]), function(d) {
    f <- B %*% t(coef_draws[d, , , drop = TRUE])        # weeks x sites
    wk[max.col(t(f), ties.method = "first")]
  }, numeric(dim(coef_draws)[2])))
  colnames(peak) <- dimnames(coef_draws)[[2]]
  at_edge <- colMeans(peak <= min(wk) + step | peak >= max(wk) - step)
  list(peak = peak, at_edge = at_edge)
}

# -----------------------------------------------------------------------------
# Year-week index: columns 1..(n_years*52), contiguous within year
# -----------------------------------------------------------------------------
make_tw_index <- function(years) {
  yrs <- sort(unique(years))
  idx <- expand.grid(woy = 1:52, year = yrs)
  idx$tw <- seq_len(nrow(idx))
  idx[, c("year", "woy", "tw")]
}

# -----------------------------------------------------------------------------
# Site x year-week matrix of anomalies (NA where the trap was not sampled)
# -----------------------------------------------------------------------------
anomaly_matrix <- function(df, site_ids, tw_index, value_col = "anom") {
  M <- matrix(NA_real_, length(site_ids), nrow(tw_index),
              dimnames = list(site_ids, tw_index$tw))
  df <- df %>% inner_join(tw_index, by = c("year", "woy"))
  r <- match(df$location2, site_ids)
  ok <- !is.na(r)
  M[cbind(r[ok], df$tw[ok])] <- df[[value_col]][ok]
  M
}

# -----------------------------------------------------------------------------
# Shift columns by k weeks WITHIN a year (k > 0 = lag, k < 0 = lead).
# Column c of the output holds column c-k of the input; NA across year bounds.
# -----------------------------------------------------------------------------
shift_within_year <- function(M, k, col_woy, fill = NA_real_) {
  out <- matrix(fill, nrow(M), ncol(M), dimnames = dimnames(M))
  valid <- which(col_woy - k >= 1 & col_woy - k <= 52)
  out[, valid] <- M[, valid - k]
  out
}

# -----------------------------------------------------------------------------
# Distances (km) and bearings (degrees, source -> sink) between site sets
# -----------------------------------------------------------------------------
pair_geometry <- function(sinks, sources) {
  # sinks / sources: data frames with location2, longitude, latitude
  D <- distm(cbind(sinks$longitude, sinks$latitude),
             cbind(sources$longitude, sources$latitude)) / 1000
  g <- expand.grid(i = seq_len(nrow(sinks)), j = seq_len(nrow(sources)))
  b <- bearing(cbind(sources$longitude[g$j], sources$latitude[g$j]),
               cbind(sinks$longitude[g$i],  sinks$latitude[g$i]))
  Bg <- matrix(b, nrow(sinks), nrow(sources))
  dimnames(D) <- dimnames(Bg) <- list(sinks$location2, sources$location2)
  list(D = D, bearing = Bg)
}

# -----------------------------------------------------------------------------
# Kernel weights: exponential distance decay x optional von Mises direction
#   lambda = Inf gives a flat (spatially non-specific) kernel -> placebo
#   theta0 = preferred source->sink bearing (deg), kappa = concentration
# -----------------------------------------------------------------------------
kernel_weights <- function(geom, lambda, kappa = 0, theta0 = 0) {
  W <- if (is.infinite(lambda)) matrix(1, nrow(geom$D), ncol(geom$D)) else exp(-geom$D / lambda)
  if (kappa > 0) W <- W * exp(kappa * (cos((geom$bearing - theta0) * pi / 180) - 1))
  dimnames(W) <- dimnames(geom$D)
  W
}

# -----------------------------------------------------------------------------
# Source pressure: kernel-weighted mean of SAMPLED source anomalies.
#   A: sources x tw anomaly matrix (NA = not sampled)
#   lag_mat (optional): sinks x sources integer lag (weeks) for the
#     travelling-wave version (lag = round(distance / speed)); pairs with
#     lag > max_lag get zero weight.
#   min_den: minimum summed kernel weight required; below it P is NA so a
#     sink isn't driven by one far-away trap when nothing nearby was sampled.
# Returns sinks x tw matrix.
# -----------------------------------------------------------------------------
source_pressure <- function(W, A, col_woy, lag_mat = NULL, max_lag = 8, min_den = 0.05) {
  O  <- !is.na(A) * 1
  A0 <- A; A0[is.na(A0)] <- 0
  if (is.null(lag_mat)) {
    num <- W %*% A0
    den <- W %*% O
  } else {
    num <- den <- matrix(0, nrow(W), ncol(A))
    for (l in 0:max_lag) {
      Wl <- W * (lag_mat == l)
      if (!any(Wl > 0)) next
      num <- num + Wl %*% shift_within_year(A0, l, col_woy, fill = 0)
      den <- den + Wl %*% shift_within_year(O,  l, col_woy, fill = 0)
    }
  }
  # relative threshold so it works for any lambda
  thr <- min_den * apply(W, 1, max)
  P <- num / den
  P[sweep(den, 1, thr, "<")] <- NA
  dimnames(P) <- list(rownames(W), colnames(A))
  P
}

# -----------------------------------------------------------------------------
# Pull lagged pressure values onto sink observations
#   obs: data frame with location2, year, woy; returns obs with P_lag{k} cols
# -----------------------------------------------------------------------------
attach_lags <- function(obs, P, tw_index, lags = 0:8, prefix = "P_lag") {
  col_woy <- tw_index$woy
  obs <- obs %>% left_join(tw_index, by = c("year", "woy"))
  r <- match(obs$location2, rownames(P))
  for (k in lags) {
    Pk <- shift_within_year(P, k, col_woy)
    nm <- if (k >= 0) paste0(prefix, k) else paste0("P_lead", -k)
    obs[[nm]] <- Pk[cbind(r, obs$tw)]
  }
  obs
}

# -----------------------------------------------------------------------------
# Constrained distributed lag: spline basis over lags.
#   Q_m = sum_k B[k, m] * P_lag_k   -> fit Q_m;  beta_k = B %*% theta
# -----------------------------------------------------------------------------
lag_basis <- function(lags = 0:8, df = 4) {
  B <- ns(lags, df = df, intercept = TRUE)
  B <- unclass(B)[, , drop = FALSE]
  attr(B, "knots") <- NULL; attr(B, "Boundary.knots") <- NULL
  rownames(B) <- lags
  colnames(B) <- paste0("Q", seq_len(ncol(B)))
  B
}

add_Q <- function(obs, B, prefix = "P_lag") {
  Pm <- as.matrix(obs[, paste0(prefix, rownames(B))])
  Q <- Pm %*% B
  obs[, colnames(B)] <- Q
  obs
}

# -----------------------------------------------------------------------------
# Fast screening: leave-one-year-out CV gain of adding pressure terms
#   gain = 1 - SSE(with P) / SSE(without P), on identical rows
# -----------------------------------------------------------------------------
loyo_gain <- function(dat, p_terms, base_terms, response = "anom") {
  vars <- unique(c(response, p_terms, base_terms, "year"))
  d <- dat[complete.cases(dat[, vars]), vars]
  f0 <- reformulate(base_terms, response)
  f1 <- reformulate(c(base_terms, p_terms), response)
  sse0 <- sse1 <- 0
  for (y in unique(d$year)) {
    tr <- d[d$year != y, ]; te <- d[d$year == y, ]
    sse0 <- sse0 + sum((te[[response]] - predict(lm(f0, tr), te))^2)
    sse1 <- sse1 + sum((te[[response]] - predict(lm(f1, tr), te))^2)
  }
  c(gain = 1 - sse1 / sse0, n = nrow(d))
}

# -----------------------------------------------------------------------------
# Year-permutation null: shuffle whole source years (keeps within-season
# structure of the source field, breaks the year link to sinks)
# -----------------------------------------------------------------------------
permute_source_years <- function(A, tw_index) {
  yrs <- unique(tw_index$year)
  perm <- sample(yrs)
  while (any(perm == yrs) && length(yrs) > 1) perm <- sample(yrs)  # derangement
  out <- A
  for (k in seq_along(yrs)) {
    to   <- tw_index$tw[tw_index$year == yrs[k]]
    from <- tw_index$tw[tw_index$year == perm[k]]
    out[, to] <- A[, from]
  }
  out
}

# -----------------------------------------------------------------------------
# Compare kernel / speed settings on the SAME rows (otherwise sharper kernels
# drop more rows via min_den and gains aren't comparable).
#   P_list: named list of sinks x tw pressure matrices
# -----------------------------------------------------------------------------
grid_screen <- function(obs, P_list, tw_index, lags = 0:8, df = 4,
                        base_terms = "zone_fac", response = "anom") {
  B <- lag_basis(lags, df)
  built <- lapply(P_list, function(P) add_Q(attach_lags(obs, P, tw_index, lags), B))
  keep <- Reduce(`&`, lapply(built, function(o) complete.cases(o[, colnames(B)])))
  out <- t(sapply(built, function(o)
    loyo_gain(o[keep, ], colnames(B), base_terms, response)))
  data.frame(setting = names(P_list), out, row.names = NULL)
}

# -----------------------------------------------------------------------------
# Out-of-sample scores
#   lpd_pointwise: log predictive density per observation from log_lik()
#     draws (draws x obs) computed on the held-out data
#   oos_scores: y, predictive draws (draws x obs), pointwise lpd ->
#     RMSE/MAE/bias of the predictive mean, central 90% interval coverage, ELPD
#   score_by: oos_scores within groups plus "all"
# -----------------------------------------------------------------------------
lpd_pointwise <- function(ll) apply(ll, 2, function(l) max(l) + log(mean(exp(l - max(l)))))

oos_scores <- function(y, pred, lpd) {
  mu <- colMeans(pred)
  q  <- apply(pred, 2, quantile, c(0.05, 0.95))
  data.frame(n = length(y), rmse = sqrt(mean((y - mu)^2)), mae = mean(abs(y - mu)),
             bias = mean(mu - y), cover90 = mean(y >= q[1, ] & y <= q[2, ]),
             elpd = sum(lpd), elpd_per_obs = mean(lpd))
}

score_by <- function(y, pred, lpd, group) {
  idx <- c(split(seq_along(y), group), list(all = seq_along(y)))
  bind_rows(lapply(idx, function(i) oos_scores(y[i], pred[, i, drop = FALSE], lpd[i])), .id = "group")
}
