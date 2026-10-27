# =============================================================================
# cew_granger_functions.R
# Helpers for step 2 of cew_peak_granger.R: regional clustering, regional
# weekly series, and Granger-causality tests from Southern to Transitional/
# Northern regions - in-sample F tests, leave-one-year-out (LOYO) out-of-
# sample gain and a year-permutation null.
# Uses fourier_matrix() from cew_source_sink_functions.R; cluster (for
# silhouettes) ships with R.
# =============================================================================

library(dplyr)

zone_code <- c(`2` = "S", `1` = "T", `0` = "N")

# -----------------------------------------------------------------------------
# Regions: Ward hierarchical clustering of each zone's sites.
# Dissimilarity = space_weight x distance on km coordinates
#               + (1 - space_weight) x distance between the sites' seasonal
#                 curves (curves: sites x weeks matrix, rownames = location2,
#                 e.g. the mean filled weekly log(CEW + 1) over the site's years),
# each scaled to a maximum of 1. Without curves: coordinates only. The
# spatial part keeps regions compact enough to map and to act as sources/
# sinks; the curve part groups sites with similar abundance and timing.
# Among cuts in k_range whose clusters all have >= min_sites sites, k per zone
# = the largest k with mean silhouette >= sil_min (0.5: "reasonable
# structure"; the maximum is usually k = 2, an east/west split too coarse to
# say where influence comes from), else the best silhouette.
# k_by_zone (e.g. c(`2` = 4)) overrides it.
# Regions are numbered west to east: S1, S2, ..., T1, ..., N1, ...
# Returns list(sites = sites + region, silhouette = zone x k table).
# -----------------------------------------------------------------------------
cluster_regions <- function(sites, curves = NULL, space_weight = 0.5, k_range = 2:6,
                            min_sites = 10, sil_min = 0.5, k_by_zone = NULL) {
  if (!is.null(curves)) sites <- sites %>% filter(as.character(location2) %in% rownames(curves))
  res <- lapply(split(sites, sites$zone_30y), function(s) {
    z  <- as.character(s$zone_30y[1])
    D  <- dist(cbind(s$x_km, s$y_km))
    if (!is.null(curves)) {
      Dc <- dist(curves[as.character(s$location2), , drop = FALSE])
      D  <- space_weight * D / max(D) + (1 - space_weight) * Dc / max(Dc)
    }
    hc <- hclust(D, "ward.D2")
    ks <- k_range[k_range < nrow(s)]
    sil <- vapply(ks, function(k) {
      cl <- cutree(hc, k)
      if (min(table(cl)) < min_sites) NA_real_ else mean(cluster::silhouette(cl, D)[, "sil_width"])
    }, numeric(1))
    good <- which(sil >= sil_min)
    k <- if (z %in% names(k_by_zone)) k_by_zone[[z]] else if (length(good)) ks[max(good)] else
      if (all(is.na(sil))) 1L else ks[which.max(sil)]
    cl <- cutree(hc, k)
    rank_we <- rank(tapply(s$x_km, cl, mean))           # west -> east
    s$region <- paste0(zone_code[[z]], rank_we[as.character(cl)])
    list(sites = s, sil = data.frame(zone_30y = s$zone_30y[1], k = ks, silhouette = sil,
                                     chosen = ks == k))
  })
  list(sites = bind_rows(lapply(res, `[[`, "sites")),
       silhouette = bind_rows(lapply(res, `[[`, "sil")))
}

# Region centroids (mean site coordinates), site counts and zone
region_info <- function(regions) {
  regions %>% group_by(region, zone_30y) %>%
    summarise(n_sites = n(), longitude = mean(longitude), latitude = mean(latitude),
              x_km = mean(x_km), y_km = mean(y_km), .groups = "drop")
}

# -----------------------------------------------------------------------------
# Regional weekly means as arrays region x year x week (1-52), NA = no site of
# the region sampled that week. d needs region, year, woy and the value and
# control columns. n = number of sites averaged per cell. Weeks outside
# `weeks` are NA.
# resp (optional): column of observed values (NA where not trapped). Y_resp =
# mean over the region's trapped sites (the Granger response) and n_obs =
# their number; Y (e.g. observed + step-1 fills) supplies the lags. Without
# resp, Y_resp = Y and n_obs = n.
# -----------------------------------------------------------------------------
region_series <- function(d, value = "log_CEW", controls = character(), years = sort(unique(d$year)),
                          weeks = 1:52, resp = NULL) {
  regs <- sort(unique(d$region))
  d <- d %>% filter(year %in% years, woy %in% weeks)
  to_array <- function(v, f) {
    a <- array(NA_real_, c(length(regs), length(years), 52), dimnames = list(regs, years, 1:52))
    s <- d %>% group_by(region, year, woy) %>% summarise(m = f(.data[[v]]), .groups = "drop")
    a[cbind(match(s$region, regs), match(s$year, years), s$woy)] <- s$m
    a
  }
  avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  cnt <- function(x) sum(!is.na(x))
  Y <- to_array(value, avg); n <- to_array(value, cnt)
  list(Y = Y, n = n,
       Y_resp = if (is.null(resp)) Y else to_array(resp, avg),
       n_obs  = if (is.null(resp)) n else to_array(resp, cnt),
       X = setNames(lapply(controls, to_array, f = avg), controls),
       years = years, regions = regs)
}

# -----------------------------------------------------------------------------
# Shift a years x 52 matrix k weeks within each year: k > 0 gives the value k
# weeks earlier (lag), k < 0 the value -k weeks later (lead); NA across the
# year boundary (no Dec -> Jan carry-over between seasons).
# -----------------------------------------------------------------------------
lag_cols <- function(M, k) {
  out <- matrix(NA_real_, nrow(M), ncol(M))
  nc <- ncol(M)
  if (k >= 0) out[, (k + 1):nc] <- M[, 1:(nc - k)] else out[, 1:(nc + k)] <- M[, (1 - k):nc]
  out
}

# -----------------------------------------------------------------------------
# Seasonal climatology of a years x 52 matrix from the rows (years) in `rows`:
# intercept + H Fourier pairs fit to the observed cells. Weeks observed in
# fewer than min_years of those years are NA (the curve would extrapolate).
# Returns a 52-vector.
# -----------------------------------------------------------------------------
seasonal_clim <- function(Y, rows, H = 3, min_years = 3) {
  Yr <- Y[rows, , drop = FALSE]
  y  <- as.vector(Yr)                                  # column-major: years vary fastest
  w  <- rep(1:52, each = nrow(Yr))
  ok <- !is.na(y)
  B  <- cbind(1, fourier_matrix(1:52, H))
  cf <- qr.coef(qr(B[w[ok], , drop = FALSE]), y[ok])
  cf[is.na(cf)] <- 0
  clim <- drop(B %*% cf)
  clim[colSums(!is.na(Yr)) < min_years] <- NA
  clim
}

# -----------------------------------------------------------------------------
# Regression rows for one sink/source pair: response = sink anomaly at week
# w; own lags, source lags (negative = leads) and the sink's controls.
# A_resp: the response series (default A_sink; e.g. the mean over trapped
# sites only, while A_sink - observed + filled - supplies the own lags).
# resp_ok (years x 52 logical): weeks allowed as responses (e.g. the sink's
# season window); the others still serve as lags. Complete cases only, so
# the restricted and full models share rows.
# -----------------------------------------------------------------------------
granger_rows <- function(A_sink, A_src, X, own, src, resp_ok = NULL, A_resp = A_sink) {
  nm_src <- ifelse(src >= 0, paste0("src_lag", src), paste0("src_lead", -src))
  y <- A_resp
  if (!is.null(resp_ok)) y[!resp_ok] <- NA
  cols <- c(list(y = y),
            setNames(lapply(own, function(k) lag_cols(A_sink, k)), paste0("own_lag", own)),
            setNames(lapply(src, function(k) lag_cols(A_src, k)), nm_src),
            X)
  M  <- vapply(cols, as.vector, numeric(length(A_sink)))
  yr <- rep(seq_len(nrow(A_sink)), times = ncol(A_sink))
  wk <- rep(seq_len(ncol(A_sink)), each = nrow(A_sink))
  ok <- complete.cases(M)
  list(M = M[ok, , drop = FALSE], year = yr[ok], woy = wk[ok], src = nm_src)
}

ols_fit <- function(X, y) { cf <- qr.coef(qr(X), y); cf[is.na(cf)] <- 0; cf }

# -----------------------------------------------------------------------------
# Granger comparison of one sink and one source series for several lag
# settings (list of list(own = own lags, src = source lags), named).
#   Restricted R: sink anomaly ~ 1 + own lags + controls
#   Full       U: R + source lags
# Anomaly = series - seasonal_clim(). In-sample: climatology and fits on all
# years -> F test of the source lags. LOYO: for each year, climatologies and
# fits use the other years only, and that year is predicted; gain =
# 1 - SSE_U / SSE_R summed over held-out years (> 0: Southern lags improve
# out-of-sample prediction). years_better = held-out years where U beat R.
#   Ys, Yk: sink / source years x 52; X: list of sink years x 52 controls
#   Ys_resp: sink response series (default Ys), e.g. the mean of trapped sites
#     only; it shares the sink's climatology
#   resp_ok: sink weeks allowed as responses (see granger_rows)
#   Ys_oos, Yk_oos: optional years x 52 series whose row t replaces row t of
#     Ys / Yk when year t is held out (step-1 fills from the fit without year t)
# Returns list(table = one row per setting, pred = held-out predictions).
# -----------------------------------------------------------------------------
granger_pair <- function(Ys, Yk, X, settings, Ys_resp = Ys, resp_ok = NULL, Ys_oos = NULL,
                         Yk_oos = NULL, H = 3, min_years = 3, keep_pred = FALSE) {
  ny <- nrow(Ys); all_y <- seq_len(ny)
  anoms <- function(rows, S = Ys, K = Yk) {        # climatology from `rows` only
    cs <- seasonal_clim(S, rows, H, min_years); ck <- seasonal_clim(K, rows, H, min_years)
    list(s = sweep(S, 2, cs), k = sweep(K, 2, ck), r = sweep(Ys_resp, 2, cs), clim_s = cs)
  }
  # in-sample F tests
  a <- anoms(all_y)
  ins <- lapply(settings, function(st) {
    r  <- granger_rows(a$s, a$k, X, st$own, st$src, resp_ok, a$r)
    y  <- r$M[, "y"]; XU <- cbind(1, r$M[, -1, drop = FALSE]); XR <- XU[, !colnames(XU) %in% r$src, drop = FALSE]
    n  <- length(y); q <- length(r$src)
    sse <- function(Xm) sum((y - Xm %*% ols_fit(Xm, y))^2)
    sR <- sse(XR); sU <- sse(XU); df2 <- n - ncol(XU)
    Fst <- if (df2 > 0) ((sR - sU) / q) / (sU / df2) else NA
    cfU <- if (n > ncol(XU)) ols_fit(XU, y)[r$src] else NA
    data.frame(n = n, F = Fst, df1 = q, df2 = df2, p_F = pf(Fst, q, df2, lower.tail = FALSE),
               src_coef_sum = sum(cfU))
  })
  # leave-one-year-out
  sse <- matrix(0, length(settings), 2 * ny, dimnames = list(names(settings), NULL))
  ntest <- matrix(0L, length(settings), ny)
  preds <- list()
  for (t in all_y) {
    S <- Ys; K <- Yk
    if (!is.null(Ys_oos)) S[t, ] <- Ys_oos[t, ]
    if (!is.null(Yk_oos)) K[t, ] <- Yk_oos[t, ]
    a <- anoms(all_y[-t], S, K)
    for (s in seq_along(settings)) {
      st <- settings[[s]]
      r  <- granger_rows(a$s, a$k, X, st$own, st$src, resp_ok, a$r)
      XU <- cbind(1, r$M[, -1, drop = FALSE]); XR <- XU[, !colnames(XU) %in% r$src, drop = FALSE]
      y  <- r$M[, "y"]; tr <- r$year != t; te <- r$year == t
      if (!any(te) || sum(tr) <= ncol(XU)) next
      pR <- XR[te, , drop = FALSE] %*% ols_fit(XR[tr, , drop = FALSE], y[tr])
      pU <- XU[te, , drop = FALSE] %*% ols_fit(XU[tr, , drop = FALSE], y[tr])
      sse[s, t] <- sum((y[te] - pR)^2); sse[s, ny + t] <- sum((y[te] - pU)^2)
      ntest[s, t] <- sum(te)
      if (keep_pred)
        preds[[length(preds) + 1]] <- data.frame(setting = names(settings)[s], year_i = t, woy = r$woy[te],
                                                 clim = a$clim_s[r$woy[te]], anom = y[te],
                                                 pred_R = drop(pR), pred_U = drop(pU))
    }
  }
  sR <- sse[, 1:ny, drop = FALSE]; sU <- sse[, ny + 1:ny, drop = FALSE]
  tab <- as_tibble(bind_rows(ins, .id = "setting")) %>%
    mutate(n_test = rowSums(ntest), loyo_gain = 1 - rowSums(sU) / rowSums(sR),
           years_better = rowSums(sU < sR & ntest > 0), years_tested = rowSums(ntest > 0))
  list(table = tab, pred = bind_rows(preds))
}

# -----------------------------------------------------------------------------
# Year-permutation null: the source's seasons are reassigned to other years
# (a derangement), keeping each season's within-year shape but breaking the
# link to the sink's year; the whole LOYO procedure is rerun. The source's
# held-out-year series (Yk_oos) is permuted with it. Returns an n_perm x
# settings matrix of null LOYO gains.
# -----------------------------------------------------------------------------
derange <- function(n) { repeat { p <- sample.int(n); if (all(p != seq_len(n))) return(p) } }

granger_perm <- function(Ys, Yk, X, settings, n_perm = 199, Yk_oos = NULL, ...) {
  t(vapply(seq_len(n_perm), function(i) {
    p <- derange(nrow(Yk))
    granger_pair(Ys, Yk[p, , drop = FALSE], X, settings,
                 Yk_oos = if (is.null(Yk_oos)) NULL else Yk_oos[p, , drop = FALSE], ...)$table$loyo_gain
  }, numeric(length(settings))))
}

# Permutation p-values: per setting, and for the best of the settings in
# max_cols (max over them), which accounts for choosing the best lag
perm_p <- function(obs, null, max_cols = rep(TRUE, length(obs))) {
  per <- (1 + colSums(sweep(null, 2, obs, ">="), na.rm = TRUE)) / (1 + nrow(null))
  mx  <- (1 + sum(apply(null[, max_cols, drop = FALSE], 1, max, na.rm = TRUE) >=
                    max(obs[max_cols], na.rm = TRUE))) / (1 + nrow(null))
  list(per_setting = per, max = mx)
}

# -----------------------------------------------------------------------------
# Peak-timing delay between regions from the step-1 peak draws.
#   P: draws x site-years peak weeks; sy: site_years with region and year.
# Region-year peak = mean over the region's site-years, per draw. A pair's
# delay = sink - source, averaged over the years both were trapped.
# -----------------------------------------------------------------------------
region_year_peaks <- function(P, sy) {
  g <- paste(sy$region, sy$year, sep = "|")
  vapply(split(seq_len(ncol(P)), g), function(j) rowMeans(P[, j, drop = FALSE]), numeric(nrow(P)))
}

peak_delay <- function(RP, source, sink) {
  key <- do.call(rbind, strsplit(colnames(RP), "|", fixed = TRUE))
  ys  <- intersect(key[key[, 1] == source, 2], key[key[, 1] == sink, 2])
  if (!length(ys)) return(data.frame(delay_med = NA, delay_lo = NA, delay_hi = NA, p_later = NA, n_years = 0))
  dl <- rowMeans(RP[, paste(sink, ys, sep = "|"), drop = FALSE] - RP[, paste(source, ys, sep = "|"), drop = FALSE])
  data.frame(delay_med = median(dl), delay_lo = unname(quantile(dl, 0.025)),
             delay_hi = unname(quantile(dl, 0.975)), p_later = mean(dl > 0), n_years = length(ys))
}
