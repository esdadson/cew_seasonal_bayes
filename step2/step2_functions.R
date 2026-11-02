# =============================================================================
# step2/step2_functions.R
# Helpers for step2/step2_granger.R: regional clustering, regional weekly
# series, and Granger-causality tests from Southern to Transitional/
# Northern regions - fit on the training years (in-sample F test and a
# year-permutation null), validated on the held-out year (as in step 1).
# Conditional tests (one series given another) compare the two DAGs,
# South -> North directly vs South -> Transitional -> North. cluster (for
# silhouettes) ships with R.
# =============================================================================

library(dplyr)

zone_code <- c(`2` = "S", `1` = "T", `0` = "N")

# Fourier basis matrix of week (period = 52), columns s1, c1, s2, c2, ...
# (as in cew_source_sink_functions.R)
fourier_matrix <- function(weeks, H = 2, period = 52) {
  B <- do.call(cbind, lapply(seq_len(H), function(h)
    cbind(sin(2 * pi * h * weeks / period), cos(2 * pi * h * weeks / period))))
  colnames(B) <- as.vector(rbind(paste0("s", 1:H), paste0("c", 1:H)))
  B
}

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
# Regression rows for one sink: response = sink anomaly at week w; the sink's
# own lags; the lags of each predictor series (lags: named list, series name
# -> lags, negative = leads; A_pred: named list of anomaly matrices); the
# sink's controls. A_resp: the response series (default A_sink; e.g. the mean
# over trapped sites only, while A_sink - observed + filled - supplies the own
# lags). resp_ok (years x 52 logical): weeks allowed as responses (e.g. the
# sink's season window); the others still serve as lags. Complete cases only,
# so nested models share rows.
# -----------------------------------------------------------------------------
lag_names <- function(s, k) paste0(s, ifelse(k >= 0, "_lag", "_lead"), abs(k))

granger_design <- function(A_sink, A_pred, lags, own, X, resp_ok = NULL, A_resp = A_sink) {
  y <- A_resp
  if (!is.null(resp_ok)) y[!resp_ok] <- NA
  pcols <- do.call(c, c(list(list()), lapply(names(lags), function(s)
    setNames(lapply(lags[[s]], function(k) lag_cols(A_pred[[s]], k)), lag_names(s, lags[[s]])))))
  cols <- c(list(y = y),
            setNames(lapply(own, function(k) lag_cols(A_sink, k)), paste0("own_lag", own)),
            pcols, X)
  M  <- vapply(cols, as.vector, numeric(length(A_sink)))
  yr <- rep(seq_len(nrow(A_sink)), times = ncol(A_sink))
  wk <- rep(seq_len(ncol(A_sink)), each = nrow(A_sink))
  ok <- complete.cases(M)
  list(M = M[ok, , drop = FALSE], year = yr[ok], woy = wk[ok])
}

ols_fit <- function(X, y) { cf <- qr.coef(qr(X), y); cf[is.na(cf)] <- 0; cf }

# A derangement of 1..n (no element stays in place)
derange <- function(n) { repeat { p <- sample.int(n); if (all(p != seq_len(n))) return(p) } }
# n_perm derangements of the training years, shared by the settings of one
# scan so that a max-over-settings p-value is valid
make_perms <- function(n_train, n_perm) t(replicate(n_perm, derange(n_train)))

# -----------------------------------------------------------------------------
# Granger test for one sink: do the lags in `test` improve on a model with
# the sink's own lags, its controls and the lags in `base`?
#   R: a_sink(w) = c + own lags + controls + base lags + e
#   U: R + test lags
# a_ = series minus its seasonal climatology (Fourier, H pairs), fit on the
# training years only. Regressions are fit on the training years (`train`,
# row indices of the years x 52 matrices); the held-out years (`holdout`)
# are predicted from that fit.
#   Ys: sink years x 52 (observed + step-1 fills: the own lags); preds: named
#   list of predictor series (years x 52); X: list of the sink's controls;
#   Ys_resp: sink response (trapped sites only; shares the sink's climatology);
#   resp_ok: sink weeks allowed as responses.
# Returns list(table = one row, null = permutation F's, pred = held-out rows):
#   F, p_F       in-sample F test of the test lags (nominal: serially
#                correlated errors make it optimistic, hence p_perm)
#   gain         1 - SSE_U / SSE_R in the training years: the share of the
#                sink's remaining variation the test lags explain (partial R^2)
#   coef_sum     sum of the test lags' coefficients (sign of the association)
#   p_perm       year-permutation p of F: the test series' seasons are
#                reassigned among the training years (perms rows), keeping
#                each season's shape but breaking its link to the sink's year
#   gain_holdout 1 - SSE_U / SSE_R in the held-out years (> 0: out of sample,
#                the test lags improve the prediction)
# -----------------------------------------------------------------------------
granger_test <- function(Ys, preds, X, own = 1:2, test, base = list(), Ys_resp = Ys, resp_ok = NULL,
                         train, holdout = integer(), H = 3, min_years = 3, perms = NULL, keep_pred = FALSE) {
  cs  <- seasonal_clim(Ys, train, H, min_years)
  A_s <- sweep(Ys, 2, cs); A_r <- sweep(Ys_resp, 2, cs)
  A_p <- lapply(preds[unique(c(names(base), names(test)))],
                function(P) sweep(P, 2, seasonal_clim(P, train, H, min_years)))
  tcols <- unlist(lapply(names(test), function(s) lag_names(s, test[[s]])))
  lags  <- c(base, test)
  if (any(duplicated(names(lags)))) stop("granger_test: a series is in both base and test")
  fit1 <- function(Ap, full = TRUE) {
    r  <- granger_design(A_s, Ap, lags, own, X, resp_ok, A_r)
    y  <- r$M[, "y"]; XU <- cbind(1, r$M[, -1, drop = FALSE])
    XR <- XU[, !colnames(XU) %in% tcols, drop = FALSE]
    tr <- r$year %in% train
    bU <- ols_fit(XU[tr, , drop = FALSE], y[tr]); bR <- ols_fit(XR[tr, , drop = FALSE], y[tr])
    sU <- sum((y[tr] - XU[tr, , drop = FALSE] %*% bU)^2); sR <- sum((y[tr] - XR[tr, , drop = FALSE] %*% bR)^2)
    n <- sum(tr); q <- length(tcols); df2 <- n - ncol(XU)
    Fst <- if (df2 > 0 && sU > 0) ((sR - sU) / q) / (sU / df2) else NA_real_
    if (!full) return(Fst)
    te <- r$year %in% holdout
    pR <- drop(XR[te, , drop = FALSE] %*% bR); pU <- drop(XU[te, , drop = FALSE] %*% bU)
    list(tab = data.frame(n = n, df1 = q, df2 = df2, F = Fst, p_F = pf(Fst, q, df2, lower.tail = FALSE),
                          gain = 1 - sU / sR, coef_sum = sum(bU[tcols]), n_holdout = sum(te),
                          gain_holdout = if (any(te)) 1 - sum((y[te] - pU)^2) / sum((y[te] - pR)^2) else NA_real_),
         pred = if (keep_pred && any(te)) data.frame(year_i = r$year[te], woy = r$woy[te], clim = cs[r$woy[te]],
                                                     anom = y[te], pred_R = pR, pred_U = pU) else NULL)
  }
  out  <- fit1(A_p)
  null <- NULL
  if (!is.null(perms)) {
    null <- apply(perms, 1, function(p) {
      Ap <- A_p
      for (s in names(test)) Ap[[s]][train, ] <- A_p[[s]][train[p], ]
      fit1(Ap, full = FALSE)
    })
    out$tab$p_perm <- (1 + sum(null >= out$tab$F, na.rm = TRUE)) / (1 + length(null))
  }
  list(table = out$tab, null = null, pred = out$pred)
}

# Permutation p of the best of several settings scanned with the same perms
# (null: n_perm x settings F's): accounts for having picked the best one
perm_max_p <- function(F_obs, null)
  (1 + sum(apply(null, 1, max, na.rm = TRUE) >= max(F_obs, na.rm = TRUE))) / (1 + nrow(null))

# -----------------------------------------------------------------------------
# Lag scan of one source series for one sink: each lag (negative = lead) on
# its own, with the sink's own lags, all sharing one set of permutations.
# Returns the per-lag table plus p_perm_max (the best lag against the null of
# the best lag) and the held-out predictions (keep_pred).
# -----------------------------------------------------------------------------
lag_scan <- function(Ys, Yk, X, lags, own = 1:2, n_perm = 0, train, holdout = integer(),
                     keep_pred = FALSE, ...) {
  perms <- if (n_perm > 0) make_perms(length(train), n_perm) else NULL
  res <- lapply(lags, function(k) granger_test(Ys, list(src = Yk), X, own, test = list(src = k),
                                               train = train, holdout = holdout, perms = perms,
                                               keep_pred = keep_pred, ...))
  setting <- ifelse(lags >= 0, paste0("lag", lags), paste0("lead", -lags))
  tab <- as_tibble(bind_rows(lapply(res, `[[`, "table"))) %>% mutate(setting = setting, lag = lags, .before = 1)
  if (n_perm > 0) {
    is_lag <- lags > 0
    null <- do.call(cbind, lapply(res, `[[`, "null"))
    tab$p_perm_max <- perm_max_p(tab$F[is_lag], null[, is_lag, drop = FALSE])
  }
  pred <- bind_rows(lapply(seq_along(res), function(i)
    if (is.null(res[[i]]$pred)) NULL else res[[i]]$pred %>% mutate(setting = setting[i])))
  list(table = tab, pred = pred)
}

# -----------------------------------------------------------------------------
# DAG tests for one Northern sink with a Southern (S) and a Transitional (T)
# series - which path carries the Southern signal north?
#   DAG direct : S -> N              (S predicts N, T adds nothing given S)
#   DAG chain  : S -> T -> N         (S predicts N only through T)
#   S -> N       S lags vs own lags only        (is there a Southern signal?)
#   T -> N       T lags vs own lags only
#   S -> N | T   S lags given T lags            (the direct path)
#   T -> N | S   T lags given S lags            (the mediated path)
# Verdict: chain = T -> N | S significant, S -> N | T not, and S -> T
# established (s_to_t, from the S -> T test); direct = S -> N | T significant,
# T -> N | S not; both = both significant.
# Also the four nested models' AIC / BIC (training years) and held-out MSE.
# -----------------------------------------------------------------------------
dag_tests <- function(Ys, S, T, X, lags_S, lags_T, own = 1:2, n_perm = 199, train, holdout = integer(),
                      Ys_resp = Ys, resp_ok = NULL, H = 3, min_years = 3, alpha = 0.05, s_to_t = NA) {
  pr <- list(S = S, T = T)
  tst <- function(test, base) granger_test(Ys, pr, X, own, test = test, base = base, Ys_resp = Ys_resp,
                                           resp_ok = resp_ok, train = train, holdout = holdout, H = H,
                                           min_years = min_years,
                                           perms = if (n_perm > 0) make_perms(length(train), n_perm))$table
  tab <- as_tibble(bind_rows(
    `S -> N`     = tst(list(S = lags_S), list()),
    `T -> N`     = tst(list(T = lags_T), list()),
    `S -> N | T` = tst(list(S = lags_S), list(T = lags_T)),
    `T -> N | S` = tst(list(T = lags_T), list(S = lags_S)), .id = "test"))
  # Nested models on common rows (the design with every lag): own lags only,
  # + S (direct), + T (chain), + S + T (both)
  mods <- list(own = list(), direct = list(S = lags_S), chain = list(T = lags_T), both = list(S = lags_S, T = lags_T))
  cs  <- seasonal_clim(Ys, train, H, min_years)
  A_p <- lapply(pr, function(P) sweep(P, 2, seasonal_clim(P, train, H, min_years)))
  r_all <- granger_design(sweep(Ys, 2, cs), A_p, mods$both, own, X, resp_ok, sweep(Ys_resp, 2, cs))
  y <- r_all$M[, "y"]; tr <- r_all$year %in% train; te <- r_all$year %in% holdout
  ic <- bind_rows(lapply(names(mods), function(m) {
    keep <- c("y", paste0("own_lag", own), names(X),
              unlist(lapply(names(mods[[m]]), function(s) lag_names(s, mods[[m]][[s]]))))
    Xm <- cbind(1, r_all$M[, setdiff(keep, "y"), drop = FALSE])
    b <- ols_fit(Xm[tr, , drop = FALSE], y[tr]); n <- sum(tr); k <- ncol(Xm) + 1
    sse <- sum((y[tr] - Xm[tr, , drop = FALSE] %*% b)^2)
    data.frame(model = m, k = k, AIC = n * log(sse / n) + 2 * k, BIC = n * log(sse / n) + log(n) * k,
               mse_holdout = if (any(te)) mean((y[te] - Xm[te, , drop = FALSE] %*% b)^2) else NA_real_)
  })) %>% mutate(dAIC = AIC - min(AIC), dBIC = BIC - min(BIC))
  p <- setNames(tab$p_perm, tab$test)
  sig <- function(t) !is.na(p[[t]]) && p[[t]] < alpha
  verdict <- if (sig("S -> N | T") && sig("T -> N | S")) "both paths" else
    if (sig("T -> N | S")) { if (isTRUE(s_to_t)) "chain (S -> T -> N)" else "T -> N, but S -> T not shown" } else
    if (sig("S -> N | T")) "direct (S -> N)" else
    if (sig("S -> N") || sig("T -> N")) "signal, path unresolved" else "no signal"
  list(tests = tab, models = ic, verdict = verdict)
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
