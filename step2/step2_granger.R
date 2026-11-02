# =============================================================================
# step2/step2_granger.R
# Step 2: do Southern populations (zone_30y == 2) drive Transitional (1) and
# Northern (0) ones, from where, how fast, and along which path?
#   1. Sites are clustered into regions within each zone (space + seasonal
#      curves of the step-1 weekly populations).
#   2. Regional weekly series: observed counts + the selected step-1 model's
#      predictions for unsampled weeks.
#   3. Zone level: Southern lags 1-8 weeks -> Transitional and Northern
#      (Granger tests), with falsification (leads, reverse direction).
#   4. Region level: which Southern regions, at which lags, predict which
#      Transitional/Northern regions; lags vs distance and step-1 peak delays.
#   5. Two DAGs for the Northern populations: South -> North directly, or
#      South -> Transitional -> North (conditional Granger tests).
# Every regression is fit on the training years (2008-2019) and predicts the
# held-out year (2020), as in step 1: significance from the training years
# (year-permutation null), out-of-sample gain from the held-out year.
# Run from the project directory, after step 1 (step1/step1_hpc.R compare)
# has written results/step1/. Writes results/step2/ for CEW_JIPM_Figures.R.
#
# Sections
#   0. Step-1 outputs
#   1. Weekly populations and regional clustering
#   2. Regional weekly series
#   3. Southern lags -> Transitional / Northern zones (3b: step-1 uncertainty)
#   4. Southern regions -> Transitional / Northern regions: where and how fast
#   5. DAGs: South -> North vs South -> Transitional -> North
#   6. Save outputs
# =============================================================================

library(dplyr); library(tidyr)
source("step2/step2_functions.R")
step1_dir <- "results/step1"
step2_dir <- "results/step2"
dir.create(step2_dir, recursive = TRUE, showWarnings = FALSE)
zl <- c(`0` = "Northern", `1` = "Transitional", `2` = "Southern")

# =============================================================================
# 0. Step-1 outputs ----------------------------------------------------------------
#    Copy results/step1/ back from Hazel (hpc/README.md step 11):
#    model_data.rds (modelling data), predictions.rds (E[log(CEW + 1)] for
#    weeks 1-52 of every monitored site-year: mean + draws) and
#    peak_draws.rds (peak week of every site-year in every draw).
# =============================================================================
need <- file.path(step1_dir, c("model_data.rds", "predictions.rds", "peak_draws.rds"))
if (!all(file.exists(need)))
  stop("Missing step-1 outputs: ", paste(need[!file.exists(need)], collapse = ", "),
       ". Run ./hpc/submit_step1.sh compare and copy results/step1/ here.")
md <- readRDS(need[1]); pr <- readRDS(need[2]); pk <- readRDS(need[3])
d_all <- md$d_all; sites <- md$sites; anom_vars <- md$anom_vars; yrs <- md$yrs; peak_win <- md$peak_win
d_obs <- d_all %>% filter(observed)
holdout <- if (is.null(md$holdout) || is.na(md$holdout)) integer() else md$holdout
train_years <- setdiff(yrs, holdout)
message("Step-1 model: ", pr$model, "; training years ", paste(range(train_years), collapse = "-"),
        "; held out: ", if (length(holdout)) holdout else "none")

# =============================================================================
# 1. Weekly populations and regional clustering ------------------------------------
#    fill_mode
#      "gaps"  : observed log(CEW + 1) where the trap was sampled, the step-1
#                prediction in unsampled weeks (default)
#      "model" : the step-1 prediction everywhere. The model has no site-year
#                terms, so its week-to-week changes within a site-year come
#                only from EVI/GDD: Granger would then test the model's
#                covariates rather than the moths. For comparison only.
#    The held-out year's unsampled weeks are filled by a model that never saw
#    that year; the Granger responses are observed counts only (section 2).
#    Clustering: single traps are sampled a few weeks a year, so sites are
#    grouped into regions within their zone. Ward clustering on km coordinates
#    and each site's mean weekly curve over its training years (season
#    weeks), mixed by space_weight; k per zone = the largest with mean
#    silhouette >= 0.5 and every region >= min_sites sites (cluster_regions()).
#    Check the silhouette table; fix k with k_by_zone (e.g. c(`2` = 4)) if wanted.
# =============================================================================
fill_mode    <- "gaps"
space_weight <- 0.5
wk_season <- min(peak_win$wk_lo):max(peak_win$wk_hi)            # zone windows
wk_range  <- (min(peak_win$wk_lo) - 8):max(peak_win$wk_hi)      # + 8 lag weeks
ctrl <- d_all %>% select(location2, year, woy, all_of(anom_vars))
# The step-1 model is Gaussian on log(CEW + 1), so its expectation can dip
# below 0 where the curve extrapolates (season edges); log(CEW + 1) cannot,
# so predictions are floored at 0 (the raw values stay in predictions.rds)
fill_at <- function(observed, log_CEW, fill) {
  fill <- pmax(fill, 0)
  if (fill_mode == "gaps") ifelse(observed, log_CEW, fill) else fill
}
d_g <- pr$pred %>%
  mutate(value = fill_at(observed, log_CEW, pred_mean),
         obs_value = ifelse(observed, log_CEW, NA_real_)) %>%
  left_join(ctrl, by = c("location2", "year", "woy"))

# Each site's mean weekly curve over its training years (season weeks); sites
# trapped only in the held-out year use that year (their fill is out of sample)
cw <- d_g %>% filter(woy %in% wk_season) %>%
  group_by(location2) %>% filter(year %in% train_years | !any(year %in% train_years)) %>%
  group_by(location2, woy) %>% summarise(v = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = woy, values_from = v)
curves <- as.matrix(cw[, -1]); rownames(curves) <- cw$location2
miss <- which(is.na(curves), arr.ind = TRUE)                     # rare grid gaps
curves[miss] <- colMeans(curves, na.rm = TRUE)[miss[, "col"]]

sites_g <- sites %>% filter(location2 %in% d_obs$location2)
cl <- cluster_regions(sites_g, curves = curves, space_weight = space_weight,
                      k_range = 2:6, min_sites = 10, sil_min = 0.5)
cl$silhouette
regions <- cl$sites
reg_info <- region_info(regions)

# =============================================================================
# 2. Regional weekly series ----------------------------------------------------------
#    Y      = mean of `value` (observed + predictions) over the region's
#             monitored site-years: the lagged predictors, complete in season
#    Y_resp = mean over the sites actually trapped that week: the response,
#             so every Granger regression is scored on observed moths
#    Controls = the region's mean local weather/EVI anomalies (anom_vars).
#    Seasonality is removed inside each Granger fit (training years only).
# =============================================================================
add_region <- function(d) d %>% inner_join(regions %>% select(location2, region), by = "location2")
series <- function(d, zone_level = FALSE) {
  if (zone_level) d <- d %>% mutate(region = zone_code[as.character(zone_30y)])
  region_series(d, value = "value", resp = "obs_value", controls = anom_vars,
                years = yrs, weeks = wk_range)
}
rs_reg  <- series(add_region(d_g))
rs_zone <- series(d_g, zone_level = TRUE)
tr_i <- match(train_years, rs_zone$years)                       # rows of the year x 52 matrices
ho_i <- match(holdout, rs_zone$years)

# Response weeks: in the sink zone's season window, >= min_obs trapped sites
min_obs <- 1
zone_of <- function(r) as.integer(names(zone_code)[zone_code == substr(r, 1, 1)])
resp_ok <- function(rs, r) {
  w  <- peak_win %>% filter(zone_30y == zone_of(r))
  ok <- !is.na(rs$Y_resp[r, , ]) & rs$n_obs[r, , ] >= min_obs
  ok[, setdiff(1:52, w$wk_lo:w$wk_hi)] <- FALSE
  ok
}

# Coverage (training years): share of in-window weeks with a trapped site,
# sites per such week
cov_tab <- function(rs) {
  bind_rows(lapply(rs$regions, function(r) {
    w  <- peak_win %>% filter(zone_30y == zone_of(r))
    nw <- rs$n_obs[r, tr_i, w$wk_lo:w$wk_hi]; nw[is.na(nw)] <- 0
    data.frame(region = r, weeks_with_data = mean(nw > 0), sites_per_week = mean(nw[nw > 0]))
  }))
}
cov_tab(rs_zone)
# Regions trapped in fewer than min_cover of their in-window weeks give too
# few responses to test; they stay on the maps but get no Granger pairs
min_cover <- 0.5
reg_info <- reg_info %>% left_join(cov_tab(rs_reg), by = "region") %>%
  mutate(tested = weeks_with_data >= min_cover)
reg_info

# One sink/source pair: lags from the filled series, responses from trapped sites only
pair_scan <- function(rs, source, sink, lags, n_perm = 0, keep_pred = FALSE) {
  s <- lag_scan(rs$Y[sink, , ], rs$Y[source, , ], lapply(rs$X, function(x) x[sink, , ]), lags,
                own = own_p, n_perm = n_perm, train = tr_i, holdout = ho_i, keep_pred = keep_pred,
                Ys_resp = rs$Y_resp[sink, , ], resp_ok = resp_ok(rs, sink))
  list(table = s$table %>% mutate(source = source, sink = sink, .before = 1),
       pred = if (keep_pred && nrow(s$pred)) s$pred %>% mutate(source = source, sink = sink, year = rs$years[year_i]))
}

# =============================================================================
# 3. Southern lags -> Transitional / Northern zones -----------------------------------
#    For sink zone j and source zone i, with a_ = log(CEW + 1) minus its
#    seasonal climatology (Fourier, H = 3, training years):
#      R: a_j(w) = c + sum_l phi_l a_j(w - l) + controls_j(w) + e
#      U: R + gamma_k a_i(w - k)
#    One Southern lag k at a time (k = 1-8 weeks; own lags 1-2): the lag with
#    the largest F is the apparent arrival delay. Lags stay within a season.
#    Scores (granger_test): gain = in-sample partial R^2 of the Southern lag;
#    p_perm = year-permutation p; p_perm_max = the best lag against the null
#    of the best lag (accounts for picking it); gain_holdout = out of sample.
#    Falsification: leads 1-3 (FUTURE Southern weeks) and the reverse
#    direction (Transitional/Northern -> Southern) should not predict.
# =============================================================================
set.seed(42)
n_perm  <- 199
own_p   <- 1:2
lags_fw <- c(1:8, -(1:3))                      # lags 1-8, leads 1-3
zone_pairs <- data.frame(source = c("S", "S", "T", "N"), sink = c("T", "N", "S", "S"),
                         direction = c("forward", "forward", "reverse", "reverse"))
gz <- lapply(seq_len(nrow(zone_pairs)), function(i)
  pair_scan(rs_zone, zone_pairs$source[i], zone_pairs$sink[i], lags_fw, n_perm = n_perm, keep_pred = TRUE))
granger_zone <- bind_rows(lapply(gz, `[[`, "table")) %>% left_join(zone_pairs, by = c("source", "sink"))
granger_pred_zone <- bind_rows(lapply(gz, `[[`, "pred"))

granger_zone %>% filter(direction == "forward") %>%
  select(source, sink, setting, n, gain, coef_sum, p_F, p_perm, gain_holdout) %>% print(n = Inf)
# Falsification: leads and the reverse direction should show no signal
granger_zone %>% filter(lag < 0 | direction == "reverse") %>%
  select(source, sink, setting, gain, p_perm, gain_holdout) %>% print(n = Inf)

# Best Southern lag per sink zone (largest F among lags 1-8)
best_zone <- granger_zone %>% filter(direction == "forward", lag > 0) %>%
  group_by(sink) %>% slice_max(F, n = 1, with_ties = FALSE) %>% ungroup()
best_zone %>% select(source, sink, setting, gain, coef_sum, p_perm, p_perm_max, gain_holdout)

# ---- 3b. Step-1 uncertainty ----------------------------------------------------
#    The tests above use the posterior mean predictions. Rerun the forward
#    zone lag scans with n_fill posterior draws of them (predictions.rds keeps
#    STEP1_PRED_DRAWS): the spread of gain is the step-1 uncertainty carried
#    into step 2. Observed weeks are the same in every draw.
n_fill <- 50
granger_fill_unc <- NULL
if (ncol(pr$draws) > 0) {
  m_use <- seq_len(min(n_fill, ncol(pr$draws)))
  unc <- lapply(m_use, function(m) {
    d_m <- d_g %>% mutate(value = fill_at(observed, log_CEW, pr$draws[, m]))
    rs_m <- series(d_m, zone_level = TRUE)
    bind_rows(lapply(c("T", "N"), function(snk) pair_scan(rs_m, "S", snk, 1:8)$table)) %>% mutate(draw = m)
  })
  granger_fill_unc <- bind_rows(unc) %>% group_by(source, sink, setting, lag) %>%
    summarise(gain_med = median(gain), gain_lo = quantile(gain, 0.05), gain_hi = quantile(gain, 0.95),
              gain_holdout_med = median(gain_holdout), P_gain_holdout_pos = mean(gain_holdout > 0),
              .groups = "drop") %>%
    left_join(granger_zone %>% select(source, sink, setting, gain_at_mean = gain), by = c("source", "sink", "setting"))
  print(granger_fill_unc, n = Inf)
}

# =============================================================================
# 4. Southern regions -> Transitional / Northern regions: where and how fast --------
#    Every tested Southern region -> every tested Transitional/Northern region
#    (coverage >= min_cover, section 2), lags 1-8: which Southern locations
#    predict each sink, and at which lag. Best lag = the largest F; a link =
#    p_perm_max (the best lag against its permutation null) < 0.05 after
#    Benjamini-Hochberg across pairs, with a positive association (more
#    Southern moths -> more sink moths). Geometry: centroid distance (km) and
#    bearing; distance / best lag = apparent northward speed. Step-1 peak
#    delay: the posterior delay of the sink's mean peak week behind the
#    source's, over common years.
# =============================================================================
tested <- reg_info$region[reg_info$tested]
reg_pairs <- expand.grid(source = tested[startsWith(tested, "S")],
                         sink = tested[!startsWith(tested, "S")],
                         stringsAsFactors = FALSE)
if (!nrow(reg_pairs))
  stop("No Southern -> Transitional/Northern pair of tested regions (coverage >= min_cover = ",
       min_cover, "): see reg_info; lower min_cover or change the clustering (section 1)")
t0 <- Sys.time()
granger_reg <- bind_rows(lapply(seq_len(nrow(reg_pairs)), function(i)
  pair_scan(rs_reg, reg_pairs$source[i], reg_pairs$sink[i], 1:8, n_perm = n_perm)$table))
message(sprintf("Region pairs: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ri <- reg_info %>% select(region, longitude, latitude)
granger_pairs <- granger_reg %>%
  group_by(source, sink) %>% slice_max(F, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(source, sink, best_lag = lag, n, best_gain = gain, coef_sum, gain_holdout, p_F, p_perm_max) %>%
  mutate(p_adj = p.adjust(p_perm_max, "BH"),
         link = p_adj < 0.05 & coef_sum > 0) %>%
  left_join(ri %>% rename_with(~ paste0(.x, "_src"), -region), by = c(source = "region")) %>%
  left_join(ri %>% rename_with(~ paste0(.x, "_snk"), -region), by = c(sink = "region")) %>%
  mutate(dist_km = geosphere::distHaversine(cbind(longitude_src, latitude_src),
                                            cbind(longitude_snk, latitude_snk)) / 1000,
         bearing = geosphere::bearing(cbind(longitude_src, latitude_src),
                                      cbind(longitude_snk, latitude_snk)),
         speed_km_wk = dist_km / best_lag,
         sink_zone = ifelse(startsWith(sink, "T"), "Transitional", "Northern"))

# Step-1 peak delay of each pair (sink - source, weeks)
sy <- pk$site_years %>% left_join(regions %>% select(location2, region), by = "location2")
RP <- region_year_peaks(pk$draws[, !is.na(sy$region), drop = FALSE], sy %>% filter(!is.na(region)))
granger_pairs <- granger_pairs %>%
  bind_cols(bind_rows(lapply(seq_len(nrow(.)), function(i)
    peak_delay(RP, granger_pairs$source[i], granger_pairs$sink[i]))))
granger_pairs %>% arrange(sink, desc(best_gain)) %>%
  select(source, sink, best_lag, best_gain, coef_sum, p_adj, link, gain_holdout, dist_km, delay_med) %>% print(n = Inf)

# Does the predictive lag grow with distance (a travelling front)? Over the
# links, or over all positive associations when there are too few links
lag_dist <- granger_pairs %>% filter(if (sum(link) >= 5) link else coef_sum > 0)
if (nrow(lag_dist) >= 3) print(summary(lm(best_lag ~ dist_km, data = lag_dist)))

# =============================================================================
# 5. DAGs: South -> North vs South -> Transitional -> North ---------------------------
#    For a Northern sink with a Southern (S) and a Transitional (T) series
#    (dag_tests()): S -> N and T -> N alone, then each given the other:
#      S -> N | T significant, T -> N | S not : direct path   (DAG 1)
#      T -> N | S significant, S -> N | T not : chain via T   (DAG 2; needs S -> T,
#                                               section 3)
#      both significant                       : both paths
#    Lags: dag_lags_S / dag_lags_T weeks (fixed in advance, so the tests don't
#    reuse the lag search of sections 3-4). The four nested models' AIC/BIC
#    (training years) and held-out MSE give the same comparison as model fit.
#    Zone level: the S, T and N zone series. Region level: each tested
#    Northern region with its strongest Southern and Transitional sources
#    (best F at lags 1-8; exploratory - choosing the sources on the same
#    data makes those p-values optimistic).
# =============================================================================
dag_lags_S <- 1:6
dag_lags_T <- 1:4
run_dag <- function(rs, s_src, t_src, sink, s_to_t) {
  dag_tests(rs$Y[sink, , ], rs$Y[s_src, , ], rs$Y[t_src, , ], lapply(rs$X, function(x) x[sink, , ]),
            lags_S = dag_lags_S, lags_T = dag_lags_T, own = own_p, n_perm = n_perm,
            train = tr_i, holdout = ho_i, Ys_resp = rs$Y_resp[sink, , ], resp_ok = resp_ok(rs, sink),
            s_to_t = s_to_t)
}
# S -> T: the best Southern lag of the Transitional zone (section 3)
s_to_t <- best_zone %>% filter(sink == "T")
dag_zone <- run_dag(rs_zone, "S", "T", "N", s_to_t = s_to_t$p_perm_max < 0.05 && s_to_t$coef_sum > 0)
message(sprintf("Zone level: S -> T best lag %d, p_perm_max = %.3f; DAG verdict: %s",
                s_to_t$lag, s_to_t$p_perm_max, dag_zone$verdict))
print(dag_zone$tests %>% select(test, df1, gain, coef_sum, p_F, p_perm, gain_holdout))
print(dag_zone$models)

# Region level: strongest Transitional source of each Northern region
t_reg <- tested[startsWith(tested, "T")]; n_reg <- tested[startsWith(tested, "N")]
granger_reg_tn <- if (length(t_reg) && length(n_reg)) bind_rows(lapply(n_reg, function(snk)
  bind_rows(lapply(t_reg, function(src) pair_scan(rs_reg, src, snk, 1:8)$table)))) else NULL
dag_reg <- NULL
if (!is.null(granger_reg_tn)) {
  best_src <- function(tab, snk) tab %>% filter(sink == snk) %>% slice_max(F, n = 1, with_ties = FALSE) %>% pull(source)
  dag_reg <- bind_rows(lapply(n_reg, function(snk) {
    s_src <- best_src(granger_reg, snk); t_src <- best_src(granger_reg_tn, snk)
    st <- granger_pairs %>% filter(source == s_src, sink == t_src)   # S -> T link (section 4)
    d <- run_dag(rs_reg, s_src, t_src, snk, s_to_t = nrow(st) == 1 && st$link)
    d$tests %>% mutate(sink = snk, S_source = s_src, T_source = t_src, S_to_T_link = nrow(st) == 1 && st$link,
                       verdict = d$verdict, .before = 1)
  }))
  dag_reg %>% select(sink, S_source, T_source, test, gain, p_perm, gain_holdout, verdict) %>% print(n = Inf)
}

# =============================================================================
# 6. Save outputs (results/step2/, read by CEW_JIPM_Figures.R) -----------------------
# =============================================================================
write.csv(regions %>% select(location2, location, zone_30y, longitude, latitude, x_km, y_km, region),
          file.path(step2_dir, "cew_regions.csv"), row.names = FALSE)
write.csv(granger_zone, file.path(step2_dir, "granger_zone.csv"), row.names = FALSE)
write.csv(granger_reg, file.path(step2_dir, "granger_region_lags.csv"), row.names = FALSE)
write.csv(granger_pairs, file.path(step2_dir, "granger_pairs.csv"), row.names = FALSE)
write.csv(bind_rows(dag_zone$tests %>% mutate(level = "zone", sink = "N", S_source = "S", T_source = "T",
                                             verdict = dag_zone$verdict),
                    if (!is.null(dag_reg)) dag_reg %>% mutate(level = "region")),
          file.path(step2_dir, "granger_dag.csv"), row.names = FALSE)
write.csv(dag_zone$models, file.path(step2_dir, "granger_dag_models.csv"), row.names = FALSE)
if (!is.null(granger_fill_unc))
  write.csv(granger_fill_unc, file.path(step2_dir, "granger_fill_uncertainty.csv"), row.names = FALSE)
# Each region's mean seasonal curve (observed + predictions, averaged over years)
region_curves <- add_region(d_g) %>% filter(woy %in% wk_season) %>%
  group_by(region, zone_30y, woy) %>% summarise(value = mean(value), .groups = "drop")
saveRDS(list(regions = regions, reg_info = reg_info, silhouette = cl$silhouette,
             region_curves = region_curves,
             granger_zone = granger_zone, granger_pred_zone = granger_pred_zone, best_zone = best_zone,
             granger_reg = granger_reg, granger_pairs = granger_pairs, granger_reg_tn = granger_reg_tn,
             dag_zone = dag_zone, dag_reg = dag_reg, granger_fill_unc = granger_fill_unc,
             settings = list(n_perm = n_perm, own_lags = own_p, train_years = train_years, holdout = holdout,
                             dag_lags_S = dag_lags_S, dag_lags_T = dag_lags_T, fill_mode = fill_mode,
                             space_weight = space_weight, step1_model = pr$model)),
        file.path(step2_dir, "granger_results.rds"))
