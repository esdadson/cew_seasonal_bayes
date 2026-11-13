# =============================================================================
# step1_inla/step1_inla.R
# Step 1 with R-INLA (laptop, Windows): weekly CEW trap catches with smooth
# seasonal curves (no Fourier terms), a spatial / space-time SPDE field, and
# posterior peak-week distributions for every monitored site-year.
# A parallel prototype to the brms models in step1/ - same data, same
# held-out year - so the two can be compared.
#
#   Rscript step1_inla/step1_inla.R                            # default model set
#   INLA_MODELS="spatial_free spatial_free_nb" Rscript step1_inla/step1_inla.R
#   INLA_QUICK=1 Rscript step1_inla/step1_inla.R               # coarse mesh, fast check
#
# Model name = structure [+ family suffix]. t year, i site, w ISO week 1-52.
#
#  Zones prescribed (zone_30y enters the model):
#   seasonal       : y = a_zone + b_trap + x' beta + s_zone(w) + u_t
#                    s_zone = one cyclic second-order random walk (RW2) per
#                    zone, weeks 1-52 wrapping to week 1, summing to 0 (the
#                    level is a_zone); smoothness estimated (one precision per
#                    zone) instead of a fixed number of harmonics.
#                    u_t ~ N(0, sd_year^2): the held-out year's u is drawn
#                    from it, as in brms.
#   spatial        : seasonal + v_i + g(s_i)
#                    v_i ~ N(0, sd_site^2), the non-spatial site effect
#                    (nugget); g = Matern (nu = 1) field on a mesh over the
#                    sites (SPDE), range and sd estimated, averaging 0 over
#                    each zone's sites (the zone terms keep each zone's level).
#   spacetime      : seasonal + v_i + g(s_i, w)
#                    g varies through the season: one field per temporal knot
#                    (every INLA_KNOT_WEEKS weeks, quadratic B-splines between
#                    knots), AR1-correlated across knots, averaging 0 over each
#                    zone's sites at every knot (s_zone stays the zone's mean
#                    curve, g the within-zone departures).
#  Zone-free (zones used only to summarise the fit afterwards):
#   seasonal_free  : y = a + b_trap + x' beta + s(w) + u_t, one shared curve
#   spatial_free   : seasonal_free + v_i + g(s_i), g unconstrained
#   spacetime_free : seasonal_free + v_i + g(s_i, w), g unconstrained: any
#                    regional difference in level or timing must come from
#                    the field. INLA_ZONE_RE=1 adds zone random intercepts
#                    (iid, 3 levels) to the *_free models.
#  Family suffix:
#   (none)  Gaussian on log(CEW + 1)
#   _nb     negative binomial on round(CEW) (1.6% of weekly values are means
#           of several trap dates, not integers), log link
#   _pois   Poisson-lognormal on round(CEW): Poisson with an observation-level
#           random effect e_itw ~ N(0, sd_obs^2) (plain Poisson cannot carry
#           this overdispersion: variance ~ 600 x the mean)
#   x = weekly EVI, weekly GDD, six crop shares (standardised, as in brms).
#
# Priors (PC priors; weakly informative):
#   fixed effects N(0, 1), zone intercepts N(0, 3^2), shared intercept N(0, 10^2)
#   season curves: P(sd > 1) = 0.05 (scaled RW2)
#   sd_year, sd_site, sd_zone, sd_obs: P(sd > 1) = 0.01
#   field: P(range < INLA_RANGE0 km) = 0.1, P(sd > 1) = 0.01
#   AR1 across knots: P(rho > 0.5) = 0.9 (favours a smoothly changing field)
#   NB size: INLA default
#
# Peak week of a site-year = the week of its largest linear predictor (E[log(CEW
# + 1)] for Gaussian, log E[CEW] for the count families: the same timing),
# per posterior sample, searched within ONE window for every site - the 5th to
# 95th percentile of all sampled weeks - so zone-free and zone-prescribed
# models are judged alike (a per-zone window would prescribe zones). Every
# summary is computed sample by sample, then summarised.
#
# Comparing models across families: the native elpd's (CPO, held-out) are on
# different scales (density of log(y + 1) vs probability of a count).
# elpd_test_count scores every model on the held-out counts: for Gaussian the
# probability of the observed count implied by the normal on log(y + 1)
# (z = log(y + 1) in [log(y + 0.5), log(y + 1.5))), for NB / Poisson-lognormal
# their probability mass. elpd_cpo_count is the same for the training rows
# (Gaussian: Jacobian-adjusted, approximate). MAE / RMSE / R2 / cov90 of the
# held-out year are on log(CEW + 1) for every family (count families: from
# predictive draws).
#
# Zone summaries (zone_contrasts.csv, every model): Transitional and Northern
# minus Southern, from the fitted predictions grouped by zone_30y AFTER
# fitting: peak delay (weeks; zone-year mean site-year peak), level difference
# (mean linear predictor over the peak window), and eta2 = share of the
# variance of site peak weeks that lies between zones.
#
# Settings (environment variables):
#   INLA_MODELS      models to fit (default: see default_models below)
#   INLA_OUT         output folder                       results/step1_inla
#   STEP1_HOLDOUT    held-out year, or none              2020
#   INLA_SAMPLES     joint posterior samples             1000
#   INLA_THREADS     INLA threads "outer:inner"          "<min(cores, 16)>:1"
#   INLA_INT         hyperparameter integration          auto (ccd); eb for spacetime*
#                    (eb fixes the hyperparameters at their mode: ~5-10x faster)
#   INLA_MAX_EDGE    mesh edge near sites, km            40
#   INLA_CUTOFF      merge sites closer than, km         15
#   INLA_HULL        dense mesh within, km of sites      100
#   INLA_RANGE0      range prior: P(range < this) = 0.1  20
#   INLA_KNOT_WEEKS  spacetime: weeks between knots      4
#   INLA_ZONE_RE=1   zone random intercepts in the *_free models
#   INLA_REFIT=1     refit even if results exist
#   INLA_QUICK=1     coarse mesh (80 / 30 km), eb, 200 samples
#   STEP1_DATA_CACHE an .rds of list(d_all, sites); built from step1_data.R
#                    on the first run if given, loaded after
#
# Outputs (INLA_OUT):
#   comparison.csv                one row per model: fit, held-out and zone summaries
#   zone_contrasts.csv/.png       zone delays, level differences and eta2 by model
#   zone_curves.png               zone mean predicted curve (over years) by model
#   peak_zone_year.png, holdout_ribbon.png, mesh.png
#   <model>/summary.txt, fixed.csv, hyper.csv, season_component.csv
#   <model>/peak_site_year.csv    per site-year: median, mean, 95% interval, sd,
#                                 P(peak within its own sampled weeks), observed peak
#   <model>/peak_zone_year.csv    zone-year mean peak week (95% interval)
#   <model>/peak_zone_delay.csv   Transitional / Northern delay behind Southern
#   <model>/peak_site.csv         site peak averaged over its monitored years
#   <model>/peak_draws.rds        samples x site-years matrix + the site-year table
#   <model>/zone_curve.csv        zone-year mean linear predictor by week (95% interval)
#   <model>/pred_rows.rds         posterior mean / sd of the linear predictor, every row of d_all
#   <model>/test_rows.csv, holdout_zone_week.csv   held-out year (log(CEW + 1) scale)
#   <model>/peak_map.png, field_map.png (spatial models)
#   <model>/fit.rds, check.rds    the INLA fit and the cached results
# =============================================================================

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(Matrix); library(ggplot2)
  library(dplyr); library(tidyr)
})

env <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
default_models <- c("seasonal", "spatial", "spacetime", "spatial_free", "spacetime_free",
                    "spatial_free_nb", "spatial_free_pois", "spacetime_free_nb")
quick    <- env("INLA_QUICK", "0") == "1"
out_dir  <- env("INLA_OUT", "results/step1_inla")
models   <- strsplit(trimws(env("INLA_MODELS", paste(default_models, collapse = " "))), "\\s+")[[1]]
n_samp   <- as.integer(env("INLA_SAMPLES", if (quick) "200" else "1000"))
threads  <- env("INLA_THREADS", paste0(min(parallel::detectCores(), 16), ":1"))
int_for  <- function(m) env("INLA_INT", if (quick || grepl("^spacetime", m)) "eb" else "auto")
max_edge <- as.numeric(env("INLA_MAX_EDGE", if (quick) "80" else "40"))
cutoff   <- as.numeric(env("INLA_CUTOFF", if (quick) "30" else "15"))
hull_km  <- as.numeric(env("INLA_HULL", "100"))
range0   <- as.numeric(env("INLA_RANGE0", "20"))
knot_wk  <- as.numeric(env("INLA_KNOT_WEEKS", "4"))
zone_re  <- env("INLA_ZONE_RE", "0") == "1"
refit    <- env("INLA_REFIT", "0") == "1"
ho       <- env("STEP1_HOLDOUT", "2020")
holdout  <- if (ho %in% c("", "none")) NA_integer_ else as.integer(ho)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

# Model name -> structure and family
bases <- c("seasonal", "spatial", "spacetime", "seasonal_free", "spatial_free", "spacetime_free")
parse_model <- function(m) {
  fam  <- if (grepl("_nb$", m)) "nbinomial" else if (grepl("_pois$", m)) "poisson" else "gaussian"
  base <- sub("_(nb|pois)$", "", m)
  if (!base %in% bases) stop("Unknown model '", m, "'. Structures: ", paste(bases, collapse = ", "),
                             "; optional suffix _nb or _pois")
  list(name = m, base = base, family = fam, zoned = !grepl("_free$", base),
       field = grepl("^spatial", base), st = grepl("^spacetime", base),
       site = !grepl("^seasonal", base), obs = fam == "poisson",
       zone_re = zone_re && grepl("_free$", base))
}
specs <- setNames(lapply(models, parse_model), models)
message(sprintf("INLA %s | models: %s | samples %d | threads %s | mesh %g/%g km%s%s",
                packageVersion("INLA"), paste(models, collapse = " "), n_samp, threads,
                max_edge, cutoff, if (zone_re) " | zone RE in *_free" else "", if (quick) " | QUICK" else ""))

# =============================================================================
# Data -----------------------------------------------------------------------------
#   d = d_all (observed weeks + the weeks 1-52 grid of every monitored
#   site-year), sorted so each site-year is 52 consecutive rows.
# =============================================================================
cache <- Sys.getenv("STEP1_DATA_CACHE")
if (nzchar(cache) && file.exists(cache)) {
  message("Loading ", cache); dc <- readRDS(cache); d_all <- dc$d_all; sites <- dc$sites
} else {
  source("step1/step1_data.R")
  if (nzchar(cache)) saveRDS(list(d_all = d_all, sites = sites), cache)
}

covs <- c("EVI_filled", "gdd", "corn", "cotton", "soybeans", "tobacco", "peanuts", "other")
d <- d_all %>% arrange(location2, year, woy)
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))
d <- d %>% mutate(fit_row  = observed & (is.na(holdout) | year != holdout),
                  test_row = observed & !is.na(holdout) & year == holdout,
                  zone_id  = as.integer(as.character(zone_30y)) + 1L,
                  count    = round(CEW))
yrs_all  <- sort(unique(d$year));      d$year_id <- match(d$year, yrs_all)
site_ids <- sort(unique(d$location2)); d$site_id <- match(d$location2, site_ids)
n_site   <- length(site_ids)

# Fixed-effect designs: zone intercepts (zoned) or one intercept (zone-free),
# trap type against BL, covariates
trap_X <- model.matrix(~ trap_type, d)[, -1, drop = FALSE]
colnames(trap_X) <- make.names(sub("^trap_type", "trap_", colnames(trap_X)))
X_zone <- cbind(z0 = as.numeric(d$zone_id == 1), z1 = as.numeric(d$zone_id == 2),
                z2 = as.numeric(d$zone_id == 3), trap_X, as.matrix(d[, covs]))
X_free <- cbind(intercept = 1, trap_X, as.matrix(d[, covs]))
Xof <- function(s) if (s$zoned) X_zone else X_free
if (anyNA(X_zone[d$observed, ])) stop("NA in the fixed-effect design of observed rows")
fit_idx  <- which(d$fit_row); test_idx <- which(d$test_row)
message("Training rows: ", length(fit_idx), "; held-out (", holdout, ") rows: ", length(test_idx),
        "; site-years (all years): ", nrow(d) / 52)

zone_lab <- c(`0` = "Northern limits", `1` = "Transitional zone", `2` = "Southern range")
peak_win <- c(lo = floor(quantile(d$woy[fit_idx], 0.05, names = FALSE)),
              hi = ceiling(quantile(d$woy[fit_idx], 0.95, names = FALSE)))
in_win   <- 1:52 >= peak_win[["lo"]] & 1:52 <= peak_win[["hi"]]
message(sprintf("Peak window (all sites): weeks %d-%d", peak_win[["lo"]], peak_win[["hi"]]))

# =============================================================================
# Mesh, SPDE and the temporal basis -----------------------------------------------
#   Mesh over every site (fitting and held-out), in CONUS Albers km: fine
#   (INLA_MAX_EDGE) within INLA_HULL km of a site, coarse beyond so the
#   field is not distorted at the edge.
# =============================================================================
site_xy <- d %>% distinct(site_id, location2, x_km, y_km, longitude, latitude, zone_30y) %>% arrange(site_id)
stopifnot(nrow(site_xy) == n_site)
loc  <- as.matrix(site_xy[, c("x_km", "y_km")])
mesh <- fm_mesh_2d(loc = loc,
                   boundary = list(fm_nonconvex_hull(loc, convex = hull_km),
                                   fm_nonconvex_hull(loc, convex = max(400, 4 * hull_km))),
                   max.edge = c(max_edge, 300), cutoff = cutoff)
spde <- inla.spde2.pcmatern(mesh, prior.range = c(range0, 0.1), prior.sigma = c(1, 0.01))
A_site <- fm_basis(mesh, loc = loc)                      # sites x mesh nodes
n_mesh <- mesh$n

knots <- unique(c(seq(1, 52, by = knot_wk), 52))
tmesh <- fm_mesh_1d(knots, degree = 2, boundary = "free")
B_wk  <- as.matrix(fm_basis(tmesh, loc = 1:52))          # 52 weeks x K knots
n_knot <- ncol(B_wk)
message(sprintf("Mesh: %d nodes (edge %g km near sites, cutoff %g km); spacetime knots: %d (every %g weeks)",
                n_mesh, max_edge, cutoff, n_knot, knot_wk))

# Row-wise Kronecker product: column (k - 1) * m + j = As[, j] * Bt[, k],
# matching INLA's ordering of a grouped effect (index fastest, then group)
row_kron <- function(As, Bt) {
  a <- as(As, "TsparseMatrix"); b <- as(as(Bt, "CsparseMatrix"), "TsparseMatrix")
  t <- merge(data.frame(i = a@i, j = a@j, x = a@x), data.frame(i = b@i, k = b@j, y = b@x), by = "i")
  sparseMatrix(i = t$i + 1, j = t$k * ncol(As) + t$j + 1, x = t$x * t$y,
               dims = c(nrow(As), ncol(As) * ncol(Bt)))
}

# =============================================================================
# Model terms ----------------------------------------------------------------------
# =============================================================================
hp_season <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
hp_sd1    <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hp_rho    <- list(rho = list(prior = "pc.cor1", param = c(0.5, 0.9)))

# Zone-prescribed models only: the field averaged over each zone's sites is 0
# (per knot in spacetime). Without it a broad field reproduces the zone
# intercepts and curves (their sd's collapsed to ~0 in a test). INLA repeats a
# grouped effect's constraint for every group (knot).
constr_field <- local({
  Z <- vapply(0:2, function(z) Matrix::colSums(A_site[site_xy$zone_30y == z, , drop = FALSE]) /
                sum(site_xy$zone_30y == z), numeric(n_mesh))
  list(A = t(Z), e = rep(0, 3))
})

rw2 <- function(v) sprintf('f(%s, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE, values = 1:52, hyper = hp_season)', v)
st_f <- 'f(st, model = spde, group = st.group, control.group = list(model = "ar1", hyper = hp_rho)%s)'
terms_for <- function(s) {
  t <- c(if (s$zoned) rw2(c("w0", "w1", "w2")) else rw2("w"),
         'f(year_id, model = "iid", values = seq_along(yrs_all), hyper = hp_sd1)')
  if (s$site)    t <- c(t, 'f(site_id, model = "iid", values = seq_len(n_site), hyper = hp_sd1)')
  if (s$field)   t <- c(t, if (s$zoned) 'f(field, model = spde, extraconstr = constr_field)' else 'f(field, model = spde)')
  if (s$st)      t <- c(t, sprintf(st_f, if (s$zoned) ", extraconstr = constr_field" else ""))
  if (s$zone_re) t <- c(t, 'f(zone_id, model = "iid", values = 1:3, hyper = hp_sd1)')
  if (s$obs)     t <- c(t, 'f(obs_id, model = "iid", hyper = hp_sd1)')
  t
}
inla_formula <- function(s)
  as.formula(paste("y ~ 0 +", paste(colnames(Xof(s)), collapse = " + "), "+",
                   paste(terms_for(s), collapse = " + ")), env = globalenv())

# Week of each row in its zone's curve (NA = not this zone: no contribution); w = shared curve
week_cols <- function(idx) {
  w <- d$woy[idx]; z <- d$zone_id[idx]
  data.frame(w = w, w0 = ifelse(z == 1, w, NA), w1 = ifelse(z == 2, w, NA), w2 = ifelse(z == 3, w, NA))
}

inla_stack <- function(s, idx) {
  A   <- list(1)
  eff <- list(cbind(as.data.frame(Xof(s)[idx, , drop = FALSE]), week_cols(idx),
                    year_id = d$year_id[idx], site_id = d$site_id[idx], zone_id = d$zone_id[idx],
                    obs_id = seq_along(idx)))
  if (s$field) {
    A   <- c(list(A_site[d$site_id[idx], , drop = FALSE]), A)
    eff <- c(list(list(field = seq_len(n_mesh))), eff)
  }
  if (s$st) {
    A   <- c(list(row_kron(A_site[d$site_id[idx], , drop = FALSE], B_wk[d$woy[idx], , drop = FALSE])), A)
    eff <- c(list(list(st = rep(seq_len(n_mesh), n_knot), st.group = rep(seq_len(n_knot), each = n_mesh))), eff)
  }
  y <- if (s$family == "gaussian") d$log_CEW[idx] else d$count[idx]
  inla.stack(data = list(y = y), A = A, effects = eff, tag = "fit", remove.unused = FALSE)
}

prec_fixed <- c(list(default = 1, intercept = 0.01), setNames(as.list(rep(1 / 9, 3)), c("z0", "z1", "z2")))

fit_inla <- function(s) {
  stk <- inla_stack(s, fit_idx)
  inla(inla_formula(s), family = s$family,
       data = inla.stack.data(stk, spde = spde),
       control.predictor = list(A = inla.stack.A(stk), compute = FALSE, link = 1),
       control.compute = list(config = TRUE, waic = TRUE, dic = TRUE, cpo = TRUE,
                              return.marginals.predictor = FALSE),
       control.fixed = list(mean = 0, prec = prec_fixed),
       control.inla = list(int.strategy = int_for(s$name)),
       num.threads = threads, verbose = FALSE)
}

# =============================================================================
# Posterior samples -> linear predictor for any rows ----------------------------------
#   One joint sample of every latent component (and the hyperparameters);
#   eta of a row = X beta + s(week) + u_year + v_site + g(site[, week]) (+ zone RE).
#   The held-out year's u and a new site's v come from their priors given the
#   sampled sd (INLA includes them in the latent field, no data). The
#   Poisson-lognormal observation effect is noise: left out of eta, drawn for
#   predictive draws only.
# =============================================================================
draw_posterior <- function(fit, s, n) {
  X <- Xof(s); sea <- if (s$zoned) c("w0", "w1", "w2") else "w"
  sel <- c(setNames(as.list(rep(1, ncol(X))), colnames(X)),
           setNames(as.list(rep(0, length(sea))), sea), list(year_id = 0))
  if (s$site)    sel$site_id <- 0
  if (s$field)   sel$field <- 0
  if (s$st)      sel$st <- 0
  if (s$zone_re) sel$zone_id <- 0
  smp <- inla.posterior.sample(n, fit, selection = sel, num.threads = threads)
  rn  <- rownames(smp[[1]]$latent)
  grab <- function(p) {
    i <- startsWith(rn, paste0(p, ":"))
    if (!any(i)) stop("No latent rows for ", p)
    matrix(vapply(smp, function(x) x$latent[i], numeric(sum(i))), sum(i))
  }
  hyp <- function(pat) vapply(smp, function(x) {
    h <- x$hyperpar; j <- grep(pat, names(h), fixed = TRUE)
    if (!length(j)) stop("No hyperparameter matching '", pat, "'"); h[[j[1]]] }, 1)
  P <- list(beta = do.call(rbind, lapply(colnames(X), grab)), season = lapply(sea, grab), year = grab("year_id"))
  if (s$family == "gaussian")  P$sigma  <- 1 / sqrt(hyp("Precision for the Gaussian observations"))
  if (s$family == "nbinomial") P$size   <- hyp("size for the nbinomial observations")
  if (s$obs)                   P$sd_obs <- 1 / sqrt(hyp("Precision for obs_id"))
  stopifnot(nrow(P$beta) == ncol(X), all(vapply(P$season, nrow, 1) == 52), nrow(P$year) == length(yrs_all))
  if (s$site)    { P$site  <- grab("site_id"); stopifnot(nrow(P$site) == n_site) }
  if (s$field)   { P$field <- grab("field");   stopifnot(nrow(P$field) == n_mesh) }
  if (s$st)      { P$st    <- grab("st");      stopifnot(nrow(P$st) == n_mesh * n_knot) }
  if (s$zone_re) { P$zone  <- grab("zone_id"); stopifnot(nrow(P$zone) == 3) }
  P
}

# linear predictor, rows idx x samples b
eta_rows <- function(P, s, idx, b) {
  e <- Xof(s)[idx, , drop = FALSE] %*% P$beta[, b, drop = FALSE]
  w <- d$woy[idx]; si <- d$site_id[idx]
  if (s$zoned) {
    z <- d$zone_id[idx]
    for (k in 1:3) { j <- which(z == k); if (length(j)) e[j, ] <- e[j, ] + P$season[[k]][w[j], b, drop = FALSE] }
  } else e <- e + P$season[[1]][w, b, drop = FALSE]
  e <- e + P$year[d$year_id[idx], b, drop = FALSE]
  if (!is.null(P$site))  e <- e + P$site[si, b, drop = FALSE]
  if (!is.null(P$zone))  e <- e + P$zone[d$zone_id[idx], b, drop = FALSE]
  if (!is.null(P$field)) e <- e + as.matrix(A_site %*% P$field[, b, drop = FALSE])[si, , drop = FALSE]
  if (!is.null(P$st)) {
    for (q in seq_along(b)) {                     # sites x 52 weeks for sample b[q]
      G <- as.matrix(A_site %*% matrix(P$st[, b[q]], n_mesh, n_knot)) %*% t(B_wk)
      e[, q] <- e[, q] + G[cbind(si, w)]
    }
  }
  e
}

# =============================================================================
# Predictions, peaks and held-out scores from the samples ------------------------------
#   All rows of d (every site-year x 52 weeks), in batches of samples.
# =============================================================================
q95 <- function(x) c(med = median(x), lo = quantile(x, 0.025, names = FALSE),
                     hi = quantile(x, 0.975, names = FALSE))
season_block <- function(woy) cut(woy, c(0, 22, 27, 32, 37, 53),
                                  labels = c("<=22", "23-27", "28-32", "33-37", ">=38"))

predict_all <- function(P, s, batch = 100) {
  S    <- ncol(P$beta)
  sy   <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "year", "zone_30y")]
  n_sy <- nrow(sy)
  zk   <- paste(d$zone_30y, d$year, d$woy)                                   # zone-year-week
  zkey <- unique(zk)
  G    <- sparseMatrix(seq_len(nrow(d)), match(zk, zkey), x = 1)
  G    <- G %*% Diagonal(x = 1 / colSums(G))                                 # row means
  PK   <- matrix(NA_integer_, S, n_sy, dimnames = list(NULL, sy$loc_year))
  ZC   <- matrix(NA_real_, length(zkey), S, dimnames = list(zkey, NULL))
  ET   <- matrix(NA_real_, length(test_idx), S)
  m1 <- m2 <- numeric(nrow(d))
  for (b in split(seq_len(S), ceiling(seq_len(S) / batch))) {
    e  <- eta_rows(P, s, seq_len(nrow(d)), b)
    m1 <- m1 + rowSums(e); m2 <- m2 + rowSums(e^2)
    ZC[, b] <- as.matrix(crossprod(G, e))
    if (length(test_idx)) ET[, b] <- e[test_idx, , drop = FALSE]
    for (q in seq_along(b)) {
      M <- matrix(e[, q], 52); M[!in_win, ] <- -Inf
      PK[b[q], ] <- max.col(t(M), ties.method = "first")
    }
    rm(e); invisible(gc())
  }
  pm <- m1 / S
  list(PK = PK, sy = sy, ZC = ZC, ET = ET,
       pred = d %>% select(location2, loc_year, zone_30y, year, woy, observed, log_CEW, count) %>%
         mutate(eta_mean = pm, eta_sd = sqrt(pmax(m2 / S - pm^2, 0))))
}

peak_summaries <- function(PK, sy) {
  obs <- d %>% filter(observed) %>% group_by(loc_year) %>%
    summarise(samp_lo = min(woy), samp_hi = max(woy), n_weeks = n(),
              obs_peak = woy[which.max(log_CEW)], .groups = "drop")
  sy <- left_join(sy, obs, by = "loc_year")
  site_year <- sy %>%
    mutate(peak_med = apply(PK, 2, median), peak_mean = colMeans(PK),
           peak_lo = apply(PK, 2, quantile, 0.025), peak_hi = apply(PK, 2, quantile, 0.975),
           peak_sd = apply(PK, 2, sd), ci_width = peak_hi - peak_lo,
           in_sampled = colMeans(sweep(PK, 2, samp_lo, ">=") & sweep(PK, 2, samp_hi, "<="))) %>%
    left_join(site_xy %>% select(location2, longitude, latitude, x_km, y_km), by = "location2")

  ZY <- vapply(split(seq_len(ncol(PK)), paste(sy$zone_30y, sy$year)),
               function(j) rowMeans(PK[, j, drop = FALSE]), numeric(nrow(PK)))
  zone_year <- data.frame(key = colnames(ZY), t(apply(ZY, 2, q95))) %>%
    separate(key, c("zone_30y", "year"), convert = TRUE) %>%
    left_join(count(sy, zone_30y, year, name = "n_sites"), by = c("zone_30y", "year")) %>%
    mutate(zone = zone_lab[as.character(zone_30y)])

  zone_delay <- bind_rows(lapply(c(1, 0), function(z) {
    ys <- intersect(sy$year[sy$zone_30y == z], sy$year[sy$zone_30y == 2])
    Dz <- sapply(ys, function(y) ZY[, paste(z, y)] - ZY[, paste(2, y)])
    colnames(Dz) <- ys
    Dz <- cbind(Dz, all_years = rowMeans(Dz))
    data.frame(zone = zone_lab[as.character(z)], year = colnames(Dz), t(apply(Dz, 2, q95)),
               P_later = colMeans(Dz > 0), row.names = NULL)
  }))

  # site peak averaged over its monitored years, per sample (S x sites)
  Ws <- sparseMatrix(seq_len(nrow(sy)), match(sy$location2, site_ids), x = 1, dims = c(nrow(sy), n_site))
  ns <- colSums(Ws); has_sy <- ns > 0
  SP <- as.matrix(PK %*% Ws)[, has_sy, drop = FALSE] / rep(ns[has_sy], each = nrow(PK))
  site <- site_xy[has_sy, ] %>% mutate(n_years = ns[has_sy]) %>%
    bind_cols(as.data.frame(t(apply(SP, 2, q95)))) %>% mutate(ci_width = hi - lo)
  list(site_year = site_year, zone_year = zone_year, zone_delay = zone_delay, site = site,
       SP = SP, site_zone = site_xy$zone_30y[has_sy])
}

# Zone summaries after fitting (zones are labels here, not model terms)
zone_summaries <- function(pr, pk, model) {
  key <- do.call(rbind, strsplit(rownames(pr$ZC), " "))                 # zone, year, week
  zw  <- paste(key[, 1], key[, 3]); zwk <- unique(zw)
  W   <- sparseMatrix(seq_along(zw), match(zw, zwk), x = 1)
  ZW  <- as.matrix(crossprod(W %*% Diagonal(x = 1 / colSums(W)), pr$ZC))  # zone-week (mean over years) x S
  zk  <- do.call(rbind, strsplit(zwk, " "))
  curve <- data.frame(model = model, zone_30y = as.integer(zk[, 1]), woy = as.integer(zk[, 2]),
                      mean = rowMeans(ZW), lo = matrixStats::rowQuantiles(ZW, probs = 0.025),
                      hi = matrixStats::rowQuantiles(ZW, probs = 0.975))
  lev <- sapply(0:2, function(z) colMeans(ZW[as.integer(zk[, 1]) == z & in_win[as.integer(zk[, 2])], , drop = FALSE]))
  lev_c <- bind_rows(lapply(c(1, 0), function(z) {
    v <- lev[, z + 1] - lev[, 3]
    data.frame(zone = zone_lab[as.character(z)], measure = "level difference vs Southern (linear predictor)",
               t(q95(v)), P_pos = mean(v > 0))
  }))
  del_c <- pk$zone_delay %>% filter(year == "all_years") %>%
    transmute(zone, measure = "peak delay vs Southern (weeks)", med, lo, hi, P_pos = P_later)
  SP <- pk$SP; zs <- pk$site_zone; mu <- rowMeans(SP)
  tot <- rowSums((SP - mu)^2)
  btw <- Reduce(`+`, lapply(unique(zs), function(z) { j <- zs == z; sum(j) * (rowMeans(SP[, j, drop = FALSE]) - mu)^2 }))
  eta2 <- ifelse(tot > 0, btw / tot, NA)
  e2 <- data.frame(zone = "all", measure = "eta2: share of site peak-week variance between zones",
                   t(if (any(!is.na(eta2))) q95(eta2[!is.na(eta2)]) else c(med = NA, lo = NA, hi = NA)), P_pos = NA)
  list(curve = curve, contrasts = bind_rows(del_c, lev_c, e2) %>% mutate(model = model, .before = 1))
}

# Held-out year. EL = E[log(CEW + 1)] per row and sample (Gaussian: eta; count
# families: mean of log(y_rep + 1) over 20 predictive draws); YR = one
# predictive draw on the log(CEW + 1) scale; lp_cnt = log probability of the
# observed (rounded) count.
holdout_scores <- function(ET, P, s, model) {
  if (!length(test_idx)) return(NULL)
  dt <- d[test_idx, ]; S <- ncol(ET); n <- nrow(ET)
  y_log <- dt$log_CEW; Y <- matrix(dt$count, n, S)
  set.seed(42)
  if (s$family == "gaussian") {
    SG <- matrix(P$sigma, n, S, byrow = TRUE)
    EL <- ET
    YR <- ET + SG * matrix(rnorm(n * S), n)
    lp_native <- matrix(dnorm(matrix(y_log, n, S), ET, SG, log = TRUE), n)
    up <- matrix(pnorm(log(Y + 1.5), ET, SG), n); lo <- matrix(ifelse(Y == 0, 0, pnorm(log(Y + 0.5), ET, SG)), n)
    lp_cnt <- log(pmax(up - lo, 1e-300))
  } else if (s$family == "nbinomial") {
    SZ <- matrix(P$size, n, S, byrow = TRUE); MU <- exp(ET)
    rdraw <- function() matrix(rnbinom(n * S, size = SZ, mu = MU), n)
    EL <- Reduce(`+`, lapply(1:20, function(k) log(rdraw() + 1))) / 20
    YR <- log(rdraw() + 1)
    lp_cnt <- lp_native <- matrix(dnbinom(Y, size = SZ, mu = MU, log = TRUE), n)
  } else {                                        # Poisson-lognormal
    SO <- matrix(P$sd_obs, n, S, byrow = TRUE)
    rdraw <- function() matrix(rpois(n * S, exp(ET + SO * matrix(rnorm(n * S), n))), n)
    EL <- Reduce(`+`, lapply(1:20, function(k) log(rdraw() + 1))) / 20
    YR <- log(rdraw() + 1)
    eps <- qnorm((1:20 - 0.5) / 20)              # quadrature over e ~ N(0, sd_obs^2)
    LP  <- vapply(eps, function(ep) matrix(dpois(Y, exp(ET + SO * ep), log = TRUE), n), matrix(0, n, S))
    lp_cnt <- lp_native <- matrix(matrixStats::rowLogSumExps(matrix(LP, n * S)) - log(length(eps)), n)
  }
  lpd   <- matrixStats::rowLogSumExps(lp_native) - log(S)
  lpd_c <- matrixStats::rowLogSumExps(lp_cnt) - log(S)
  pi90 <- matrixStats::rowQuantiles(YR, probs = c(0.05, 0.95))
  ci   <- matrixStats::rowQuantiles(EL, probs = c(0.025, 0.975))
  mt   <- rowMeans(EL)
  metrics <- data.frame(n_test = n, elpd_test = sum(lpd), elpd_test_count = sum(lpd_c),
                        elpd_test_count_se = sd(lpd_c) * sqrt(n),
                        mae_test = mean(abs(y_log - mt)), rmse_test = sqrt(mean((y_log - mt)^2)),
                        R2_test = 1 - sum((y_log - mt)^2) / sum((y_log - mean(y_log))^2),
                        cov90_test = mean(y_log >= pi90[, 1] & y_log <= pi90[, 2]))
  test_rows <- dt %>% select(location2, loc_year, zone_30y, year, woy, log_CEW, count) %>%
    mutate(model = model, pred_mean = mt, ci_lo = ci[, 1], ci_hi = ci[, 2],
           pi_lo = pi90[, 1], pi_hi = pi90[, 2], lpd_count = lpd_c, block = season_block(woy))
  g  <- droplevels(interaction(dt$zone_30y, dt$woy, lex.order = TRUE))
  A  <- model.matrix(~ 0 + g); A <- sweep(A, 2, colSums(A), "/")
  em <- crossprod(A, EL); ym <- crossprod(A, YR)
  zone_week <- dt %>% mutate(g = g) %>% group_by(g) %>%
    summarise(zone_30y = first(zone_30y), woy = first(woy), n = n(), obs = mean(log_CEW), .groups = "drop") %>%
    arrange(g) %>% select(-g) %>%
    mutate(model = model, year = holdout, pred_mean = rowMeans(em),
           ci_lo = matrixStats::rowQuantiles(em, probs = 0.025), ci_hi = matrixStats::rowQuantiles(em, probs = 0.975),
           pi_lo = matrixStats::rowQuantiles(ym, probs = 0.025), pi_hi = matrixStats::rowQuantiles(ym, probs = 0.975),
           level = 0.95)
  list(metrics = metrics, test_rows = test_rows, zone_week = zone_week)
}

# Hyperparameters on interpretable scales (sd instead of precision)
hyper_table <- function(fit) {
  h <- fit$summary.hyperpar
  out <- data.frame(parameter = rownames(h), mean = h$mean, lo = h$`0.025quant`, med = h$`0.5quant`,
                    hi = h$`0.975quant`, row.names = NULL)
  pr <- grepl("^Precision", out$parameter)
  out[pr, c("mean", "lo", "med", "hi")] <- 1 / sqrt(out[pr, c("mean", "hi", "med", "lo")])
  out$parameter[pr] <- sub("^Precision for", "SD for", out$parameter[pr])
  out$parameter <- sub("Range for (field|st)", "Range (km) for \\1", out$parameter)
  out$parameter <- sub("Stdev for (field|st)", "SD for \\1", out$parameter)
  out
}

# =============================================================================
# Maps ------------------------------------------------------------------------------
# =============================================================================
outline <- tryCatch({
  v  <- terra::unwrap(readRDS("gadm/gadm41_USA_0_pk.rds"))
  ol <- sf::st_transform(sf::st_as_sf(v), 5070)
  bb <- sf::st_bbox(c(xmin = (min(loc[, 1]) - 300) * 1000, ymin = (min(loc[, 2]) - 300) * 1000,
                      xmax = (max(loc[, 1]) + 300) * 1000, ymax = (max(loc[, 2]) + 300) * 1000),
                    crs = sf::st_crs(5070))
  ol <- sf::st_simplify(suppressWarnings(sf::st_crop(ol, bb)), dTolerance = 3000)
  xy <- as.data.frame(sf::st_coordinates(ol))
  xy %>% mutate(x_km = X / 1000, y_km = Y / 1000, grp = paste(L1, L2, L3))
}, error = function(e) { message("No map outline (", conditionMessage(e), ")"); NULL })

map_base <- function() {
  p <- ggplot() + coord_equal(xlim = range(loc[, 1]) + c(-60, 60), ylim = range(loc[, 2]) + c(-60, 60)) +
    theme_bw(base_size = 9) + theme(axis.title = element_blank(), axis.text = element_blank(),
                                    axis.ticks = element_blank(), panel.grid = element_blank())
  if (!is.null(outline)) p <- p + geom_polygon(data = outline, aes(x_km, y_km, group = grp),
                                               fill = "grey97", colour = "grey55", linewidth = 0.25)
  p
}

# field mean on a grid, masked to within 75 km of a site
field_grid <- function(values) {
  pr <- fm_evaluator(mesh, xlim = range(loc[, 1]) + c(-60, 60), ylim = range(loc[, 2]) + c(-60, 60),
                     dims = c(260, 200))
  g  <- expand.grid(x_km = pr$x, y_km = pr$y)
  g$value <- as.vector(fm_evaluate(pr, field = values))
  near <- vapply(seq_len(nrow(g)), function(i) min((loc[, 1] - g$x_km[i])^2 + (loc[, 2] - g$y_km[i])^2), 1)
  g[sqrt(near) <= 75, ]
}

plot_model_maps <- function(s, fit, pk, mdir) {
  p1 <- map_base() + geom_point(data = pk$site, aes(x_km, y_km, colour = med), size = 1.3) +
    scale_colour_viridis_c(name = "Peak week\n(median)", option = "D") +
    labs(subtitle = paste0("A  ", s$name, ": site peak week, mean over monitored years"))
  p2 <- map_base() + geom_point(data = pk$site, aes(x_km, y_km, colour = ci_width), size = 1.3) +
    scale_colour_viridis_c(name = "95% CrI\nwidth (wk)", option = "C") + labs(subtitle = "B  Uncertainty")
  ggsave(file.path(mdir, "peak_map.png"), patchwork::wrap_plots(p1, p2, ncol = 1), width = 6.5, height = 8, dpi = 300)

  if (s$field) {
    g <- field_grid(fit$summary.random$field$mean)
    p <- map_base() + geom_tile(data = g, aes(x_km, y_km, fill = value), alpha = 0.9) +
      geom_point(data = site_xy, aes(x_km, y_km), size = 0.3) +
      scale_fill_gradient2(name = "Field") + labs(subtitle = paste(s$name, "- posterior mean of the spatial field"))
    ggsave(file.path(mdir, "field_map.png"), p, width = 6.5, height = 4.5, dpi = 300)
  }
  if (s$st) {
    Wm <- matrix(fit$summary.random$st$mean, n_mesh, n_knot)
    wks <- c(18, 24, 30, 36)
    g <- bind_rows(lapply(wks, function(w) field_grid(as.vector(Wm %*% B_wk[w, ])) %>% mutate(week = paste("Week", w))))
    g$week <- factor(g$week, paste("Week", wks))
    p <- map_base() + geom_tile(data = g, aes(x_km, y_km, fill = value), alpha = 0.9) +
      facet_wrap(~ week, ncol = 2) + scale_fill_gradient2(name = "Field") +
      labs(subtitle = paste(s$name, "- posterior mean of the field at four weeks"))
    ggsave(file.path(mdir, "field_map.png"), p, width = 7, height = 6, dpi = 300)
  }
}

# =============================================================================
# Run -------------------------------------------------------------------------------
# =============================================================================
png(file.path(out_dir, "mesh.png"), width = 2000, height = 1500, res = 200)
plot(mesh, main = sprintf("Mesh: %d nodes", n_mesh)); points(loc, col = "red", pch = 20, cex = 0.4)
invisible(dev.off())

res <- list()
for (m in models) {
  s <- specs[[m]]
  mdir <- file.path(out_dir, m); dir.create(mdir, showWarnings = FALSE)
  ck_file <- file.path(mdir, "check.rds")
  if (!refit && file.exists(ck_file)) {
    ck <- readRDS(ck_file)
    if (!is.null(ck$zone$contrasts) && !is.null(ck$test_lpd)) { message("Reusing ", ck_file); res[[m]] <- ck; next }
  }

  message("\n==== ", m, " (", s$family, if (s$zoned) ", zones prescribed" else ", zone-free", ") ====")
  t0  <- Sys.time()
  fit <- fit_inla(s)
  fit_min <- mins(t0)
  if (!is.null(fit$mode$mode.status) && fit$mode$mode.status != 0)
    warning(m, ": INLA mode.status = ", fit$mode$mode.status, " (optimiser did not converge cleanly)")
  message(sprintf("Fitted in %.1f min", fit_min))
  saveRDS(fit, file.path(mdir, "fit.rds"))
  hyp <- hyper_table(fit)
  print(hyp, digits = 3)

  t1 <- Sys.time()
  P  <- draw_posterior(fit, s, n_samp)
  pr <- predict_all(P, s)
  pk <- peak_summaries(pr$PK, pr$sy)
  zs <- zone_summaries(pr, pk, m)
  ho_res <- holdout_scores(pr$ET, P, s, m)
  samp_min <- mins(t1)
  message(sprintf("%d samples, predictions, peaks and scores in %.1f min", n_samp, samp_min))

  cpo <- fit$cpo$cpo; ok <- !is.na(cpo) & cpo > 0
  jac <- if (s$family == "gaussian") sum(log(d$CEW[fit_idx][ok] + 1)) else 0
  zc  <- zs$contrasts
  zval <- function(z, meas) { r <- zc[zc$zone == z & startsWith(zc$measure, meas), ]; if (nrow(r)) r$med else NA }
  # CPO is unreliable where INLA flags failures (e.g. the Poisson-lognormal's
  # observation-level effect: leaving out y_i leaves its own effect unidentified)
  n_fail <- sum(fit$cpo$failure > 0, na.rm = TRUE) + sum(!ok)
  cpo_ok <- n_fail <= 0.01 * length(cpo)
  stats <- data.frame(model = m, family = s$family, zones = if (s$zoned) "prescribed" else "free",
                      n_train = length(fit_idx), mesh_nodes = n_mesh,
                      mlik = unname(fit$mlik[1, 1]),
                      elpd_cpo = if (cpo_ok) sum(log(cpo[ok])) else NA,
                      elpd_cpo_count = if (cpo_ok) sum(log(cpo[ok])) - jac else NA,
                      cpo_failures = n_fail,
                      waic = fit$waic$waic, p_waic = fit$waic$p.eff, dic = fit$dic$dic,
                      delay_T = zval("Transitional zone", "peak delay"), delay_N = zval("Northern limits", "peak delay"),
                      level_T = zval("Transitional zone", "level"), level_N = zval("Northern limits", "level"),
                      eta2 = zval("all", "eta2"),
                      fit_min = fit_min, sample_min = samp_min, samples = n_samp, int_strategy = int_for(m))
  rng <- hyp$med[grepl("^Range", hyp$parameter)]
  stats$range_km <- if (length(rng)) rng[1] else NA
  if (!is.null(ho_res)) stats <- cbind(stats, ho_res$metrics)

  # season component(s): intercept + RW2 curve (mean covariates, BL traps, other effects 0)
  sc <- bind_rows(lapply(seq_along(P$season), function(k) {
    C <- P$season[[k]] + matrix(P$beta[k, ], 52, ncol(P$beta), byrow = TRUE)   # z_k or intercept is row k
    data.frame(model = m, curve = if (s$zoned) zone_lab[k] else "shared", woy = 1:52, mean = rowMeans(C),
               lo = matrixStats::rowQuantiles(C, probs = 0.025), hi = matrixStats::rowQuantiles(C, probs = 0.975))
  }))
  zcv <- data.frame(key = rownames(pr$ZC), mean = rowMeans(pr$ZC),
                    lo = matrixStats::rowQuantiles(pr$ZC, probs = 0.025),
                    hi = matrixStats::rowQuantiles(pr$ZC, probs = 0.975)) %>%
    separate(key, c("zone_30y", "year", "woy"), sep = " ", convert = TRUE) %>% mutate(model = m)

  sink(file.path(mdir, "summary.txt")); print(summary(fit)); cat("\nInterpretable hyperparameters\n"); print(hyp, digits = 3); sink()
  write.csv(data.frame(term = rownames(fit$summary.fixed), fit$summary.fixed[, 1:5], row.names = NULL),
            file.path(mdir, "fixed.csv"), row.names = FALSE)
  write.csv(hyp, file.path(mdir, "hyper.csv"), row.names = FALSE)
  write.csv(sc, file.path(mdir, "season_component.csv"), row.names = FALSE)
  write.csv(zcv, file.path(mdir, "zone_curve.csv"), row.names = FALSE)
  write.csv(pk$site_year, file.path(mdir, "peak_site_year.csv"), row.names = FALSE)
  write.csv(pk$zone_year, file.path(mdir, "peak_zone_year.csv"), row.names = FALSE)
  write.csv(pk$zone_delay, file.path(mdir, "peak_zone_delay.csv"), row.names = FALSE)
  write.csv(pk$site, file.path(mdir, "peak_site.csv"), row.names = FALSE)
  saveRDS(list(draws = pr$PK, site_years = pr$sy), file.path(mdir, "peak_draws.rds"))
  saveRDS(pr$pred, file.path(mdir, "pred_rows.rds"))
  if (!is.null(ho_res)) {
    write.csv(ho_res$test_rows, file.path(mdir, "test_rows.csv"), row.names = FALSE)
    write.csv(ho_res$zone_week, file.path(mdir, "holdout_zone_week.csv"), row.names = FALSE)
  }
  plot_model_maps(s, fit, pk, mdir)

  res[[m]] <- list(stats = stats, hyper = hyp, season = sc, zone = zs,
                   peaks = pk[c("site_year", "zone_year", "zone_delay", "site")], zone_week = ho_res$zone_week,
                   test_lpd = if (!is.null(ho_res)) ho_res$test_rows$lpd_count else numeric(0))
  saveRDS(res[[m]], ck_file)
  print(stats, digits = 3)
  rm(fit, P, pr); invisible(gc())
}

# =============================================================================
# Across models -----------------------------------------------------------------------
# =============================================================================
mods <- names(res)
tab  <- bind_rows(lapply(res, `[[`, "stats"))
# Held-out count-scale score relative to the best model, paired over the same
# trap-weeks (se of the summed per-row difference)
if (length(test_idx)) {
  L    <- sapply(res, `[[`, "test_lpd")
  best <- which.max(colSums(L))
  tab$elpd_test_diff    <- colSums(L) - sum(L[, best])
  tab$elpd_test_diff_se <- apply(L, 2, function(l) sd(l - L[, best]) * sqrt(nrow(L)))
}
write.csv(tab, file.path(out_dir, "comparison.csv"), row.names = FALSE)
show <- intersect(c("model", "elpd_cpo_count", "elpd_test_count", "elpd_test_diff", "elpd_test_diff_se",
                    "R2_test", "mae_test", "cov90_test", "delay_T", "delay_N", "eta2", "range_km", "fit_min"), names(tab))
message("\nComparison (elpd_*_count: count scale, comparable across families; higher is better):")
print(tab[show], digits = 4)

zcon <- bind_rows(lapply(res, function(r) r$zone$contrasts))
write.csv(zcon, file.path(out_dir, "zone_contrasts.csv"), row.names = FALSE)

lab_zone <- function(z) factor(zone_lab[as.character(z)], zone_lab)
zcp <- zcon %>% mutate(model = factor(model, rev(mods)),
                       zones = setNames(tab$zones, tab$model)[as.character(model)],
                       panel = ifelse(zone == "all", "eta2 (site peaks)", paste(sub(" vs Southern.*", "", measure), "-", sub(" zone| limits", "", zone))))
p <- ggplot(zcp, aes(med, model, colour = zones)) +
  geom_vline(data = data.frame(panel = unique(zcp$panel[zcp$zone != "all"])), aes(xintercept = 0), colour = "grey60", linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25) +
  facet_wrap(~ panel, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = c(prescribed = "#1B9E77", free = "#D95F02"), name = "Zones") +
  labs(x = "Posterior median (95% CrI)", y = NULL,
       subtitle = "Zone differences after fitting: vs Southern range (peak delay in weeks; level on the linear-predictor scale)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "zone_contrasts.png"), p, width = 8, height = 1.8 + 0.25 * length(mods) * 2, dpi = 300)

zcu <- bind_rows(lapply(res, function(r) r$zone$curve)) %>% mutate(zone = lab_zone(zone_30y), model = factor(model, mods))
p <- ggplot(zcu, aes(woy, mean, colour = zone, fill = zone)) +
  annotate("rect", xmin = peak_win[["lo"]], xmax = peak_win[["hi"]], ymin = -Inf, ymax = Inf, alpha = 0.06) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.5) +
  facet_wrap(~ model, scales = "free_y", ncol = 3) + scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Week of year", y = "Zone mean linear predictor (mean over years)", colour = NULL, fill = NULL,
       subtitle = "Zone mean predicted curves (95% CrI); shaded = peak window. Gaussian: E[log(CEW + 1)]; NB / Poisson: log E[CEW]") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "zone_curves.png"), p, width = 8, height = 1 + 2.2 * ceiling(length(mods) / 3), dpi = 300)

zy <- bind_rows(lapply(mods, function(m) res[[m]]$peaks$zone_year %>% mutate(model = m))) %>%
  mutate(zone = factor(zone, zone_lab), model = factor(model, mods))
p <- ggplot(zy, aes(year, med, colour = zone)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.15) +
  facet_wrap(~ model, ncol = 2) + scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = yrs_all) +
  labs(x = NULL, y = "Mean peak week (95% CrI)", colour = NULL,
       subtitle = paste0("Zone-year peak week (mean over the zone's monitored site-years)",
                         if (!is.na(holdout)) paste0("; ", holdout, " held out") else "")) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(out_dir, "peak_zone_year.png"), p, width = 8, height = 1 + 2 * ceiling(length(mods) / 2), dpi = 300)

if (!is.na(holdout)) {
  zw <- bind_rows(lapply(res, `[[`, "zone_week")) %>% mutate(zone = lab_zone(zone_30y), model = factor(model, mods))
  p <- ggplot(zw, aes(woy)) +
    geom_ribbon(aes(ymin = pi_lo, ymax = pi_hi, fill = "pred"), alpha = 0.35) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = "cred"), alpha = 0.7) +
    geom_line(aes(y = pred_mean), colour = "#08306B", linewidth = 0.5) +
    geom_point(aes(y = obs), size = 0.5, colour = "grey15") +
    facet_grid(zone ~ model) +
    scale_fill_manual(values = c(cred = "#6BAED6", pred = "#C6DBEF"), breaks = c("cred", "pred"),
                      labels = c("95% credible interval (mean)", "95% predictive interval (observed mean)"), name = NULL) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(x = "Week of year", y = "Mean log(CEW + 1)",
         subtitle = paste0(holdout, " (held out): zone weekly means, observed (points) vs predicted")) +
    theme_bw(base_size = 8) + theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "holdout_ribbon.png"), p, width = max(5, 1.9 * length(mods) + 1), height = 6, dpi = 300)
}
message("Done. Outputs in ", out_dir)
