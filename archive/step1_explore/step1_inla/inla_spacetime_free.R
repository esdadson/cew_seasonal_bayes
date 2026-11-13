# =============================================================================
# step1_inla/inla_spacetime_free.R
# The zone-free space-time model in R-INLA (Gaussian on log(CEW + 1)) and the
# simpler models nested inside it, fitted one step at a time on a laptop:
#
#   fixed     : y = a + b_trap + x' beta
#   seasonal  : fixed + s(w) + u_t
#   spatial   : seasonal + v_i + g(s_i)
#   spacetime : seasonal + v_i + g(s_i, w)
#
#   y        log(CEW + 1) of site i, year t, ISO week w (1-52)
#   x        weekly EVI, weekly GDD, six crop shares (standardised)
#   s(w)     ONE seasonal curve shared by every site: a cyclic second-order
#            random walk over weeks 1-52 (week 52 wraps to week 1). Its
#            smoothness is estimated, so no harmonics have to be chosen.
#   u_t      year random intercept; the held-out year has no data, so its u
#            is drawn from N(0, sd_year^2)
#   v_i      site random intercept (non-spatial site differences)
#   g(s_i)   spatial Matern field (SPDE on a triangular mesh): sites closer
#            than the estimated range have similar levels
#   g(s_i,w) the same field, but allowed to change through the season: one
#            field every knot_wk weeks, AR1-correlated from knot to knot,
#            smoothly interpolated between knots. A field that is higher
#            early and lower late at a site moves that site's peak earlier,
#            so regional differences in TIMING come from here.
#
# Zones (zone_30y) are not in any model. They are used only after fitting,
# to summarise the predicted peak weeks by zone.
#
# Each section can be run on its own after the ones above it. Fits are saved
# in out_dir and reloaded on the next run (set refit <- TRUE to refit).
# Time on the default mesh (~1,150 nodes): fixed, seasonal, spatial take
# seconds; spacetime ~30 min. The finer mesh (max_edge 40, cutoff 15, ~3,900
# nodes) takes several hours for spacetime.
#
# INLA for R 4.4 on Windows (not on CRAN; see README.md):
#   utils::install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable",
#                           type = "binary", lib = .libPaths()[1], dependencies = FALSE)
#   renv::install(c("fmesher", "sn"))
# =============================================================================

library(INLA); library(fmesher); library(Matrix)
library(dplyr); library(tidyr); library(ggplot2)

# =============================================================================
# 0. Settings ------------------------------------------------------------------
# =============================================================================
out_dir    <- "results/inla_spacetime_free"
data_cache <- file.path(out_dir, "step1_data.rds")  # d_all saved after the first run
yr_test    <- 2020        # held out of the fit, predicted afterwards
n_samp     <- 1000         # joint posterior samples (1000 for final results)
n_threads  <- "8:1"       # INLA threads, "outer:inner"
max_edge   <- 80          # mesh triangle edge near the sites (km); 40 = fine mesh
cutoff     <- 30          # sites closer than this share a mesh node (km); 15 = fine mesh
knot_wk    <- 4           # spacetime: weeks between the field's knots
range0     <- 20          # field range prior: P(range < range0 km) = 0.1
refit      <- FALSE       # TRUE refits models that are already saved
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_lab <- c(`0` = "Northern limits", `1` = "Transitional zone", `2` = "Southern range")

# =============================================================================
# 1. Data ----------------------------------------------------------------------
#    d_all from step1/step1_data.R: the observed weeks plus every other week
#    (1-52) of each monitored site-year, with covariates for all of them, so
#    full seasonal curves and peak weeks can be predicted. Only observed,
#    non-held-out rows enter the fit.
# =============================================================================
if (file.exists(data_cache)) {
  d_all <- readRDS(data_cache)
} else {
  source("step1/step1_data.R")
  saveRDS(d_all, data_cache)
}

covs <- c("EVI_filled", "gdd", "corn", "cotton", "soybeans", "tobacco", "peanuts", "other")

# Sorted so every site-year is 52 consecutive rows (weeks 1-52):
# matrix(any_row_vector, 52) is then weeks x site-years
d <- d_all %>%
  arrange(location2, year, woy) %>%
  mutate(fit_row  = observed & year != yr_test,
         test_row = observed & year == yr_test,
         year_id  = match(year, sort(unique(year))),
         site_id  = match(location2, sort(unique(location2))))
stopifnot(all(table(d$loc_year) == 52), all(d$woy == rep(1:52, length.out = nrow(d))))

n_year <- max(d$year_id)
n_site <- max(d$site_id)
fit_idx  <- which(d$fit_row)
test_idx <- which(d$test_row)
message(sprintf("Training rows: %d; held-out %d rows: %d; sites: %d; site-years: %d",
                length(fit_idx), yr_test, length(test_idx), n_site, nrow(d) / 52))

# Fixed-effect design: intercept, trap type (BL_PH and PH against BL), covariates.
# INLA gets the columns directly, so every model uses the same X.
trap_X <- model.matrix(~ trap_type, d)[, -1, drop = FALSE]
colnames(trap_X) <- sub("^trap_type", "trap_", colnames(trap_X))
X <- cbind(intercept = 1, trap_X, as.matrix(d[, covs]))
colnames(X)

# One row per site, with its coordinates in CONUS Albers km (equal-distance
# units, unlike degrees of longitude)
site_xy <- d %>% distinct(site_id, location2, longitude, latitude, x_km, y_km, zone_30y) %>% arrange(site_id)
stopifnot(nrow(site_xy) == n_site)
loc <- as.matrix(site_xy[, c("x_km", "y_km")])

# =============================================================================
# 2. Mesh and SPDE (spatial and spacetime) ---------------------------------------
#    The field is defined at the mesh nodes, not at the sites. A_site
#    (sites x nodes) maps node values to sites: each row holds the 3 weights
#    of the triangle the site falls in. Field at the sites = A_site %*% g.
#    Fine triangles near the sites, coarse ones in a buffer so the field is
#    not distorted at the edge.
# =============================================================================
mesh <- fm_mesh_2d(loc = loc,
                   boundary = list(fm_nonconvex_hull(loc, convex = 100),
                                   fm_nonconvex_hull(loc, convex = 400)),
                   max.edge = c(max_edge, 300), cutoff = cutoff)
n_mesh <- mesh$n
A_site <- fm_basis(mesh, loc = loc)
message("Mesh nodes: ", n_mesh)

png(file.path(out_dir, "mesh.png"), width = 2000, height = 1500, res = 200)
plot(mesh, main = sprintf("Mesh: %d nodes", n_mesh)); points(loc, col = "red", pch = 20, cex = 0.4)
invisible(dev.off())

# Matern field with PC priors: P(range < range0 km) = 0.1, P(sd > 1) = 0.01
spde <- inla.spde2.pcmatern(mesh, prior.range = c(range0, 0.1), prior.sigma = c(1, 0.01))

# =============================================================================
# 3. Temporal basis (spacetime) ---------------------------------------------------
#    The space-time field has one copy of the spatial field per knot (every
#    knot_wk weeks). B_wk (52 weeks x knots) holds quadratic B-spline weights,
#    so the field in week w is a smooth mix of the neighbouring knots' fields:
#      g(site, w) = sum_k A_site[site, ] %*% g_k * B_wk[w, k]
#    A_st below combines both maps for every data row (row-wise Kronecker
#    product of the row's A_site row and its B_wk row).
# =============================================================================
knots  <- unique(c(seq(1, 52, by = knot_wk), 52))
tmesh  <- fm_mesh_1d(knots, degree = 2, boundary = "free")
B_wk   <- as.matrix(fm_basis(tmesh, loc = 1:52))
n_knot <- ncol(B_wk)
message("Space-time knots: ", n_knot, " (weeks ", paste(knots, collapse = ", "), ")")

png(file.path(out_dir, "temporal_basis.png"), width = 1600, height = 900, res = 200)
matplot(1:52, B_wk, type = "l", lty = 1, xlab = "Week", ylab = "Weight",
        main = "Temporal basis of the space-time field (one curve per knot)")
invisible(dev.off())

# Column (k - 1) * n_mesh + j of the result = As[, j] * Bt[, k]: INLA's order
# for a grouped effect (mesh node fastest, then knot)
row_kron <- function(As, Bt) {
  a <- as(As, "TsparseMatrix"); b <- as(as(Bt, "CsparseMatrix"), "TsparseMatrix")
  t <- merge(data.frame(i = a@i, j = a@j, x = a@x), data.frame(i = b@i, k = b@j, y = b@x), by = "i")
  sparseMatrix(i = t$i + 1, j = t$k * ncol(As) + t$j + 1, x = t$x * t$y,
               dims = c(nrow(As), ncol(As) * ncol(Bt)))
}

# =============================================================================
# 4. Priors and formulas ------------------------------------------------------------
#    Penalised-complexity (PC) priors shrink each term towards "absent":
#      seasonal curve: P(sd > 1) = 0.05 (scale.model: sd of the curve itself)
#      sd_year, sd_site: P(sd > 1) = 0.01
#      AR1 across knots: P(rho > 0.5) = 0.9, i.e. favours a field that
#      changes slowly through the season
#      fixed effects N(0, 1) (standardised covariates), intercept N(0, 10^2)
#    f(week, ...) is the curve s(w); constr = TRUE makes it sum to 0 so the
#    intercept keeps the level. values = 1:n_year includes the held-out year.
# =============================================================================
hp_season <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
hp_sd     <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
hp_rho    <- list(rho  = list(prior = "pc.cor1", param = c(0.5, 0.9)))

f_fixed <- as.formula(paste("y ~ 0 +", paste(colnames(X), collapse = " + ")))

f_seasonal <- update(f_fixed, . ~ . +
  f(week, model = "rw2", cyclic = TRUE, scale.model = TRUE, constr = TRUE,
    values = 1:52, hyper = hp_season) +
  f(year_id, model = "iid", values = 1:n_year, hyper = hp_sd))

f_spatial <- update(f_seasonal, . ~ . +
  f(site_id, model = "iid", values = 1:n_site, hyper = hp_sd) +
  f(field, model = spde))

f_spacetime <- update(f_seasonal, . ~ . +
  f(site_id, model = "iid", values = 1:n_site, hyper = hp_sd) +
  f(st, model = spde, group = st.group,
    control.group = list(model = "ar1", hyper = hp_rho)))

f_spacetime

# =============================================================================
# 5. Data for INLA ------------------------------------------------------------------
#    fixed and seasonal: an ordinary data frame (one row per observation).
#    spatial and spacetime: an inla.stack, because the field lives on the
#    mesh nodes and reaches the data through A. A stack pairs each block of
#    effects with its A matrix:
#      eta = A_field %*% g + 1 %*% (X beta + s(week) + u_year + v_site)
# =============================================================================
df_fit <- data.frame(y = d$log_CEW[fit_idx], X[fit_idx, ],
                     week = d$woy[fit_idx], year_id = d$year_id[fit_idx], site_id = d$site_id[fit_idx])

stk_spatial <- inla.stack(
  data    = list(y = df_fit$y),
  A       = list(A_site[d$site_id[fit_idx], ], 1),
  effects = list(list(field = 1:n_mesh), df_fit[, -1]),
  tag = "fit", remove.unused = FALSE)

A_st <- row_kron(A_site[d$site_id[fit_idx], ], B_wk[d$woy[fit_idx], ])
stk_spacetime <- inla.stack(
  data    = list(y = df_fit$y),
  A       = list(A_st, 1),
  effects = list(list(st = rep(1:n_mesh, n_knot), st.group = rep(1:n_knot, each = n_mesh)),
                 df_fit[, -1]),
  tag = "fit", remove.unused = FALSE)
dim(A_st)   # rows x (mesh nodes * knots)

# =============================================================================
# 6. Fit -------------------------------------------------------------------------
#    config = TRUE keeps what inla.posterior.sample() needs (section 8);
#    dic / waic / cpo for model comparison. int.strategy "eb" fixes the
#    hyperparameters at their posterior mode (much faster; used for
#    spacetime), "auto" integrates over them.
# =============================================================================
fit_model <- function(name, formula, data, A = NULL, int = "auto") {
  f_rds <- file.path(out_dir, paste0("fit_", name, ".rds"))
  if (!refit && file.exists(f_rds)) { message("Loading ", f_rds); return(readRDS(f_rds)) }
  t0  <- Sys.time()
  fit <- inla(formula, family = "gaussian", data = data,
              control.predictor = if (is.null(A)) list(compute = FALSE) else list(A = A, compute = FALSE),
              control.compute   = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE,
                                       return.marginals.predictor = FALSE),
              control.fixed     = list(mean = 0, prec = list(default = 1, intercept = 0.01)),
              control.inla      = list(int.strategy = int),
              num.threads = n_threads)
  fit$fit_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("%s fitted in %.1f min", name, fit$fit_min))
  saveRDS(fit, f_rds)
  fit
}

fits <- list()
fits$fixed     <- fit_model("fixed",    f_fixed,    df_fit)
fits$seasonal  <- fit_model("seasonal", f_seasonal, df_fit)
fits$spatial   <- fit_model("spatial",  f_spatial,  inla.stack.data(stk_spatial, spde = spde),
                            A = inla.stack.A(stk_spatial))
fits$spacetime <- fit_model("spacetime", f_spacetime, inla.stack.data(stk_spacetime, spde = spde),
                            A = inla.stack.A(stk_spacetime), int = "eb")

# =============================================================================
# 7. What each model estimated --------------------------------------------------------
# =============================================================================
summary(fits$spacetime)

# Fixed effects side by side (posterior means)
sapply(fits, function(f) round(f$summary.fixed[colnames(X), "mean"], 3)) %>%
  `rownames<-`(colnames(X))

# Hyperparameters on the sd scale (INLA reports precisions = 1 / sd^2)
hyper_sd <- function(fit) {
  h <- fit$summary.hyperpar[, c("0.5quant", "0.025quant", "0.975quant")]
  colnames(h) <- c("median", "lo", "hi")
  is_prec <- grepl("^Precision", rownames(h))
  h[is_prec, ] <- 1 / sqrt(h[is_prec, c(1, 3, 2)])   # sd = 1 / sqrt(precision): the bounds swap
  rownames(h) <- sub("^Precision for (the )?", "SD for ", rownames(h))
  round(h, 3)
}
lapply(fits, hyper_sd)
# Range for field / st: distance (km) at which the field's correlation drops
# to ~0.13; Stdev for field / st: its sd on the log(CEW + 1) scale;
# GroupRho for st: correlation of the field between neighbouring knots.

# The shared seasonal curve s(w) in each model
curves <- bind_rows(lapply(names(fits)[-1], function(m)
  fits[[m]]$summary.random$week %>%
    transmute(model = m, week = ID, mean, lo = `0.025quant`, hi = `0.975quant`)))
p <- ggplot(curves, aes(week, mean, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) + geom_line() +
  labs(x = "Week of year", y = "s(w) (log(CEW + 1) scale)", colour = NULL, fill = NULL,
       subtitle = "Shared seasonal curve (95% CrI)") +
  theme_bw(base_size = 9)
ggsave(file.path(out_dir, "season_curves.png"), p, width = 6, height = 3.5, dpi = 300)

# Maps of the fields (posterior means) on a grid, masked to within 75 km of a site
outline <- tryCatch({
  ol <- sf::st_transform(sf::st_as_sf(terra::unwrap(readRDS("gadm/gadm41_USA_0_pk.rds"))), 5070)
  ol <- sf::st_simplify(suppressWarnings(sf::st_crop(ol, sf::st_bbox(c(
          xmin = (min(loc[, 1]) - 300) * 1000, ymin = (min(loc[, 2]) - 300) * 1000,
          xmax = (max(loc[, 1]) + 300) * 1000, ymax = (max(loc[, 2]) + 300) * 1000), crs = sf::st_crs(5070)))),
        dTolerance = 3000)
  as.data.frame(sf::st_coordinates(ol)) %>% mutate(x_km = X / 1000, y_km = Y / 1000, grp = paste(L1, L2, L3))
}, error = function(e) NULL)

map_base <- function() {
  p <- ggplot() + coord_equal(xlim = range(loc[, 1]) + c(-60, 60), ylim = range(loc[, 2]) + c(-60, 60)) +
    theme_bw(base_size = 9) + theme(axis.title = element_blank(), axis.text = element_blank(),
                                    axis.ticks = element_blank(), panel.grid = element_blank())
  if (!is.null(outline)) p <- p + geom_polygon(data = outline, aes(x_km, y_km, group = grp),
                                               fill = "grey97", colour = "grey55", linewidth = 0.25)
  p
}
grid_proj <- fm_evaluator(mesh, xlim = range(loc[, 1]) + c(-60, 60), ylim = range(loc[, 2]) + c(-60, 60),
                          dims = c(260, 200))
grid_xy <- expand.grid(x_km = grid_proj$x, y_km = grid_proj$y)
near    <- apply(as.matrix(grid_xy), 1, function(g) min(sqrt((loc[, 1] - g[1])^2 + (loc[, 2] - g[2])^2))) <= 75
field_on_grid <- function(node_values) mutate(grid_xy, value = as.vector(fm_evaluate(grid_proj, field = node_values)))[near, ]

p <- map_base() + geom_tile(data = field_on_grid(fits$spatial$summary.random$field$mean),
                            aes(x_km, y_km, fill = value)) +
  geom_point(data = site_xy, aes(x_km, y_km), size = 0.3) + scale_fill_gradient2(name = "g(s)") +
  labs(subtitle = "spatial: posterior mean of the field (time-constant)")
ggsave(file.path(out_dir, "field_spatial.png"), p, width = 6.5, height = 4.5, dpi = 300)

# Space-time field at four weeks: node values per knot (n_mesh x n_knot) times
# the week's B-spline weights
st_mean <- matrix(fits$spacetime$summary.random$st$mean, n_mesh, n_knot)
wks <- c(20, 22, 24, 26, 28, 30, 32, 34)
g <- bind_rows(lapply(wks, function(w) field_on_grid(as.vector(st_mean %*% B_wk[w, ])) %>% mutate(week = w)))
p <- map_base() + geom_tile(data = g, aes(x_km, y_km, fill = value)) +
  facet_wrap(~ week, ncol = 2, labeller = label_both) + scale_fill_gradient2(name = "g(s, w)") +
  labs(subtitle = "spacetime: posterior mean of the field through the season")
ggsave(file.path(out_dir, "field_spacetime.png"), p, width = 7, height = 6, dpi = 300)

# =============================================================================
# 8. Posterior samples -> predictions and peak weeks --------------------------------
#    inla.posterior.sample() draws the latent effects (beta, s, u, v, field)
#    jointly, with the hyperparameters. From each draw b, the linear
#    predictor of EVERY row of d (all weeks of every site-year, held-out
#    year included) is rebuilt by hand:
#      eta = X beta + s[week] + u[year] + v[site] + (A_site g)[site]
#      spacetime: + (A_site G B_wk')[site, week], G = nodes x knots
#    The peak week of a site-year in draw b is the week of its largest eta,
#    searched within one window for all sites (5th-95th percentile of the
#    sampled weeks), so no zone information enters.
# =============================================================================
peak_win <- c(20, 40)
in_win   <- 1:52 >= floor(peak_win[1]) & 1:52 <= ceiling(peak_win[2])
message("Peak window: weeks ", paste(range(which(in_win)), collapse = "-"))

draw_samples <- function(fit, n) {
  rnd <- names(fit$summary.random)
  sel <- c(setNames(as.list(rep(1, ncol(X))), colnames(X)), setNames(as.list(rep(0, length(rnd))), rnd))
  smp <- inla.posterior.sample(n, fit, selection = sel, num.threads = n_threads)
  rn  <- rownames(smp[[1]]$latent)
  grab <- function(p) sapply(smp, function(s) s$latent[startsWith(rn, paste0(p, ":"))])
  P <- list(beta = do.call(rbind, lapply(colnames(X), function(v) grab(v))))   # ncol(X) x n
  for (r in rnd) P[[r]] <- matrix(grab(r), ncol = n)                           # levels x n
  P$sigma <- sapply(smp, function(s) 1 / sqrt(s$hyperpar[["Precision for the Gaussian observations"]]))
  P
}

eta_draw <- function(P, b) {
  e <- as.vector(X %*% P$beta[, b])
  if (!is.null(P$week))    e <- e + P$week[d$woy, b]
  if (!is.null(P$year_id)) e <- e + P$year_id[d$year_id, b]
  if (!is.null(P$site_id)) e <- e + P$site_id[d$site_id, b]
  if (!is.null(P$field))   e <- e + as.vector(A_site %*% P$field[, b])[d$site_id]
  if (!is.null(P$st)) {
    G <- as.matrix(A_site %*% matrix(P$st[, b], n_mesh, n_knot) %*% t(B_wk))  # sites x 52 weeks
    e <- e + G[cbind(d$site_id, d$woy)]
  }
  e
}

sy <- d[seq(1, nrow(d), by = 52), c("loc_year", "location2", "site_id", "year", "zone_30y")]
post <- list()
for (m in names(fits)) {
  t0 <- Sys.time()
  P  <- draw_samples(fits[[m]], n_samp)
  PK <- matrix(NA_integer_, n_samp, nrow(sy))        # peak week: draws x site-years
  ET <- matrix(NA_real_, length(test_idx), n_samp)   # eta of the held-out rows
  e_sum <- e_sq <- numeric(nrow(d))
  for (b in 1:n_samp) {
    e <- eta_draw(P, b)
    e_sum <- e_sum + e; e_sq <- e_sq + e^2
    ET[, b] <- e[test_idx]
    M <- matrix(e, 52); M[!in_win, ] <- -Inf          # weeks x site-years
    PK[b, ] <- max.col(t(M), ties.method = "first")
  }
  eta_mean <- e_sum / n_samp
  post[[m]] <- list(PK = PK, ET = ET, sigma = P$sigma,
                    eta = data.frame(eta_mean = eta_mean, eta_sd = sqrt(pmax(e_sq / n_samp - eta_mean^2, 0))))
  message(sprintf("%s: %d samples in %.1f min", m, n_samp, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(list(post = post, sy = sy, peak_window = range(which(in_win))), file.path(out_dir, "posterior_draws.rds"))

# =============================================================================
# 9. Model comparison ------------------------------------------------------------------
#    In-sample: DIC, WAIC (lower is better), LPML = sum of log CPO (leave-one-
#    out; higher is better; CPO failures are flagged by INLA).
#    Held-out year: elpd = sum over its rows of log p(y | training data),
#    p averaged over the draws (y ~ N(eta, sigma^2)); elpd_diff and its se
#    are paired against the best model. R2 / MAE / coverage of the 90%
#    predictive interval on log(CEW + 1).
# =============================================================================
y_test <- d$log_CEW[test_idx]
lpd <- sapply(post, function(p) {
  LP <- matrix(dnorm(y_test, p$ET, rep(p$sigma, each = length(y_test)), log = TRUE),
               length(y_test))                                                  # rows x draws
  matrixStats::rowLogSumExps(LP) - log(n_samp)
})
best <- which.max(colSums(lpd))

comparison <- bind_rows(lapply(names(fits), function(m) {
  f <- fits[[m]]; p <- post[[m]]
  mu   <- rowMeans(p$ET)
  yrep <- p$ET + matrix(rnorm(length(p$ET)), nrow(p$ET)) * rep(p$sigma, each = nrow(p$ET))
  q    <- matrixStats::rowQuantiles(yrep, probs = c(0.05, 0.95))
  dlp  <- lpd[, m] - lpd[, best]
  h    <- f$summary.hyperpar
  data.frame(model = m, fit_min = round(f$fit_min, 2),
             dic = f$dic$dic, waic = f$waic$waic,
             lpml = sum(log(f$cpo$cpo), na.rm = TRUE), cpo_failures = sum(f$cpo$failure > 0, na.rm = TRUE),
             elpd_test = sum(lpd[, m]), elpd_diff = sum(dlp), elpd_diff_se = sqrt(length(dlp)) * sd(dlp),
             R2_test = 1 - sum((y_test - mu)^2) / sum((y_test - mean(y_test))^2),
             mae_test = mean(abs(y_test - mu)),
             cov90_test = mean(y_test >= q[, 1] & y_test <= q[, 2]),
             range_km = if (any(grepl("^Range", rownames(h)))) h[grep("^Range", rownames(h)), "mean"] else NA)
}))
print(comparison, digits = 4)
write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)

# Held-out year: predicted vs observed for a few sites per zone (spacetime)
set.seed(1)
ex_sites <- d %>% filter(test_row) %>% distinct(site_id, zone_30y) %>%
  group_by(zone_30y) %>% slice_sample(n = 3) %>% pull(site_id)
ex <- bind_rows(lapply(names(fits)[-1], function(m)
  d %>% mutate(model = m, post[[m]]$eta) %>% filter(year == yr_test, site_id %in% ex_sites))) %>%
  mutate(site = paste0(location2, " (", sub(" .*", "", zone_lab[as.character(zone_30y)]), ")"),
         model = factor(model, names(fits)))
p <- ggplot(ex, aes(woy)) +
  geom_ribbon(aes(ymin = eta_mean - 2 * eta_sd, ymax = eta_mean + 2 * eta_sd, fill = model), alpha = 0.2) +
  geom_line(aes(y = eta_mean, colour = model)) +
  geom_point(data = filter(ex, observed, model == "spacetime"), aes(y = log_CEW), size = 0.8) +
  facet_wrap(~ site, ncol = 3) +
  labs(x = "Week of year", y = "log(CEW + 1)", colour = NULL, fill = NULL,
       subtitle = paste0(yr_test, " (held out): mean +/- 2 sd of the linear predictor; points = observed")) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "holdout_sites.png"), p, width = 8, height = 6, dpi = 300)

# =============================================================================
# 10. Peak weeks and zones (after fitting) ---------------------------------------------
#    Everything per draw, then summarised:
#      site-year: posterior median and 95% interval of the peak week
#      zone delay: mean peak over the zone's site-years minus Southern's
#      eta2: share of the variance of site peak weeks (averaged over a
#            site's years) that lies between zones
#    fixed has no seasonal curve, so its "peak" only follows the covariates.
# =============================================================================
zone_of_sy <- sy$zone_30y
Ws <- sparseMatrix(seq_len(nrow(sy)), sy$site_id, x = 1)          # site-years -> sites
Ws <- Ws %*% Diagonal(x = 1 / colSums(Ws))                        # column means
zone_of_site <- site_xy$zone_30y

zone_stats <- bind_rows(lapply(names(fits), function(m) {
  PK <- post[[m]]$PK
  zmean <- sapply(c(`2` = 2, `1` = 1, `0` = 0), function(z) rowMeans(PK[, zone_of_sy == z, drop = FALSE]))
  SP <- as.matrix(PK %*% Ws)                                       # draws x sites
  eta2 <- apply(SP, 1, function(x) {
    zm <- tapply(x, zone_of_site, mean)[as.character(zone_of_site)]
    sum((zm - mean(x))^2) / sum((x - mean(x))^2)
  })
  q <- function(v) sprintf("%.2f (%.2f, %.2f)", median(v), quantile(v, 0.025), quantile(v, 0.975))
  data.frame(model = m,
             peak_Southern = q(zmean[, "2"]), peak_Transitional = q(zmean[, "1"]), peak_Northern = q(zmean[, "0"]),
             delay_Transitional = q(zmean[, "1"] - zmean[, "2"]), delay_Northern = q(zmean[, "0"] - zmean[, "2"]),
             eta2 = q(eta2))
}))
print(zone_stats)
write.csv(zone_stats, file.path(out_dir, "zone_stats.csv"), row.names = FALSE)

# Site-year peak summaries (spacetime)
peak_site_year <- sy %>%
  mutate(peak_med = apply(post$spacetime$PK, 2, median),
         peak_lo  = apply(post$spacetime$PK, 2, quantile, 0.025),
         peak_hi  = apply(post$spacetime$PK, 2, quantile, 0.975))
write.csv(peak_site_year, file.path(out_dir, "peak_site_year_spacetime.csv"), row.names = FALSE)

# Histogram of all peak-week draws by zone (density within zone, so zones
# with more site-years don't look taller)
pk_long <- bind_rows(lapply(names(fits)[-1], function(m)
  data.frame(model = m, zone_30y = rep(zone_of_sy, each = n_samp), week = as.vector(post[[m]]$PK)))) %>%
  mutate(zone = factor(zone_lab[as.character(zone_30y)], zone_lab), model = factor(model, names(fits)))
p <- ggplot(pk_long, aes(week, fill = zone, colour = zone)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                 position = "identity", alpha = 0.35, linewidth = 0.3) +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_brewer(palette = "Dark2") + scale_colour_brewer(palette = "Dark2") +
  labs(x = "Peak week", y = "Density", fill = NULL, colour = NULL,
       subtitle = "Posterior peak-week draws, all site-years, by zone (zones assigned after fitting)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "peak_hist.png"), p, width = 7, height = 7, dpi = 300)

# Map of each site's peak week (mean over its years; spacetime)
SP <- as.matrix(post$spacetime$PK %*% Ws)
site_pk <- site_xy %>% mutate(med = apply(SP, 2, median),
                              width = apply(SP, 2, quantile, 0.975) - apply(SP, 2, quantile, 0.025))
p <- map_base() + geom_point(data = site_pk, aes(x_km, y_km, colour = med), size = 1.3) +
  scale_colour_viridis_c(name = "Peak week\n(median)") +
  labs(subtitle = "spacetime: site peak week, mean over monitored years")
ggsave(file.path(out_dir, "peak_map_spacetime.png"), p, width = 6.5, height = 4.5, dpi = 300)

message("Done. Outputs in ", out_dir)
