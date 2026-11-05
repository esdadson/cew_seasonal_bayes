# =============================================================================
# step1/step1_krige_test.R
# Trial, before changing step1_models.R: the spatial term as an exact GP over
# the monitored sites with a FIXED range, fit as one site random effect, and
# kriging of that effect to unmonitored locations.
#
#   u = site effects ~ N(0, sd^2 K),  K = rho C + (1 - rho) I,
#   C_ij = exp(-d_ij / phi), d = CONUS Albers km between sites
#   rho (spatial share) and phi (range) are estimated once, from the
#   correlation of site-mean residuals of the fixed-effects trend by distance
#   (section 1), and held fixed: the model adds one parameter (sd), the same
#   as an independent site intercept, and samples like one - no GP length-
#   scale or basis functions. New site s0: u0 | u ~ N(rho c0' K^-1 u,
#   sd^2 (1 - rho^2 c0' K^-1 c0)), c0 = its correlations with the sites.
#
# Models (step1_models.R's fixed terms + ...):
#   fixed     : as in step1_models.R (no site effect)
#   site_iid  : + (1 | location2)                independent site intercepts
#   spatialK  : + (1 | gr(location2, cov = K))   the fixed-range GP
#   fourierK  : spatialK + Fourier curve + zone deviations (as step1 'fourier')
#   spatialK_phi<x> : spatialK with phi = x km (KT_PHI_SENS, range sensitivity)
#
# Checks, per model: sampler efficiency, PSIS-LOO, R^2/MAE and the 2020
# hold-out (step1_metrics), the residual correlogram after the fit, and a
# leave-blocks-of-sites-out check of kriging (section 5); a map of the kriged
# spatial effect (section 6). Writes to results/step1_krige_test/ only.
#
# Run from the project directory, on the laptop or inside a Slurm job:
#   Rscript step1/step1_krige_test.R
# Environment (defaults):
#   KT_MODELS   "fixed site_iid spatialK fourierK"
#   KT_CHAINS 4, KT_ITER 2000, KT_WARMUP 1000, KT_THIN 2 (2000 draws), KT_THREADS 1
#   KT_RHO, KT_PHI   override the estimated spatial share and range (km)
#   KT_PHI_SENS      extra spatialK fits at these ranges, e.g. "25 90"
#   KT_BLOCK_KM 100  block size of the leave-sites-out check
#   KT_CHECK_DRAWS 200   posterior draws used by the kriging checks and map
#   KT_OUT      results/step1_krige_test
# =============================================================================

source("step1/step1_models.R")
library(ggplot2)

kt_dir <- Sys.getenv("KT_OUT", "results/step1_krige_test")
dir.create(kt_dir, recursive = TRUE, showWarnings = FALSE)
env_int <- function(name, default) as.integer(Sys.getenv(name, default))
env_num <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) as.numeric(v) else default }
env_vec <- function(name, default = "") { v <- Sys.getenv(name, default); if (nzchar(v)) strsplit(trimws(v), " +")[[1]] else character() }
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

kt <- list(chains = env_int("KT_CHAINS", 4), iter = env_int("KT_ITER", 2000), warmup = env_int("KT_WARMUP", 1000),
           thin = env_int("KT_THIN", 2), threads = env_int("KT_THREADS", 1))
phi_sens  <- as.numeric(env_vec("KT_PHI_SENS"))
sens_names <- if (length(phi_sens)) paste0("spatialK_phi", phi_sens) else character()
kt_models <- c(env_vec("KT_MODELS", "fixed site_iid spatialK fourierK"), sens_names)
block_km  <- env_num("KT_BLOCK_KM", 100)
n_check   <- env_int("KT_CHECK_DRAWS", 200)
message(sprintf("Training: %d rows, %d sites; held out (%s): %d rows", nrow(d_fit), n_distinct(d_fit$location2),
                holdout, nrow(d_test)))

# =============================================================================
# 1. Spatial share and range ----------------------------------------------------
#    Trend = the fixed model's terms with years as fixed effects and a
#    seasonal curve per zone, by OLS on the training years. Site-mean
#    residuals, standardised; their products averaged by distance bin give the
#    correlogram; rho exp(-h / phi) fit by weighted least squares (weights =
#    pairs per bin).
# =============================================================================
site_xy <- d_fit %>% distinct(location2, x_km, y_km) %>% arrange(location2)
XY <- as.matrix(site_xy[, c("x_km", "y_km")])
dist_km  <- function(a, b) sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
corr_exp <- function(D, phi) exp(-D / phi)

correlogram <- function(r, d, breaks = c(0, 5, 10, 15, 20, 30, 40, 50, 65, 80, 100, 130, 160, 200, 300)) {
  sm <- d %>% mutate(r = r) %>% group_by(location2) %>% summarise(r = mean(r), x = first(x_km), y = first(y_km))
  D  <- dist_km(as.matrix(sm[, c("x", "y")]), as.matrix(sm[, c("x", "y")]))
  z  <- (sm$r - mean(sm$r)) / sd(sm$r)
  ut <- upper.tri(D)
  tibble(dist = D[ut], prod = outer(z, z)[ut]) %>% filter(dist > 0, dist <= max(breaks)) %>%
    mutate(bin = cut(dist, breaks)) %>% group_by(bin) %>%
    summarise(h = mean(dist), pairs = n(), cor = mean(prod), .groups = "drop")
}

trend <- lm(log_CEW ~ 0 + zone_fac + trap_type + year_f + EVI_filled + gdd + corn + cotton + soybeans +
              tobacco + peanuts + other + zone_fac:(f1 + f2 + f3 + f4), data = d_fit)
cg_trend <- correlogram(resid(trend), d_fit)
cfit <- nls(cor ~ rho * exp(-h / phi), data = cg_trend, weights = pairs, start = list(rho = 0.8, phi = 40))
rho <- env_num("KT_RHO", coef(cfit)[["rho"]])
phi <- env_num("KT_PHI", coef(cfit)[["phi"]])
message(sprintf("Trend R2 %.3f. Correlogram fit: rho = %.2f, phi = %.0f km%s", summary(trend)$r.squared,
                coef(cfit)[["rho"]], coef(cfit)[["phi"]],
                if (nzchar(Sys.getenv("KT_RHO")) || nzchar(Sys.getenv("KT_PHI")))
                  sprintf(" (using rho = %.2f, phi = %.0f km)", rho, phi) else ""))

make_K <- function(phi_km) {
  K <- rho * corr_exp(dist_km(XY, XY), phi_km) + (1 - rho) * diag(nrow(XY))
  dimnames(K) <- list(site_xy$location2, site_xy$location2)
  K
}

# =============================================================================
# 2. Models ---------------------------------------------------------------------
#    Added to step1_terms in this session only, so step1_formula(),
#    step1_priors(), re_no_year() and step1_metrics() work on them unchanged.
#    Sampler: step1_control + a dense mass matrix (6x the effective draws per
#    gradient for site intercepts in the earlier benchmark).
# =============================================================================
site_term  <- "(1 | gr(location2, cov = K))"
four_terms <- c("f1 + f2 + f3 + f4", "(0 + f1 + f2 + f3 + f4 || zone_fac)")
kt_terms <- list(fixed    = terms_fixed,
                 site_iid = c(terms_fixed, "(1 | location2)"),
                 spatialK = c(terms_fixed, site_term),
                 fourierK = c(terms_fixed, site_term, four_terms))
for (s in sens_names) kt_terms[[s]] <- kt_terms$spatialK
kt_phi <- setNames(rep(phi, length(kt_terms)), names(kt_terms))
if (length(phi_sens)) kt_phi[sens_names] <- phi_sens
stopifnot(all(kt_models %in% names(kt_terms)))
step1_terms[names(kt_terms)] <- kt_terms
has_site_K <- function(model) site_term %in% step1_terms[[model]]
kt_control <- c(step1_control, list(metric = "dense_e"))

kt_fit <- function(model) {
  brm(step1_formula(model), data = droplevels(d_fit), data2 = list(K = make_K(kt_phi[[model]])),
      family = gaussian(), prior = step1_priors(model),
      chains = kt$chains, cores = kt$chains, iter = kt$iter, warmup = kt$warmup, thin = kt$thin,
      threads = if (kt$threads > 1) threading(kt$threads) else NULL,
      seed = 42, init = 0, control = kt_control, refresh = max(kt$iter %/% 10, 1),
      file = file.path(kt_dir, paste0("fit_", model)), file_refit = "on_change")
}

# Sampling time and efficiency (post-warmup); with thinning, the kept draws'
# step counts stand in for all sampling iterations
sampler_stats <- function(fit, model, wall_min) {
  np <- nuts_params(fit)
  v  <- function(p) np$Value[np$Parameter == p]
  et <- rstan::get_elapsed_time(fit$fit)                      # chains x (warmup, sample), seconds
  steps <- mean(v("n_leapfrog__"))
  sm <- summarise_draws(as_draws_df(fit), "rhat", "ess_bulk", "ess_tail")
  data.frame(model = model, wall_min = wall_min,
             warmup_min = mean(et[, 1]) / 60, sampling_min = mean(et[, 2]) / 60,
             steps_per_iter = steps, pct_max_treedepth = 100 * mean(v("treedepth__") >= step1_control$max_treedepth),
             divergent = sum(v("divergent__")),
             ms_per_grad = 1000 * mean(et[, 2]) / ((kt$iter - kt$warmup) * steps),
             max_rhat = max(sm$rhat, na.rm = TRUE), min_ess_bulk = min(sm$ess_bulk, na.rm = TRUE),
             min_ess_tail = min(sm$ess_tail, na.rm = TRUE), worst = sm$variable[which.min(sm$ess_bulk)])
}

# =============================================================================
# 3. Kriging the site effect --------------------------------------------------------
#    u: draws x sites (columns in the order of xy_obs); sdu: draws of the site
#    sd. Returns draws x new sites, and the conditional variance share v
#    (0 = the effect is known, 1 = nothing nearby: N(0, sd^2)).
# =============================================================================
site_draws <- function(fit, draw_ids) {
  u  <- as_draws_matrix(fit, variable = "^r_location2\\[.*,Intercept\\]$", regex = TRUE)
  id <- as.integer(sub("^r_location2\\[(\\d+),.*", "\\1", colnames(u)))
  u  <- unclass(u)[draw_ids, match(site_xy$location2, id), drop = FALSE]
  colnames(u) <- site_xy$location2
  list(u = u, sd = as_draws_df(fit, variable = "sd_location2__Intercept")$sd_location2__Intercept[draw_ids])
}
krige_site <- function(u, sdu, xy_obs, xy_new, phi_km, seed = 42) {
  Kinv <- chol2inv(chol(rho * corr_exp(dist_km(xy_obs, xy_obs), phi_km) + (1 - rho) * diag(nrow(xy_obs))))
  c0 <- corr_exp(dist_km(xy_new, xy_obs), phi_km)
  A  <- rho * c0 %*% Kinv                                     # kriging weights
  v  <- pmax(1 - rowSums(A * (rho * c0)), 0)
  mu <- u %*% t(A)                                            # conditional mean, per draw
  set.seed(seed)
  list(draws = mu + outer(sdu, sqrt(v)) * matrix(rnorm(nrow(u) * nrow(xy_new)), nrow(u)), v = v,
       mean = colMeans(mu), sd = sqrt(matrixStats::colVars(mu) + mean(sdu^2) * v))   # exact moments
}

# E[log(CEW + 1)] (draws x rows) at sites whose effect is supplied (u_site:
# draws x sites, colnames = location2): every group-level term but the site
# effect from the fit, plus u_site; unseen years get a new-year draw
# (u_year[s] added to draw s), as epred_rows() does
predict_new_sites <- function(fit, model, nd, u_site, draw_ids, u_year = new_year_draws(fit)[draw_ids]) {
  re   <- setdiff(grep("|", step1_terms[[model]], fixed = TRUE, value = TRUE), c(site_term, "(1 | location2)"))
  re_f <- function(t) if (length(t)) as.formula(paste("~", paste(t, collapse = " + "))) else NA
  seen <- nd$year_f %in% levels(fit$data$year_f)
  out  <- matrix(NA_real_, length(draw_ids), nrow(nd))
  if (any(seen)) out[, seen] <- posterior_epred(fit, newdata = droplevels(nd[seen, ]), re_formula = re_f(re),
                                                draw_ids = draw_ids)
  if (any(!seen)) out[, !seen] <- posterior_epred(fit, newdata = droplevels(nd[!seen, ]),
                                                  re_formula = re_f(setdiff(re, year_term)), draw_ids = draw_ids,
                                                  allow_new_levels = TRUE) + u_year
  out + u_site[, match(as.character(nd$location2), colnames(u_site)), drop = FALSE]
}

score <- function(y, ep, sig) {
  m  <- colMeans(ep)
  pi <- matrixStats::colQuantiles(ep + sig * matrix(rnorm(length(ep)), nrow(ep)), probs = c(0.05, 0.95))
  data.frame(mae = mean(abs(y - m)), rmse = sqrt(mean((y - m)^2)),
             R2 = 1 - sum((y - m)^2) / sum((y - mean(y))^2), cov90 = mean(y >= pi[, 1] & y <= pi[, 2]))
}

# =============================================================================
# 4. Fit, LOO, metrics ------------------------------------------------------------
# =============================================================================
# Each model's results are saved to check_<model>.rds and reused while newer
# than its fit, so an interrupted run resumes where it stopped. One fit is in
# memory at a time: LOO and the metrics hold several draws x 36,000-row
# matrices, and an R session that grows too large starves the next model's
# chains (they then stop with no error message)
fit_file   <- function(m) file.path(kt_dir, paste0("fit_", m, ".rds"))
check_file <- function(m) file.path(kt_dir, paste0("check_", m, ".rds"))
stats <- list(); mets <- list(); los <- list(); cgs <- list(trend = cg_trend)
for (m in kt_models) {
  message("\n==== ", m, " ====")
  if (file.exists(check_file(m)) && file.exists(fit_file(m)) && file.mtime(check_file(m)) > file.mtime(fit_file(m))) {
    ck <- readRDS(check_file(m))
    message("Reusing ", check_file(m))
  } else {
    print(step1_formula(m))
    t0  <- Sys.time()
    fit <- kt_fit(m)
    st  <- sampler_stats(fit, m, mins(t0))
    print(st, digits = 3, row.names = FALSE)
    lo <- loo(fit, save_psis = TRUE)
    set.seed(42)
    mt <- step1_metrics(fit, m, lo, d_test)
    lo$psis_object <- NULL
    saveRDS(mt, file.path(kt_dir, paste0("metrics_", m, ".rds")))
    # Residual correlogram after the fit (posterior mean over n_check draws)
    ids <- sort(sample(ndraws(fit), min(n_check, ndraws(fit))))
    cg  <- correlogram(d_fit$log_CEW - colMeans(posterior_epred(fit, draw_ids = ids)), d_fit)
    ck  <- list(stats = st, loo = lo, metrics = mt$metrics, correlogram = cg)
    saveRDS(ck, check_file(m))
    rm(fit, lo, mt); invisible(gc())
  }
  stats[[m]] <- ck$stats; los[[m]] <- ck$loo; mets[[m]] <- ck$metrics; cgs[[m]] <- ck$correlogram
}

cmp <- bind_rows(mets)
if (length(los) >= 2) {
  lc  <- loo::loo_compare(los)
  cmp <- cmp %>% left_join(data.frame(model = rownames(lc), elpd_diff = lc[, "elpd_diff"], se_diff = lc[, "se_diff"]),
                           by = "model")
}
cmp <- cmp %>% left_join(bind_rows(stats), by = "model") %>%
  mutate(rho = ifelse(sapply(model, has_site_K), rho, NA), phi_km = ifelse(sapply(model, has_site_K), kt_phi[model], NA))
write.csv(cmp, file.path(kt_dir, "comparison.csv"), row.names = FALSE)
message("\n==== Comparison (training years: LOO; held-out ", holdout, ": *_test) ====")
print(cmp %>% select(model, any_of(c("elpd_diff", "se_diff")), p_loo, loo_R2, mae_loo,
                     any_of(c("elpd_test", "mae_test", "R2_test", "cov90_test")),
                     sampling_min, steps_per_iter, pct_max_treedepth, min_ess_bulk), digits = 3, row.names = FALSE)

cg_all <- bind_rows(cgs, .id = "fit") %>% mutate(fit = factor(fit, c("trend", kt_models)))
write.csv(cg_all, file.path(kt_dir, "correlogram.csv"), row.names = FALSE)
message("\nCorrelation of site-mean residuals by distance (trend = before any site effect):")
print(cg_all %>% filter(h < 160) %>% mutate(cor = round(cor, 2), h = round(h)) %>% select(fit, h, cor) %>%
        pivot_wider(names_from = fit, values_from = cor), n = Inf)
ggsave(file.path(kt_dir, "krige_correlogram.png"),
       ggplot(cg_all, aes(h, cor, colour = fit)) +
         geom_hline(yintercept = 0, colour = "grey60") +
         geom_function(fun = function(h) coef(cfit)[["rho"]] * exp(-h / coef(cfit)[["phi"]]),
                       colour = "black", linetype = "dashed") +
         geom_line() + geom_point(aes(size = pairs)) + scale_size_area(max_size = 2.5, guide = "none") +
         labs(x = "Distance between sites (km)", y = "Correlation of site-mean residuals", colour = NULL,
              subtitle = sprintf("Dashed: fitted %.2f exp(-d / %.0f km) to the trend residuals",
                                 coef(cfit)[["rho"]], coef(cfit)[["phi"]])) +
         theme_bw(), width = 7, height = 4.5, dpi = 200)

# =============================================================================
# 5. Leave-blocks-of-sites-out check of kriging -----------------------------------
#    Sites grouped into block_km x block_km blocks. Each block's site effects
#    are kriged from the posterior effects of the sites outside it, as for a
#    location never trapped, and its training-year weeks predicted:
#      own     = the site's own fitted effect (what the 2020 hold-out uses)
#      kriged  = kriged from the other blocks (an unmonitored location)
#      none    = no site information: u = 0 with N(0, sd^2) uncertainty
#    Optimistic: the other effects (and the covariate effects) were fit with
#    these sites in; a strict version refits without each block.
# =============================================================================
blk <- with(site_xy, paste(floor(x_km / block_km), floor(y_km / block_km)))
message(sprintf("\nLeave-sites-out blocks: %d blocks of %d km (median %d sites per block)", n_distinct(blk),
                block_km, as.integer(median(table(blk)))))
chk <- list(); site_chk <- list()
for (m in kt_models[sapply(kt_models, has_site_K)]) {
  fit <- readRDS(fit_file(m))
  set.seed(42); ids <- sort(sample(ndraws(fit), min(n_check, ndraws(fit))))
  sd_draws <- site_draws(fit, ids)
  U <- sd_draws$u
  U_cv <- U; v_cv <- setNames(numeric(ncol(U)), colnames(U))
  for (b in unique(blk)) {
    i <- blk == b
    k <- krige_site(U[, !i, drop = FALSE], sd_draws$sd, XY[!i, , drop = FALSE], XY[i, , drop = FALSE], kt_phi[[m]])
    U_cv[, i] <- k$draws; v_cv[i] <- k$v
  }
  U_none <- matrix(rnorm(length(U)), nrow(U)) * sd_draws$sd
  colnames(U_none) <- colnames(U)
  sig <- as_draws_df(fit, variable = "sigma")$sigma[ids]
  chk[[m]] <- bind_rows(lapply(list(own = U, kriged = U_cv, none = U_none), function(Ux)
    score(d_fit$log_CEW, predict_new_sites(fit, m, d_fit, Ux, ids), sig)), .id = "site_effect") %>%
    mutate(model = m, .before = 1)
  site_chk[[m]] <- site_xy %>% mutate(model = m, block = blk, u_own = colMeans(U), u_kriged = colMeans(U_cv),
                                      v_kriged = v_cv)
  message(sprintf("%s: corr(own, kriged site effect) = %.2f over %d sites", m,
                  cor(colMeans(U), colMeans(U_cv)), ncol(U)))
}
if (length(chk)) {
  chk <- bind_rows(chk); write.csv(chk, file.path(kt_dir, "block_check.csv"), row.names = FALSE)
  write.csv(bind_rows(site_chk), file.path(kt_dir, "block_check_sites.csv"), row.names = FALSE)
  print(chk, digits = 3, row.names = FALSE)
}

# =============================================================================
# 6. Map of the kriged spatial effect ----------------------------------------------
#    25 km grid within 200 km of a monitored site (farther, the effect is its
#    prior, N(0, sd^2)); posterior mean and sd of u at each cell, from the
#    first spatial model fitted. The prediction at a new location adds the
#    fixed effects, year and seasonal terms at that location's covariates.
# =============================================================================
m_map <- kt_models[sapply(kt_models, has_site_K)][1]
if (!is.na(m_map)) {
  fit <- readRDS(fit_file(m_map))
  set.seed(42); ids <- sort(sample(ndraws(fit), min(n_check, ndraws(fit))))
  sd_draws <- site_draws(fit, ids)
  g <- expand.grid(x_km = seq(min(XY[, 1]) - 100, max(XY[, 1]) + 100, by = 25),
                   y_km = seq(min(XY[, 2]) - 100, max(XY[, 2]) + 100, by = 25))
  g <- g[apply(dist_km(as.matrix(g), XY), 1, min) <= 200, ]
  k <- krige_site(sd_draws$u, sd_draws$sd, XY, as.matrix(g), kt_phi[[m_map]])
  g$u_mean <- k$mean; g$u_sd <- k$sd; g$v <- k$v
  write.csv(g, file.path(kt_dir, "krige_grid.csv"), row.names = FALSE)
  gl <- bind_rows(g %>% transmute(x_km, y_km, value = u_mean, what = "Posterior mean (log scale)"),
                  g %>% transmute(x_km, y_km, value = u_sd, what = "Posterior sd"))
  pm <- lapply(split(gl, gl$what), function(d)
    ggplot(d, aes(x_km, y_km, fill = value)) + geom_raster() +
      geom_point(data = site_xy, aes(x_km, y_km), inherit.aes = FALSE, size = 0.3) +
      (if (grepl("mean", d$what[1])) scale_fill_gradient2(name = NULL) else scale_fill_viridis_c(name = NULL)) +
      coord_equal() + labs(title = d$what[1], x = "Albers x (km)", y = "Albers y (km)") + theme_bw())
  ggsave(file.path(kt_dir, "krige_map.png"), gridExtra::arrangeGrob(grobs = pm, ncol = 1),
         width = 8, height = 10, dpi = 200)
  message(sprintf("Kriged map (%s): %d cells; saved krige_map.png", m_map, nrow(g)))
}
saveRDS(list(rho = rho, phi = phi, kt_phi = kt_phi, correlogram_fit = coef(cfit), site_xy = site_xy,
             block_km = block_km, settings = kt), file.path(kt_dir, "krige_settings.rds"))
message("\nSaved to ", kt_dir)
