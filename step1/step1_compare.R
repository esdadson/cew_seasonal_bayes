# =============================================================================
# step1/step1_compare.R
# Fits the step 1 candidate models (step1_models.R) to 2008-2019, scores them
# with PSIS-LOO and on the held-out year 2020, and writes a comparison table
# and two diagnostic figures. Runs on the laptop (no site terms yet).
#
#   Rscript step1/step1_compare.R
#   STEP1_MODELS="fixed seasonal"   models to run (default: all in step1_specs)
#   STEP1_CHAINS / STEP1_ITER / STEP1_WARMUP / STEP1_THIN / STEP1_THREADS
#
# One model is in memory at a time; each model's results are cached in
# results/step1/check_<model>.rds and reused while the fit is unchanged.
#
# Outputs (results/step1/):
#   comparison.csv          loo_compare + metrics + sampler stats, one row per model
#   test_rows.csv           2020 predictions per trap-week and model
#   coverage_by_season.csv  2020 90% interval coverage and MAE by season block
#   holdout_zone_week.png   2020 zone weekly means: observed vs each model
#   coverage_by_season.png  2020 coverage by season block (target 0.90)
#   spatial_range.csv       correlogram of the trend's site-mean residuals + the
#                           rho / phi used by the spatial model
#   time_estimates.csv      projected time of a final fit (all years) per model
#                           and sampler setting, from this run's timings
# =============================================================================

source("step1/step1_models.R")
library(ggplot2)

models <- strsplit(Sys.getenv("STEP1_MODELS", paste(names(step1_specs), collapse = " ")), "\\s+")[[1]]
stopifnot(all(models %in% names(step1_specs)))
check_file <- function(m) file.path(step1_dir, paste0("check_", m, ".rds"))
mins <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "mins"))

message("Training rows: ", nrow(d_fit), "; held-out (", holdout, ") rows: ", nrow(d_test))

res <- list()
for (m in models) {
  ff <- paste0(fit_file(m), ".rds")
  if (file.exists(check_file(m)) && file.exists(ff) && file.mtime(check_file(m)) > file.mtime(ff)) {
    ck <- readRDS(check_file(m))
    if ("post_min" %in% names(ck$stats)) {         # older check files lack the timings
      message("Reusing ", check_file(m)); res[[m]] <- ck; next
    }
  }
  message("\n==== ", m, " ====\n", deparse1(step1_formula(m)$formula))
  t0  <- Sys.time()
  fit <- fit_step1(m)
  st  <- cbind(sampler_stats(fit, m), fit_wall_min = mins(t0))
  print(st)
  t1  <- Sys.time()
  lo  <- loo(fit, save_psis = TRUE)
  mt  <- step1_metrics(fit, m, lo)
  lo$psis_object <- NULL
  st$post_min <- mins(t1)                      # loo + metrics
  st$draws    <- ndraws(fit)
  st$iter     <- step1_sampler$iter; st$warmup <- step1_sampler$warmup; st$chains <- step1_sampler$chains
  res[[m]] <- list(stats = st, loo = lo, metrics = mt$metrics, test_rows = mt$test_rows,
                   fixef = fixef(fit), vc = VarCorr(fit))
  saveRDS(res[[m]], check_file(m))
  print(summary(fit))
  rm(fit, lo, mt); invisible(gc())
}

# ---- Comparison table -------------------------------------------------------------
tab <- bind_rows(lapply(res, function(r) left_join(r$metrics, r$stats, by = "model")))
if (length(res) > 1) {
  lc  <- loo_compare(lapply(res, `[[`, "loo"))
  tab <- left_join(tab, data.frame(model = rownames(lc), elpd_diff = lc[, "elpd_diff"],
                                   se_diff = lc[, "se_diff"]), by = "model")
}
write.csv(tab, file.path(step1_dir, "comparison.csv"), row.names = FALSE)
show <- intersect(c("model", "elpd_diff", "se_diff", "elpd_loo", "p_loo", "loo_R2", "mae_loo",
                    "elpd_test", "mae_test", "R2_test", "cov90_test",
                    "fit_wall_min", "divergent", "max_rhat", "min_ess_bulk"), names(tab))
print(tab[show], digits = 3)

# ---- Spatial range used by the spatial model -----------------------------------------
write.csv(cg_trend %>% mutate(bin = as.character(bin),
                              fitted = sp_rho * exp(-h / sp_phi),
                              rho_est = sp_est[["rho"]], phi_est = sp_est[["phi"]],
                              rho_used = sp_rho, phi_used = sp_phi),
          file.path(step1_dir, "spatial_range.csv"), row.names = FALSE)

# ---- Projected time of the final fit ----------------------------------------------------
# Per model, from this run: minutes per warmup and per sampling iteration
# (mean over chains, which run in parallel), scaled to each target setting and
# by rows (the final fit uses every year: gradient cost ~ rows). loo + metrics
# scale with draws x rows. Measured on THIS machine: a Hazel core can be
# slower or faster - compare one model's fit_wall_min between the two.
# Not included: queue wait, compilation (~1 min), 52-week predictions.
targets <- data.frame(iter = c(2000, 4000, 10000), warmup = c(1000, 2000, 5000), draws = 1000)
row_scale <- nrow(d_obs) / nrow(d_fit)
est <- bind_rows(lapply(res, function(r) {
  s <- r$stats
  w_per <- s$warmup_min / s$warmup
  s_per <- s$sampling_min / (s$iter - s$warmup)
  targets %>% mutate(model = s$model,
                     fit_hours  = (w_per * warmup + s_per * (iter - warmup)) * row_scale / 60,
                     post_hours = s$post_min * (draws / s$draws) * row_scale / 60,
                     total_hours = fit_hours + post_hours,
                     cpus = s$chains * step1_sampler$threads,
                     measured_iter = s$iter, measured_steps_per_iter = s$steps_per_iter,
                     measured_ms_per_grad = s$ms_per_grad)
})) %>% relocate(model)
write.csv(est, file.path(step1_dir, "time_estimates.csv"), row.names = FALSE)
message("\nProjected final-fit time (all years, this machine; hours):")
print(est %>% select(model, iter, warmup, fit_hours, post_hours, total_hours), digits = 2)

# ---- Held-out year diagnostics ------------------------------------------------------
tr <- bind_rows(lapply(res, `[[`, "test_rows"))
if (nrow(tr)) {
  write.csv(tr, file.path(step1_dir, "test_rows.csv"), row.names = FALSE)
  zone_lab <- c(`0` = "Northern limits", `1` = "Transitional zone", `2` = "Southern range")

  cov <- tr %>% group_by(model, block) %>%
    summarise(n = n(), cov90 = mean(log_CEW >= pi_lo & log_CEW <= pi_hi),
              mae = mean(abs(log_CEW - pred_mean)), sd_resid = sd(log_CEW - pred_mean),
              .groups = "drop")
  write.csv(cov, file.path(step1_dir, "coverage_by_season.csv"), row.names = FALSE)
  print(as.data.frame(cov), digits = 3)

  zw <- tr %>% group_by(model, zone_30y, woy) %>%
    summarise(obs = mean(log_CEW), pred = mean(pred_mean), n = n(), .groups = "drop") %>%
    mutate(zone = factor(zone_lab[as.character(zone_30y)], zone_lab))
  p1 <- ggplot(zw, aes(woy)) +
    geom_point(aes(y = obs, size = n), colour = "grey40", alpha = 0.6) +
    geom_line(aes(y = pred, colour = model), linewidth = 0.6) +
    facet_wrap(~ zone, ncol = 1) +
    scale_size_area(max_size = 2, name = "Traps") +
    labs(x = "Week of year", y = "Mean log(CEW + 1)", colour = "Model",
         subtitle = paste0(holdout, " (held out): observed zone weekly means vs predictions")) +
    theme_bw(base_size = 9)
  ggsave(file.path(step1_dir, "holdout_zone_week.png"), p1, width = 6, height = 6, dpi = 200)

  p2 <- ggplot(cov, aes(block, cov90, colour = model, group = model)) +
    geom_hline(yintercept = 0.9, linetype = 2, colour = "grey50") +
    geom_point() + geom_line() +
    labs(x = "Week of year", y = "90% interval coverage", colour = "Model",
         subtitle = paste0(holdout, " (held out): coverage by season block")) +
    theme_bw(base_size = 9)
  ggsave(file.path(step1_dir, "coverage_by_season.png"), p2, width = 5, height = 3, dpi = 200)
}
message("Done. Outputs in ", step1_dir)
