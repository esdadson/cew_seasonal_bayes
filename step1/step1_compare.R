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
    message("Reusing ", check_file(m)); res[[m]] <- readRDS(check_file(m)); next
  }
  message("\n==== ", m, " ====\n", deparse1(step1_formula(m)$formula))
  t0  <- Sys.time()
  fit <- fit_step1(m)
  st  <- cbind(sampler_stats(fit, m), wall_min = mins(t0))
  print(st)
  lo  <- loo(fit, save_psis = TRUE)
  mt  <- step1_metrics(fit, m, lo)
  lo$psis_object <- NULL
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
                    "wall_min", "divergent", "max_rhat", "min_ess_bulk"), names(tab))
print(tab[show], digits = 3)

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
