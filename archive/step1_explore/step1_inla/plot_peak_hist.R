# Histograms of the full posterior peak-week draws by zone, from the saved
# step1_inla.R results (<model>/peak_draws.rds: samples x site-years).
# Every draw for every site-year in a zone goes into that zone's histogram,
# so it combines the posterior uncertainty with the spread between sites and
# years. The y axis is density (each zone's histogram integrates to 1), so
# zones with more site-years don't look taller. Dashed lines are zone medians.
#
#   Rscript step1_inla/plot_peak_hist.R
#   INLA_OUT=results/step1_inla INLA_MODELS="spacetime spacetime_free" Rscript step1_inla/plot_peak_hist.R
#
# Writes peak_hist.png (all years) and peak_hist_year.png (by year) to INLA_OUT.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

out_dir <- Sys.getenv("INLA_OUT", "results/step1_inla")
mods    <- strsplit(Sys.getenv("INLA_MODELS", "spacetime spacetime_free spacetime_free_nb"), "[ ,]+")[[1]]
mods    <- mods[file.exists(file.path(out_dir, mods, "peak_draws.rds"))]
if (!length(mods)) stop("no peak_draws.rds found for the requested models in ", out_dir)

zone_lab <- c(`0` = "Northern limits", `1` = "Transitional zone", `2` = "Southern range")

# one row per draw x site-year
dr <- bind_rows(lapply(mods, function(m) {
  x  <- readRDS(file.path(out_dir, m, "peak_draws.rds"))
  sy <- x$site_years
  data.frame(model = m,
             zone_30y = rep(sy$zone_30y, each = nrow(x$draws)),
             year     = rep(sy$year,     each = nrow(x$draws)),
             week     = as.vector(x$draws))
})) %>% mutate(zone = factor(zone_lab[as.character(zone_30y)], zone_lab),
               model = factor(model, mods))
cat(format(nrow(dr), big.mark = ","), "draws over", length(mods), "models\n")

med <- dr %>% group_by(model, zone) %>%
  summarise(n_draws = n(), med = median(week), mean = mean(week),
            lo = quantile(week, 0.025), hi = quantile(week, 0.975), .groups = "drop")
print(as.data.frame(med), digits = 3)

hist_layer <- function(lw) geom_histogram(aes(y = after_stat(density)), binwidth = 1, boundary = 0.5,
                                          position = "identity", alpha = 0.35, linewidth = lw)
scales <- list(scale_fill_brewer(palette = "Dark2"), scale_colour_brewer(palette = "Dark2"))

p <- ggplot(dr, aes(week, fill = zone, colour = zone)) +
  hist_layer(0.3) +
  geom_vline(data = med, aes(xintercept = med, colour = zone), linetype = 2, linewidth = 0.5) +
  facet_wrap(~ model, ncol = 1) + scales +
  scale_x_continuous(breaks = seq(0, 52, 2)) +
  labs(x = "Peak week (week of year)", y = "Density", fill = NULL, colour = NULL,
       subtitle = "Posterior peak-week draws, all site-years, by zone (density within zone; dashed: zone median)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "peak_hist.png"), p, width = 7, height = 1 + 2.2 * length(mods), dpi = 300)

p <- ggplot(dr, aes(week, fill = zone, colour = zone)) +
  hist_layer(0.2) +
  facet_grid(year ~ model) + scales +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(x = "Peak week (week of year)", y = "Density", fill = NULL, colour = NULL,
       subtitle = "Posterior peak-week draws by zone and year (density within zone-year)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
ggsave(file.path(out_dir, "peak_hist_year.png"), p, width = 2.5 + 2.5 * length(mods),
       height = 1.5 + 0.9 * length(unique(dr$year)), dpi = 300)
cat("Wrote peak_hist.png and peak_hist_year.png to", out_dir, "\n")
