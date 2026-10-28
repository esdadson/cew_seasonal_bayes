# Publication figures for JIPM
# 2. Get map boundaries for countries and Great Lakes

# Figure 1
us_states <- map_data("state")
canada_prov <- map_data("world", region = c("Canada", "Mexico"))
lakes <- map_data("lakes")

# 3. Create the usmap-styled plot
p = ggplot() +
  # Draw Canada base layer (light aesthetic)
  geom_polygon(data = canada_prov, aes(x = long, y = lat, group = group), 
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3) +
  # Draw US States layer 
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group), 
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3) +
  # Overlay the Great Lakes (using a light blue fill to punch them out)
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), 
               fill = "#cbdbee", color = "#9c9c9c", linewidth = 0.3) +
  # Plot points directly 
  geom_point(
    data = loc_year_min5[loc_year_min5$nyr < 3,], 
    aes(x = longitude, y = latitude), shape = 4, size = 1
  ) + geom_point(data = loc_year_min5[loc_year_min5$nyr >= 3,], 
                 aes(x = longitude, y = latitude, color = nyr), size = 2.5, 
                 alpha = 0.8
  ) + scale_color_viridis_c(name = "Years \nSampled", limits = c(1,31)) + 
  # Set the usmap-style projection (Lambert Conformal Conic) 
  # and frame it cleanly to exclude Alaska/Hawaii/Northern Canada
  coord_sf(crs = st_crs(4269), xlim = c(-125, -66), ylim = c(24, 52)) + 
  # Apply clean minimalist theme settings
  theme_void() +
  theme(
    panel.background = element_rect(fill = NA, color = NA), # Blue ocean background
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=10)),
    plot.margin = margin(20, 20, 20, 20)
  )

# Figure 2
ggplot() +
  # Use tidyterra's specialized geometry for categorical rasters
  geom_spatvector(data = us_states, fill = "#f9f9f9") +
  geom_spatraster(data = ag_land_rast, aes(fill = class)) +
  
  # Apply the custom color palette
  scale_fill_manual(
    values = nass_colors_alt,
    na.value = "transparent",
    name = "Land Cover Type"
  ) +
  geom_spatvector(data = us_states, fill = NA, color = 'black', linewidth = 0.7) +
  # Format coordinate system labels professionally
  # scale_x_continuous(labels = function(x) paste0(abs(x), "°W")) +
  # scale_y_continuous(labels = function(y) paste0(y, "°N")) +
  
  # Apply strict, clean publication themes
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    # Grid lines and panel modifications
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # Legend styling for layout compactness
    legend.position        = "right",
    legend.title           = element_text(face = "bold", size = 10),
    legend.text            = element_text(size = 9),
    legend.background      = element_blank(),
    legend.box.background  = element_blank(),
    
    # Text formatting
    axis.title  = element_blank(), # Lat/Long are self-explanatory
    axis.text   = element_text(color = "black", size = 9),
    plot.title  = element_text(face = "bold", size = 13, hjust = 0),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

library(cowplot)
library(ggpubr)

legend_plot = ggpubr::get_legend(p)
legend_plot = ggpubr::as_ggplot(legend_plot)

p_no_legend = p + theme(legend.position = "none")

ggsave(filename = "C:/Users/edadson/Desktop/CEW/ag_rast_legend.png", legend_plot, width = 3,
       height = 3, units = 'in', dpi = 300)
ggsave(filename = "C:/Users/edadson/Desktop/CEW/ag_rast.png", p_no_legend, width = 14,
       height = 7, units = 'in', dpi = 300)

# peak week plot
# peak_sf from step 1 (results/step1/peak_site_year.csv: posterior peak week per site-year), when present
if (file.exists("C:/Users/edadson/Desktop/cew_seasonal_bayes/results/step1/peak_site_year.csv")) {
  peak_sf = read.csv("C:/Users/edadson/Desktop/cew_seasonal_bayes/results/step1/peak_site_year.csv") %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
}
p = ggplot() +
  geom_polygon(data = canada_prov, aes(x = long, y = lat, group = group), 
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3) +
  # Draw US States layer 
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group), 
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3) +
  # Overlay the Great Lakes (using a light blue fill to punch them out)
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), 
               fill = "#cbdbee", color = "#9c9c9c", linewidth = 0.3) +
  # trap points colored by peak week
  geom_sf(
    data  = peak_sf,
    aes(color = peak_week_sd, shape = factor(zone_30y)),
    size  = 2, alpha = 0.85
  ) +
  # one panel per year
  facet_wrap(~ year, ncol = 3) +        # ncol = 7 gives 6 rows for 41 years
  # diverging palette centered on median peak week
  scale_shape_manual(labels = c("0" = "Northern limits", "1" = "Transitional zone",
                                "2" = "Southern range"),
                     values = c("0" = 16, "1" = 17,
                                "2" = 15),
                     name = "Overwintering \nzone") + 
  scale_color_viridis_c(
    name   = "Peak Week \nCI Width", direction = -1,
    option = "plasma",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  coord_sf(crs = st_crs(4269), xlim = c(-125, -66), ylim = c(24, 52)) + 
  theme_void() +
  theme(
    panel.background = element_rect(fill = NA, color = NA), # Blue ocean background
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=10)),
    plot.margin = margin(20, 20, 20, 20)
  )

legend_plot = ggpubr::get_legend(p)
legend_plot = ggpubr::as_ggplot(legend_plot)

p_no_legend = p + theme(legend.position = "none")

ggsave(filename = "C:/Users/edadson/Desktop/CEW/peak_cew_mean_legend.png", legend_plot, width = 3,
       height = 3, units = 'in', dpi = 300)
ggsave(filename = "C:/Users/edadson/Desktop/CEW/peak_cew_mean.png", p_no_legend, width = 7,
       height = 14, units = 'in', dpi = 300)


# Plot model averaged curve for overwintering regions by year
p = ggplot(data = dat6_brms, aes(x = woy, y = exp(fitted)-1, 
                        color = factor(zone_30y), 
                        fill = factor(zone_30y))) + 
  geom_smooth(method = "gam", formula = y ~ s(x)) + 
  scale_color_manual(name  = "Overwintering \nzone", 
                     values = c("0" = "gray",
                                "1" = "yellow3",
                                "2" = "deepskyblue4"), 
                     labels = c("0" = "Northern limits",
                                "1" = "Transitional zone", 
                                "2" = "Southern range")) + 
  scale_fill_manual(name  = "Overwintering \nzone", 
                     values = c("0" = "gray",
                                "1" = "yellow3",
                                "2" = "deepskyblue4"), 
                     labels = c("0" = "Northern limits",
                                "1" = "Transitional zone", 
                                "2" = "Southern range")) + 
  scale_alpha_manual(values = c(0.5, 0.5, 0.5)) + 
  ylab("") + xlab("Week of Year") + 
  facet_wrap(~year, ncol = 3) + 
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA), # Blue ocean background
    plot.title = element_text(face = "bold", size = 14, margin = margin(b=10)),
    plot.margin = margin(20, 20, 20, 20)
  )
legend_plot = ggpubr::get_legend(p)
legend_plot = ggpubr::as_ggplot(legend_plot)

p_no_legend = p + theme(legend.position = "none")

ggsave(filename = "C:/Users/edadson/Desktop/CEW/peak_cew_smooth_legend.png", legend_plot, width = 3,
       height = 3, units = 'in', dpi = 300)
ggsave(filename = "C:/Users/edadson/Desktop/CEW/peak_cew_smooth.png", p_no_legend, width = 7,
       height = 7, units = 'in', dpi = 300)


# site-years plot
site_years = dat6 %>% 
  group_by(year) %>% 
  summarize(sites = length(unique(location2))) %>%
  arrange(year)

sites_weeks = dat6 %>% 
  group_by(woy) %>% 
  summarize(observations = sum(!is.na(CEW_sum))) %>% 
  arrange(woy)

cew_weeks = dat6 %>% 
  group_by(woy) %>%
  summarize(average_cew = mean(CEW_sum, na.rm = T)) %>% 
  arrange(woy)

cew_years = dat6 %>% 
  group_by(year) %>%
  summarize(sum_cew = sum(CEW_sum, na.rm = T)) %>% 
  arrange(year)

a = ggplot(data = site_years, aes(x = year , y = sites)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_line(linetype = 2) + 
  scale_x_continuous(breaks = seq(2008, 2020, by = 1), minor_breaks = NULL) + 
  xlab("") + ylab("Sites") + 
  theme_minimal()

b = ggplot(data = sites_weeks, aes(x = woy , y = observations)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_line(linetype = 2) + 
  scale_x_continuous(breaks = seq(1, 52, by = 4), minor_breaks = NULL) + 
  xlab("") + ylab("Observations") + 
  theme_minimal()

d = ggplot(data = dat6, aes(x = woy , y = CEW_sum)) + 
  geom_point(size = 2, alpha = 0.8) + 
  geom_smooth(color = "lightgray", fill = "gray") + 
  scale_x_continuous(breaks = seq(1, 52, by = 4), minor_breaks = NULL) + 
  xlab("Week") + ylab("Moth count") + 
  theme_minimal()

c = ggplot(data = cew_years, aes(x = year, y = sum_cew)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_line(linetype = 2) + 
  scale_x_continuous(breaks = seq(2008, 2020, by = 1), minor_breaks = NULL) + 
  xlab("Year") + ylab("Cumulative count") + 
  theme_minimal()

final_plot = cowplot::plot_grid(a, b, c, d, labels = c("a", "b", "c", "d"), ncol = 2, nrow = 2)

ggsave(filename = "C:/Users/edadson/Desktop/CEW/sites_over_time.png", final_plot, width = 10,
       height = 10, units = 'in', dpi = 300)

# landscape coefficient estimates
beta_summary = summary(samples[,grep("^beta", colnames(posterior))])

plot_data <- data.frame(
  Parameter = rownames(beta_summary),
  Mean      = beta_summary$Estimate,
  Lower     = beta_summary$`l-95% CI`,
  Upper     = beta_summary$`u-95% CI`
)

plot_data$Significant <- ifelse(plot_data$Lower > 0, 
                                "Significant (+)", ifelse(plot_data$Upper < 0,
                                                          "Significant (-)",
                                                          "Not Significant"))

plot_data1 = plot_data %>% filter(Parameter != "Intercept",
                                 Parameter != "EVI_filled",
                                 Parameter != "gdd")

p1 = ggplot(plot_data1, aes(x = Parameter, y = Mean, fill = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_crossbar(aes(ymin = Lower, ymax = Upper), width = 0.5, size = 0.6) + 
  scale_x_discrete(labels = c("Corn", "Cotton", "Other Crops", "Peanuts", "Soybeans", "Tobacco")) + 
  scale_fill_manual(values = c("Significant (+)" = "forestgreen",
                               "Significant (-)" = "darkred",
                               "Not Significant" = "lightgray")) + 
  ylab("") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

legend_plot = ggpubr::get_legend(p1)
legend_plot = ggpubr::as_ggplot(legend_plot)

p1 = p1 + theme(legend.position = "none")

plot_data2 = plot_data %>% filter(Parameter != "Intercept",
                                  (Parameter == "EVI_filled") | (Parameter == "gdd"))
p2 = ggplot(plot_data2, aes(x = Parameter, y = Mean, fill = Significant)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_crossbar(aes(ymin = Lower, ymax = Upper), width = 0.5, size = 0.6) +
  scale_x_discrete(labels = c("EVI", "GDD")) + 
  scale_fill_manual(values = c("Significant (+)" = "forestgreen",
                               "Significant (-)" = "darkred", 
                               "Not Significant" = "lightgray")) + 
  ylab("") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

p2 = p2 + theme(legend.position = "none")

final_plot = cowplot::plot_grid(p1, p2, labels = c("a", "b"), ncol = 2, nrow = 1)
ggsave(filename = "C:/Users/edadson/Desktop/CEW/posterior_credible_intervals.png", final_plot,
       width = 10,
       height = 10, units = 'in', dpi = 300)
ggsave(filename = "C:/Users/edadson/Desktop/CEW/posterior_credible_intervals_legend.png", legend_plot,
       width = 4,
       height = 6, units = 'in', dpi = 300)


# =============================================================================
# Step 1 and step 2 figures
#   Step 1 (results/step1/, from step1/step1_hpc.R compare): model comparison,
#     predicted weekly curves (weeks 1-52) and peak week.
#   Step 2 (results/step2/, from step2/step2_granger.R): regions and
#     Southern -> Transitional/Northern Granger causality and lag times.
# The base map layers (us_states, canada_prov, lakes) are those of Figure 1 above.
# =============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)

proj_dir = "C:/Users/edadson/Desktop/cew_seasonal_bayes"
step1_dir = file.path(proj_dir, "results/step1")
step2_dir = file.path(proj_dir, "results/step2")
fig_dir = "C:/Users/edadson/Desktop/CEW"

zone_labels = c("0" = "Northern limits", "1" = "Transitional zone", "2" = "Southern range")
zone_colors = c("0" = "gray40", "1" = "yellow3", "2" = "deepskyblue4")
zone_shapes = c("0" = 16, "1" = 17, "2" = 15)

# Base layers and frame for the maps (zoomed to the trapped area)
base_layers = list(
  geom_polygon(data = canada_prov, aes(x = long, y = lat, group = group),
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3),
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group),
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.3),
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group),
               fill = "#cbdbee", color = "#9c9c9c", linewidth = 0.3)
)
map_frame = list(
  coord_sf(crs = st_crs(4269), xlim = c(-106, -66), ylim = c(24, 50)),
  theme_void(),
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 10)),
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(20, 20, 20, 20)
  )
)

# Step-2 regions, when step 2 has run (outlined on the step-1 peak map too)
g2 = NULL
hulls = NULL
if (file.exists(file.path(step2_dir, "granger_results.rds"))) {
  g2 = readRDS(file.path(step2_dir, "granger_results.rds"))
  regions = g2$regions
  reg_info = g2$reg_info
  pairs = g2$granger_pairs %>%
    mutate(sink_zone = factor(sink_zone, levels = c("Transitional", "Northern")))
  # Convex hull of each region's sites (tested = enough weekly data for Granger tests)
  hulls = regions %>%
    group_by(region, zone_30y) %>%
    slice(chull(longitude, latitude)) %>%
    ungroup() %>%
    left_join(reg_info %>% select(region, tested), by = "region")
}

# -----------------------------------------------------------------------------
# Step 1
# -----------------------------------------------------------------------------
if (file.exists(file.path(step1_dir, "loo_compare.csv"))) {
  # Model comparison: ELPD difference from the best model (PSIS-LOO, +/- 2 SE)
  lc = read.csv(file.path(step1_dir, "loo_compare.csv"))
  if ("elpd_diff" %in% names(lc)) {
    lc = lc %>%
      mutate(model = factor(model, levels = c("fixed", "spatial", "fourier"),
                            labels = c("Fixed effects", "+ Spatial GP", "+ Fourier seasonality")))
    p_loo = ggplot(lc, aes(x = elpd_diff, y = model)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_pointrange(aes(xmin = elpd_diff - 2 * se_diff, xmax = elpd_diff + 2 * se_diff),
                      color = "deepskyblue4") +
      labs(x = "ELPD difference from the best model (PSIS-LOO, \u00b1 2 SE)", y = NULL) +
      theme_minimal()

    ggsave(filename = file.path(fig_dir, "step1_loo_compare.png"), p_loo, width = 7,
           height = 3, units = 'in', dpi = 300)
  }
}

if (file.exists(file.path(step1_dir, "pred_zone_week.csv"))) {
  # Predicted weekly curves, weeks 1-52: zone mean E[log(CEW + 1)] over the
  # zone's monitored site-years, by year (95% CrI)
  pzw = read.csv(file.path(step1_dir, "pred_zone_week.csv"))
  p_pred = ggplot(pzw, aes(x = woy, y = med, ymin = lo, ymax = hi,
                           color = factor(zone_30y), fill = factor(zone_30y))) +
    geom_ribbon(alpha = 0.25, color = NA) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    scale_fill_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    facet_wrap(~ year, ncol = 3) +
    labs(x = "Week of Year", y = "Predicted log(CEW + 1) (95% CrI)") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "step1_pred_zone_week.png"), p_pred, width = 9,
         height = 11, units = 'in', dpi = 300)
}

if (file.exists(file.path(step1_dir, "peak_site.csv"))) {
  peak_site = read.csv(file.path(step1_dir, "peak_site.csv"))
  peak_zone_year = read.csv(file.path(step1_dir, "peak_zone_year.csv"))

  # Site peak week, averaged over the site's monitored years (step-2 region
  # outlines when available)
  p_peak_site = ggplot() + base_layers +
    (if (!is.null(hulls)) geom_polygon(data = hulls, aes(x = longitude, y = latitude, group = region),
                                       fill = NA, color = "gray30", linewidth = 0.3, linetype = 2)) +
    geom_point(data = peak_site, aes(x = longitude, y = latitude, color = med,
                                     shape = factor(zone_30y)),
               size = 2.2, alpha = 0.9) +
    scale_color_viridis_c(name = "Peak week \n(posterior \nmedian)", option = "plasma") +
    scale_shape_manual(values = zone_shapes, labels = zone_labels, name = "Overwintering \nzone") +
    map_frame

  ggsave(filename = file.path(fig_dir, "peak_week_site.png"), p_peak_site, width = 9,
         height = 6, units = 'in', dpi = 300)

  # Zone mean peak week by year (95% credible intervals)
  p_peak_zone = ggplot(peak_zone_year, aes(x = year, y = med, ymin = lo, ymax = hi,
                                           color = factor(zone_30y))) +
    geom_pointrange(position = position_dodge(0.5), size = 0.3) +
    scale_color_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    scale_x_continuous(breaks = seq(2008, 2020, by = 1), minor_breaks = NULL) +
    labs(x = "Year", y = "Mean peak week (95% CrI)") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "peak_week_zone_year.png"), p_peak_zone, width = 9,
         height = 5, units = 'in', dpi = 300)
}

# -----------------------------------------------------------------------------
# Step 2
# -----------------------------------------------------------------------------
if (!is.null(g2)) {
  # Regions map: clusters within each overwintering zone; dashed = too sparse to test
  p_regions = ggplot() + base_layers +
    geom_polygon(data = hulls, aes(x = longitude, y = latitude, group = region,
                                   fill = factor(zone_30y), linetype = tested),
                 alpha = 0.2, color = "black", linewidth = 0.3) +
    scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                          labels = c("TRUE" = "Yes", "FALSE" = "No (sparse)"),
                          name = "Granger \ntested") +
    geom_point(data = regions, aes(x = longitude, y = latitude,
                                   color = factor(zone_30y), shape = factor(zone_30y)),
               size = 1.6, alpha = 0.85) +
    geom_label(data = reg_info, aes(x = longitude, y = latitude, label = region),
               size = 2.8, label.padding = unit(0.12, "lines"), fill = "white", alpha = 0.8) +
    scale_fill_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    scale_color_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    scale_shape_manual(values = zone_shapes, labels = zone_labels, name = "Overwintering \nzone") +
    map_frame

  ggsave(filename = file.path(fig_dir, "granger_regions.png"), p_regions, width = 9,
         height = 6, units = 'in', dpi = 300)

  # Region seasonal curves: mean weekly log(CEW + 1) (observed weeks + step-1
  # predictions for unsampled weeks, averaged over years) - what the clustering grouped
  p_curves = ggplot(g2$region_curves, aes(x = woy, y = value, color = region)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ factor(zone_30y, levels = c("2", "1", "0"), labels = zone_labels[c("2", "1", "0")]),
               ncol = 3) +
    labs(x = "Week of Year", y = "Mean log(CEW + 1), observed + predicted", color = "Region") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_region_curves.png"), p_curves, width = 11,
         height = 4, units = 'in', dpi = 300)

  # Granger link map: arrows from each Southern region to the sink regions it
  # improves out of sample (gain > 0); bold = permutation test p < 0.05 (BH)
  links = pairs %>%
    filter(best_gain > 0) %>%
    mutate(test = ifelse(link, "p < 0.05 (BH)", "n.s."),
           # stop arrows short of the sink centroid so the label stays readable
           xend = longitude_src + 0.85 * (longitude_snk - longitude_src),
           yend = latitude_src + 0.85 * (latitude_snk - latitude_src)) %>%
    arrange(best_gain)

  p_links = ggplot() + base_layers +
    geom_polygon(data = hulls, aes(x = longitude, y = latitude, group = region,
                                   fill = factor(zone_30y)),
                 alpha = 0.15, color = "gray50", linewidth = 0.2) +
    geom_segment(data = links,
                 aes(x = longitude_src, y = latitude_src, xend = xend, yend = yend,
                     color = 100 * best_gain, linewidth = test, alpha = test),
                 arrow = arrow(length = unit(0.1, "inches"), type = "closed")) +
    geom_label(data = reg_info, aes(x = longitude, y = latitude, label = region),
               size = 2.4, label.padding = unit(0.1, "lines"), fill = "white", alpha = 0.8) +
    scale_color_viridis_c(name = "Out-of-sample \ngain (%)", option = "viridis", end = 0.85) +
    scale_linewidth_manual(values = c("p < 0.05 (BH)" = 1.1, "n.s." = 0.4), name = "Permutation \ntest") +
    scale_alpha_manual(values = c("p < 0.05 (BH)" = 0.95, "n.s." = 0.6), name = "Permutation \ntest") +
    scale_fill_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    facet_wrap(~ sink_zone) +
    map_frame

  ggsave(filename = file.path(fig_dir, "granger_links.png"), p_links, width = 14,
         height = 6, units = 'in', dpi = 300)

  # Lag map: each sink region filled by the lag (weeks) at which its strongest
  # Southern source best predicts it; faded = not significant; unfilled = no
  # Southern lag improves out-of-sample prediction
  best_src = pairs %>%
    filter(best_gain > 0) %>%
    group_by(sink) %>%
    slice_max(best_gain, n = 1, with_ties = FALSE) %>%
    ungroup()
  lag_hulls = hulls %>%
    inner_join(best_src %>% select(region = sink, best_lag, link), by = "region")

  p_lag = ggplot() + base_layers +
    geom_polygon(data = hulls %>% filter(zone_30y == 2),
                 aes(x = longitude, y = latitude, group = region),
                 fill = "gray85", color = "gray40", linewidth = 0.3) +
    geom_polygon(data = lag_hulls, aes(x = longitude, y = latitude, group = region,
                                       fill = best_lag, alpha = link),
                 color = "black", linewidth = 0.3) +
    geom_label(data = reg_info %>% inner_join(best_src, by = c(region = "sink")),
               aes(x = longitude, y = latitude,
                   label = paste0(region, " (", source, ", ", best_lag, " wk)")),
               size = 2.3, label.padding = unit(0.1, "lines"), fill = "white", alpha = 0.8) +
    geom_label(data = reg_info %>% filter(zone_30y == 2),
               aes(x = longitude, y = latitude, label = region),
               size = 2.3, label.padding = unit(0.1, "lines"), fill = "gray95") +
    scale_fill_viridis_c(name = "Most predictive \nSouthern lag \n(weeks)", limits = c(1, 8),
                         breaks = 1:8) +
    scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.3),
                       labels = c("TRUE" = "p < 0.05 (BH)", "FALSE" = "n.s."),
                       name = "Permutation \ntest") +
    map_frame

  ggsave(filename = file.path(fig_dir, "granger_best_lag.png"), p_lag, width = 9,
         height = 6, units = 'in', dpi = 300)

  # Lag profiles: out-of-sample gain by Southern lag, one panel per sink region
  lag_prof = g2$granger_reg %>%
    mutate(lag = as.integer(sub("lag", "", setting)),
           sig = !is.na(p_perm) & p_perm < 0.05)

  p_prof = ggplot(lag_prof, aes(x = lag, y = 100 * loyo_gain, color = source)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line() +
    geom_point(aes(shape = sig), size = 1.8) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Permutation \np < 0.05") +
    scale_x_continuous(breaks = 1:8) +
    facet_wrap(~ sink) +
    labs(x = "Southern lag (weeks)", y = "Out-of-sample gain (%)", color = "Southern \nregion") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_lag_profiles.png"), p_prof, width = 10,
         height = 7, units = 'in', dpi = 300)

  # Zone level: lags (forward and reverse) and leads (falsification)
  gz = g2$granger_zone %>%
    filter(grepl("^(lag|lead)", setting)) %>%
    mutate(k = as.integer(sub("lag|lead", "", setting)) * ifelse(grepl("^lag", setting), 1, -1),
           pair = factor(paste(source, "->", sink),
                         levels = c("S -> T", "S -> N", "T -> S", "N -> S")),
           sig = !is.na(p_perm) & p_perm < 0.05)

  p_zone = ggplot(gz, aes(x = k, y = 100 * loyo_gain, color = direction)) +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line() +
    geom_point(aes(shape = sig), size = 2) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Permutation \np < 0.05") +
    scale_color_manual(values = c(forward = "deepskyblue4", reverse = "gray50"), name = "") +
    scale_x_continuous(breaks = -3:8) +
    facet_wrap(~ pair, ncol = 2) +
    labs(x = "Lag of source zone (weeks; negative = lead)", y = "Out-of-sample gain (%)") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_zone_lags.png"), p_zone, width = 8,
         height = 6, units = 'in', dpi = 300)

  # Leave-one-year-out predictions at each sink zone's best Southern lag:
  # observed vs own-history-only vs + Southern lag (one week ahead, held-out year)
  best_set = g2$granger_zone %>%
    filter(direction == "forward", grepl("^lag", setting)) %>%
    group_by(sink) %>%
    slice_max(loyo_gain, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(source, sink, setting)
  oos = g2$granger_pred_zone %>%
    semi_join(best_set, by = c("source", "sink", "setting")) %>%
    transmute(sink, setting, year, woy,
              `Observed` = clim + anom,
              `Own history` = clim + pred_R,
              `+ Southern lag` = clim + pred_U) %>%
    pivot_longer(c(`Observed`, `Own history`, `+ Southern lag`), names_to = "series") %>%
    mutate(series = factor(series, levels = c("Observed", "Own history", "+ Southern lag")))

  for (s in unique(oos$sink)) {
    d = oos %>% filter(sink == s)
    p_oos = ggplot(d, aes(x = woy, y = value, color = series)) +
      geom_line() +
      scale_color_manual(values = c("Observed" = "black", "Own history" = "gray60",
                                    "+ Southern lag" = "deepskyblue4"), name = "") +
      facet_wrap(~ year, ncol = 3) +
      labs(x = "Week of Year", y = "Regional mean log(CEW + 1)",
           title = paste0(c(T = "Transitional", N = "Northern")[s], " zone, held-out years (",
                          sub("lag", "Southern lag ", d$setting[1]), " wk)")) +
      theme_minimal()
    ggsave(filename = file.path(fig_dir, paste0("granger_loyo_", s, ".png")), p_oos, width = 8,
           height = 9, units = 'in', dpi = 300)
  }

  # Region pairs: step-1 peak delay vs step-2 most predictive lag, by distance
  p_delay = ggplot(pairs, aes(x = dist_km)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(y = delay_med, ymin = delay_lo, ymax = delay_hi,
                        color = "Peak delay (step 1)"), size = 0.3, alpha = 0.7) +
    geom_point(aes(y = best_lag, color = "Most predictive lag (step 2)", shape = link),
               size = 2.2) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "p < 0.05 (BH)", "FALSE" = "n.s."),
                       name = "Granger link") +
    scale_color_manual(values = c("Peak delay (step 1)" = "gray40",
                                  "Most predictive lag (step 2)" = "deepskyblue4"), name = "") +
    facet_wrap(~ sink_zone) +
    labs(x = "Distance between region centroids (km)", y = "Weeks behind Southern region") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_delay_distance.png"), p_delay, width = 10,
         height = 5, units = 'in', dpi = 300)
}
