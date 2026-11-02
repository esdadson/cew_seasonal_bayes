# Publication figures for JIPM
# 2. Get map boundaries for countries and Great Lakes

# Figure 1
library(dplyr)
library(ggplot2)
library(sf)
library(cowplot)

us_states <- map_data("state")
canada_prov <- map_data("world", region = c("Canada", "Mexico"))
lakes <- map_data("lakes")

# Sites of the step-1 analysis, from the same code that builds the model data
# (step1/step1_data.R: 2008-2020, longitude >= -105, site-years with >= 5
# sampled weeks at sites with >= 3 such years). Creates counts_raw, sites,
# keep_ly (site-years with >= 5 weeks), d_obs and root.
if (!exists("d_obs")) local({
  owd = setwd("C:/Users/edadson/Desktop/cew_seasonal_bayes"); on.exit(setwd(owd))
  source("step1/step1_data.R")
})

# Full Lawton et al. (2022) dataset (incl. sites west of 105W)
lawton = read.csv(file.path(root, "pest-occurrence/corn_earworm/doi_10_5061_dryad_m0cfxpp5x__v20221206/Lawton_etal_data.csv"))

# Years sampled = years 2008-2020 with >= 5 sampled weeks; selected = in the fit
site_nyr = keep_ly %>% count(location2, name = "nyr")
fig1_sites = lawton %>%
  distinct(longitude, latitude) %>%
  left_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude")) %>%
  left_join(site_nyr, by = "location2") %>%
  mutate(nyr = coalesce(nyr, 0L),
         selected = location2 %in% d_obs$location2)

# JIPM print size: measure the double-column width on a recent JIPM PDF.
# Text >= 7 pt and lines 0.25-1 pt at this size (OUP artwork guidance);
# ggplot linewidth is in mm, so 0.25-1 pt = linewidth 0.12-0.47.
fig_width = 6.8                                   # in, double column
fig_font = "Arial"
yrs_scale = scale_color_viridis_c(name = "Years\nsampled", limits = c(3, 13), breaks = c(3, 5, 7, 9, 11, 13),
                                  end = 0.9, aesthetics = c("color", "fill"),   # no pale yellow
                                  na.value = "gray70")

fig1_map = ggplot() +
  geom_polygon(data = canada_prov, aes(x = long, y = lat, group = group),
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.15) +
  geom_polygon(data = us_states, aes(x = long, y = lat, group = group),
               fill = "#f9f9f9", color = "#9c9c9c", linewidth = 0.15) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group),
               fill = "#cbdbee", color = "#9c9c9c", linewidth = 0.15) +
  # Unselected sites gray crosses (shape as well as color), selected on top
  geom_point(data = filter(fig1_sites, !selected), aes(x = longitude, y = latitude, shape = "Not selected"),
             color = "gray55", size = 0.6, stroke = 0.3) +
  geom_point(data = filter(fig1_sites, selected) %>% arrange(nyr),
             aes(x = longitude, y = latitude, color = nyr, shape = "Selected"), size = 1.2, alpha = 0.9) +
  scale_shape_manual(values = c("Selected" = 16, "Not selected" = 4), name = NULL) +
  yrs_scale +
  guides(shape = guide_legend(order = 2, override.aes = list(size = 1.5, color = c("gray55", "black"))),
         color = guide_colorbar(order = 1, barheight = unit(25, "mm"), barwidth = unit(3, "mm"))) +
  coord_sf(crs = st_crs(4269), xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_void(base_size = 8, base_family = fig_font) +
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0))

# Years sampled per site (sites with >= 1 qualifying year); bars < 3 years gray
fig1_hist = ggplot(fig1_sites %>% filter(nyr > 0) %>% count(nyr),
                   aes(x = nyr, y = n, fill = ifelse(nyr >= 3, nyr, NA))) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "gray93") +
  annotate("text", x = 1.5, y = Inf, label = "Not selected", vjust = 1.3,
           size = 7 / .pt, family = fig_font, color = "gray30") +
  geom_col(width = 0.85) +
  geom_text(aes(label = n), vjust = -0.4, size = 7 / .pt, family = fig_font) +
  yrs_scale +
  scale_x_continuous(breaks = 1:13, expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Years sampled (2008\u20132020)", y = "Sites") +
  theme_classic(base_size = 8, base_family = fig_font) +
  theme(legend.position = "none", axis.text = element_text(size = 7, color = "black"),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
        plot.margin = margin(12, 5.5, 5.5, 12))

# Data-selection flowchart; numbers computed from the data so they match the fit
fmt = function(x) format(x, big.mark = ",")
fig1_flow_txt = c(
  sprintf("Full Lawton et al. (2022) dataset\nn = %s, sites = %s",
          fmt(nrow(lawton)), fmt(nrow(distinct(lawton, longitude, latitude)))),
  sprintf("Analysis sites: 2008\u20132020, east of 105\u00b0W,\n\u2265 5 weeks sampled per year, \u2265 3 such years\nn = %s, site-years = %s, sites = %s",
          fmt(nrow(d_obs)), fmt(n_distinct(d_obs$loc_year)), fmt(n_distinct(d_obs$location2))))
fig1_flow = ggplot() +
  annotate("label", x = 0.5, y = c(0.8, 0.2), label = fig1_flow_txt, size = 7 / .pt,
           family = fig_font, lineheight = 0.95, linewidth = 0.4,
           label.padding = unit(2, "mm"), label.r = unit(1.5, "mm")) +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.66, yend = 0.37, linewidth = 0.4,
           arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_void()

fig1 = plot_grid(
  plot_grid(fig1_flow, fig1_hist, ncol = 2, rel_widths = c(1, 1.15), labels = c("a", "b"),
            label_size = 9, label_fontfamily = fig_font),
  fig1_map, ncol = 1, rel_heights = c(1, 1.45), labels = c("", "c"),
  label_size = 9, label_fontfamily = fig_font)

# Vector PDF (fonts embedded by cairo) + 600 dpi TIFF (combination art)
ggsave("C:/Users/edadson/Desktop/CEW/fig1.pdf", fig1, width = fig_width, height = 5.4,
       units = "in", device = cairo_pdf, bg = "white")
ggsave("C:/Users/edadson/Desktop/CEW/fig1.tif", fig1, width = fig_width, height = 5.4,
       units = "in", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")

# Figure 2
p = ggplot() +
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
#   Step 1 (results/step1/, from step1/step1_hpc.R compare): model comparison
#     (LOO, R^2, MAE), out-of-sample prediction of the held-out year 2020,
#     predicted weekly curves (weeks 1-52) and peak week.
#   Step 2 (results/step2/, from step2/step2_granger.R): regions, Southern ->
#     Transitional/Northern Granger causality and lag times, the held-out
#     year, and the two DAGs (South -> North vs South -> Transitional -> North).
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
model_levels = c("fixed", "spatial", "fourier", "fourier_site")
model_labels = c("Fixed effects", "+ Spatial GP", "+ Seasonal curves (zone)", "+ Site curves")

if (file.exists(file.path(step1_dir, "model_comparison.csv"))) {
  # Model comparison: fit to 2008-2019 (PSIS-LOO) and prediction of the
  # held-out year. ELPD difference from the best model (+/- 2 SE); LOO R^2
  # (95% CrI); MAE of the LOO and of the held-out predictions; held-out R^2.
  # Filled = the selected model (simplest within 2 SE of the best ELPD).
  cmp = read.csv(file.path(step1_dir, "model_comparison.csv")) %>%
    mutate(model = factor(model, levels = model_levels, labels = model_labels))
  ho_lab = function(x) paste0(x, " (", if ("n_test" %in% names(cmp)) "2020" else "n/a", ")")
  cmp_long = bind_rows(
    cmp %>% transmute(model, selected, metric = "ELPD difference (LOO)", est = elpd_diff,
                      lo = elpd_diff - 2 * se_diff, hi = elpd_diff + 2 * se_diff),
    cmp %>% transmute(model, selected, metric = "R\u00b2 (LOO)", est = loo_R2, lo = loo_R2_lo, hi = loo_R2_hi),
    cmp %>% transmute(model, selected, metric = "MAE (LOO)", est = mae_loo, lo = NA, hi = NA),
    if ("R2_test" %in% names(cmp)) cmp %>% transmute(model, selected, metric = ho_lab("R\u00b2"), est = R2_test,
                                                     lo = NA, hi = NA),
    if ("mae_test" %in% names(cmp)) cmp %>% transmute(model, selected, metric = ho_lab("MAE"), est = mae_test,
                                                      lo = NA, hi = NA)) %>%
    mutate(metric = factor(metric, levels = unique(metric)))

  p_cmp = ggplot(cmp_long, aes(x = est, y = model)) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", color = "deepskyblue4",
                  na.rm = TRUE) +
    geom_point(aes(shape = selected), size = 2.4, color = "deepskyblue4") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       name = "Selected") +
    facet_wrap(~ metric, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))

  ggsave(filename = file.path(fig_dir, "step1_model_comparison.png"), p_cmp, width = 13,
         height = 3.2, units = 'in', dpi = 300)
}

if (file.exists(file.path(step1_dir, "holdout_zone_week.csv"))) {
  # Out-of-sample validation: the held-out year's zone weekly mean of the
  # trapped sites, observed (points) vs each model's prediction from
  # 2008-2019 (line, 95% CrI of the mean; the year effect is drawn, unknown)
  hzw = read.csv(file.path(step1_dir, "holdout_zone_week.csv")) %>%
    mutate(model = factor(model, levels = model_levels, labels = model_labels),
           zone = factor(zone_30y, levels = c(2, 1, 0), labels = zone_labels[c("2", "1", "0")]))
  ho_txt = cmp %>% transmute(model, label = sprintf("MAE %.2f\nR\u00b2 %.2f", mae_test, R2_test),
                             zone = factor(levels(hzw$zone)[1], levels = levels(hzw$zone)))   # top row

  p_ho = ggplot(hzw, aes(x = woy)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = zone), alpha = 0.25) +
    geom_line(aes(y = med, color = zone), linewidth = 0.7) +
    geom_point(aes(y = obs_mean, size = n_sites), shape = 21, fill = "white", color = "black", stroke = 0.4) +
    geom_text(data = ho_txt, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.2, size = 2.8) +
    scale_color_manual(values = setNames(zone_colors[c("2", "1", "0")], zone_labels[c("2", "1", "0")]),
                       name = "Overwintering \nzone") +
    scale_fill_manual(values = setNames(zone_colors[c("2", "1", "0")], zone_labels[c("2", "1", "0")]),
                      name = "Overwintering \nzone") +
    scale_size_area(max_size = 2.5, name = "Sites \ntrapped") +
    facet_grid(zone ~ model) +
    labs(x = "Week of Year", y = paste0("Mean log(CEW + 1), ", hzw$year[1]),
         title = paste0("Out-of-sample prediction of ", hzw$year[1], " (models fit to the other years)")) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))

  ggsave(filename = file.path(fig_dir, "step1_holdout_prediction.png"), p_ho, width = 13,
         height = 7.5, units = 'in', dpi = 300)
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
  # is positively associated with at its best lag (color = share of the
  # sink's remaining variation explained, 2008-2019); bold = permutation
  # test p < 0.05 (BH)
  links = pairs %>%
    filter(coef_sum > 0) %>%
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
    scale_color_viridis_c(name = "Variation \nexplained (%)", option = "viridis", end = 0.85) +
    scale_linewidth_manual(values = c("p < 0.05 (BH)" = 1.1, "n.s." = 0.4), name = "Permutation \ntest") +
    scale_alpha_manual(values = c("p < 0.05 (BH)" = 0.95, "n.s." = 0.6), name = "Permutation \ntest") +
    scale_fill_manual(values = zone_colors, labels = zone_labels, name = "Overwintering \nzone") +
    facet_wrap(~ sink_zone) +
    map_frame

  ggsave(filename = file.path(fig_dir, "granger_links.png"), p_links, width = 14,
         height = 6, units = 'in', dpi = 300)

  # Lag map: each sink region filled by the lag (weeks) at which its strongest
  # Southern source best predicts it; faded = not significant; unfilled = no
  # positive Southern association
  best_src = pairs %>%
    filter(coef_sum > 0) %>%
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

  # Lag profiles: variation explained by each Southern lag (2008-2019), one
  # panel per sink region
  lag_prof = g2$granger_reg %>%
    mutate(sig = !is.na(p_perm) & p_perm < 0.05)

  p_prof = ggplot(lag_prof, aes(x = lag, y = 100 * gain, color = source)) +
    geom_line() +
    geom_point(aes(shape = sig), size = 1.8) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Permutation \np < 0.05") +
    scale_x_continuous(breaks = 1:8) +
    facet_wrap(~ sink) +
    labs(x = "Southern lag (weeks)", y = "Variation explained (%)", color = "Southern \nregion") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_lag_profiles.png"), p_prof, width = 10,
         height = 7, units = 'in', dpi = 300)

  # Zone level: lags (forward and reverse) and leads (falsification)
  gz = g2$granger_zone %>%
    mutate(pair = factor(paste(source, "->", sink),
                         levels = c("S -> T", "S -> N", "T -> S", "N -> S")),
           sig = !is.na(p_perm) & p_perm < 0.05)

  p_zone = ggplot(gz, aes(x = lag, y = 100 * gain, color = direction)) +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_line() +
    geom_point(aes(shape = sig), size = 2) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Permutation \np < 0.05") +
    scale_color_manual(values = c(forward = "deepskyblue4", reverse = "gray50"), name = "") +
    scale_x_continuous(breaks = -3:8) +
    facet_wrap(~ pair, ncol = 2) +
    labs(x = "Lag of source zone (weeks; negative = lead)", y = "Variation explained (%), 2008-2019") +
    theme_minimal()

  ggsave(filename = file.path(fig_dir, "granger_zone_lags.png"), p_zone, width = 8,
         height = 6, units = 'in', dpi = 300)

  # Held-out year at each sink zone's best Southern lag: observed vs own
  # history only vs + Southern lag (one week ahead), fit to 2008-2019
  oos = g2$granger_pred_zone %>%
    semi_join(g2$best_zone %>% select(source, sink, setting), by = c("source", "sink", "setting")) %>%
    left_join(g2$best_zone %>% select(sink, gain_holdout), by = "sink") %>%
    transmute(sink = factor(sink, levels = c("T", "N"), labels = c("Transitional zone", "Northern zone")),
              panel = sprintf("%s (Southern lag %s wk; out-of-sample gain %.0f%%)", sink,
                              sub("lag", "", setting), 100 * gain_holdout),
              year, woy,
              `Observed` = clim + anom,
              `Own history` = clim + pred_R,
              `+ Southern lag` = clim + pred_U) %>%
    pivot_longer(c(`Observed`, `Own history`, `+ Southern lag`), names_to = "series") %>%
    mutate(series = factor(series, levels = c("Observed", "Own history", "+ Southern lag")))

  if (nrow(oos)) {
    p_oos = ggplot(oos, aes(x = woy, y = value, color = series)) +
      geom_line() +
      geom_point(data = oos %>% filter(series == "Observed"), size = 1) +
      scale_color_manual(values = c("Observed" = "black", "Own history" = "gray60",
                                    "+ Southern lag" = "deepskyblue4"), name = "") +
      facet_wrap(~ panel, ncol = 1) +
      labs(x = "Week of Year", y = "Zone mean log(CEW + 1), trapped sites",
           title = paste0("Held-out year ", oos$year[1], ": one-week-ahead predictions")) +
      theme_minimal()
    ggsave(filename = file.path(fig_dir, "granger_holdout.png"), p_oos, width = 8,
           height = 6, units = 'in', dpi = 300)
  }

  # DAGs: South -> North directly vs South -> Transitional -> North. Each
  # test's variation explained (2008-2019); filled = permutation p < 0.05.
  # S -> N | T = the direct path, T -> N | S = the path through the
  # Transitional zone. Zone level + each Northern region with its strongest
  # Southern and Transitional sources (exploratory).
  dag = read.csv(file.path(step2_dir, "granger_dag.csv")) %>%
    mutate(sig = !is.na(p_perm) & p_perm < 0.05,
           test = factor(test, levels = c("S -> N", "T -> N", "S -> N | T", "T -> N | S")),
           panel = ifelse(level == "zone", paste0("Zones: ", verdict),
                          paste0(sink, " (", S_source, ", ", T_source, "): ", verdict)),
           panel = factor(panel, levels = unique(panel[order(level != "zone", sink)])))   # zones first
  p_dag = ggplot(dag, aes(x = test, y = 100 * gain)) +
    geom_col(aes(fill = sig), width = 0.6, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("TRUE" = "deepskyblue4", "FALSE" = "gray85"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."), name = "Permutation \ntest") +
    facet_wrap(~ panel, ncol = 2) +
    labs(x = NULL, y = "Variation explained (%), 2008-2019",
         title = "Direct (S -> N | T) vs through the Transitional zone (T -> N | S)") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 9))

  ggsave(filename = file.path(fig_dir, "granger_dag.png"), p_dag, width = 9,
         height = 2.2 + 2 * ceiling(n_distinct(dag$panel) / 2), units = 'in', dpi = 300)

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
