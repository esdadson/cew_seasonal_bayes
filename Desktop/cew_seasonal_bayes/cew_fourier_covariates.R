library(dplyr)
library(lubridate)
library(ISOweek)
library(tidyverse)
library(zoo)
library(brms)
library(tidyverse)
library(cmdstanr)
library(rstan)
library(tidybayes)

options(brms.backend = "cmdstanr")

dat <- read.csv("Z:/Data/Original/pest-occurrence/corn_earworm/cew_landscape_metrics_full_timeseries.csv")
dat = dat %>% group_by(longitude, latitude) %>% mutate(location2 = cur_group_id())
dat$location3 = paste(dat$location, dat$location2, sep = "_")

dat5 <- dat %>%
  group_by(location2, year, woy) %>%
  summarise(
    location3 = first(location3),
    location = first(location),
    longitude = first(longitude), latitude = first(latitude),
    CEW_sum = as.integer(round(mean(CEW_sum, na.rm = TRUE))),
    date = first(date),
    unix_date = first(unix_date),
    trap_type = first(trap_type),
    zone_30y = first(zone_30y),
    doy = first(doy),
    agdd = first(agdd), gdd = first(gdd), tmax = first(tmax), tmin = first(tmin), 
    Corn_pland = first(Corn_pland), Corn_ed = first(Corn_ed), Corn_lpi = first(Corn_lpi), 
    Corn_pd = first(Corn_pd), Cotton_pland = first(Cotton_pland), 
    Cotton_ed = first(Cotton_ed), Cotton_lpi = first(Cotton_lpi), Cotton_pd = first(Cotton_pd),
    Soybeans_pland = first(Soybeans_pland), Soybeans_ed = first(Soybeans_ed), 
    Soybeans_lpi = first(Soybeans_lpi), Soybeans_pd = first(Soybeans_pd), 
    Tobacco_pland = first(Tobacco_pland), Tobacco_ed = first(Tobacco_ed), 
    Tobacco_lpi = first(Tobacco_lpi), Tobacco_pd = first(Tobacco_pd), 
    Peanuts_pland = first(Peanuts_pland), Peanuts_ed = first(Peanuts_ed), 
    Peanuts_lpi = first(Peanuts_lpi), Peanuts_pd = first(Peanuts_pd), 
    Other_Crops_pland = first(Other_Crops_pland), Other_Crops_ed = first(Other_Crops_ed),
    Other_Crops_lpi = first(Other_Crops_lpi), Other_Crops_pd = first(Other_Crops_pd), 
    Non_Crop_pland = first(Non_Crop_pland), Non_Crop_ed = first(Non_Crop_ed), 
    Non_Crop_lpi = first(Non_Crop_lpi), Non_Crop_pd = first(Non_Crop_pd),
    .groups = "drop"
  )

dat6 <- dat5 %>% filter(year >= 1990 & year < 2021) %>%
  mutate(date = ISOweek::ISOweek2date(paste0(year, "-W", sprintf("%02d", woy), "-3")),
         unix_date = as.Date(date), doy = yday(date),
         s1 = sin(2*pi*woy/52), 
         s2 = sin(2*pi*woy*2/52), 
         s3 = sin(2*pi*woy*3/52),
         c1 = cos(2*pi*woy/52), 
         c2 = cos(2*pi*woy*2/52), 
         c3 = cos(2*pi*woy*3/52),
         log_CEW = log(CEW_sum + 1))

loc_year_min5 <- dat6 %>%
  group_by(location2, year) %>%
  summarise(weeks_sampled = sum(!is.na(CEW_sum)),
            cew_sampled = sum(CEW_sum, na.rm=T), .groups = "drop") %>%
  filter(weeks_sampled >= 5,
         cew_sampled >= 50) %>%
  group_by(location2) %>%
  mutate(years_sampled = length(unique(year))) %>%
  filter(years_sampled >= 3)

dat6 <- dat6 %>%
  semi_join(loc_year_min5, by = c("location2", "year"))

# Create numeric indices for random effects
dat6 <- dat6 %>%
  mutate(
    zone_id     = as.integer(factor(zone_30y)),
    trap_id     = as.integer(factor(trap_type))
  )

dat6$loc_year_id = as.numeric(as.factor(paste(dat6$location2, dat6$year, sep = '-')))

y_obs_flag <- as.integer(!is.na(dat6$CEW_sum))

# Read in EVI data and take mean per location x week 
evi_df = read.csv("Z:/Data/Original/pest-occurrence/corn_earworm/MCD43_EVI_points.csv")
evi_df = evi_df %>% mutate(year = year(date),
                           woy = isoweek(date))
evi_df = evi_df %>% 
  rename(location2 = id) %>%
  group_by(location2, year, woy) %>%
  summarize(EVI = mean(EVI, na.rm = T))

dat6 = merge(dat6, evi_df, by = c("location2", "year", "woy"))

dat6 = dat6 %>% 
  arrange(location2, year, woy) %>% 
  group_by(location2) %>% 
  mutate(EVI_filled = na.approx(EVI, na.rm = FALSE))

dat6 = dat6 %>% 
  filter(woy <= 52) %>% 
  group_by(location2, woy) %>% 
  mutate(EVI_mean = mean(EVI_filled, na.rm = T),
         EVI_filled = coalesce(EVI_filled, EVI_mean))

# mu_it needs a loc_year index (already have loc_year_id)
jags_data <- list(
  y = dat6$log_CEW,
  s1 = dat6$s1, c1 = dat6$c1,
  s2 = dat6$s2, c2 = dat6$c2,
  s3 = dat6$s3, c3 = dat6$c3,
  loc_year  = dat6$loc_year_id,
  N = nrow(dat6),
  N_loc_year = n_distinct(dat6$loc_year_id), 
  corn = as.numeric(scale(dat6$Corn_pland)),
  cotton = as.numeric(scale(dat6$Cotton_pland)),
  soybeans = as.numeric(scale(dat6$Soybeans_pland)),
  tobacco = as.numeric(scale(dat6$Tobacco_pland)),
  peanuts = as.numeric(scale(dat6$Peanuts_pland)),
  other = as.numeric(scale(dat6$Other_Crops_pland)),
  non_crop = as.numeric(scale(dat6$Non_Crop_pland)),
  gdd = as.numeric(scale(dat6$gdd)),
  EVI_filled = as.numeric(scale(dat6$EVI_filled))
)

dat6_brms <- dat6 %>% ungroup() %>%
  mutate(
    location2 = as.factor(location2),
    year_id = factor(year),
    zone_fac = factor(zone_30y),
    log_CEW      = log(CEW_sum + 1),
    lon_scaled   = as.numeric(scale(longitude)),
    lat_scaled   = as.numeric(scale(latitude)),
    year_scaled = as.numeric(scale(year)),
    corn = as.numeric(scale(Corn_pland)),
    cotton = as.numeric(scale(Cotton_pland)),
    soybeans = as.numeric(scale(Soybeans_pland)),
    tobacco = as.numeric(scale(Tobacco_pland)),
    peanuts = as.numeric(scale(Peanuts_pland)),
    other = as.numeric(scale(Other_Crops_pland)),
    non_crop = as.numeric(scale(Non_Crop_pland)),
    gdd = as.numeric(scale(gdd)),
    EVI_filled = as.numeric(scale(EVI_filled))
  )

formula_gp <- bf(
  log_CEW ~ 
    EVI_filled + gdd + 
    corn + cotton + soybeans + tobacco + peanuts + other + 
    (0 + s1 + c1 + s2 + c2 + s3 + c3 | zone_fac) + 
    gp(lon_scaled, lat_scaled, cov = "exponential") + 
    (1 | location2) + 
    (1 | year),
  family = gaussian()
)

priors_gp <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(exponential(1), class = "sdgp"),
  prior(inv_gamma(5, 5), class = "lscale")
)

fit_gp3 <- brm(
  formula = formula_gp,
  data = dat6_brms,
  prior = priors_gp,
  chains = 2,
  iter = 4000,
  warmup = 2000,
  control = list(
    adapt_delta = 0.95
  ),
  cores = 2,
  seed = 42,
  save_pars = save_pars(all = TRUE),
  file = "C:/Users/edadson/Desktop/cew/brms_fit/fit_gp3"
)

summary(fit_gp3)

# Stan version of the code 
stan_data = standata(dat6_brms)

model_code = "
// generated with brms 2.22.13
functions {
  /* compute a latent Gaussian process with exponential kernel
   * also known as the Matern 1/2 kernel
   * Args:
   *   x: array of continuous predictor values
   *   sdgp: marginal SD parameter
   *   lscale: length-scale parameter
   *   zgp: vector of independent standard normal variables
   * Returns:
   *   a vector to be added to the linear predictor
   */
  vector gp_exponential(data array[] vector x, real sdgp, vector lscale, vector zgp) {
    int Dls = rows(lscale);
    int N = size(x);
    matrix[N, N] cov;
    if (Dls == 1) {
      // one dimensional or isotropic GP
      cov = gp_exponential_cov(x, sdgp, lscale[1]);
    } else {
      // multi-dimensional non-isotropic GP
      cov = gp_exponential_cov(x[, 1], sdgp, lscale[1]);
      for (d in 2:Dls) {
        cov = cov .* gp_exponential_cov(x[, d], 1, lscale[d]);
      }
    }
    for (n in 1:N) {
      // deal with numerical non-positive-definiteness
      cov[n, n] += 1e-12;
    }
    return cholesky_decompose(cov) * zgp;
  }
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data related to GPs
  int<lower=1> Kgp_1;  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Dgp_1;  // GP dimension
  // number of latent GP groups
  int<lower=1> Nsubgp_1;
  // indices of latent GP groups per observation
  array[N] int<lower=1> Jgp_1;
  array[Nsubgp_1] vector[Dgp_1] Xgp_1;  // covariates of the GP
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  vector[N] Z_3_3;
  vector[N] Z_3_4;
  vector[N] Z_3_5;
  vector[N] Z_3_6;
  int<lower=1> NC_3;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector<lower=0>[Kgp_1] sdgp_1;  // GP standard deviation parameters
  array[Kgp_1] vector<lower=0>[1] lscale_1;  // GP length-scale parameters
  vector[Nsubgp_1] zgp_1;  // latent variables of the GP
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // standardized group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_1;
  vector[N_3] r_3_2;
  vector[N_3] r_3_3;
  vector[N_3] r_3_4;
  vector[N_3] r_3_5;
  vector[N_3] r_3_6;
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  // compute actual group-level effects
  r_3 = scale_r_cor(z_3, sd_3, L_3);
  r_3_1 = r_3[, 1];
  r_3_2 = r_3[, 2];
  r_3_3 = r_3[, 3];
  r_3_4 = r_3[, 4];
  r_3_5 = r_3[, 5];
  r_3_6 = r_3[, 6];
  lprior += normal_lpdf(b | 0, 1);
  lprior += normal_lpdf(Intercept | 0, 1);
  lprior += exponential_lpdf(sdgp_1 | 1);
  lprior += inv_gamma_lpdf(lscale_1[1][1] | 1.494197, 0.056607);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_3 | 3, 0, 2.5)
    - 6 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_3 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[Nsubgp_1] gp_pred_1 = gp_exponential(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + gp_pred_1[Jgp_1];
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_3_2[J_3[n]] * Z_3_2[n] + r_3_3[J_3[n]] * Z_3_3[n] + r_3_4[J_3[n]] * Z_3_4[n] + r_3_5[J_3[n]] * Z_3_5[n] + r_3_6[J_3[n]] * Z_3_6[n];
    }
    target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zgp_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
  target += std_normal_lpdf(to_vector(z_3));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
}
"

stan1.1 <- stan(
  data = stan_data,
  model_code = model_code,
  cores = 2, seed = 42)

print(stan1.1, probs = c(0.025, 0.975))