// =============================================================================
// step1_stan/zone_harmonic.stan
// Zone harmonic seasonal curves + short-range spatial site effects (Gaussian):
//
//   y = X beta + u_t + w_i + svc * (alpha_i cos1(w) + beta_i sin1(w)) + e
//
//   X beta    zone intercepts, zone x harmonic seasonal curves, trap type,
//             covariates (columns built in zone_harmonic_common.R)
//   u_t       year intercept (held-out year: drawn from its prior)
//   w         site effect over the RAW site locations:
//               w ~ MVN(0, sigma_field^2 C(range) + sigma_site^2 I)
//             = short-range Matern field + site nugget, as one covariance.
//             The nugget keeps it well conditioned even for sites metres
//             apart, so no merging of sites is needed.
//   alpha, beta (svc = 1) site-level shifts of the first harmonic: together a
//             shift in each site's peak timing and height. Short-range
//             Matern fields with one shared sd (so the prior does not depend
//             on where the year starts); non-centered; tiny jitter.
//   Every site-level effect is centred within zone (sums to zero over the
//   sites of each zone), so zone differences belong to the zone curves and
//   cannot be absorbed by the fields.
//   One Matern 3/2 range for all fields. range = 2 * length-scale = distance
//   where the correlation is ~0.14 (INLA's definition of range).
// =============================================================================
functions {
  // subtract each zone's mean: x sums to zero over the sites of every zone
  vector zone_center(vector x, array[] int zone, int Z) {
    vector[Z] tot = rep_vector(0, Z);
    vector[Z] n = rep_vector(0, Z);
    for (i in 1:size(zone)) {
      tot[zone[i]] += x[i];
      n[zone[i]] += 1;
    }
    return x - (tot ./ n)[zone];
  }
}
data {
  int<lower=1> N;                          // training rows
  int<lower=1> P;                          // fixed-effect columns
  int<lower=1> T;                          // years (held-out year included)
  int<lower=1> S;                          // sites (raw locations)
  int<lower=1> Z;                          // zones
  vector[N] y;
  matrix[N, P] X;
  array[N] int<lower=1, upper=52> week;
  array[N] int<lower=1, upper=T> year;
  array[N] int<lower=1, upper=S> site;
  array[S] vector[2] xy;                   // site coordinates (km)
  array[S] int<lower=1, upper=Z> zone_site;

  int<lower=0, upper=1> svc;               // site-level first-harmonic fields
  int<lower=0, upper=1> centered_site;     // parameterisation of w
  real<lower=0> jitter;                    // diagonal added to C for the svc fields

  // priors
  vector<lower=0>[P] beta_sd;              // beta ~ N(0, beta_sd^2)
  real range_mu;                           // range_km ~ lognormal(range_mu, range_sd)
  real<lower=0> range_sd;
  real<lower=0> lambda_sd;                 // every random-effect sd ~ exponential(lambda_sd)
  real<lower=0> sigma_scale;               // residual sd ~ half-Student-t(3, 0, sigma_scale)

  // predictions: all rows of every site-year, 52 per site-year in week order
  int<lower=0> N_all;
  matrix[N_all, P] X_all;
  array[N_all] int<lower=1, upper=52> week_all;
  array[N_all] int<lower=1, upper=T> year_all;
  array[N_all] int<lower=1, upper=S> site_all;
  int<lower=1, upper=52> win_lo;           // peak-week search window
  int<lower=win_lo, upper=52> win_hi;
  int<lower=0> N_test;
  array[N_test] int<lower=1, upper=N_all> test_idx;
  vector[N_test] y_test;
}
transformed data {
  int N_sy = N_all %/% 52;
  vector[52] cos1;
  vector[52] sin1;
  for (k in 1:52) {
    cos1[k] = cos(2 * pi() * k / 52);
    sin1[k] = sin(2 * pi() * k / 52);
  }
}
parameters {
  vector[P] beta;
  vector[T] u_z;
  real<lower=0> sigma_year;
  vector[S] w_par;                         // w before zone centring (centered) or innovations
  real<lower=0> sigma_field;
  real<lower=0> sigma_site;
  vector[svc ? S : 0] a_z;                 // innovations of alpha, beta
  vector[svc ? S : 0] b_z;
  array[svc] real<lower=0> sigma_svc;
  real<lower=0> range_km;
  real<lower=0> sigma;
}
transformed parameters {
  vector[T] u = sigma_year * u_z;
}
model {
  matrix[S, S] C = gp_matern32_cov(xy, 1.0, range_km / 2);
  matrix[S, S] L_w = cholesky_decompose(add_diag(square(sigma_field) * C, square(sigma_site)));
  vector[S] w = zone_center(centered_site ? w_par : L_w * w_par, zone_site, Z);
  vector[N] off = u[year] + w[site];
  if (svc) {
    matrix[S, S] L_c = cholesky_decompose(add_diag(C, jitter));
    vector[S] alpha  = zone_center(sigma_svc[1] * (L_c * a_z), zone_site, Z);
    vector[S] beta_s = zone_center(sigma_svc[1] * (L_c * b_z), zone_site, Z);
    off += alpha[site] .* cos1[week] + beta_s[site] .* sin1[week];
  }
  target += normal_id_glm_lpdf(y | X, off, beta, sigma);

  beta ~ normal(0, beta_sd);
  u_z ~ std_normal();
  if (centered_site) w_par ~ multi_normal_cholesky(rep_vector(0, S), L_w);
  else w_par ~ std_normal();
  if (svc) { a_z ~ std_normal(); b_z ~ std_normal(); sigma_svc ~ exponential(lambda_sd); }
  sigma_year ~ exponential(lambda_sd);
  sigma_field ~ exponential(lambda_sd);
  sigma_site ~ exponential(lambda_sd);
  range_km ~ lognormal(range_mu, range_sd);
  sigma ~ student_t(3, 0, sigma_scale);
}
generated quantities {
  vector[S] w;                             // site effects (zone-centred)
  vector[S] alpha = rep_vector(0, S);      // site shifts of the first harmonic (cos, sin)
  vector[S] beta_s = rep_vector(0, S);
  array[N_sy] int peak;                    // peak week of every site-year (linear predictor)
  vector[N_test] eta_test;
  vector[N_test] log_lik_test;
  {
    matrix[S, S] C = gp_matern32_cov(xy, 1.0, range_km / 2);
    if (centered_site) w = zone_center(w_par, zone_site, Z);
    else w = zone_center(cholesky_decompose(add_diag(square(sigma_field) * C, square(sigma_site))) * w_par,
                         zone_site, Z);
    if (svc) {
      matrix[S, S] L_c = cholesky_decompose(add_diag(C, jitter));
      alpha  = zone_center(sigma_svc[1] * (L_c * a_z), zone_site, Z);
      beta_s = zone_center(sigma_svc[1] * (L_c * b_z), zone_site, Z);
    }
    vector[N_all] eta = X_all * beta + u[year_all] + w[site_all];
    if (svc) eta += alpha[site_all] .* cos1[week_all] + beta_s[site_all] .* sin1[week_all];
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (k in (win_lo + 1):win_hi) if (eta[r0 + k] > eta[r0 + best]) best = k;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
    for (n in 1:N_test) log_lik_test[n] = normal_lpdf(y_test[n] | eta_test[n], sigma);
  }
}
