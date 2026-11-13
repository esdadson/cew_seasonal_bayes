// =============================================================================
// step1_stan/spacetime_simple.stan
// The zone-free space-time model from standard parts (Gaussian):
//
//   y_n  = x_n' beta + u_t + g_i + v_i + a_i cos1(w) + b_i sin1(w) + e_n
//
//   x_n' beta   intercept, shared seasonal curve (harmonics), trap type, covariates
//   u_t         year random intercept, N(0, sigma_year^2)
//   g_i         spatial random intercept: g ~ MVN(0, sigma_g^2 C)
//   v_i         independent site random intercept, N(0, sigma_v^2)
//   a_i, b_i    site random slopes on cos1 = cos(2 pi w / 52) and
//               sin1 = sin(2 pi w / 52), spatially correlated:
//               a ~ MVN(0, sigma_a^2 C), b ~ MVN(0, sigma_b^2 C)
//               (each site's seasonal peak shifted earlier / later)
//   e_n         AR(1) within site-year: corr(e_n, e_m) = phi^(weeks apart),
//               marginal sd sigma
//   C           Matern 3/2 correlation between raw site locations, one range
//               for all three spatial effects (one Cholesky per gradient)
//
// Spatial effects are non-centered (g = sigma_g * L * z) with a 1e-6 jitter
// on C: sites metres apart get the same g, a, b; they still differ by v.
// =============================================================================
data {
  int<lower=1> N;                          // training rows
  int<lower=1> P;                          // fixed-effect columns
  int<lower=1> T;                          // years (held-out year included)
  int<lower=1> S;                          // sites
  vector[N] y;
  matrix[N, P] X;
  array[N] int<lower=1, upper=52> week;
  array[N] int<lower=1, upper=T> year;
  array[N] int<lower=1, upper=S> site;
  array[S] vector[2] xy;                   // site coordinates (km)

  int<lower=0, upper=1> slopes;            // site random slopes on cos1, sin1
  int<lower=0, upper=1> ar;                // AR(1) residual
  int<lower=0, upper=1> centered_site;     // parameterisation of v
  real<lower=0> jitter;

  // AR(1): first row of each site-year, and every other row with the row
  // before it and the gap in weeks
  int<lower=0> N_first;
  array[N_first] int<lower=1, upper=N> idx_first;
  int<lower=0> N_next;
  array[N_next] int<lower=1, upper=N> idx_next;
  array[N_next] int<lower=1, upper=N> idx_prev;
  vector<lower=1>[N_next] gap;

  // priors
  vector<lower=0>[P] beta_sd;              // beta ~ N(0, beta_sd^2)
  real<lower=0> lambda_sd;                 // random-effect sds ~ exponential(lambda_sd)
  real<lower=0> lambda_range;              // PC prior: range ~ inv_gamma(1, lambda_range)
  real<lower=0> sigma_scale;               // sigma ~ half-Student-t(3, 0, sigma_scale)

  // predictions: all rows of every site-year, 52 per site-year in week order
  int<lower=0> N_all;
  matrix[N_all, P] X_all;
  array[N_all] int<lower=1, upper=52> week_all;
  array[N_all] int<lower=1, upper=T> year_all;
  array[N_all] int<lower=1, upper=S> site_all;
  int<lower=1, upper=52> win_lo;
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
  vector[T] u_z;                           // year effects / sigma_year
  real<lower=0> sigma_year;
  vector[S] z_g;                           // spatial intercept / sigma_g, before L
  real<lower=0> sigma_g;
  vector[S] v_par;                         // independent site intercept (centered) or v / sigma_v
  real<lower=0> sigma_v;
  vector[slopes ? S : 0] z_a;              // spatial slopes, before L
  vector[slopes ? S : 0] z_b;
  array[slopes] real<lower=0> sigma_a;
  array[slopes] real<lower=0> sigma_b;
  real<lower=0> range_km;
  array[ar] real<lower=-1, upper=1> phi;
  real<lower=0> sigma;                     // marginal residual sd
}
transformed parameters {
  vector[T] u = sigma_year * u_z;
}
model {
  matrix[S, S] L = cholesky_decompose(add_diag(gp_matern32_cov(xy, 1.0, range_km / 2), jitter));
  vector[S] g = sigma_g * (L * z_g);
  vector[S] v = centered_site ? v_par : sigma_v * v_par;
  vector[N] mu = X * beta + u[year] + g[site] + v[site];
  if (slopes) {
    vector[S] a = sigma_a[1] * (L * z_a);
    vector[S] b = sigma_b[1] * (L * z_b);
    mu += a[site] .* cos1[week] + b[site] .* sin1[week];
  }

  if (ar) {
    // exact AR(1) likelihood with gaps: e_n | e_prev ~ N(phi^gap e_prev, sigma^2 (1 - phi^(2 gap)))
    vector[N] e = y - mu;
    vector[N_next] phig = pow(phi[1], gap);
    target += normal_lpdf(e[idx_first] | 0, sigma);
    target += normal_lpdf(e[idx_next] | phig .* e[idx_prev], sigma * sqrt(1 - square(phig)));
  } else {
    target += normal_lpdf(y | mu, sigma);
  }

  beta ~ normal(0, beta_sd);
  u_z ~ std_normal();
  z_g ~ std_normal();
  if (centered_site) v_par ~ normal(0, sigma_v); else v_par ~ std_normal();
  if (slopes) {
    z_a ~ std_normal(); z_b ~ std_normal();
    sigma_a ~ exponential(lambda_sd); sigma_b ~ exponential(lambda_sd);
  }
  sigma_year ~ exponential(lambda_sd);
  sigma_g ~ exponential(lambda_sd);
  sigma_v ~ exponential(lambda_sd);
  range_km ~ inv_gamma(1, lambda_range);
  sigma ~ student_t(3, 0, sigma_scale);   // phi: uniform on (-1, 1)
}
generated quantities {
  vector[S] g;                             // spatial site intercept
  vector[S] v = centered_site ? v_par : sigma_v * v_par;   // independent site intercept
  vector[S] a = rep_vector(0, S);          // site slopes on cos1, sin1
  vector[S] b = rep_vector(0, S);
  real sigma_innov = ar ? sigma * sqrt(1 - square(phi[1])) : sigma;   // = brms sigma
  array[N_sy] int peak;                    // peak week of every site-year (mean, no residual)
  vector[N_test] eta_test;
  {
    matrix[S, S] L = cholesky_decompose(add_diag(gp_matern32_cov(xy, 1.0, range_km / 2), jitter));
    g = sigma_g * (L * z_g);
    if (slopes) { a = sigma_a[1] * (L * z_a); b = sigma_b[1] * (L * z_b); }
    vector[N_all] eta = X_all * beta + u[year_all] + g[site_all] + v[site_all];
    if (slopes) eta += a[site_all] .* cos1[week_all] + b[site_all] .* sin1[week_all];
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (k in (win_lo + 1):win_hi) if (eta[r0 + k] > eta[r0 + best]) best = k;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
  }
}
