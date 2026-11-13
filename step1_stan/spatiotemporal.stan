// =============================================================================
// step1_stan/spatiotemporal.stan
// Fourier spatial regression of y = log(CEW + 1) with an AR(1) residual:
//
//   y_n = x_n' beta + u_t + g_i + v_i + a_i cos1(w) + b_i sin1(w) + e_n
//
//   x_n' beta   intercept, shared seasonal curve (H <= 2 harmonics), trap type,
//               weekly EVI and GDD, six crop shares
//   u_t         year random intercept, N(0, sigma_year^2)
//   g_i         spatial site intercept:   g ~ MVN(0, sigma_g^2 C)
//   v_i         independent site intercept, N(0, sigma_v^2)
//   a_i, b_i    spatial site slopes on cos1 = cos(2 pi w/52), sin1 = sin(2 pi w/52):
//               a ~ MVN(0, sigma_a^2 C), b ~ MVN(0, sigma_b^2 C). They shift and
//               reshape each site's seasonal curve, so peak timing varies in space.
//   e_n         AR(1) over weeks within each site-year: corr(e_n, e_m) =
//               phi^(weeks apart), marginal sd sigma
//   C           Matern 3/2 correlation between the raw site locations, one range
//               for g, a and b (one Cholesky factor per gradient)
//
// All site effects are non-centered (g = sigma_g * L * z_g). C gets a 1e-6
// jitter: sites metres apart share g, a, b and still differ by v.
// =============================================================================
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> T;
  int<lower=1> S;                          // sites
  vector[N] y;
  matrix[N, P] X;
  array[N] int<lower=1, upper=52> week;
  array[N] int<lower=1, upper=T> year;
  array[N] int<lower=1, upper=S> site;
  array[S] vector[2] xy;                   // site coordinates (km)
  real<lower=0> jitter;

  // AR(1): the first row of each site-year, and every other row with the row
  // before it in the same site-year and the gap between them in weeks
  int<lower=0> N_first;
  array[N_first] int<lower=1, upper=N> idx_first;
  int<lower=0> N_next;
  array[N_next] int<lower=1, upper=N> idx_next;
  array[N_next] int<lower=1, upper=N> idx_prev;
  vector<lower=1>[N_next] gap;

  vector<lower=0>[P] beta_sd;              // beta ~ N(0, beta_sd^2)
  real<lower=0> lambda_sd;                 // random-effect sds ~ exponential(lambda_sd)
  real<lower=0> lambda_range;              // PC prior: range ~ inv_gamma(1, lambda_range)
  real<lower=0> sigma_scale;               // sigma ~ half-Student-t(3, 0, sigma_scale)

  int<lower=0> N_all;
  matrix[N_all, P] X_all;
  array[N_all] int<lower=1, upper=52> week_all;
  array[N_all] int<lower=1, upper=T> year_all;
  array[N_all] int<lower=1, upper=S> site_all;
  int<lower=1, upper=52> win_lo;
  int<lower=win_lo, upper=52> win_hi;
  int<lower=0> N_test;
  array[N_test] int<lower=1, upper=N_all> test_idx;
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
  vector[S] z_g;                           // spatial effects before scaling: g = sigma_g * L * z_g
  vector[S] z_a;
  vector[S] z_b;
  vector[S] v_z;                           // v / sigma_v
  real<lower=0> sigma_g;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  real<lower=0> sigma_v;
  real<lower=0> range_km;                  // Matern 3/2 range (km); length-scale = range / 2
  real<lower=-1, upper=1> phi;             // AR(1) correlation between consecutive weeks
  real<lower=0> sigma;                     // marginal residual sd
}
transformed parameters {
  vector[T] u = sigma_year * u_z;
}
model {
  matrix[S, S] L = cholesky_decompose(add_diag(gp_matern32_cov(xy, 1.0, range_km / 2), jitter));
  vector[S] g = sigma_g * (L * z_g);
  vector[S] a = sigma_a * (L * z_a);
  vector[S] b = sigma_b * (L * z_b);
  vector[S] v = sigma_v * v_z;
  vector[N] mu = X * beta + u[year] + g[site] + v[site]
                 + a[site] .* cos1[week] + b[site] .* sin1[week];

  // exact AR(1) likelihood with gaps: e_n | e_prev ~ N(phi^gap e_prev, sigma^2 (1 - phi^(2 gap)))
  vector[N] e = y - mu;
  vector[N_next] phig = pow(phi, gap);
  target += normal_lpdf(e[idx_first] | 0, sigma);
  target += normal_lpdf(e[idx_next] | phig .* e[idx_prev], sigma * sqrt(1 - square(phig)));

  beta ~ normal(0, beta_sd);
  u_z ~ std_normal();
  z_g ~ std_normal();
  z_a ~ std_normal();
  z_b ~ std_normal();
  v_z ~ std_normal();
  sigma_year ~ exponential(lambda_sd);
  sigma_g ~ exponential(lambda_sd);
  sigma_a ~ exponential(lambda_sd);
  sigma_b ~ exponential(lambda_sd);
  sigma_v ~ exponential(lambda_sd);
  range_km ~ inv_gamma(1, lambda_range);
  sigma ~ student_t(3, 0, sigma_scale);    // phi: uniform on (-1, 1)
}
generated quantities {
  vector[S] g;
  vector[S] a;
  vector[S] b;
  vector[S] v = sigma_v * v_z;
  real sigma_innov = sigma * sqrt(1 - square(phi));   // sd of the week-to-week innovation
  vector[N] log_lik;                       // conditional on the previous week (AR(1))
  array[N_sy] int peak;                    // peak week of every site-year (mean, no residual)
  vector[N_test] eta_test;
  {
    matrix[S, S] L = cholesky_decompose(add_diag(gp_matern32_cov(xy, 1.0, range_km / 2), jitter));
    g = sigma_g * (L * z_g);
    a = sigma_a * (L * z_a);
    b = sigma_b * (L * z_b);

    vector[N] e = y - (X * beta + u[year] + g[site] + v[site]
                       + a[site] .* cos1[week] + b[site] .* sin1[week]);
    for (n in idx_first) log_lik[n] = normal_lpdf(e[n] | 0, sigma);
    for (j in 1:N_next) {
      real pg = pow(phi, gap[j]);
      log_lik[idx_next[j]] = normal_lpdf(e[idx_next[j]] | pg * e[idx_prev[j]], sigma * sqrt(1 - square(pg)));
    }

    vector[N_all] eta = X_all * beta + u[year_all] + g[site_all] + v[site_all]
                        + a[site_all] .* cos1[week_all] + b[site_all] .* sin1[week_all];
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (k in (win_lo + 1):win_hi) if (eta[r0 + k] > eta[r0 + best]) best = k;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
  }
}
