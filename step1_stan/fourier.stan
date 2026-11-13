// =============================================================================
// step1_stan/fourier.stan
// Fourier (harmonic) regression of y = log(CEW + 1) with a year intercept:
//
//   y_n ~ N(x_n' beta + u_t, sigma^2)
//
//   x_n' beta  intercept, shared seasonal curve
//              sum_h [c_h cos(2 pi h w/52) + s_h sin(2 pi h w/52)], h = 1..H (H <= 2),
//              trap type, weekly EVI and GDD, six crop shares
//   u_t        year random intercept, u_t ~ N(0, sigma_year^2); the held-out
//              year has no data, so its u is drawn from that distribution
// =============================================================================
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> T;                          // years (held-out year included)
  vector[N] y;
  matrix[N, P] X;
  array[N] int<lower=1, upper=T> year;
  vector<lower=0>[P] beta_sd;              // beta ~ N(0, beta_sd^2)
  real<lower=0> lambda_sd;                 // sigma_year ~ exponential(lambda_sd)
  real<lower=0> sigma_scale;               // sigma ~ half-Student-t(3, 0, sigma_scale)

  int<lower=0> N_all;
  matrix[N_all, P] X_all;
  array[N_all] int<lower=1, upper=T> year_all;
  int<lower=1, upper=52> win_lo;
  int<lower=win_lo, upper=52> win_hi;
  int<lower=0> N_test;
  array[N_test] int<lower=1, upper=N_all> test_idx;
}
transformed data {
  int N_sy = N_all %/% 52;
}
parameters {
  vector[P] beta;
  vector[T] u_z;                           // u / sigma_year (non-centered)
  real<lower=0> sigma_year;
  real<lower=0> sigma;
}
transformed parameters {
  vector[T] u = sigma_year * u_z;
}
model {
  y ~ normal(X * beta + u[year], sigma);
  beta ~ normal(0, beta_sd);
  u_z ~ std_normal();
  sigma_year ~ exponential(lambda_sd);
  sigma ~ student_t(3, 0, sigma_scale);
}
generated quantities {
  vector[N] log_lik;
  array[N_sy] int peak;
  vector[N_test] eta_test;
  {
    vector[N] mu = X * beta + u[year];
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    vector[N_all] eta = X_all * beta + u[year_all];
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (k in (win_lo + 1):win_hi) if (eta[r0 + k] > eta[r0 + best]) best = k;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
  }
}
