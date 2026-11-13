// =============================================================================
// step1_stan/fixed.stan
// Fixed-effects regression of y = log(CEW + 1):
//
//   y_n ~ N(x_n' beta, sigma^2)
//
//   x_n   intercept, zone intercept (Transitional, Northern against Southern),
//         trap type, weekly EVI and GDD, six crop shares
//
// No seasonal curve: a peak week can only come from the covariates.
// =============================================================================
data {
  int<lower=1> N;                          // training rows
  int<lower=1> P;                          // fixed-effect columns
  vector[N] y;
  matrix[N, P] X;
  vector<lower=0>[P] beta_sd;              // beta ~ N(0, beta_sd^2)
  real<lower=0> sigma_scale;               // sigma ~ half-Student-t(3, 0, sigma_scale)

  // predictions: weeks 1-52 of every site-year, in week order
  int<lower=0> N_all;
  matrix[N_all, P] X_all;
  int<lower=1, upper=52> win_lo;           // peak searched in weeks win_lo..win_hi
  int<lower=win_lo, upper=52> win_hi;
  int<lower=0> N_test;                     // held-out rows (indices into N_all)
  array[N_test] int<lower=1, upper=N_all> test_idx;
}
transformed data {
  int N_sy = N_all %/% 52;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(X * beta, sigma);
  beta ~ normal(0, beta_sd);
  sigma ~ student_t(3, 0, sigma_scale);
}
generated quantities {
  vector[N] log_lik;                       // pointwise log-likelihood (PSIS-LOO)
  array[N_sy] int peak;                    // peak week of every site-year
  vector[N_test] eta_test;                 // mean of the held-out rows
  {
    vector[N] mu = X * beta;
    for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | mu[n], sigma);
    vector[N_all] eta = X_all * beta;
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (k in (win_lo + 1):win_hi) if (eta[r0 + k] > eta[r0 + best]) best = k;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
  }
}
