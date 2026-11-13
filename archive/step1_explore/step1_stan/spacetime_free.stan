// =============================================================================
// step1_stan/spacetime_free.stan
// Stan mirror of the R-INLA zone-free space-time model (Gaussian):
//
//   y = X beta + s(w) + u_t + v_i + g(c_i, w) + e,     e ~ N(0, sigma^2)
//
//   s(w)      cyclic RW2 seasonal curve over weeks 1-52, scaled as INLA's
//             scale.model = TRUE, summing to ~0 (soft constraint)
//   u_t       year random intercept (held-out year: drawn from its prior)
//   v_i       site random intercept
//   g(c, w)   space-time field: g(c, w) = sum_k B[w, k] G[c, k], with
//               G[, k] = one Matern field over the field locations c
//                        (range_km and sd ESTIMATED, as INLA's pcmatern)
//               G[, k] = rho G[, k-1] + sqrt(1 - rho^2) * new field
//                        (stationary AR1 across knots, as INLA's group ar1)
//             Field locations = sites merged within a cutoff (as INLA's mesh
//             cutoff): every site in a cluster shares the cluster's field.
//
// Matern: nu = 1 (what INLA's SPDE with alpha = 2 approximates) or nu = 3/2
// (built-in, faster). Range is defined as INLA does: sqrt(8 nu) / kappa,
// where the correlation is ~0.13.
// Priors = INLA's PC priors (see step1_stan.R for their parameters).
// =============================================================================
functions {
  // Matern correlation between field locations
  matrix matern_cor(array[] vector xy, real range_km, int nu_code) {
    int M = size(xy);
    matrix[M, M] C;
    if (nu_code == 2) {
      C = gp_matern32_cov(xy, 1.0, range_km / 2);         // kappa = sqrt(3) / ell = sqrt(12) / range_km
    } else {
      real kappa = sqrt(8.0) / range_km;                  // nu = 1: C(d) = kappa d K1(kappa d)
      for (i in 1:M) {
        C[i, i] = 1;
        for (j in (i + 1):M) {
          real kd = kappa * distance(xy[i], xy[j]);
          C[i, j] = kd * modified_bessel_second_kind(1, kd);
          C[j, i] = C[i, j];
        }
      }
    }
    return C;
  }
}
data {
  // training rows
  int<lower=1> N;
  int<lower=1> P;                          // fixed-effect columns (intercept included)
  int<lower=1> T;                          // years (held-out year included)
  int<lower=1> S;                          // sites
  int<lower=1> M;                          // field locations (site clusters)
  int<lower=2> K;                          // temporal basis functions (knots)
  vector[N] y;
  matrix[N, P] X;
  array[N] int<lower=1, upper=52> week;
  array[N] int<lower=1, upper=T> year;
  array[N] int<lower=1, upper=S> site;
  array[S] int<lower=1, upper=M> cluster;  // field location of each site
  array[M] vector[2] xy;                   // field locations (km)
  matrix[52, K] B;                         // temporal basis (quadratic B-splines)

  // model choices
  int<lower=1, upper=2> nu_code;           // 1: Matern nu = 1; 2: nu = 3/2
  int<lower=0, upper=1> centered_site;     // parameterisation of v
  int<lower=0, upper=1> centered_field;    // parameterisation of G
  real<lower=0> jitter;                    // added to the correlation diagonal

  // priors
  real<lower=0> rw2_scale;                 // scale.model factor of the cyclic RW2
  vector<lower=0>[P] beta_sd;              // N(0, beta_sd^2)
  real<lower=0> lambda_season;             // sd ~ exponential(lambda): PC prior
  real<lower=0> lambda_year;
  real<lower=0> lambda_site;
  real<lower=0> lambda_field_sd;
  real<lower=0> lambda_range;              // PC prior: 1 / range_km ~ exponential(lambda_range)
  real<lower=0> theta_rho;                 // PC prior on rho (base model rho = 1)
  real<lower=0> sigma_scale;               // half-Student-t(3) scale for the residual sd

  // predictions: all rows of every monitored site-year, 52 per site-year in
  // week order, so peak weeks and held-out scores are computed per draw
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
  // linear index of (field location, week) in the M x 52 field matrix
  array[N] int cw;
  array[N_all] int cw_all;
  for (n in 1:N) cw[n] = cluster[site[n]] + (week[n] - 1) * M;
  for (n in 1:N_all) cw_all[n] = cluster[site_all[n]] + (week_all[n] - 1) * M;
  array[52] int wprev;
  array[52] int wnext;
  for (w in 1:52) { wprev[w] = w == 1 ? 52 : w - 1; wnext[w] = w == 52 ? 1 : w + 1; }
}
parameters {
  vector[P] beta;
  vector[52] s;
  real<lower=0> sigma_season;
  vector[T] u_z;
  real<lower=0> sigma_year;
  vector[S] v_par;                         // v (centered) or v / sigma_site
  real<lower=0> sigma_site;
  matrix[M, K] G_par;                      // G (centered) or standard-normal innovations
  real<lower=0> range_km;                     // km
  real<lower=0> sigma_field;
  real<lower=-1, upper=1> rho;
  real<lower=0> sigma;
}
transformed parameters {
  vector[T] u = sigma_year * u_z;
  vector[S] v = centered_site ? v_par : sigma_site * v_par;
}
model {
  real sd_innov = sigma_field * sqrt(1 - square(rho));
  matrix[M, M] L = cholesky_decompose(add_diag(matern_cor(xy, range_km, nu_code), jitter));
  matrix[M, K] G;
  if (centered_field) {
    // G[, 1] ~ N(0, sigma_field^2 C); G[, k] | G[, k-1] ~ N(rho G[, k-1], sd_innov^2 C)
    matrix[M, K] D = G_par;
    D[, 2:K] = G_par[, 2:K] - rho * G_par[, 1:(K - 1)];
    matrix[M, K] W = mdivide_left_tri_low(L, D);
    target += -0.5 * dot_self(W[, 1]) / square(sigma_field)
              - 0.5 * sum(columns_dot_self(W[, 2:K])) / square(sd_innov)
              - K * sum(log(diagonal(L)))
              - M * log(sigma_field) - M * (K - 1) * log(sd_innov);
    G = G_par;
  } else {
    matrix[M, K] LZ = L * G_par;
    to_vector(G_par) ~ std_normal();
    G[, 1] = sigma_field * LZ[, 1];
    for (k in 2:K) G[, k] = rho * G[, k - 1] + sd_innov * LZ[, k];
  }

  // linear predictor of the training rows
  vector[M * 52] gw = to_vector(G * B');   // field at every location x week
  vector[N] alpha = s[week] + u[year] + v[site] + gw[cw];
  target += normal_id_glm_lpdf(y | X, alpha, beta, sigma);

  // priors
  beta ~ normal(0, beta_sd);
  // cyclic RW2, scaled: precision rw2_scale / sigma_season^2 on the second
  // differences; rank 51 (a constant is free), removed by the soft sum-to-zero
  target += -0.5 * rw2_scale * dot_self(s[wprev] - 2 * s + s[wnext]) / square(sigma_season)
            - 51 * log(sigma_season);
  sum(s) ~ normal(0, 0.001 * 52);
  sigma_season ~ exponential(lambda_season);
  u_z ~ std_normal();
  sigma_year ~ exponential(lambda_year);
  if (centered_site) v_par ~ normal(0, sigma_site); else v_par ~ std_normal();
  sigma_site ~ exponential(lambda_site);
  sigma_field ~ exponential(lambda_field_sd);
  target += log(lambda_range) - 2 * log(range_km) - lambda_range / range_km;   // PC prior on range_km (2-d)
  target += -theta_rho * sqrt(1 - rho) - 0.5 * log1m(rho);              // PC prior on rho
  sigma ~ student_t(3, 0, sigma_scale);
}
generated quantities {
  matrix[M, K] G;                          // the space-time field at the knots
  array[N_sy] int peak;                    // peak week of every site-year
  vector[N_test] eta_test;                 // held-out rows: linear predictor
  vector[N_test] log_lik_test;             //   and log density of the observation
  {
    if (centered_field) {
      G = G_par;
    } else {
      matrix[M, M] L = cholesky_decompose(add_diag(matern_cor(xy, range_km, nu_code), jitter));
      matrix[M, K] LZ = L * G_par;
      real sd_innov = sigma_field * sqrt(1 - square(rho));
      G[, 1] = sigma_field * LZ[, 1];
      for (k in 2:K) G[, k] = rho * G[, k - 1] + sd_innov * LZ[, k];
    }
    vector[M * 52] gw = to_vector(G * B');
    vector[N_all] eta = X_all * beta + s[week_all] + u[year_all] + v[site_all] + gw[cw_all];
    for (j in 1:N_sy) {
      int r0 = (j - 1) * 52;
      int best = win_lo;
      for (w in (win_lo + 1):win_hi) if (eta[r0 + w] > eta[r0 + best]) best = w;
      peak[j] = best;
    }
    eta_test = eta[test_idx];
    for (n in 1:N_test) log_lik_test[n] = normal_lpdf(y_test[n] | eta_test[n], sigma);
  }
}
