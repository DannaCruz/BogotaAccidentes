data {
  int<lower=1> N;
  int<lower=1> K;

  array[N] int<lower=0> Y;

  vector[N] log_offset;
  matrix[N, K] X;

  matrix[N, N] W;
  vector<lower=0>[N] degree;
}

parameters {
  real alpha;
  vector[K] beta;

  vector[N] u_raw;
  real<lower=0> sigma_u;
  real<lower=0, upper=0.99> rho;

  real<lower=0> phi;
}

transformed parameters {
  matrix[N, N] Q;
  matrix[N, N] L_Q;
  vector[N] u;
  vector[N] eta;

  Q = diag_matrix(degree) - rho * W;

  for (n in 1:N) {
    Q[n, n] = Q[n, n] + 1e-6;
  }

  L_Q = cholesky_decompose(Q);

  u = sigma_u * mdivide_left_tri_low(L_Q, u_raw);

  u = u - mean(u);

  eta = alpha + X * beta + u + log_offset;
}

model {
  alpha ~ normal(0, 5);
  beta ~ normal(0, 2);

  u_raw ~ normal(0, 1);
  sigma_u ~ exponential(1);

  rho ~ beta(2, 2);
  phi ~ exponential(1);

  Y ~ neg_binomial_2_log(eta, phi);
}

generated quantities {
  vector[N] lambda;
  vector[N] log_lik;
  array[N] int Y_rep;

  for (n in 1:N) {
    lambda[n] = exp(eta[n]);
    log_lik[n] = neg_binomial_2_log_lpmf(Y[n] | eta[n], phi);
    Y_rep[n] = neg_binomial_2_log_rng(eta[n], phi);
  }
}