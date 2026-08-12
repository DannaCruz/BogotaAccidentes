
data {
  int<lower=1> n;
  matrix[n, n] Q;                 // SPD
  int<lower=0> Y[n];
  vector[n] log_E;
  int<lower=0> p_x;
  matrix[n, p_x] X;
  real<lower=0> a_tau, b_tau;
}
parameters {
  real alpha;
  vector[p_x] beta;
  vector[n] eps_std;
  real<lower=0> tau;
}
transformed parameters {
  matrix[n, n] L = cholesky_decompose(tau * Q);      // L L' = tau Q
  vector[n] theta = mdivide_left_tri_low(L, eps_std); // ~ N(0,(tau Q)^-1)
  vector[n] eta   = alpha + X * beta + log_E + theta;
}
model {
  alpha   ~ normal(0, 1000);
  beta    ~ normal(0, 1000);
  tau     ~ gamma(a_tau, b_tau);
  eps_std ~ normal(0, 1);

  Y ~ poisson_log(eta);
}
generated quantities {
  vector[n] mu_hat;
  int Y_rep[n];
  vector[n] log_lik;
  real rmse = 0;
  real dev  = 0;
  for (i in 1:n) {
    mu_hat[i]  = exp(eta[i]);
    Y_rep[i]   = poisson_log_rng(eta[i]);
    log_lik[i] = poisson_log_lpmf(Y[i] | eta[i]);
    rmse += square(Y[i] - mu_hat[i]);
    dev  += -2 * log_lik[i];
  }
  rmse = sqrt(rmse / n);
}

