data {
  int<lower=1> n;                 // areas
  int<lower=0> k;                 // # bases para omega (0 si no usas)
  int<lower=1> p;                 // # covariables en X (p=1 si solo IDHM)
  matrix[n, k] B;                 // base en nodos para omega (n x k; si k=0, matriz con 0 cols)
  matrix[n, n] Q;                 // precision CAR (SPD)
  matrix[n, p] X;                 // diseño (IDHM,...)
  int<lower=0> Y[n];              // conteos
  vector[n] log_E;                // offset log(E)
  real<lower=0> a_tau;            // hiperparas de tau
  real<lower=0> b_tau;
}

parameters {
  real alpha;                     // intercepto
  vector[p] beta_x;               // betas de X (no informativas)
  real<lower=0> Z;                // HALF-normal para inducir asimetria
  // --- Identificabilidad opcional para omega ---
  vector[k] beta_raw;             // direccion cruda
  real<lower=0> kappa;            // magnitud de omega
  // --------------------------------------------
  vector[n] eps;                  // efecto CAR
  real<lower=0> tau;              // escala CAR
}

transformed parameters {
  vector[n] u;        // direccion unitaria en nodos
  vector[n] omega;    // sesgo proyectado en nodos
  vector[n] theta;    // efecto espacial total
  vector[n] eta;      // predictor lineal

  if (k > 0) {
    vector[k] beta_unit = beta_raw / fmax(1e-12, sqrt(dot_self(beta_raw)));
    u = B * beta_unit;                       // n x 1
    // normaliza tambien en nodos por si B no es ortonormal:
    u = u / fmax(1e-12, sqrt(dot_self(u)));  // ||u||=1 en R^n
    omega = kappa * u;                       // separa magnitud/direccion
  } else {
    u = rep_vector(0, n);
    omega = rep_vector(0, n);
  }

  theta = omega * Z + eps;
  eta   = alpha + X * beta_x + log_E + theta;
}

model {
  // Priors debiles / no informativas
  alpha   ~ normal(0, 1000);
  beta_x  ~ normal(0, 1000);
  Z       ~ normal(0, 1);           // truncado a [0,inf) por la restriccion
  beta_raw ~ normal(0, 1);          // direccion
  kappa   ~ normal(0, 5);           // magnitud (half por <0)
  tau     ~ gamma(a_tau, b_tau);

  // CAR con precision tau * Q (SPD)
  eps ~ multi_normal_prec(rep_vector(0, n), tau * Q);

  // Poisson-log con offset
  Y ~ poisson_log(eta);
}

generated quantities {
  vector[n] mu_hat;               // E[Y|par] = exp(eta)
  int Y_rep[n];                   // replicas
  vector[n] log_lik;              // por observacion
  real rmse;                      // RMSE por draw
  real dev;                       // deviance por draw (para DIC)

  rmse = 0;
  dev  = 0;
  for (i in 1:n) {
    mu_hat[i] = exp(eta[i]);
    Y_rep[i]  = poisson_log_rng(eta[i]);
    log_lik[i] = poisson_log_lpmf(Y[i] | eta[i]);
    rmse += square(Y[i] - mu_hat[i]);
    dev  += -2 * log_lik[i];
  }
  rmse = sqrt(rmse / n);
}
