

source("funciones.R")
library(Metrics)
library(scales)
roundMio<-function(x, d){
  format(round(x, d), nsmall = 2)
}

theme_map <- function(base_size=9, base_family="") {
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.line=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          panel.spacing=unit(0, "lines"),
          plot.background=element_blank(),
          legend.justification = c(0,0),
          legend.position = c(0,0)
    )
}


library(sf)
shp <- st_read("BRMIE250GC_SIR.shp")


shp$NM_MICRO<-as.character(shp$NM_MICRO)


#Encoding(shp$NM_MICRO) <- "UTF-8"

esta <- substr(x = shp$CD_GEOCMI, 1, 2)

shp$esta<-esta

shp_sul1 <- shp[shp$esta == 43, ]
shp_sul2 <- shp[shp$esta ==42,]
shp_sul3 <- shp[shp$esta ==41,]
shp_sul4 <- shp[shp$esta ==35,]

shp_sul <- rbind(shp_sul1,shp_sul2,shp_sul3,shp_sul4)
names(shp_sul)[2]<-"NM_MICRO"



shp_sul <- rbind(shp_sul1, shp_sul2, shp_sul3, shp_sul4)


adj_mat<-nb2mat(poly2nb(shp_sul,queen=F,row.names=shp_sul$NM_MICRO.x), style="B")
A<-adj_mat
g <- network(A, ignore.eval=FALSE, names.eval="a", directed =FALSE, label=TRUE)
C<-as.matrix(g,matrix.type="incidence")
p=length(C[1,])
n=length(C[,1])
A<-(C)%*%t(C)-diag(diag((C)%*%t(C)))
Ar<-t(C)%*%C-2*diag(rep(1,p))
######################################################






# ============================================================
# ESCENARIOS DE SIMULACIÓN PARA RENeGe-sk
# ============================================================

library(MASS)
library(Matrix)
library(rstan)
library(loo)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# Objetos base
# ------------------------------------------------------------

n <- nrow(C)
p <- ncol(C)

E <- rep(1, n)

p_x <- 0
X_mat <- matrix(0, nrow = n, ncol = 0)

b_skew <- sqrt(2 / pi)

# ------------------------------------------------------------
# Función para simular epsilon ~ N(0, (tau Q)^-1)
# ------------------------------------------------------------

sim_eps_precision <- function(Q, tau = 1) {
  Sigma <- solve(tau * Q)
  as.numeric(MASS::mvrnorm(n = 1, mu = rep(0, nrow(Q)), Sigma = Sigma))
}

# ------------------------------------------------------------
# Función para construir eta verdadero según escenario
# ------------------------------------------------------------

make_eta_scenario <- function(scenario, p) {
  
  eta <- rep(0, p)
  
  if (scenario == "S0_symmetric") {
    
    eta <- rep(0, p)
    
  } else if (scenario == "S1_dense_weak") {
    
    eta <- rep(0.25, p)
    
  } else if (scenario == "S2_dense_moderate") {
    
    eta <- rep(0.75, p)
    
  } else if (scenario == "S3_sparse_local") {
    
    idx <- round(seq(0.35 * p, 0.50 * p))
    eta[idx] <- 1.00
    
  } else if (scenario == "S4_sign_changing") {
    
    idx_pos <- round(seq(0.20 * p, 0.35 * p))
    idx_neg <- round(seq(0.60 * p, 0.75 * p))
    
    eta[idx_pos] <-  1.00
    eta[idx_neg] <- -1.00
  }
  
  return(eta)
}

# ------------------------------------------------------------
# Función para simular un conjunto de datos desde RENeGe-sk
# ------------------------------------------------------------

simulate_renege_sk <- function(
    scenario,
    C,
    Q_edge,
    alpha_true = log(5),
    tau_true = 1,
    E = rep(1, nrow(C)),
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(C)
  p <- ncol(C)
  
  eta_true <- make_eta_scenario(scenario, p)
  
  # U half-normal
  U_true <- abs(rnorm(1, mean = 0, sd = 1))
  
  # epsilon edge
  eps_true <- sim_eps_precision(Q = Q_edge, tau = tau_true)
  
  # rho skew-normal centrado
  rho_true <- -b_skew * eta_true + eta_true * U_true + eps_true
  
  # theta nodal
  theta_true <- as.numeric(C %*% rho_true)
  
  # centrar theta para evitar explosiones numéricas en Poisson
  theta_true <- theta_true - mean(theta_true)
  
  linpred_true <- alpha_true + log(E) + theta_true
  
  # protección contra medias demasiado grandes
  linpred_true <- pmin(linpred_true, log(200))
  linpred_true <- pmax(linpred_true, log(0.01))
  
  mu_true <- exp(linpred_true)
  
  Y <- rpois(n, lambda = mu_true)
  
  list(
    Y = as.integer(Y),
    E = E,
    eta_true = eta_true,
    U_true = U_true,
    eps_true = eps_true,
    rho_true = rho_true,
    theta_true = theta_true,
    mu_true = mu_true,
    alpha_true = alpha_true,
    tau_true = tau_true,
    scenario = scenario
  )
}



# ============================================================
# FUNCIÓN PARA AJUSTAR LOS TRES MODELOS
# ============================================================

fit_models_one_dataset <- function(
    sim_data,
    C,
    Q_edge_spd,
    Q_node_spd,
    X_mat,
    sm_RENeGe_SK,
    sm_RENeGe,
    sm_CAR,
    iter = 2000,
    warmup = 500,
    chains = 1,
    seed = 123
) {
  
  n <- nrow(C)
  p <- ncol(C)
  p_x <- ncol(X_mat)
  
  data_RENeGeNorm <- list(
    n = n,
    p = p,
    C = as.matrix(C),
    Q = as.matrix(Q_edge_spd),
    Y = sim_data$Y,
    log_E = log(sim_data$E),
    p_x = p_x,
    X = X_mat,
    a_tau = 1.0,
    b_tau = 1.0
  )
  
  data_CARNorm <- list(
    n = n,
    Q = as.matrix(Q_node_spd),
    Y = sim_data$Y,
    log_E = log(sim_data$E),
    p_x = p_x,
    X = X_mat,
    a_tau = 1.0,
    b_tau = 1.0
  )
  
  # -------------------------------
  # RENeGe-sk
  # -------------------------------
  
  t0 <- proc.time()
  
  fit_sk <- sampling(
    sm_RENeGe_SK,
    data = data_RENeGeNorm,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 0
  )
  
  time_sk <- as.numeric((proc.time() - t0)["elapsed"])
  
  # -------------------------------
  # RENeGe normal
  # -------------------------------
  
  t0 <- proc.time()
  
  fit_re <- sampling(
    sm_RENeGe,
    data = data_RENeGeNorm,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 0
  )
  
  time_re <- as.numeric((proc.time() - t0)["elapsed"])
  
  # -------------------------------
  # CAR normal
  # -------------------------------
  
  t0 <- proc.time()
  
  fit_car <- sampling(
    sm_CAR,
    data = data_CARNorm,
    iter = iter,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    refresh = 0
  )
  
  time_car <- as.numeric((proc.time() - t0)["elapsed"])
  
  list(
    fit_sk = fit_sk,
    fit_re = fit_re,
    fit_car = fit_car,
    time_sk = time_sk,
    time_re = time_re,
    time_car = time_car
  )
}


# ============================================================
# MÉTRICAS DE SIMULACIÓN PARA RENeGe-sk
# ============================================================

coverage_interval <- function(lower, upper, truth) {
  mean(lower <= truth & truth <= upper)
}

extract_metrics_sk <- function(fit_sk, sim_data, time_sk, iter, warmup, chains) {
  
  post <- rstan::extract(fit_sk)
  
  # ----------------------------------------------------------
  # eta
  # ----------------------------------------------------------
  
  eta_post <- post$eta
  eta_mean <- apply(eta_post, 2, mean)
  eta_q025 <- apply(eta_post, 2, quantile, probs = 0.025)
  eta_q975 <- apply(eta_post, 2, quantile, probs = 0.975)
  
  eta_true <- sim_data$eta_true
  
  bias_eta <- mean(eta_mean - eta_true)
  mse_eta  <- mean((eta_mean - eta_true)^2)
  cov_eta  <- coverage_interval(eta_q025, eta_q975, eta_true)
  
  # ----------------------------------------------------------
  # theta
  # ----------------------------------------------------------
  
  theta_post <- post$theta
  theta_mean <- apply(theta_post, 2, mean)
  theta_q025 <- apply(theta_post, 2, quantile, probs = 0.025)
  theta_q975 <- apply(theta_post, 2, quantile, probs = 0.975)
  
  theta_true <- sim_data$theta_true
  
  bias_theta <- mean(theta_mean - theta_true)
  mse_theta  <- mean((theta_mean - theta_true)^2)
  cov_theta  <- coverage_interval(theta_q025, theta_q975, theta_true)
  
  # ----------------------------------------------------------
  # sigma_eta
  # ----------------------------------------------------------
  
  sigma_eta_post <- post$sigma_eta
  sigma_eta_mean <- mean(sigma_eta_post)
  
  # Como sigma_eta no tiene un "true" directo en esta simulación,
  # reportamos su media posterior como indicador de intensidad global.
  
  # ----------------------------------------------------------
  # tiempo por iteración
  # ----------------------------------------------------------
  
  total_iter <- iter * chains
  time_per_iter <- time_sk / total_iter
  
  tibble(
    scenario = sim_data$scenario,
    bias_eta = bias_eta,
    mse_eta = mse_eta,
    coverage_eta_95 = cov_eta,
    bias_theta = bias_theta,
    mse_theta = mse_theta,
    coverage_theta_95 = cov_theta,
    sigma_eta_mean = sigma_eta_mean,
    time_total_sec = time_sk,
    time_per_iter_sec = time_per_iter
  )
}


# ============================================================
# MÉTRICAS PARA MODELOS SIN ETA
# ============================================================

extract_metrics_theta_only <- function(fit, sim_data, model_name, time_model, iter, warmup, chains) {
  
  post <- rstan::extract(fit)
  
  theta_post <- post$theta
  theta_mean <- apply(theta_post, 2, mean)
  theta_q025 <- apply(theta_post, 2, quantile, probs = 0.025)
  theta_q975 <- apply(theta_post, 2, quantile, probs = 0.975)
  
  theta_true <- sim_data$theta_true
  
  bias_theta <- mean(theta_mean - theta_true)
  mse_theta  <- mean((theta_mean - theta_true)^2)
  cov_theta  <- coverage_interval(theta_q025, theta_q975, theta_true)
  
  total_iter <- iter * chains
  time_per_iter <- time_model / total_iter
  
  tibble(
    scenario = sim_data$scenario,
    model = model_name,
    bias_theta = bias_theta,
    mse_theta = mse_theta,
    coverage_theta_95 = cov_theta,
    time_total_sec = time_model,
    time_per_iter_sec = time_per_iter
  )
}










# ============================================================
# LOOP MONTE CARLO POR ESCENARIOS
# ============================================================

scenarios <- c(
  "S0_symmetric",
  "S1_dense_weak",
  "S2_dense_moderate",
  "S3_sparse_local",
  "S4_sign_changing"
)

R_sim <- 5   # primero prueba con 5; luego subir a 50 o 100

iter <- 2000
warmup <- 500
chains <- 1

results_sk <- list()
results_theta <- list()

counter <- 1

for (sc in scenarios) {
  
  cat("\n=====================================\n")
  cat("Scenario:", sc, "\n")
  cat("=====================================\n")
  
  for (r in 1:R_sim) {
    
    cat("Replicate:", r, "of", R_sim, "\n")
    
    seed_r <- 1000 + 100 * which(scenarios == sc) + r
    
    # -------------------------------
    # Simular datos
    # -------------------------------
    
    sim_data <- simulate_renege_sk(
      scenario = sc,
      C = C,
      Q_edge = Q_edge_spd,
      alpha_true = log(5),
      tau_true = 1,
      E = E,
      seed = seed_r
    )
    
    # -------------------------------
    # Ajustar modelos
    # -------------------------------
    
    fits <- fit_models_one_dataset(
      sim_data = sim_data,
      C = C,
      Q_edge_spd = Q_edge_spd,
      Q_node_spd = Q_node_spd,
      X_mat = X_mat,
      sm_RENeGe_SK = sm_RENeGe_SK,
      sm_RENeGe = sm_RENeGe,
      sm_CAR = sm_CAR,
      iter = iter,
      warmup = warmup,
      chains = chains,
      seed = seed_r
    )
    
    # -------------------------------
    # Métricas RENeGe-sk
    # -------------------------------
    
    met_sk <- extract_metrics_sk(
      fit_sk = fits$fit_sk,
      sim_data = sim_data,
      time_sk = fits$time_sk,
      iter = iter,
      warmup = warmup,
      chains = chains
    ) %>%
      mutate(rep = r)
    
    results_sk[[counter]] <- met_sk
    
    # -------------------------------
    # Métricas theta para los 3 modelos
    # -------------------------------
    
    met_theta_sk <- extract_metrics_theta_only(
      fit = fits$fit_sk,
      sim_data = sim_data,
      model_name = "RENeGe-sk",
      time_model = fits$time_sk,
      iter = iter,
      warmup = warmup,
      chains = chains
    ) %>%
      mutate(rep = r)
    
    met_theta_re <- extract_metrics_theta_only(
      fit = fits$fit_re,
      sim_data = sim_data,
      model_name = "RENeGe-normal",
      time_model = fits$time_re,
      iter = iter,
      warmup = warmup,
      chains = chains
    ) %>%
      mutate(rep = r)
    
    met_theta_car <- extract_metrics_theta_only(
      fit = fits$fit_car,
      sim_data = sim_data,
      model_name = "CAR",
      time_model = fits$time_car,
      iter = iter,
      warmup = warmup,
      chains = chains
    ) %>%
      mutate(rep = r)
    
    results_theta[[counter]] <- bind_rows(
      met_theta_sk,
      met_theta_re,
      met_theta_car
    )
    
    counter <- counter + 1
  }
}

results_sk_df <- bind_rows(results_sk)
results_theta_df <- bind_rows(results_theta)

save(
  results_sk_df,
  results_theta_df,
  file = "simulation_results_RENeGe_sk.RData"
)



