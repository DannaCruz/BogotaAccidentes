# --- RENeGe–Skew-Normal prior simulation -------------------------------
# Simula desde θ = C ρ con ρ ~ SN_p(-b η, σ^2 (M_e - γ A_e)^{-1}; η)
# Implementación via auxiliar U = |Z|, Z ~ N(0,1).
# Aquí usamos K_e = M_e + γ A_e en la covarianza de θ (según tu derivación).

rRENeGeSN <- function(C, eta, Me, Ae, gamma, sigma2, nsim = 1, seed = NULL) {
  # C:     (n x p) matrix
  # eta:   (p) skewness direction on edges
  # Me,Ae: (p x p) matrices
  # gamma: scalar multiplying Ae in the edge precision (your η in the precision term)
  # sigma2: scalar dispersion for θ
  # nsim: number of draws
  stopifnot(is.matrix(C), is.numeric(eta), is.matrix(Me), is.matrix(Ae))
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(C); p <- ncol(C)
  if (length(eta) != p) stop("Length of 'eta' must match ncol(C).")
  
  b <- sqrt(2/pi)
  
  # K_e = M_e + gamma A_e   (edge precision combination)
  # We need Sigma_theta = sigma2 * C %*% solve(K_e, t(C))
  # Use linear solve instead of explicit inverse
  Ke <- Me + gamma * Ae
  
  # Numerical checks (optional but helpful)
  # Symmetrize small numerical noise
  Ke <- 0.5 * (Ke + t(Ke))
  
  # Compute Sigma_theta
  # Solve Ke X = t(C) for X, then C %*% X
  X  <- solve(Ke, t(C))                      # p x n
  Sigma_theta <- sigma2 * (C %*% X)          # n x n
  # Symmetrize
  Sigma_theta <- 0.5 * (Sigma_theta + t(Sigma_theta))
  
  # Cholesky (base::chol expects SPD)
  L <- chol(Sigma_theta)                     # upper-triangular
  
  # Mean and directional term
  mu_theta <- as.vector(-b * (C %*% eta))    # n-vector
  dir_theta <- as.vector(C %*% eta)          # n-vector
  
  # Draw nsim half-normal via |Z|
  U <- abs(rnorm(nsim))                      # length nsim
  
  # Draw eps ~ MN(0, Sigma_theta) in batch: eps = L %*% N(0,I)
  Zmat <- matrix(rnorm(n * nsim), nrow = n, ncol = nsim)
  eps  <- t(L %*% Zmat)                      # nsim x n
  
  # Assemble theta draws: each row is a draw
  Theta <- sweep(eps, 2, mu_theta, FUN = "+")               # add mean
  Theta <- Theta + tcrossprod(U, dir_theta)                 # add C eta * U
  
  list(
    theta = Theta,                 # nsim x n matrix (rows are draws)
    U = U,                         # auxiliary absolute Gaussian
    mu = mu_theta,                 # mean vector of θ|U=0
    dir = dir_theta,               # C %*% eta
    Sigma = Sigma_theta,           # covariance of ε
    chol = L                       # Cholesky factor (upper)
  )
}

# ejemplo sencillo


#sim <- rRENeGeSN(C, eta, Me, Ae, gamma, sigma2, nsim = 5, seed = 2024)



