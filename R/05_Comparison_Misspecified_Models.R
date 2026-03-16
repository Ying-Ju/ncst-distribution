pacman::p_load(e1071, MASS, metR, patchwork, sn, tidyverse)
source("00_ncst_core_functions.R")

set.seed(2025)

# Set parameters
mu <- c(1, 2)
Omega <- matrix(c(4, 0, 0, 1), 2)
alpha <- c(3, 0)
r      <- 1  # degrees of freedom
n      <- 100

# Generate X ~ Skew-Normal
X <- rmsn(n = n, xi = mu, Omega = Omega, alpha = alpha)

# Generate Y ~ Chi-squared
Y <- rchisq(n, df = r)

# Construct Z = X / sqrt(Y / df)
Z <- X / sqrt(Y / r)

# Store in data frame
df_data <- as.data.frame(Z)
names(df_data) <- c("Z1", "Z2")
k <- ncol(df_data)

# Convert to matrix
Z_mat <- as.matrix(df_data)

# 1. Fit Multivariate Normal (sample estimates)
mu_norm <- colMeans(Z_mat)
Sigma_norm <- cov(Z_mat)
loglik_norm <- sum(mvtnorm::dmvnorm(Z_mat, mean = mu_norm, sigma = Sigma_norm, log = TRUE))
aic_norm <- -2 * loglik_norm + 2 * (k+k*(k+1)/2)
bic_norm <- -2 * loglik_norm + log(n) * (k+k*(k+1)/2)

# 2. Fit Multivariate Skew-Normal
fit_msn <- msn.mle(y = Z_mat)
loglik_msn <- fit_msn$logL
aic_msn <- -2 * loglik_msn + 2 * (2*k+k*(k+1)/2)
bic_msn <- -2 * loglik_msn + log(n) * (2*k+k*(k+1)/2)

# 3. Fit Multivariate Skew-t 
fit_mst <- try(mst.mple(y = Z_mat), silent = TRUE)
if (inherits(fit_mst, "try-error")) {
  loglik_mst <- NA
  aic_mst <- NA
} else {
  loglik_mst <- fit_mst$logL
  aic_mst <- -2 * loglik_mst + 2*(2*k + k*(k+1)/2 + 1)
  bic_mst <- -2 * loglik_mst + log(n) * (2*k + k*(k+1)/2 + 1)
}


################ distribution fitting
# Approximate log-density of NCST
p <- ncol(Z_mat)

# Starting values
start_mu <- as.vector(fit_mst$dp$beta) #c(2,1) #colMeans(Z_mat)
start_Omega <- fit_mst$dp$Omega #diag(2) #cov(Z_mat)  
start_L <- chol(start_Omega)[lower.tri(start_Omega, diag = TRUE)]
start_alpha <- as.vector(fit_mst$dp$alpha)
start_r <- max(round(fit_mst$dp$nu), 1)

start_par <- c(start_mu, start_L, start_alpha, start_r)

fit_ncst <- ncst.mple(y = Z_mat, start = start_par, M = 5000, opt.method = "BFGS", 
                      control = list(maxit = 1000))

fit_ncst$opt$convergence

loglik_ncst <- fit_ncst$logL
aic_ncst <- -2 * loglik_ncst + 2*(2*k + k*(k+1)/2 + 1)
bic_ncst <- -2 * loglik_ncst + log(n)*(2*k + k*(k+1)/2 + 1)

# Summary table
model_fit_summary <- data.frame(
  Model = c("Multivariate Normal", "Skew-Normal", "Azzalini Skew-t", "Skew-t"),
  LogLikelihood = c(loglik_norm, loglik_msn, loglik_mst, loglik_ncst),
  AIC = c(aic_norm, aic_msn, aic_mst, aic_ncst),
  BIC = c(bic_norm, bic_msn, bic_mst, bic_ncst)
)

print(model_fit_summary)
