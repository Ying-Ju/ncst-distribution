# ------------------------------------------------------------
# Script: 01_ncst_core_functions.R
# Purpose: Core functions for NCST likelihood evaluation and fitting
# Paper: Flexible Modeling of Multivariate Skewed and Heavy-Tailed Data
# ------------------------------------------------------------

loglik_ncst_point <- function(z, mu, Omega, alpha, r, M = 1000) {
  d <- length(z)
  z <- as.vector(z)
  
  u_samples <- rchisq(M, df = r)
  scale_factors <- sqrt(u_samples / r)
  z_mat <- matrix(rep(z, each = M), ncol = d, byrow = FALSE)
  scaled_z <- z_mat * scale_factors
  
  log_sn_densities <- sn::dmsn(scaled_z, xi = mu, Omega = Omega, alpha = alpha, log = TRUE)
  log_jacobians <- d * log(scale_factors)  # equivalent to d/2 * log(u/r)
  log_terms <- log_sn_densities + log_jacobians
  
  max_log <- max(log_terms)
  log_density <- max_log + log(mean(exp(log_terms - max_log)))
  return(log_density)
}


loglik_ncst_full <- function(data, mu, Omega, alpha, r, M = 1000) {
  sum(apply(data, 1, function(x) loglik_ncst_point(x, mu = mu, Omega = Omega, alpha = alpha, r = r, M = M)))
}

ncst.mple <- function(y, start = NULL, M = 1000, opt.method = c("nlminb", "BFGS", "Nelder-Mead"),
                      control = list()) {
  if (missing(y)) stop("Argument 'y' is required.")
  y <- data.matrix(y)
  n <- nrow(y)
  k <- ncol(y)
  
  # Initialize parameters if no custom start
  if (is.null(start)) {
    mu0 <- colMeans(y)
    Omega0 <- cov(y)
    L0 <- chol(Omega0)[lower.tri(Omega0, diag = TRUE)]
    alpha0 <- rep(0, k)
    r0 <- 5
    start <- c(mu0, L0, alpha0, r0)
  }
  
  # Negative log-likelihood wrapper
  neg_loglik_ncst_optim <- function(data, par) {
    # Unpack parameters
    mu <- par[1:k]
    L_vals <- par[(k + 1):(k + k * (k + 1) / 2)]
    alpha <- par[(k + k * (k + 1) / 2 + 1):(2 * k + k * (k + 1) / 2)]
    r <- par[length(par)]
    
    # Reconstruct Omega
    L <- matrix(0, k, k)
    L[lower.tri(L, diag = TRUE)] <- L_vals
    Omega <- L %*% t(L)
    
    # Check validity
    if (!is.finite(r) || r <= 0 || any(!is.finite(mu)) || any(!is.finite(alpha)) ||
        any(!is.finite(Omega)) || any(diag(Omega) <= 0)) {
      return(1e10)
    }
    
    # Compute log-likelihood
    loglik <- tryCatch({
      val <- loglik_ncst_full(data = data, mu = mu, Omega = Omega, alpha = alpha, r = r, M = M)
      if (!is.finite(val)) stop("Non-finite likelihood")
      val
    }, error = function(e) {
      return(-1e10)
    })
    
    return(-loglik)
  }
  
  # Choose optimization method
  opt.method <- match.arg(opt.method)
  if (opt.method == "nlminb") {
    opt <- nlminb(start, objective = neg_loglik_ncst_optim, data =y, control = control)
    par_hat <- opt$par
    logL <- -opt$objective
  } else {
    opt <- optim(start, fn = neg_loglik_ncst_optim, data = y, control = control, method = opt.method)
    par_hat <- opt$par
    logL <- -opt$value
  }
  
  # Unpack estimated parameters
  mu_hat <- par_hat[1:k]
  L_vals <- par_hat[(k + 1):(k + k * (k + 1) / 2)]
  alpha_hat <- par_hat[(k + k * (k + 1) / 2 + 1):(2 * k + k * (k + 1) / 2)]
  r_hat <- par_hat[length(par_hat)]
  
  L <- matrix(0, k, k)
  L[lower.tri(L, diag = TRUE)] <- L_vals
  Omega_hat <- L %*% t(L)
  
  dp <- list(mu = mu_hat, Omega = Omega_hat, alpha = alpha_hat, r = r_hat)
  
  cat("MLE complete.\n")
  cat("Estimated log-likelihood:", logL, "\n")
  
  list(call = match.call(), dp = dp, logL = logL, opt = opt)
}