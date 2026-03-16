pacman::p_load(e1071, MASS, metR, patchwork, sn, tidyverse)
source("00_ncst_core_functions.R")

df_data <- read_csv("../data/wdbc.data", col_names = F)

glimpse(df_data)

df_data <- df_data %>% 
  na.omit()


df_data <- df_data %>% 
  select(X16, X18, X19)

k <- ncol(df_data)
n <- nrow(df_data)


# Convert to matrix
Z_mat <- as.matrix(df_data)
n <- nrow(df_data)

# 1. Fit Multivariate Normal (sample estimates)
mu_norm <- colMeans(Z_mat)
Sigma_norm <- cov(Z_mat)
loglik_norm <- sum(mvtnorm::dmvnorm(Z_mat, mean = mu_norm, sigma = Sigma_norm, log = TRUE))
aic_norm <- -2 * loglik_norm + 2 * (k+k*(k+1)/2)
bic_norm <- -2 * loglik_norm + log(n) * (k+k*(k+1)/2)

aic_norm
bic_norm


# 2. Fit Multivariate Skew-Normal
fit_msn <- msn.mle(y = Z_mat)
loglik_msn <- fit_msn$logL
aic_msn <- -2 * loglik_msn + 2 * (2*k+k*(k+1)/2)
bic_msn <- -2 * loglik_msn + log(n) * (2*k+k*(k+1)/2)

aic_msn
bic_msn

# 3. Fit Multivariate Skew-t 
fit_mst <- try(mst.mple(y = Z_mat), silent = TRUE)
if (inherits(fit_mst, "try-error")) {
  loglik_mst <- NA
  aic_mst <- NA
  bic_mst <- NA
} else {
  loglik_mst <- fit_mst$logL
  aic_mst <- -2 * loglik_mst + 2*(2*k + k*(k+1)/2 + 1)
  bic_mst <- -2 * loglik_mst + log(n) * (2*k + k*(k+1)/2 + 1)
}

aic_mst
bic_mst

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

# Rename columns for interpretability
colnames(df_data) <- c("concavity_se", "symmetry_se", "fractal_dimension_se")

### Marginal Plots

set.seed(6292025)
n <- 1000000
mu <- c(7.64628148, 0.01893632, 0.02088567)
Omega <- matrix(c(677.51184277, -0.04478431, 0.01048751,
                  -0.04478431, 0.00003476395, 0.00003767294,
                  0.01048751, 0.00003767294, 0.0001028847), nrow = 3, byrow = TRUE)
alpha <- c(222392.924, 4046.867, 7799.454)
r <- 1.994889

X <- rmsn(n = n, xi = mu, Omega = Omega, alpha = alpha)
Y <- rchisq(n, df = r)
Z <- X / sqrt(Y / r)
df_sim <- as.data.frame(Z)
colnames(df_sim)  <- c("concavity_se", "symmetry_se", "fractal_dimension_se")


# 1. Simulated NCST sample
Z_sim <- df_sim  # simulated NCST sample with 3 columns

# 2. Observed data
Z_obs <- df_data  # your observed data

Z_sim_clipped <- Z_sim %>%
  filter(
    concavity_se > 0 & concavity_se < 600,
    symmetry_se > 0 & symmetry_se < 0.15,
    fractal_dimension_se >=0 & fractal_dimension_se < 0.4
  )

# Define consistent theme
my_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Plot for concavity_se
p1 <- ggplot(Z_obs, aes(x = concavity_se)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, 
                 color = "black", fill = "grey70", alpha = 0.6) +
  geom_density(data = Z_sim_clipped, aes(x = concavity_se), 
               color = "#648FFF", linewidth = 0.8) +
  labs(x = "concavity_se", y = "Density") +
  my_theme

# Plot for symmetry_se
p2 <- ggplot(Z_obs, aes(x = symmetry_se)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, 
                 color = "black", fill = "grey70", alpha = 0.6) +
  geom_density(data = Z_sim_clipped, aes(x = symmetry_se), 
               color = "#648FFF", linewidth = 0.8) +
  labs(x = "symmetry_se", y = NULL) +
  my_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Plot for fractal_dimension_se
p3 <- ggplot(Z_obs, aes(x = fractal_dimension_se)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, 
                 color = "black", fill = "grey70", alpha = 0.6) +
  geom_density(data = Z_sim_clipped, aes(x = fractal_dimension_se), 
               color = "#648FFF", linewidth = 0.8) +
  labs(x = "fractal_dimension_se", y = NULL) +
  my_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Combine all plots into one wide figure
(p1 | p2 | p3) 

ggsave("../figures/Figure8.pdf", width = 7.5, height = 3.5, units = "in")


### Contour Plots

df_sim_clipped <- df_sim %>%
  filter(
    concavity_se > 0 & concavity_se < 600,
    symmetry_se > 0 & symmetry_se < 0.15,
    fractal_dimension_se >=0 & fractal_dimension_se < 0.4
  )

# Kernel density estimates on simulated data
dens_xy <- kde2d(df_sim_clipped$concavity_se, df_sim_clipped$symmetry_se, n = 200, h = c(20, 0.008))
dens_xz <- kde2d(df_sim_clipped$concavity_se, df_sim_clipped$fractal_dimension_se, n = 200, h = c(20, 0.007))
dens_yz <- kde2d(df_sim_clipped$symmetry_se, df_sim_clipped$fractal_dimension_se, n = 200, h = c(0.007, 0.007))

# Convert to data frames
df_dens_xy <- expand.grid(x = dens_xy$x, y = dens_xy$y)
df_dens_xy$z <- as.vector(dens_xy$z)

df_dens_xz <- expand.grid(x = dens_xz$x, y = dens_xz$y)
df_dens_xz$z <- as.vector(dens_xz$z)

df_dens_yz <- expand.grid(x = dens_yz$x, y = dens_yz$y)
df_dens_yz$z <- as.vector(dens_yz$z)

# Quantile limits
xlim_range <- quantile(df_data$concavity_se, c(0, 0.95))
ylim_range <- quantile(df_data$symmetry_se, c(0, 0.95))
zlim_range <- quantile(df_data$fractal_dimension_se, c(0, 0.95))

# Plot 1: concavity_se vs symmetry_se
p1 <- ggplot() +
  geom_point(data = df_data, aes(x = concavity_se, y = symmetry_se), alpha = 0.2) +
  geom_contour(data = df_dens_xy, aes(x = x, y = y, z = z), color = "#648FFF", bins = 10, 
               linewidth = 0.6) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
  labs(x = "concavity_se", y = "symmetry_se") +
  theme_minimal()

# Plot 2: concavity_se vs fractal_dimension_se
p2 <- ggplot() +
  geom_point(data = df_data, aes(x = concavity_se, y = fractal_dimension_se), alpha = 0.2) +
  geom_contour(data = df_dens_xz, aes(x = x, y = y, z = z), color = "#648FFF", bins = 10,
               linewidth = 0.6) +
  coord_cartesian(xlim = xlim_range, ylim = zlim_range) +
  labs(x = "concavity_se", y = "fractal_dimension_se") +
  theme_minimal()

# Plot 3: symmetry_se vs fractal_dimension_se
p3 <- ggplot() +
  geom_point(data = df_data, aes(x = symmetry_se, y = fractal_dimension_se), alpha = 0.2) +
  geom_contour(data = df_dens_yz, aes(x = x, y = y, z = z), color = "#648FFF", bins = 10,
               linewidth = 0.6) +
  coord_cartesian(xlim = ylim_range, ylim = zlim_range) +
  labs(x = "symmetry_se", y = "fractal_dimension_se") +
  theme_minimal()

# Combine
(p1 | p2 | p3) 

ggsave("../figures/Figure9.pdf", width = 7.5, height = 3.5, units = "in")

