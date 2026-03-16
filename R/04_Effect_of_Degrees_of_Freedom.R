pacman::p_load(e1071, MASS, metR, patchwork, sn, tidyverse)

# Settings
set.seed(2025)
n <- 100000
mu <- c(1, 2)
sigma <- c(2, 1)
Omega <- diag(sigma^2)
alpha <- c(3, 3)
r_list <- c(3, 5, 10, 30)

# Generate simulation data for each r
sim_df_r <- purrr::map_dfr(r_list, function(r) {
  X <- rmsn(n = n, xi = mu, Omega = Omega, alpha = alpha)
  Y <- rchisq(n, df = r)
  Z <- X / sqrt(Y / r)
  df <- as.data.frame(Z)
  colnames(df) <- c("Z1", "Z2")
  df$r_label <- paste0("r == ", r)
  df
})

# Convert to long format for faceting
sim_long_r <- sim_df_r %>%
  pivot_longer(cols = c("Z1", "Z2"),
               names_to = "Variable",
               values_to = "Value") %>%
  mutate(
    Variable = recode(Variable,
                      "Z1" = "T[1]",
                      "Z2" = "T[2]")
  )

sim_long_r$r_label <- factor(sim_long_r$r_label, 
                             levels = c("r == 3", "r == 5", "r == 10", "r == 30"))

# Faceted plot by r and variable
ggplot(sim_long_r, aes(x = Value)) +
  geom_density(fill = "#648FFF", alpha = 0.7) +
  ylim(c(0, 0.7)) +
  facet_grid(r_label ~ Variable, labeller = label_parsed, scales = "free") +
  labs(
    x = NULL, y = "Density") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))


ggsave("../figures/Figure6.pdf", width = 7.5, height = 4.5, units = "in")

sim_df_r$r_label <- factor(sim_df_r$r_label, 
                             levels = c("r == 3", "r == 5", "r == 10", "r == 30"))

sim_df_r %>%
  group_by(r_label) %>%
  summarise(
    skew_Z1 = round(skewness(Z1), 2),
    skew_Z2 = round(skewness(Z2), 2),
    kurt_Z1 = round(kurtosis(Z1), 2),
    kurt_Z2 = round(kurtosis(Z2), 2),
    q95_Z1 = round(quantile(Z1, 0.95), 2),
    q95_Z2 = round(quantile(Z2, 0.95), 2)
  )



# Settings
set.seed(2025)

n_grid <- 300
lims <- c(-6, 6, -1, 6)
smooth_mult <- 1.5
n_levels <- 7   # same across panels

# Collect density estimates for each r
contour_df <- map_dfr(r_list, function(r) {
  X <- rmsn(n = n, xi = mu, Omega = Omega, alpha = alpha)
  Y <- rchisq(n, df = r)
  Z <- X / sqrt(Y / r)
  hx <- MASS::bandwidth.nrd(Z[,1])
  hy <- MASS::bandwidth.nrd(Z[,2])
  dens <- MASS::kde2d(
    Z[,1], Z[,2],
    n = n_grid,
    lims = lims,
    h = smooth_mult * c(hx, hy)
  )
  df <- with(dens, expand.grid(x = x, y = y))
  df$z <- as.vector(dens$z)
  df$r_label <- paste0("r == ", r)
  df
})

contour_df$r_label <- factor(contour_df$r_label, 
                             levels = c("r == 3", "r == 5", "r == 10", "r == 30"))

# Plot 2x2 contour facets
ggplot(contour_df, aes(x = x, y = y, z = z)) +
  geom_contour(color = "#648FFF", linewidth = 0.8, bins=n_levels) +
  facet_wrap(~r_label, labeller = label_parsed, ncol = 2) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = c(-2, 6), ylim = c(0, 5)) +
  labs(x = expression(T[1]), y = expression(T[2])) +
  theme(strip.text = element_text(face = "bold"))


ggsave("../figures/Figure7.pdf", width = 7.5, height = 5, units = "in")

