pacman::p_load(e1071, MASS, metR, patchwork, sn, tidyverse)

# Set simulation parameters
set.seed(2025)
n <- 100000
mu <- c(1, 2)
sigma <- c(2, 1)
Omega <- diag(sigma^2)
r <- 3
alpha_list <- list(
  "alpha = (0, 0)"   = c(0, 0),
  "alpha = (3, 3)"   = c(3, 3),
  "alpha = (15, 15)" = c(15, 15)
)

# Generate data for each alpha setting
sim_data <- map2_dfr(alpha_list, names(alpha_list), function(alpha, label) {
  X <- rmsn(n = n, xi = mu, Omega = Omega, alpha = alpha)
  Y <- rchisq(n, df = r)
  Z <- X / sqrt(Y / r)
  df <- as.data.frame(Z)
  colnames(df) <- c("Z1", "Z2")
  df$Skewness <- label
  return(df)
})

# Label expressions for Z1, Z2 and Skewness values
sim_long <- sim_data %>%
  pivot_longer(cols = c("Z1", "Z2"),
               names_to = "Variable", values_to = "Value") %>%
  mutate(
    Variable = recode(Variable,
                      "Z1" = "T[1]",
                      "Z2" = "T[2]"),
    Skewness = recode(Skewness,
                      "alpha = (0, 0)"   = "alpha^T == '(' * 0 * ',' * 0 * ')'",
                      "alpha = (3, 3)"   = "alpha^T == '(' * 3 * ',' * 3 * ')'",
                      "alpha = (15, 15)" = "alpha^T == '(' * 15 * ',' * 15 * ')'")
  )

sim_long$Skewness <- factor(sim_long$Skewness,
                            levels = c("alpha^T == '(' * 0 * ',' * 0 * ')'",
                                       "alpha^T == '(' * 3 * ',' * 3 * ')'",
                                       "alpha^T == '(' * 15 * ',' * 15 * ')'"))

# Faceted density plot with parsed math labels
ggplot(sim_long, aes(x = Value)) +
  geom_density(fill = "#648FFF", alpha = 0.7) +
  facet_grid(Skewness ~ Variable, labeller = label_parsed, scales = "free") +
  labs(
    x = NULL, y = "Density") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

ggsave("../figures/Figure4.pdf", width = 7.5, height = 4, units = "in")


sim_data %>%
  group_by(Skewness) %>%
  summarise(skew_Z1 = skewness(Z1),
            skew_Z2 = skewness(Z2),
            q95_Z1 = quantile(Z1, 0.95),
            q95_Z2 = quantile(Z2, 0.95))

# Generate contour data for each alpha setting
set.seed(2025)

lims <- c(-6, 6, -1, 6)
n_grid <- 300
smooth_mult <- 1.5
n_levels <- 7   # same across panels


contour_data <- map2_dfr(alpha_list, names(alpha_list), function(alpha, label) {
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
  df$alpha_label <- label
  df
})

# Label expressions for Z1, Z2 and Skewness values
contour_data <- contour_data %>%
  mutate(
    alpha_label = recode(alpha_label,
                         "alpha = (0, 0)"   = "alpha^T == '(' * 0 * ',' * 0 * ')'",
                         "alpha = (3, 3)"   = "alpha^T == '(' * 3 * ',' * 3 * ')'",
                         "alpha = (15, 15)" = "alpha^T == '(' * 15 * ',' * 15 * ')'")
  )

contour_data$alpha_label <- factor(contour_data$alpha_label,
                                   levels = c("alpha^T == '(' * 0 * ',' * 0 * ')'",
                                              "alpha^T == '(' * 3 * ',' * 3 * ')'",
                                              "alpha^T == '(' * 15 * ',' * 15 * ')'"))
# Plot
ggplot(contour_data, aes(x = x, y = y, z = z)) +
  geom_contour(color = "#648FFF", linewidth = 1, bins=n_levels) +
  facet_wrap(~alpha_label, labeller = label_parsed) +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(
    x = expression(T[1]), y = expression(T[2]))

ggsave("../figures/Figure5.pdf", width = 10, height = 6, units = "in")

