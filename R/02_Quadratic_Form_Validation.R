pacman::p_load(e1071, MASS, metR, patchwork, sn, tidyverse)

# Step 1: Set parameters
set.seed(2025)
n <- 100000

# Create a rank-1 projection matrix W
a <- c(1, 1)
W <- outer(a, a) / sum(a^2)  # W is 2x2, rank 1, symmetric, idempotent

# Choose B = I
B <- diag(2)

# Confirm BWB' = W is idempotent
M <- B %*% W %*% t(B)
all.equal(M %*% M, M)  # should be TRUE

# Eigen decomposition to get P
eig <- eigen(M)
P <- eig$vectors
m <- sum(eig$values > 1e-8)  # m = rank of M
P1 <- P[, 1:m]
P2 <- P[, (m+1):ncol(P)]

# For theoretical expressions
mu <- c(1, 2)
alpha <- c(3, 3)
r <- 3

nu <- t(P1) %*% B %*% W %*% mu
c_val <- sqrt(1 + t(alpha) %*% P1 %*% t(P1) %*% alpha)
alpha_star <- sqrt((1 / c_val^2) * t(alpha) %*% P1 %*% t(P1) %*% alpha)
lambda <- t(mu) %*% W %*% mu

# Step 2: Generate X ~ MSN(mu, B, alpha)
X <- rmsn(n = n, xi = mu, Omega = B, alpha = alpha)

# Step 3: Generate U ~ Chi-squared(r)
U <- rchisq(n, df = r)

# Step 4: Construct Y = X / sqrt(U / r)
Y <- X / sqrt(U / r)

# Step 5: Compute Q = Y^T W Y
Q <- rowSums((Y %*% W) * Y)

Q_df <- data.frame(Q = Q)


Q_df$logQ <- log(Q_df$Q)

ggplot(Q_df, aes(x = logQ)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 fill = "#648FFF", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  labs(title = "Empirical Distribution of log(Q)",
       x = expression(log(Q)), y = "Density") +
  theme_minimal()



# --- Parameters ---
df1 <- 1
df2 <- 3

# Step 1: Skewed non-central chi-squared approximation
# Use skew-normal Z, then square to get skewed chi-squared

Z <- rsn(n, xi = sqrt(lambda), omega = 1, alpha = alpha_star)
X1 <- Z^2

# Step 2: Central chi-squared
Y1 <- rchisq(n, df = df2)

# Step 3: Non-central skew F variable
Q_skewF <- (X1/df1) / (Y1 / df2)

Q_df$skewF <- Q_skewF

q_seq <- seq(0.001, 0.999, by = 0.001)  # avoid extremes
df_quantiles <- data.frame(
  Q_quantiles = quantile(Q, probs = q_seq),
  F_quantiles = quantile(Q_skewF, probs = q_seq)
)

colnames(df_quantiles) <- c("Q_quantiles", "F_quantiles")


p1 <- ggplot(df_quantiles, aes(x = Q_quantiles, y = F_quantiles)) +
  geom_point(color = "#2A9D8F", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "#E76F51", linewidth = 0.8) +
  labs(
    x = "Empirical Q Quantiles",
    y = "Simulated Skew-F Quantiles") +
  theme_minimal()

# --- Density Plot (zoomed) ---
Q_df <- data.frame(Q = Q, skewF = Q_skewF)
Q_df_long <- pivot_longer(Q_df, cols = everything(), names_to = "Type", values_to = "value")
Q_df_long$Type <- factor(Q_df_long$Type, levels = c("skewF", "Q"))

upper_limit <- quantile(c(Q, Q_skewF), 0.999)

density_colors <- c("Q" = "#2A9D8F", "skewF" = "#E76F51")

p2 <- ggplot(Q_df_long, aes(x = value, fill = Type)) +
  geom_density(alpha = 0.4, color = "black") +
  scale_fill_manual(
    values = density_colors,
    labels = c("Q" = "Empirical Q", "skewF" = "Empirical Skew-F")
  ) +
  coord_cartesian(xlim = c(0, upper_limit)) +
  labs(x = expression(italic(Q) * " or Skew-"*italic(F)), 
       y = "Density", fill = "Distribution") +
  theme_minimal() +
  theme(
    legend.position= c(0.95, 0.95),           # (x, y) from bottom-left in [0, 1]
    legend.justification = c("right", "top"),  # anchor legend box to this point
    legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
    legend.box.background = element_blank()
  )


# Combine with patchwork
combined_plot <- p1 + p2 + plot_layout(ncol = 2)
combined_plot

ggsave("../figures/Figure3.pdf", width = 7.5, height = 4, units = "in")
