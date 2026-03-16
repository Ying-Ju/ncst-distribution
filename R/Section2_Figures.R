pacman::p_load(Cairo, gridExtra, MASS, metR, processx, sn, tidyverse)

## Figure 1: Marginal Density
# Set parameters
set.seed(123)
n <- 100000
mu <- c(1, 2)
sigma <- c(2, 1)
alpha <- c(3, 3)
df <- 3

# Simulate skew-normal samples
x <- rmsn(n = n, xi = mu, Omega = diag(sigma^2), alpha = alpha)

# Simulate independent chi-squared samples
y <- rchisq(n, df = df)

# Construct skew-t samples: T = X / sqrt(Y / df)
sample_T <- x / sqrt(y / df)

# Convert to data frame
df_T <- as.data.frame(sample_T)
colnames(df_T) <- c("Z1", "Z2")

dens <- MASS::kde2d(df_T[,1], df_T[,2], n = 100, lims = c(-2, 7, 0, 6))

# Marginal density plots
p1 <- ggplot(df_T, aes(x = Z1)) +
  geom_density(fill = "#648FFF", alpha = 0.7) +
  theme_minimal() +
  ylim(c(0, 0.4)) +
  labs(title = expression(paste("Marginal Density of ", T[1])),
       x = expression(T[1]))

# Marginal density plots
p2 <- ggplot(df_T, aes(x = Z2)) +
  geom_density(fill = "#648FFF", alpha = 0.7) +
  theme_minimal() +
  ylim(c(0, 0.4)) +
  labs(title = expression(paste("Marginal Density of ", T[2])),
       x = expression(T[2]))

m_plot <- grid.arrange(p1, p2, ncol = 2)

ggsave("figures/Figure1.pdf", plot = m_plot, width = 7.5, height = 3.5, units = "in")

## Figure 2: Contour and 3D surface plots
# ---- KDE bandwidths ----
hx <- MASS::bandwidth.nrd(df_T$Z1)
hy <- MASS::bandwidth.nrd(df_T$Z2)

# Mild smoothing to reduce wiggles
h <- 1.5 * c(hx, hy)

# KDE grid
lims <- c(-2, 7, 0, 6)
n_grid <- 300

dens <- MASS::kde2d(
  df_T$Z1, df_T$Z2,
  n = n_grid,
  lims = lims,
  h = h
)


# Convert KDE to long df
df_dens <- with(dens, expand.grid(x = x, y = y))
df_dens$z <- as.vector(dens$z)


# Compute HDR density thresholds for given probability masses
# (approx integration on grid)
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
w  <- dx * dy


# Normalize KDE so total mass over the grid is 1
total_mass <- sum(df_dens$z) * w
df_dens$z_norm <- df_dens$z / total_mass
z_sorted <- sort(df_dens$z_norm, decreasing = TRUE)
mass_cum <- cumsum(z_sorted) * w

hdr_level <- function(p) z_sorted[which(mass_cum >= p)[1]]

probs <- c(0.1, 0.25, 0.50, 0.75, 0.90)
levs  <- sapply(probs, hdr_level)
labs_levs <- paste0(probs * 100, "%")

# Plot: contours correspond to HDR regions
p3 <- ggplot(df_dens, aes(x = x, y = y, z = z_norm)) +
  geom_contour(color = "#648FFF", breaks = levs, linewidth = 0.7) +
  geom_text_contour(breaks = levs, stroke = 0.2, check_overlap = TRUE, size = 3,
                    label.placer = metR::label_placer_fraction(frac = 0.75),
                    aes(label = after_stat(
                      labs_levs[match(round(level, 12), round(levs, 12))]
                    ))) +
  theme_minimal() +
  labs(x = expression(T[1]), y = expression(T[2]))

p3


ggsave("figures/Figure2_Contour.pdf", plot = p3, width = 4, height = 3, units = "in")


library(plotly)

# Convert z to matrix (100 x 100)
z_matrix <- matrix(df_dens$z, nrow = length(dens$x), byrow = FALSE)

# 3D surface plot
plot_ly(x = dens$x, y = dens$y, z = z_matrix) %>%
  add_surface(    colorscale = list(
    list(0, "#DCE5FF"),  # very light blue for low values
    list(1, "#648FFF")   # your vivid blue for high values
  ), showscale = FALSE) %>%  # disables color bar
  layout(
    scene = list(
      xaxis = list(title = "T<sub>1</sub>", showspikes = FALSE, showbackground = FALSE),
      yaxis = list(title = "T<sub>2</sub>", showspikes = FALSE, showbackground = FALSE),
      zaxis = list(title = "f(T<sub>1</sub>, T<sub>2</sub>)", showspikes = FALSE, showbackground = FALSE)
    )
  ) -> p

# Save HTML
htmlwidgets::saveWidget(p, "surface.html")





