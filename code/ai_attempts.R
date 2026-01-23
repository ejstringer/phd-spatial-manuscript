
# ---- Packages ----
# If needed, uncomment:
# install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---- Parameters ----
sigmas <- c(0.1, 0.25,0.4, 0.5)    # dispersal SD (σ) values
r_max  <- 5                  # max geographic distance
n_pts  <- 400                  # resolution for the curve

# ---- Build data ----
r <- seq(0, r_max, length.out = n_pts)

dat <- expand.grid(r = r, sigma = sigmas) %>%
  mutate(cov = exp(-r / sigma),
         sigma_lab = paste0("\u03C3 = ", sigma))  # σ = ...

# ---- Plot ----
p <- ggplot(dat, aes(x = r, y = cov, color = factor(sigma_lab))) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Set1", name = "Dispersal SD") +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  labs(
    title = "Decline of genetic covariance with distance",
    subtitle = "Heuristic: Cov(r) \u221D exp(-r/\u03C3)",
    x = "Distance r",
    y = "Relative covariance (arbitrary units)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p)

# ---- (Optional) Save ----
# ggsave("wright_ibd_left_panel.png", p, width = 7, height = 5, dpi = 300)
``


# ---- Packages ----
# install.packages(c("ggplot2", "dplyr", "tidyr", "scales", "patchwork")) # if needed
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)

# -----------------------------
# Panel A: Exponential heuristic
# -----------------------------
sigmas <- c(0.5, 1.0, 2.0)  # dispersal SDs (per generation)
r_max  <- 10
n_pts  <- 400
r <- seq(0, r_max, length.out = n_pts)

dat_exp <- expand.grid(r = r, sigma = sigmas) %>%
  mutate(cov = exp(-r / sigma),
         sigma_lab = paste0("\u03C3 = ", sigma)) # σ = ...

p_exp <- ggplot(dat_exp, aes(x = r, y = cov, color = sigma_lab)) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Set1", name = "Dispersal SD") +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  labs(
    title = "Decline of genetic covariance with distance",
    subtitle = "Heuristic shape: Cov(r) \u221D exp(-r/\u03C3) (density does not enter)",
    x = "Distance r",
    y = "Relative covariance (arbitrary units)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# -----------------------------------------
# Panel B: Rousset (1997) linearized model
# F/(1-F) = r / (4*pi*D*sigma^2) + constant
# Only the slope matters; set intercept = 0
# -----------------------------------------
sigmas <- c(0.01, 0.05, 0.1)  # dispersal SDs (per generation)
r_max  <- 10
n_pts  <- 400
r <- seq(0, r_max, length.out = n_pts)

dat_exp <- expand.grid(r = r, sigma = sigmas) %>%
  mutate(cov = exp(-r / sigma),
         sigma_lab = paste0("\u03C3 = ", sigma)) # σ = ...


densities <- c(10, 50)  # individuals per unit area (choose units consistent with r)
# We'll show combinations of (D, sigma) to see both effects

dat_lin <- expand.grid(r = r, sigma = sigmas, D = densities) %>%
  mutate(
    slope = 1 / (4 * pi * D * sigma^2),
    y = slope * r,    # intercept set to 0 for visualization
    D_lab = paste0("D = ", D),
    sigma_lab = paste0("\u03C3 = ", sigma),
    NeighSize = 4 * pi * D * sigma^2 # Wright's neighborhood size
  )

# A helper label that shows both D and neighborhood size (Ne = 4πDσ²)
dat_lin <- dat_lin %>%
  mutate(group_lab = paste0(D_lab, ", N\u2091 = ", round(NeighSize, 2),
                            ", \u03C3 = ", sigma))

p_lin <- ggplot(dat_lin, aes(x = r, y = y, color = interaction(D_lab, sigma_lab, sep = " | "))) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Dark2", name = "Density | Dispersal SD") +
  labs(
    title = "Linearized IBD (Rousset 1997): F/(1-F) vs distance",
    subtitle = "Slope b = 1 / (4 \u03C0 D \u03C3\u00B2)\nLower D or smaller \u03C3 \u2192 steeper slope (stronger IBD)",
    x = "Distance r",
    y = "F(r) / (1 - F(r))  (offset set to 0 for plotting)"
  ) +
  theme_minimal(base_size = 13) +
  facet_wrap(~sigma_lab, scale = 'free')+
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_lin
# -----------------------
# Combine the two panels
# -----------------------
p <- p_exp + p_lin + plot_layout(ncol = 2, widths = c(1, 1.15))
print(p)

# ---- (Optional) Save ----
# ggsave("wright_ibd_with_density.png", p, width = 12, height = 5.5, dpi = 300)



# ---- Packages ----
# install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))  # if needed
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---- Parameters ----
sigmas    <- c(0.01, .1, 1)              # dispersal SDs
densities <- c(100, 1000)*0.2           # individuals per unit area (Dickman et al 1999)
mu        <- 1e-3                          # per-generation mutation rate (choose for your system)
r_max     <- 5
n_pts     <- 400
eps       <- 0.05                          # small-scale cutoff (sampling radius)
r         <- seq(0, r_max, length.out = n_pts)

# ---- Build data ----
dat <- expand.grid(r = r, sigma = sigmas, D = densities) %>%
  mutate(
    alpha   = sigma / sqrt(2 * mu),                       # correlation length
    r_eff   = sqrt(r^2 + eps^2),                          # regularization near 0
    # Raw K0 curve with density scaling (amplitude depends on D)
    cov_raw = (1 / D) * besselK(r_eff / alpha, nu = 0),
    # Normalized version to compare shapes (set value at r=eps to 1)
    cov_norm = cov_raw / ((1 / D) * besselK(eps / alpha, nu = 0)),
    sigma_lab = paste0("\u03C3 = ", sigma),
    D_lab     = paste0("D = ", D)
  )

# ---- Choose which y to plot ----
# 1) To show density's effect on amplitude, use y = cov_raw.
# 2) To compare decay shapes independent of amplitude, use y = cov_norm.
plot_y <- "cov_raw"  # or "cov_norm"
#plot_y <- "cov_norm"

p_bessel <- ggplot(dat, aes(x = r, y = .data[[plot_y]],
                            color = D_lab)) +
  geom_line(linewidth = 1.1, linetype = 2) +
  scale_color_brewer(palette = "Set2", name = "Density") +
 # scale_linetype_discrete(name = "Dispersal SD") +
  labs(
    title = "Isolation by distance in 2D (Bessel K0 form)",
    subtitle = "Cov(r; D, \u03C3, \u03BC) \u221D (1/D)·K0(r/\u03B1),   \u03B1 = \u03C3/\u221A(2\u03BC)",
    x = "Distance r", y = "Relative covariance (arbitrary units)"
  ) +
  facet_grid(~sigma_lab)+
  theme_minimal(base_size = 13) +
 # geom_smooth(data = filter(dat, r < 1), method = 'lm')+
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p_bessel)

# ---- (Optional) Save ----
# ggsave("ibd_leftpanel_bessel.png", p_bessel, width = 7.2, height = 5.2, dpi = 300)


# ---- Packages ----
# install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))  # if needed
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---- Parameters ----
sigmas    <- c(0.1, 0.2)        # dispersal SDs (distance units per generation)
densities <- c(20, 500)     # density D (breeders per unit area)
r_max     <- 10
n_pts     <- 300
r         <- seq(1, r_max, length.out = n_pts)

# ---- Build theoretical lines ----
dat_th <- expand.grid(r = r, sigma = sigmas, D = densities) %>%
  mutate(
    slope = 1 / (4 * pi * D * sigma^2),        # b = 1/(4π D σ^2)
    y = slope * r,                              # intercept set to 0 for visualization
    sigma_lab = paste0("\u03C3 = ", sigma),     # label: σ = ...
    D_lab     = paste0("D = ", D)
  )

p_theory <- ggplot(dat_th, aes(x = r, y = y,
                               color = sigma_lab, linetype = D_lab)) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Set1", name = "Dispersal SD") +
  scale_linetype_discrete(name = "Density") +
  labs(
    title = expression("Linearized IBD: " ~ F[ST]/(1 - F[ST]) ~ " vs distance"),
    subtitle = "Slope b = 1 / (4 \u03C0 D \u03C3\u00B2); lower D or smaller \u03C3 \u2192 steeper slope",
    x = "Distance r",
    y = expression(F[ST]/(1 - F[ST]))
  ) +
  facet_grid(~sigma_lab, scale = 'free')+
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p_theory)

y <- seq(0, 10, 0.1)
y/0.5/y
1/0.5
# # Optional: save
# ggsave("ibd_theory_fst_over_one_minus_fst.png", p_theory, width = 7.2, height = 5.0, dpi = 300)


ggplot()
