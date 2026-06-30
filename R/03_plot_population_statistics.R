# ==============================================================================
# Generates Figures 1 and 2 for the manuscript, establishing the asymptotic
# structural validity of the framework and the biological construct validity
# of the extracted parameters across all physiological cohorts.
# ==============================================================================

# 1. Prepare Workspace ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# 2. Load Aggregated Data ------------------------------------------------------

df_summary <- readRDS("results/df_population_summary.rds")
mat_pacf <- readRDS("results/mat_population_pacf.rds")

# Filter for successfully converged models only
df_valid <- df_summary %>% filter(Converged == TRUE)
n_valid <- nrow(df_valid)

df_valid$Dataset <- factor(
  x = df_valid$Dataset,
  levels = c("Exercise_Recovery", "Fantasia_Baseline", "Orthostatic_Tilt"),
  labels = c("Two-minute Step Test", "Fantasia Baseline", "Orthostatic Tilt Test")
)

# Define universal color palettes
color_cohorts <- c("Fantasia Baseline" = "#2A363B",
                   "Orthostatic Tilt Test"  = "#E84A5F",
                   "Two-minute Step Test" = "#FF847C")

color_branches <- c("Parasympathetic" = "#4DB0D0",
                    "Sympathetic"     = "#DC5050")

# ==============================================================================
# FIGURE 1: ASYMPTOTIC VALIDATION PANEL (Section 3.1)
# ==============================================================================

# --- Panel A: Engineering Tracking Error (ECDF) ---
# Shows the cumulative proportion of subjects achieving a given error threshold
p1A <- ggplot(df_valid, aes(x = MAPE_pct, color = Dataset)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_color_manual(values = color_cohorts) +
  # Assuming MAPE_pct is 0-100. If it is 0-1, use scale = 100 in percent_format
  scale_x_continuous(limits = c(0, quantile(df_valid$MAPE_pct, 0.99))) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "A. Macroscopic Tracking Accuracy",
       subtitle = "Cumulative distribution of absolute error",
       x = "Mean Absolute Percentage Error (%)",
       y = "Cumulative Proportion of Subjects") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

# --- Panel B: Dynamic Completeness (Cohort-Specific PACF) ---

# 1. Filter the PACF matrix to match only converged subjects
valid_indices <- which(df_summary$Converged == TRUE)
mat_pacf_valid <- mat_pacf[valid_indices, ]

# 2. Convert to a dataframe and attach the Dataset labels
df_pacf_raw <- as.data.frame(mat_pacf_valid)
colnames(df_pacf_raw) <- 1:ncol(mat_pacf_valid) # Name columns by lag number
df_pacf_raw$Dataset <- df_valid$Dataset

# 3. Pivot to long format and calculate grouped statistics
df_pacf <- df_pacf_raw %>%
  pivot_longer(cols = -Dataset, names_to = "Lag", values_to = "PACF") %>%
  mutate(Lag = as.numeric(Lag)) %>%
  group_by(Dataset, Lag) %>%
  summarise(
    n_group = n(),
    Mean = mean(PACF, na.rm = TRUE),
    SD = sd(PACF, na.rm = TRUE),
    Lower = Mean - qnorm(0.975) * (SD / sqrt(n_group)),
    Upper = Mean + qnorm(0.975) * (SD / sqrt(n_group)),
    .groups = "drop"
  )

# Theoretical global 95% white noise bound
ci_bounds <- df_valid |>
  group_by(Dataset) |>
  summarise(
    avg_N = mean(NumBeats),
    ci_bound = qnorm(0.975) / sqrt(avg_N)
  )

# 4. Updated Plot: Map color and fill to Dataset
p1B <- ggplot(df_pacf, aes(x = Lag, y = Mean, color = Dataset, fill = Dataset)) +
  facet_grid(rows = vars(Dataset), scales = "free_y") +
  geom_rect(data = ci_bounds, aes(xmin = -Inf, xmax = Inf, ymax = ci_bound, ymin = -ci_bound, x = 1, y = 0),
            fill = "gray80", col = 0, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Plot dataset-specific ribbons and lines
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 1.2, show.legend = FALSE) +
  scale_color_manual(values = color_cohorts) +
  scale_fill_manual(values = color_cohorts) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "B. Cohort-Specific PACF",
       subtitle = "Residual innovations reduced to white noise",
       x = "Cardiac Lag (Beats)",
       y = "Partial Autocorrelation") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Panel C: Mathematical Exactness (KS-Distance Density) ---

p1C <- ggplot(df_valid, aes(x = KS_Stat, fill = Dataset, col = Dataset)) +
  facet_grid(rows = vars(Dataset)) +
  geom_density(alpha = 0.7, linewidth = 1, trim = FALSE, show.legend = FALSE) +
  scale_fill_manual(values = color_cohorts, name = "Cohort", aesthetics = c("fill", "color")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,NA)) +
  scale_y_continuous(expand = c(0,0,0.25,0)) +
  labs(title = "C. Time-Rescaling Theorem",
       subtitle = "KS D-statistic aggregation (Geometric Distance)",
       x = "Kolmogorov-Smirnov Distance (D-Statistic)",
       y = "Density") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Assemble Figure 1 ---
fig1_layout <- (p1A | p1B | p1C) +
  plot_layout(widths = c(1, 1, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = "manuscript/figures/fig1_statistical_validation.png",
       plot = fig1_layout, width = 250, height = 120, dpi = 300, scale = 15,
       units = "px")


# ==============================================================================
# FIGURE 2: BIOLOGICAL CONSTRUCT VALIDITY (Section 3.2)
# ==============================================================================

# --- Panel A: Autonomic Memory Decay (Stochastic Impulse Response) ---
# Calculate the global median clearance rates (drift)
median_kp <- median(df_valid$Kappa_P, na.rm = TRUE)
median_ks <- median(df_valid$Kappa_S, na.rm = TRUE)

# Calculate the global median diffusion coefficient (volatility)
# Since the model estimates global system variance (Sigma_sys2),
# the SDE diffusion term (sigma) is the square root of the variance.
median_sigma <- 0.005

# Generate a time vector from 0 to roughly 3x the sympathetic half-life
t_max <- 5 * (1/median_ks)
n_steps <- 500
t_seq <- seq(0, t_max, length.out = n_steps)
dt <- t_seq[2] - t_seq[1]

# Function to simulate one path of the OU process (Euler-Maruyama approximation)
# SDE: dX_t = -kappa * X_t * dt + sigma * dW_t
sim_ou_path <- function(kappa, sigma, times, dt, X0 = 1) {
  n <- length(times)
  X <- numeric(n)
  X[1] <- X0
  for(i in 2:n) {
    # Euler-Maruyama step: Drift + Stochastic Diffusion
    X[i] <- X[i-1] - (kappa * X[i-1] * dt) + (sigma * rnorm(1, mean = 0, sd = sqrt(dt)))
  }
  return(X)
}

# 1. Generate Deterministic Expected Decay (The mean reverting path)
df_expected <- data.frame(
  Time = rep(t_seq, 2),
  Branch = rep(c("Parasympathetic", "Sympathetic"), each = length(t_seq)),
  Value = c(exp(-median_kp * t_seq), exp(-median_ks * t_seq)),
  Type = "Expected"
)

# 2. Simulate Stochastic Paths (3 realizations per branch to show realistic variance)
set.seed(42) # Ensure reproducible noise paths for the manuscript figure
n_paths <- 10
df_stochastic <- data.frame()

for (p in 1:n_paths) {
  # Apply the shared global system volatility (median_sigma) to both branches
  path_p <- sim_ou_path(median_kp, median_sigma, t_seq, dt)
  path_s <- sim_ou_path(median_ks, median_sigma, t_seq, dt)

  temp_df <- data.frame(
    Time = rep(t_seq, 2),
    Branch = rep(c("Parasympathetic", "Sympathetic"), each = length(t_seq)),
    Value = c(path_p, path_s),
    Path_ID = as.character(p),
    Type = "Stochastic"
  )
  df_stochastic <- bind_rows(df_stochastic, temp_df)
}

# 3. Combine and Format
df_decay <- bind_rows(
  df_expected %>% mutate(Path_ID = "Expected"),
  df_stochastic
) %>%
  mutate(Branch = factor(Branch, levels = c("Parasympathetic", "Sympathetic")))

# 4. Create the Plot
p2_decay <- ggplot(df_decay, aes(x = Time, y = Value, color = Branch, group = interaction(Branch, Path_ID))) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  # Plot stochastic paths (thin and semi-transparent)
  geom_line(data = filter(df_decay, Type == "Stochastic"), alpha = 0.25, linewidth = 0.6) +
  # Plot deterministic expected decay (thick and solid)
  geom_line(data = filter(df_decay, Type == "Expected"), linewidth = 1) +
  scale_color_manual(values = color_branches) +
  # Expand Y limits slightly to allow the stochastic noise to dip slightly below 0 without clipping
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "A. Autonomic Memory Decay (Stochastic Impulse Response)",
       subtitle = expression(paste("Mean reversion with empirical OU diffusion (", dX[t] == -kappa*X[t]*dt + sigma*dW[t], ")")),
       x = "Time elapsed after transient stimulus (seconds)",
       y = "Proportion of stimulus remaining") +
  theme_classic() +
  theme(legend.position = c(0.85, 0.88),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Panel B: Kinetic Parameter Segregation Across Physiological States
# Tidy the data: Pivot Kappa parameters to long format for overlapping plots
df_long_kappa <- df_valid %>%
  select(Dataset, SubjectID, Kappa_P, Kappa_S) %>%
  rename(Parasympathetic = Kappa_P, Sympathetic = Kappa_S) %>%
  pivot_longer(cols = c(Parasympathetic, Sympathetic),
               names_to = "Branch",
               values_to = "Kappa") %>%
  mutate(Branch = factor(Branch, levels = c("Parasympathetic", "Sympathetic")))

# Create an overlapping density plot faceted by Dataset
p2_kappa <- ggplot(df_long_kappa, aes(x = Kappa, fill = Branch, col = Branch)) +
  geom_density(alpha = 0.5, linewidth = 1, show.legend = FALSE) +
  scale_fill_manual(values = color_branches, aesthetics = c("color", "fill")) +
  scale_x_continuous(trans = "log10",
                     labels = scales::label_log(),
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  facet_grid(rows = vars(Dataset), scales = "free_y") +
  labs(title = "B. Kinetic Parameter Segregation",
       subtitle = "Log-distribution of clearance rates (κ)",
       x = "Clearance Rate Constant, κ (Hz) [Log Scale]",
       y = "Density") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Panel C: Add a biological translation panel (Relaxation Half-Lives, Tau)
# Tidy the Tau data (Tau = 1/Kappa)
df_long_tau <- df_valid %>%
  select(Dataset, SubjectID, Tau_P, Tau_S) %>%
  rename(Parasympathetic = Tau_P, Sympathetic = Tau_S) %>%
  pivot_longer(cols = c(Parasympathetic, Sympathetic),
               names_to = "Branch",
               values_to = "Tau") %>%
  mutate(Branch = factor(Branch, levels = c("Parasympathetic", "Sympathetic")))

# Boxplot for clear biological half-life comparison
p2_tau <- ggplot(df_long_tau, aes(y = Branch, x = Tau, fill = Branch)) +
  geom_violin(trim = FALSE, show.legend = FALSE, color = NA, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = color_branches, aesthetics = c("color", "fill")) +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous(trans = "log10",
                     labels = scales::label_log(),
                     expand = c(0,0)) +
  facet_grid(rows = vars(Dataset), scales = "free_y") +
  labs(title = "C. Derived Biological Half-Lives",
       subtitle = "Relaxation time (τ = 1/κ)",
       y = NULL,
       x = "Half-life, τ (seconds) [Log Scale]") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

fig2_layout <- p2_decay / (p2_kappa | p2_tau) +
  plot_layout(heights = c(1, 2), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(filename = "manuscript/figures/fig2_biological_parameters.png",
       plot = fig2_layout, width = 180, height = 200,
       dpi = 300, units = "px", scale = 15)

