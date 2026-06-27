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

# --- Panel A: Engineering Tracking Error (RMSE vs MAPE) ---
# Log-log scatter plot to show tightly bounded tracking performance
p1A <- ggplot(df_valid, aes(x = MAPE_pct, y = RMSE_ms, fill = Dataset)) +
  facet_grid(rows = vars(Dataset), scales = "free_y") +
  geom_point(shape = 21, color = "white", size = 3, alpha = 0.8, stroke = 0.5) +
  scale_fill_manual(values = color_cohorts, name = "Cohort", aesthetics = c("fill", "color")) +
  scale_x_log10(label = scales::label_log()) +
  scale_y_log10(label = scales::label_log(digits = 1)) +
  geom_smooth(method = "glm", aes(color = Dataset, fill = Dataset)) +
  labs(title = "A. Macroscopic Tracking Accuracy",
       subtitle = "Bounded predictive error across cohorts",
       x = "Mean Absolute Percentage Error (%)",
       y = "Root Mean Square Error (ms)") +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

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
  facet_grid(rows = vars(Dataset)) +
  geom_rect(data = ci_bounds, aes(xmin = -Inf, xmax = Inf, ymax = ci_bound, ymin = -ci_bound, x = NA, y = 0),
            fill = "gray80", col = 0, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Plot dataset-specific ribbons and lines
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_manual(values = color_cohorts) +
  scale_fill_manual(values = color_cohorts) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "B. Cohort-Specific PACF",
       subtitle = "Residual innovations reduced to white noise",
       x = "Cardiac Lag (Beats)",
       y = "Partial Autocorrelation") +
  theme_classic() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Panel C: Mathematical Exactness (KS-Distance Density) ---

p1C <- ggplot(df_valid, aes(x = KS_Stat, fill = Dataset, col = Dataset)) +
  facet_grid(rows = vars(Dataset)) +
  geom_density(alpha = 0.7, linewidth = 1, trim = FALSE) +
  scale_fill_manual(values = color_cohorts, name = "Cohort", aesthetics = c("fill", "color")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,NA)) +
  scale_y_continuous(expand = c(0,0,0.25,0)) +
  labs(title = "C. Time-Rescaling Theorem Proof",
       subtitle = "KS D-statistic aggregation (Geometric Distance)",
       x = "Kolmogorov-Smirnov Distance (D-Statistic)",
       y = "Density") +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Assemble Figure 1 ---
fig1_layout <- (p1A | p1B | p1C) +
  plot_layout(widths = c(1, 1, 1))

ggsave(filename = "manuscript/figures/fig1_statistical_validation.png",
       plot = fig1_layout, width = 250, height = 100, dpi = 300, scale = 15,
       units = "px")


# ==============================================================================
# FIGURE 2: BIOLOGICAL CONSTRUCT VALIDITY (Section 3.2)
# ==============================================================================

# --- Panel A: Autonomic Memory Decay (Impulse Response) ---
# Calculate the global median clearance rates to draw the decay curves
median_kp <- median(df_valid$Kappa_P, na.rm = TRUE)
median_ks <- median(df_valid$Kappa_S, na.rm = TRUE)

# Generate a time vector from 0 to roughly 3x the sympathetic half-life (e.g., 0 to 30 seconds)
t_max <- 8 * (1/median_ks)
t_seq <- seq(0, t_max, length.out = 300)

# Calculate the Ornstein-Uhlenbeck exponential decay: exp(-kappa * t)
df_decay <- data.frame(
  Time = rep(t_seq, 2),
  Branch = rep(c("Parasympathetic", "Sympathetic"), each = length(t_seq)),
  Decay = c(exp(-median_kp * t_seq), exp(-median_ks * t_seq))
) %>%
  mutate(Branch = factor(Branch, levels = c("Parasympathetic", "Sympathetic")))

# Create the decay plot
p2_decay <- ggplot(df_decay, aes(x = Time, y = Decay, color = Branch)) +
  # Add a dashed line at 50% to visually cross-reference the half-life
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = color_branches) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.05)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "A. Autonomic Memory Decay (Ornstein-Uhlenbeck Impulse Response)",
       subtitle = expression(paste("Mean reversion of transient stimuli (", e^{-kappa*t}, ") based on population median half-lives")),
       x = "Time elapsed after transient stimulus (seconds)",
       y = "Proportion of stimulus remaining") +
  theme_classic() +
  theme(legend.position = "none",
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
  geom_density(alpha = 0.5, linewidth = 1) +
  scale_fill_manual(values = color_branches, aesthetics = c("color", "fill")) +
  scale_x_continuous(trans = "log10",
                     labels = scales::label_log(),
                     expand = c(0,0)) +
  facet_grid(rows = vars(Dataset), scales = "free_y") +
  labs(title = "B. Kinetic Parameter Segregation Across Physiological States",
       subtitle = "Log-distribution of clearance rates (κ)",
       x = "Clearance Rate Constant, κ (Hz) [Log Scale]",
       y = "Density") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_blank(),
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
  geom_violin(trim = FALSE) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.5, width = 0.1, fill = "white") +
  scale_fill_manual(values = color_branches, aesthetics = c("color", "fill")) +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous(trans = "log10",
                     labels = scales::label_log()) +
  facet_grid(rows = vars(Dataset), scales = "free_y",
             labeller = labeller(Dataset = c("Fantasia_Baseline" = "Resting Baseline (Fantasia)",
                                             "Orthostatic_Tilt" = "Orthostatic Tilt-Table",
                                             "Exercise_Recovery" = "Exercise & Recovery"))) +
  labs(title = "C. Derived Biological Half-Lives",
       subtitle = "Relaxation time (τ = 1/κ)",
       y = "",
       x = "Half-life, τ (seconds) [Log Scale]") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

fig2_layout <- p2_decay / (p2_kappa | p2_tau) +
  plot_layout(heights = c(1, 2))

ggsave(filename = "manuscript/figures/fig2_biological_parameters.png",
       plot = fig2_layout, width = 240, height = 150,
       dpi = 300, units = "px", scale = 15)

