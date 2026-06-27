# ==============================================================================
# Generates detailed Phase-Domain reconstruction panels for two representative
# subjects (Orthostatic Tilt and Exercise Recovery). It maps 95% CIs, exact
# thermodynamic phase-space energy, and compares the empirical vs. filtered PSD.
# ==============================================================================

# 1. Prepare Workspace ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(stats)

# Define universal branch colors
color_branches <- c("Parasympathetic" = "#4DB0D0",
                    "Sympathetic"     = "#DC5050")

# 2. Load Processed Models and Summaries ---------------------------------------
df_summary <- readRDS("results/df_population_summary.rds")
raw_models <- readRDS("results/list_raw_models.rds")

# Dynamically select representative subjects (longest valid recording for each protocol)
target_tilt <- df_summary %>%
  filter(Dataset == "Orthostatic_Tilt", Converged == TRUE) %>%
  arrange(desc(NumBeats)) %>%
  slice(1) %>% pull(SubjectID)

target_exercise <- df_summary %>%
  filter(Dataset == "Exercise_Recovery", Converged == TRUE) %>%
  arrange(desc(NumBeats)) %>%
  slice(1) %>% pull(SubjectID)

# ==============================================================================
# 3. CORE PLOTTING FUNCTION (Upgraded for Paper Aesthetics)
# ==============================================================================
generate_topology_figure <- function(fit_obj, protocol_name = "Protocol", n = NULL) {

  if(is.null(n) | !is.numeric(n)) {
    n <- length(fit_obj$dy)
  }
  dy <- fit_obj$dy[seq_len(length.out = n)]
  N <- length(dy)
  time_cum <- fit_obj$time[seq_len(length.out = n)]
  time_min <- time_cum / 60 # Convert to minutes for cleaner X-axes
  p <- fit_obj$params

  # DIRECT STATE EXTRACTION
  X_p <- fit_obj$states[, "X_p"][seq_len(length.out = n)]
  X_s <- fit_obj$states[, "X_s"][seq_len(length.out = n)]
  X_p_se  <- sqrt(pmax(fit_obj$states[, "X_p_var"][seq_len(length.out = n)], 0))
  X_s_se  <- sqrt(pmax(fit_obj$states[, "X_s_var"][seq_len(length.out = n)], 0))
  rr_se   <- sqrt(pmax(fit_obj$states[, "rr_var"][seq_len(length.out = n)], 0))
  Phase_p <- fit_obj$states[, "Phase_p"][seq_len(length.out = n)]
  Phase_s <- fit_obj$states[, "Phase_s"][seq_len(length.out = n)]

  # ----------------------------------------------------------------------------
  # PANEL A: PHASE-DOMAIN RECONSTRUCTION WITH 95% CI
  # ----------------------------------------------------------------------------
  implied_rr <- (1.0 + Phase_p - Phase_s) / p$nu0
  df_fit <- data.frame(
    Time_Min = time_min, Observed = dy, Implied = implied_rr,
    Lower = implied_rr - qnorm(0.975) * rr_se,
    Upper = implied_rr + qnorm(0.975) * rr_se
  )

  pA <- ggplot(df_fit, aes(x = Time_Min)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.15) +
    geom_line(aes(y = Observed, color = "Observed RR"), linewidth = 0.5, alpha = 0.5) +
    geom_line(aes(y = Implied, color = "Filtered RR"), linewidth = 0.5) +
    scale_color_manual(values = c("Observed RR" = "gray60", "Filtered RR" = "darkred")) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "I. Phase-Domain Filtering",
         subtitle = "Observed (gray) and Filtered (red)",
         y = "RR Interval (s)", x = "Time (minutes)") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"))

  # ----------------------------------------------------------------------------
  # PANEL B: NATIVE AUTONOMIC DRIVERS WITH 95% CI
  # ----------------------------------------------------------------------------
  df_X <- data.frame(
    Time_Min = rep(time_min, 2),
    Value = c(X_p, X_s),
    SE    = c(X_p_se, X_s_se),
    State = factor(rep(c("Parasympathetic", "Sympathetic"), each = N),
                   levels = c("Parasympathetic", "Sympathetic"))
  )
  df_X$Lower <- df_X$Value - qnorm(0.975) * df_X$SE
  df_X$Upper <- df_X$Value + qnorm(0.975) * df_X$SE

  pB <- ggplot(df_X, aes(x = Time_Min, y = Value, color = State, fill = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.15, color = NA) +
    geom_line() +
    scale_color_manual(values = color_branches) +
    scale_fill_manual(values = color_branches) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "II. Latent Autonomic Drivers",
         subtitle = "Sympathetic (red) and Parasympathetic (blue)",
         y = "Amplitude (Hz)", x = "Time (minutes)") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"))

  # ----------------------------------------------------------------------------
  # PANEL C: DYNAMIC ENERGY MAP FRAMING (Topology)
  # Exact Mahalanobis energy: U = 0.5 * x^T * Sigma^-1 * x
  # ----------------------------------------------------------------------------
  Sigma_stat <- diag(c(p$sig_p^2 / (2 * p$kp), p$sig_s^2 / (2 * p$ks)))
  inv_Sigma  <- solve(Sigma_stat)

  p_range <- range(X_p); s_range <- range(X_s)
  margin_p <- diff(p_range) * 0.2; if(margin_p == 0) margin_p <- 0.01
  margin_s <- diff(s_range) * 0.2; if(margin_s == 0) margin_s <- 0.01

  grid_p <- seq(p_range[1] - margin_p, p_range[2] + margin_p, length.out = 150)
  grid_s <- seq(s_range[1] - margin_s, s_range[2] + margin_s, length.out = 150)
  grid   <- expand.grid(X_p = grid_p, X_s = grid_s)

  grid$Energy <- apply(grid, 1, function(v) { 0.5 * as.numeric(t(v) %*% inv_Sigma %*% v) })
  max_E <- quantile(grid$Energy, 0.98) # Cap extreme outliers for smoother gradient
  grid$Energy <- pmin(grid$Energy, max_E)

  df_traj <- data.frame(X_p = X_p, X_s = X_s, Time = time_cum)

  pC <- ggplot() +
    geom_raster(data = grid, aes(x = X_p, y = X_s, fill = Energy)) +
    geom_contour(data = grid, aes(x = X_p, y = X_s, z = Energy),
                 color = "white", alpha = 0.2, bins = 15) +
    geom_hline(yintercept = 0, color = "white", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "white", linetype = "dashed", alpha = 0.5) +
    geom_path(data = df_traj, aes(x = X_p, y = X_s, color = Time),
              arrow = arrow(type = "closed", length = unit(0.08, "inches")),
              linewidth = 0.7, show.legend = FALSE) +
    scale_fill_viridis_c(option = "mako", name = "Thermodynamic\nPotential (U)", direction = -1) +
    scale_color_viridis_c(option = "plasma", guide = "none") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(title = "IV. Phase Space Topology & Trajectory",
         subtitle = "Exact Potential Autonomic Energy Basin",
         x = expression("Parasympathetic Drive " * (X[p])),
         y = expression("Sympathetic Drive " * (X[s]))) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  # ----------------------------------------------------------------------------
  # PANEL D: EMPIRICAL VS RECONSTRUCTED PSD
  # ----------------------------------------------------------------------------
  fs <- 4
  t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)

  hr_obs <- 1 / dy
  hr_obs_interp <- spline(time_cum, hr_obs, xout = t_grid)$y
  hr_obs_ts <- ts(hr_obs_interp - mean(hr_obs_interp), frequency = fs)
  spec_obs <- spectrum(hr_obs_ts, spans = c(3, 5), plot = FALSE)

  hr_filt <- 1 / implied_rr
  hr_filt_interp <- spline(time_cum, hr_filt, xout = t_grid)$y
  hr_filt_ts <- ts(hr_filt_interp - mean(hr_filt_interp), frequency = fs)
  spec_filt <- spectrum(hr_filt_ts, spans = c(3, 5), plot = FALSE)

  df_psd <- data.frame(
    Frequency = rep(spec_obs$freq, 2),
    Power     = c(spec_obs$spec, spec_filt$spec),
    Component = factor(rep(c("Observed PSD", "Filtered PSD"), each = length(spec_obs$freq)),
                       levels = c("Observed PSD", "Filtered PSD"))
  )

  df_psd <- df_psd[df_psd$Frequency >= 0.01 & df_psd$Frequency <= 0.50, ]
  y_floor <- min(df_psd$Power, na.rm = TRUE) * 0.5

  pD <- ggplot(df_psd, aes(x = Frequency, y = Power, color = Component)) +
    annotate("rect", xmin = 0.04, xmax = 0.15, ymin = y_floor, ymax = Inf, fill = "gray50", alpha = 0.15) +
    annotate("rect", xmin = 0.15, xmax = 0.40, ymin = y_floor, ymax = Inf, fill = "gray50", alpha = 0.08) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("Observed PSD" = "gray60", "Filtered PSD" = "darkred")) +
    scale_x_continuous(expand = c(0,0), breaks = c(0.04, 0.15, 0.40)) +
    scale_y_continuous(transform = "log10", labels = scales::label_log(), expand = c(0,0)) +
    labs(title = "III. Spectral Conservation",
         subtitle = "Observed (gray) vs Filtered PSD (red)",
         x = "Frequency (Hz)", y = "Power (log scale)") +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.title = element_text(face = "bold"))

  # ----------------------------------------------------------------------------
  # ASSEMBLE COMPOSITE (Patchwork)
  # ----------------------------------------------------------------------------
  # Left column: A, B, D stacked. Right column: C spanning the height.
  composite <- (pA / pB / pD / pC) +
    plot_layout(heights = c(2,2,1,3)) +
    plot_annotation(
      title = protocol_name,
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )

  return(composite)
}

# ==============================================================================
# 4. GENERATE AND EXPORT FIGURES (Side-by-Side Comparison)
# ==============================================================================

## Generating Topology Figure for Tilt Protocol...
fit_tilt <- raw_models[["12726"]]
fig_tilt <- generate_topology_figure(fit_tilt, protocol_name = "A. Orthostatic Tilt Transition", n = 1100)

## Generating Topology Figure for Exercise Protocol...
fit_exercise <- raw_models[["007"]]
fig_exercise <- generate_topology_figure(fit_exercise, protocol_name = "B. Exercise & Metabolic Recovery")

## Combining into side-by-side composite...
fig_combined <- wrap_elements(fig_tilt) | wrap_elements(fig_exercise)

ggsave(filename = "manuscript/figures/fig3_dynamic_topology_comparison.png",
       plot = fig_combined,
       width = 200, height = 300, dpi = 300, scale = 15, units = "px")
