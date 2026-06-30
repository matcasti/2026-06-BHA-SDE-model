# ==============================================================================
# Generates Figure 4 (Section 3.4). This script isolates a rapid physiological
# transition (Orthostatic Tilt) and compares the instantaneous beat-to-beat
# state-space tracking against a continuous Morlet wavelet transform (CWT),
# explicitly demonstrating the resolution limitations (Heisenberg smearing)
# inherent to classical time-frequency analysis.
# ==============================================================================

# 1. Prepare Workspace ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)
library(biwavelet)

# Define universal branch colors
color_branches <- c("Parasympathetic" = "#4DB0D0",
                    "Sympathetic"     = "#DC5050")
color_wavelet  <- c("HF Power (Wavelet)" = "#85C1E9",
                    "LF Power (Wavelet)"  = "#E59866")

# 2. Load Representative Subject Data ------------------------------------------
raw_models <- readRDS("results/list_raw_models.rds")

# We use the Orthostatic Tilt subject as it contains a sharp, mechanical transition
target_subj <- "12726"
fit_case <- raw_models[[target_subj]]

time_cum <- fit_case$time
time_min <- time_cum / 60
dy <- fit_case$dy

tilt_time_min <- 6 # Center of the transition vector

zoom_start <- 0
zoom_end   <- 15

# Extract SDE States
df_sde <- data.frame(
  Time_Min = time_min,
  X_p = fit_case$states[, "X_p"],
  X_s = fit_case$states[, "X_s"]
) %>% filter(Time_Min >= zoom_start & Time_Min <= zoom_end)


# ==============================================================================
# 3. COMPUTE CONTINUOUS WAVELET TRANSFORM (CWT Moving PSD)
# ==============================================================================

# Interpolate RR intervals to a regular 4 Hz time grid for wavelet analysis
fs <- 4
t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)
hr_hz <- 1 / dy
hr_interp <- spline(time_cum, hr_hz, xout = t_grid)$y
hr_centered <- hr_interp - mean(hr_interp)

# Compute Continuous Wavelet Transform (CWT) matching the initialization engine
data_mat <- cbind(t_grid, hr_centered)
invisible(capture.output(wt_res <- biwavelet::wt(data_mat, do.sig = FALSE)))

# Extract frequency scales
freqs <- 1 / wt_res$period
df <- abs(c(diff(freqs), 0))

# Define physiological boundaries
lf_idx <- which(freqs >= 0.04 & freqs < 0.15)
hf_idx <- which(freqs >= 0.15 & freqs <= 0.40)

# Integrate time-varying energy across bands (columns represent time steps)
wavelet_lf <- apply(wt_res$power[lf_idx, , drop = FALSE], 2, function(x) sum(x * df[lf_idx]))
wavelet_hf <- apply(wt_res$power[hf_idx, , drop = FALSE], 2, function(x) sum(x * df[hf_idx]))

# Build tidy dataframe and normalize scale amplitudes for direct comparison
df_wavelet <- data.frame(
  Time_Min = t_grid / 60,
  LF = wavelet_lf / max(wavelet_lf, na.rm = TRUE),
  HF = wavelet_hf / max(wavelet_hf, na.rm = TRUE)
) %>% filter(Time_Min >= zoom_start & Time_Min <= zoom_end)


# ==============================================================================
# 4. GENERATE BENCHMARKING PANELS
# ==============================================================================

# --- PANEL A: Proposed Generative SDE Architecture ---
pA <- ggplot(df_sde, aes(x = Time_Min)) +
  # Highlight the exact physiological event
  annotate("rect", xmin = tilt_time_min, xmax = tilt_time_min + 3,
           ymin = -Inf, ymax = Inf, fill = "darkorange", alpha = 0.1) +
  geom_vline(xintercept = tilt_time_min, color = "black", linewidth = 0.8) +
  annotate("text", x = tilt_time_min - 0.3, y = max(df_sde$X_p)*3, vjust = 0,
           label = "Exact Stress Onset", angle = 90, fontface = "italic") +

  geom_hline(yintercept = 0, color = "black", alpha = 0.3) +
  geom_line(aes(y = X_p, color = "Parasympathetic"), linewidth = 0.8) +
  geom_line(aes(y = X_s, color = "Sympathetic"), linewidth = 0.8) +
  scale_color_manual(values = color_branches) +
  scale_y_continuous(expand = c(0,0,0.5,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(zoom_start, zoom_end)) +
  labs(title = "A. Proposed Continuous-Time SDE Framework",
       subtitle = "Instantaneous, beat-to-beat tracking of latent autonomic drivers",
       y = "State Amplitude (Hz)", x = "") +
  theme_classic() +
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_blank())

# --- PANEL B: Classical Wavelet Time-Frequency PSD ---
pB <- ggplot(df_wavelet, aes(x = Time_Min)) +
  # Highlight the Heisenberg uncertainty smearing region around the step change
  # LF support requires wider windows, blending pre- and post-tilt spectra
  annotate("rect", xmin = tilt_time_min, xmax = tilt_time_min + 3,
           ymin = -Inf, ymax = Inf, fill = "darkorange", alpha = 0.1) +
  annotate("text", x = tilt_time_min + 0.3, y = 0.3, hjust = 0,
           label = "Heisenberg\nResolution\nSmear Zone") +

  geom_vline(xintercept = tilt_time_min, color = "black", linewidth = 0.8) +
  geom_line(aes(y = HF, color = "HF Power (Wavelet)"), linewidth = 1) +
  geom_line(aes(y = LF, color = "LF Power (Wavelet)"), linewidth = 1) +
  scale_color_manual(values = color_wavelet) +
  scale_y_continuous(expand = c(0,0,0.5,0)) +
  scale_x_continuous(expand = c(0,0), limits = c(zoom_start, zoom_end)) +
  labs(title = "B. Continuous Morlet Wavelet Transform (CWT Moving PSD)",
       subtitle = "Bidirectional temporal blurring and loss of high-frequency precision during sharp changes",
       y = "Normalized Spectral Power", x = "Time (Minutes)") +
  theme_classic() +
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))


# ==============================================================================
# 5. ASSEMBLE AND EXPORT
# ==============================================================================
# Stack panels vertically to ensure aligned time axis grids
fig_benchmark <- (pA / pB)

ggsave("manuscript/figures/fig4_benchmarking_comparison.png", plot = fig_benchmark,
       width = 180, height = 180, dpi = 300, scale = 15, units = "px")
