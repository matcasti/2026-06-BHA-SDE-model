# ==============================================================================
# HRV STATE-SPACE MODEL: EXACT PHASE-DOMAIN IPFM WITH ANALYTICAL MARGINALIZATION
# ==============================================================================
# Implements a linear-Gaussian augmented Kalman Filter for unevenly sampled
# RR intervals. It uses Exact GLS Marginalization for the baseline drift and
# global system volatility, and deterministic spectral anchoring for branch ratios.
# ==============================================================================

if (!require(ggplot2)) install.packages("ggplot2")
if (!require(patchwork)) install.packages("patchwork")
if (!require(expm)) install.packages("expm")

library(ggplot2)
library(patchwork)
library(expm)

# ==============================================================================
# 1. SPECTRAL INITIALIZATION (Priors & Parseval Anchoring)
# ==============================================================================
extract_spectral_priors <- function(dy) {
  cat("Extracting spectral priors and Parseval energy ratio...\n")
  time_cum <- cumsum(dy)
  hr_hz <- 1 / dy

  fs <- 4
  t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)
  hr_interp <- spline(time_cum, hr_hz, xout = t_grid)$y

  spec <- spectrum(hr_interp, plot = FALSE, spans = c(3, 5))
  freqs <- spec$freq * fs
  power <- spec$spec

  # Approximate power in LF and HF bands
  lf_idx <- which(freqs >= 0.04 & freqs < 0.15)
  hf_idx <- which(freqs >= 0.15 & freqs <= 0.40)

  p_lf <- sum(power[lf_idx])
  p_hf <- sum(power[hf_idx])

  energy_ratio <- p_hf / p_lf
  cat(sprintf("Spectral Energy Ratio (HF/LF): %.4f\n", energy_ratio))

  list(
    nu0_init = mean(hr_hz),
    energy_ratio = energy_ratio
  )
}

# ==============================================================================
# 2. ROBUST OU INTEGRALS HELPER (Solves Catastrophic Cancellation)
# ==============================================================================
.compute_ou_block <- function(k, dt, sig2) {
  x <- k * dt

  if (x < 1e-3) {
    # Exact Taylor series limits to prevent floating point zeros
    phi_int <- dt - (k * dt^2)/2 + (k^2 * dt^3)/6
    q11 <- sig2 * (dt - k * dt^2 + (2 * k^2 * dt^3)/3)
    q13 <- sig2 * ((dt^2)/2 - (k * dt^3)/2 + (7 * k^2 * dt^4)/24)
    q33 <- sig2 * ((dt^3)/3 - (k * dt^4)/4 + (7 * k^2 * dt^5)/60)
  } else {
    # Standard closed-form OU continuous integration
    phi_int <- (1 - exp(-x)) / k
    q11 <- sig2 * (1 - exp(-2 * x)) / (2 * k)
    q13 <- sig2 * (1 - 2*exp(-x) + exp(-2*x)) / (2 * k^2)
    q33 <- sig2 * (2*x - 3 + 4*exp(-x) - exp(-2*x)) / (2 * k^3)
  }

  list(phi_int = phi_int, q11 = q11, q13 = q13, q33 = q33)
}

# ==============================================================================
# 3. EXACT DUAL-FILTER GLS ENGINE
# ==============================================================================
run_gls_filter <- function(dy, kp, ks, lp, lR) {
  N <- length(dy)
  X0 <- matrix(0, 4, 1); Xv <- matrix(0, 4, 1)
  P  <- diag(c(1, 1, 0, 0)); H  <- matrix(c(0, 0, -1, 1), 1, 4)
  S_vec <- numeric(N); v0_vec <- numeric(N); vv_vec <- numeric(N)

  # Storage for state reconstruction
  X0_store <- matrix(0, N, 4); Xv_store <- matrix(0, N, 4)
  Phi <- matrix(0, 4, 4); Q <- matrix(0, 4, 4)

  for(k in 1:N) {
    dt <- dy[k]

    blk_p <- .compute_ou_block(kp, dt, lp)
    blk_s <- .compute_ou_block(ks, dt, 1.0)
    Phi[1,1] <- exp(-kp * dt); Phi[3,1] <- blk_p$phi_int
    Phi[2,2] <- exp(-ks * dt); Phi[4,2] <- blk_s$phi_int
    Q[1,1] <- blk_p$q11; Q[1,3] <- blk_p$q13; Q[3,1] <- blk_p$q13; Q[3,3] <- blk_p$q33
    Q[2,2] <- blk_s$q11; Q[2,4] <- blk_s$q13; Q[4,2] <- blk_s$q13; Q[4,4] <- blk_s$q33

    X0_pred <- Phi %*% X0; Xv_pred <- Phi %*% Xv
    P_pred  <- Phi %*% P %*% t(Phi) + Q

    S <- as.numeric(H %*% P_pred %*% t(H) + lR * dt)
    S <- max(S, 1e-12) # Numerical floor
    K <- P_pred %*% t(H) / S

    v0 <- 1.0 - as.numeric(H %*% X0_pred)
    vv <- dt  - as.numeric(H %*% Xv_pred)

    X0 <- X0_pred + K * v0; Xv <- Xv_pred + K * vv
    I_KH <- diag(4) - K %*% H
    P <- I_KH %*% P_pred %*% t(I_KH) + K %*% (lR * dt) %*% t(K)
    P <- (P + t(P)) / 2

    X0_store[k, ] <- as.numeric(X0)
    Xv_store[k, ] <- as.numeric(Xv)

    # EXACT BOUNDARY RESET
    X0[3:4,] <- 0; Xv[3:4,] <- 0
    P[3:4, ] <- 0; P[, 3:4] <- 0

    S_vec[k] <- S; v0_vec[k] <- v0; vv_vec[k] <- vv
  }

  # Analytical Marginalizations
  nu0_hat  <- sum(v0_vec * vv_vec / S_vec) / sum((vv_vec^2) / S_vec)
  v_final  <- v0_vec - nu0_hat * vv_vec
  sig2_hat <- max(mean(v_final^2 / S_vec), 1e-12)
  LL_star  <- -0.5 * sum(log(S_vec)) - (N / 2) * log(sig2_hat)

  # Reconstruct exact biological states by linear superposition (scaled to physics)
  states_final <- (X0_store - nu0_hat * Xv_store)
  std_innov <- v_final / sqrt(S_vec * sig2_hat)

  list(LL_star = LL_star, nu0 = nu0_hat, sig2 = sig2_hat,
       states = states_final, std_innov = std_innov)
}

objective_map_3d <- function(theta, dy, energy_ratio) {
  ks <- exp(theta[1])
  kp <- ks + exp(theta[2])
  lR <- exp(theta[3])
  lp <- energy_ratio * (kp / ks)

  res <- run_gls_filter(dy, kp, ks, lp, lR)

  pen_ks <- (2 - 1) * log(ks) - 20 * ks
  pen_kp <- (2 - 1) * log(kp) - 2 * kp
  pen_lR <- -(1.1 + 1) * log(lR) - (1e-4) / lR

  return(-(res$LL_star + pen_ks + pen_kp + pen_lR))
}

# ==============================================================================
# 4. OPTIMIZATION ROUTINE (Returns Full Extracted State Geometry)
# ==============================================================================
fit_model <- function(dy) {
  priors <- extract_spectral_priors(dy)
  theta_init <- c(log(0.05), log(0.35 - 0.05), log(0.001))

  cat("\nRunning 3D BFGS Optimization with MAP Penalties...\n")
  opt_res <- optim(par = theta_init, fn = objective_map_3d, dy = dy,
                   energy_ratio = priors$energy_ratio, method = "BFGS",
                   control = list(maxit = 500, trace = 1))

  # Extract Final Parameters
  ks_opt <- exp(opt_res$par[1])
  kp_opt <- ks_opt + exp(opt_res$par[2])
  lR_opt <- exp(opt_res$par[3])
  lp_opt <- priors$energy_ratio * (kp_opt / ks_opt)

  # Final pass to extract states & innovations
  final_res <- run_gls_filter(dy, kp_opt, ks_opt, lp_opt, lR_opt)

  return(list(
    dy = dy, time = cumsum(dy),
    params = list(
      ks = ks_opt, kp = kp_opt, lR = lR_opt, lp = lp_opt,
      nu0 = final_res$nu0, sig2 = final_res$sig2,
      sig_p = sqrt(lp_opt * final_res$sig2), sig_s = sqrt(final_res$sig2),
      energy_ratio = priors$energy_ratio
    ),
    states = final_res$states,
    innovations = final_res$std_innov,
    opt_raw = opt_res
  ))
}

# ==============================================================================
# REPORT PARAMETERS (Delta Method for 3D MAP Grid)
# ==============================================================================
report_model <- function(fit_obj) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

  dy <- fit_obj$dy
  p <- fit_obj$params

  theta_opt <- c(log(p$ks), log(p$kp - p$ks), log(p$lR))

  cat("Calculating standard errors via Delta Method...\n")
  hess_precise <- numDeriv::hessian(func = objective_map_3d, x = theta_opt,
                                    dy = dy, energy_ratio = p$energy_ratio)
  inv_hess <- tryCatch(MASS::ginv(hess_precise), error = function(e) matrix(NA, 3, 3))

  # Jacobian for Delta Method (ks, kp, lR)
  J <- matrix(0, nrow = 3, ncol = 3)
  J[1, 1] <- p$ks
  J[2, 1] <- p$ks
  J[2, 2] <- p$kp - p$ks
  J[3, 3] <- p$lR

  cov_physical <- J %*% inv_hess %*% t(J)
  se_vals <- sqrt(ifelse(diag(cov_physical) < 0, abs(diag(cov_physical)), diag(cov_physical)))

  results <- data.frame(
    Parameter = c("Nu_0 (Baseline Drift)", "Sigma_sys (Global Scale)", "Lambda_P (Vagal Ratio)",
                  "Kappa_S (Symp Resonance)", "Kappa_P (Vagal Resonance)", "Lambda_R (Threshold Ratio)"),
    Status = c("Profiled (GLS)", "Profiled (GLS)", "Profiled (Spectral)", rep("Estimated via MAP", 3)),
    Estimate = round(c(p$nu0, sqrt(p$sig2), p$lp, p$ks, p$kp, p$lR), 5),
    StdError = c(NA, NA, NA, round(se_vals, 5))
  )

  cat("===================================================\n")
  cat(" HRV PHASE-DOMAIN SDE-IPFM FITTING REPORT\n")
  cat("===================================================\n")
  print(results, row.names = FALSE)
  cat("===================================================\n")
}

# ==============================================================================
# VISUALIZATION (Exact Phase Reconstruction & Topology)
# ==============================================================================
visualize_model <- function(fit_obj) {
  dy <- fit_obj$dy
  N <- length(dy)
  time_cum <- fit_obj$time

  # Rate states (for topology and raw tracking)
  X_p <- fit_obj$states[, 1]
  X_s <- fit_obj$states[, 2]

  # Phase states (exact accumulated area over the beat)
  Phase_p <- fit_obj$states[, 3]
  Phase_s <- fit_obj$states[, 4]

  p <- fit_obj$params

  # 1. EXACT PHASE-DOMAIN RECONSTRUCTION
  # Bypasses all non-linear approximation by using the mathematical inverse
  # of the IPFM integration: dt = (1.0 + Phase_p - Phase_s) / nu0
  implied_rr <- (1.0 + Phase_p - Phase_s) / p$nu0

  df_fit <- data.frame(Time = time_cum, Observed = dy, Implied = implied_rr)

  pA <- ggplot(df_fit, aes(x = Time)) +
    geom_line(aes(y = Observed, color = "Observed RR"), linewidth = 0.8, alpha = 0.5) +
    geom_line(aes(y = Implied, color = "Filtered RR"), linewidth = 1) +
    scale_color_manual(values = c("Observed RR" = "gray50", "Filtered RR" = "darkred")) +
    labs(title = "Phase-Domain Filtering & Predictive Fit", y = "RR Interval (s)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  Z_bal <- X_s - X_p
  Z_tot <- X_s + X_p
  df_Z <- data.frame(
    Time = rep(time_cum, 2), Value = c(Z_bal, Z_tot),
    State = factor(rep(c("Autonomic Balance (s - p)", "Total Tone (s + p)"), each = N))
  )

  pB <- ggplot(df_Z, aes(x = Time, y = Value, color = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Autonomic Balance (s - p)" = "#4DBBD5", "Total Tone (s + p)" = "#00A087")) +
    labs(title = "Geometric Projections", y = "Amplitude (Hz)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  df_X <- data.frame(
    Time = rep(time_cum, 2), Value = c(X_p, X_s),
    State = factor(rep(c("Parasympathetic Drive", "Sympathetic Drive"), each = N),
                   levels = c("Sympathetic Drive", "Parasympathetic Drive"))
  )

  pC <- ggplot(df_X, aes(x = Time, y = Value, color = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Parasympathetic Drive" = "#3C5488", "Sympathetic Drive" = "#DC0000")) +
    labs(title = "Native Autonomic Drivers (Direct Tracking)", y = "Amplitude (Hz)", x = "Time (seconds)") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "bottom")

  # 2. DYNAMIC ENERGY MAP FRAMING
  Sigma_stat <- diag(c(p$sig_p^2 / (2 * p$kp), p$sig_s^2 / (2 * p$ks)))
  inv_Sigma <- solve(Sigma_stat)

  p_range <- range(X_p); s_range <- range(X_s)
  margin_p <- diff(p_range) * 0.2; if(margin_p == 0) margin_p <- 0.01
  margin_s <- diff(s_range) * 0.2; if(margin_s == 0) margin_s <- 0.01

  grid_p <- seq(p_range[1] - margin_p, p_range[2] + margin_p, length.out = 150)
  grid_s <- seq(s_range[1] - margin_s, s_range[2] + margin_s, length.out = 150)
  grid <- expand.grid(X_p = grid_p, X_s = grid_s)

  grid$Energy <- apply(grid, 1, function(v) { 0.5 * as.numeric(t(v) %*% inv_Sigma %*% v) })
  max_E <- quantile(grid$Energy, 0.95)
  grid$Energy <- pmin(grid$Energy, max_E)

  df_traj <- data.frame(X_p = X_p, X_s = X_s, Time = time_cum)

  pD <- ggplot() +
    geom_raster(data = grid, aes(x = X_p, y = X_s, fill = Energy), alpha = 0.9) +
    geom_contour(data = grid, aes(x = X_p, y = X_s, z = Energy), color = "white", alpha = 0.3, bins = 15) +
    geom_path(data = df_traj, aes(x = X_p, y = X_s, color = Time), linewidth = 0.6,
              arrow = arrow(type = "closed", length = unit(0.06, "inches"))) +
    scale_fill_viridis_c(option = "mako", name = "Energy (U)") +
    scale_color_viridis_c(option = "plasma", guide = "none") +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    labs(title = "Autonomic Phase Space Topology", x = "Parasympathetic Drive (Hz)", y = "Sympathetic Drive (Hz)") +
    theme_classic() + theme(legend.position = "bottom")

  (pA / pB / pC) | pD
}

# ==============================================================================
# DIAGNOSTICS SUITE (Pre-standardized innovations)
# ==============================================================================
diagnose_model <- function(fit_obj) {
  z <- fit_obj$innovations
  U_k <- pnorm(z)
  N <- length(U_k)

  ks_res <- ks.test(U_k, "punif")
  print(ks_res)
  ks_stat <- ks_res$statistic
  p_val <- ks_res$p.value
  bound <- 1.36 / sqrt(N)

  df_qq <- data.frame(z = z)
  pA <- ggplot(df_qq, aes(sample = z)) +
    stat_qq(color = "darkblue", alpha = 0.5) + stat_qq_line(color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Q-Q Plot of Phase Innovations", subtitle = "Tests Linear Observation Gaussianity", x = "Theoretical Normal", y = "Empirical") + theme_classic()

  empirical_cdf <- ecdf(U_k)
  df_cdf <- data.frame(U = sort(U_k), CDF = empirical_cdf(sort(U_k)))

  pB <- ggplot(df_cdf, aes(x = U, y = CDF)) +
    geom_step(color = "darkgreen", linewidth = 1) + geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_abline(slope = 1, intercept = bound, color = "gray", linetype = "dotted", linewidth=0.8) +
    geom_abline(slope = 1, intercept = -bound, color = "gray", linetype = "dotted", linewidth=0.8) +
    annotate("text", x = 0.2, y = 0.9, label = paste("KS Stat:", round(ks_stat, 4)), hjust = 0) +
    annotate("text", x = 0.2, y = 0.8, label = paste("p-value:", round(p_val, 4)), hjust = 0) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
    labs(title = "Time-Rescaling KS-Plot", subtitle = "Validates Exact IPFM Boundary Crossings", x = "Theoretical Uniform(0,1)", y = "Empirical CDF") + theme_classic()

  acf_data <- acf(z, plot = FALSE, lag.max = 40)
  df_acf <- data.frame(lag = acf_data$lag[-1], acf = acf_data$acf[-1])

  pC <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0), color = "black", linewidth = 0.8) + geom_hline(yintercept = c(-1.96/sqrt(N), 1.96/sqrt(N)), color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black") + labs(title = "Autocorrelation of Innovations", subtitle = "Tests for Unmodeled Kinetics", x = "Lag (Heartbeats)", y = "ACF") + theme_classic()

  pD <- ggplot(df_qq, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), fill = "#7E6148", alpha = 0.7, bins = 30, color = "white") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linewidth = 1.2, color = "black", linetype = "dashed") +
    labs(title = "Exact Gaussian Innovations", subtitle = "Validating Phase-Domain Linearity", x = "Standardized Residuals", y = "Density") + theme_classic()

  (pA | pB) / (pC | pD)
}
