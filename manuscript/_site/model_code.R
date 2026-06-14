
# Prepare workspace -------------------------------------------------------

## Import libraries
library(biwavelet)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. SPECTRAL INITIALIZATION (Priors & Parseval Anchoring via Morlet CWT)
# ==============================================================================
extract_spectral_priors <- function(dy) {
  cat("Extracting Parseval energy ratio via Morlet Wavelet integration...\n")
  time_cum <- cumsum(dy)
  hr_hz <- 1 / dy

  # Standardize grid to 4Hz (0.25s) for continuous integration
  fs <- 4
  t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)
  hr_interp <- spline(time_cum, hr_hz, xout = t_grid)$y

  # Mean-center the signal to prevent DC artifact dominating the wavelet scales
  hr_centered <- hr_interp - mean(hr_interp)

  # 1. Continuous Wavelet Transform (Morlet)
  # biwavelet expects a 2-column matrix: [Time, Value]
  data_mat <- cbind(t_grid, hr_centered)

  # Run CWT (suppress noisy console outputs from the package)
  invisible(capture.output(
    wt_res <- biwavelet::wt(data_mat, do.sig = FALSE)
  ))

  # 2. Extract the Global Wavelet Spectrum
  # Average the time-frequency power surface across the time dimension
  global_power <- rowMeans(wt_res$power)
  freqs <- 1 / wt_res$period

  # 3. Integrate Power Over Physiological Bands
  lf_idx <- which(freqs >= 0.04 & freqs < 0.15)
  hf_idx <- which(freqs >= 0.15 & freqs <= 0.40)

  # Because wavelet scales are logarithmic, we multiply the global power
  # by the frequency differential (df) for a true mathematical integral.
  df <- abs(c(diff(freqs), 0))

  p_lf <- sum(global_power[lf_idx] * df[lf_idx])
  p_hf <- sum(global_power[hf_idx] * df[hf_idx])

  energy_ratio <- p_hf / p_lf
  cat(sprintf("Wavelet Energy Ratio (HF/LF): %.4f\n", energy_ratio))

  # Return the intrinsic baseline drift and the robust spectral ratio
  list(
    nu0_init = mean(hr_hz),
    energy_ratio = energy_ratio
  )
}

# ==============================================================================
# 2. OU INTEGRALS HELPER
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
# 3. DUAL-FILTER GLS
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

    # BOUNDARY RESET
    X0[3:4,] <- 0; Xv[3:4,] <- 0
    P[3:4, ] <- 0; P[, 3:4] <- 0

    S_vec[k] <- S; v0_vec[k] <- v0; vv_vec[k] <- vv
  }

  # Analytical Marginalizations
  nu0_hat  <- sum(v0_vec * vv_vec / S_vec) / sum((vv_vec^2) / S_vec)
  v_final  <- v0_vec - nu0_hat * vv_vec
  sig2_hat <- max(mean(v_final^2 / S_vec), 1e-12)
  LL_star  <- -0.5 * sum(log(S_vec)) - (N / 2) * log(sig2_hat)

  # Reconstruct biological states by linear superposition
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
# 4. OPTIMIZATION ROUTINE
# ==============================================================================
fit_model <- function(dy) {
  priors <- extract_spectral_priors(dy)
  theta_init <- c(log(0.05), log(0.35 - 0.05), log(0.001))

  cat("\nRunning 3D BFGS Optimization with MAP Penalties...\n")
  opt_res <- optim(par = theta_init, fn = objective_map_3d, dy = dy,
                   energy_ratio = priors$energy_ratio, method = "BFGS",
                   control = list(maxit = 1000, trace = 1))

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
# REPORT PARAMETERS
# ==============================================================================
report_model <- function(fit_obj) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

  dy <- fit_obj$dy
  p <- fit_obj$params

  # Transformed optimization space
  theta_opt <- c(log(p$ks), log(p$kp - p$ks), log(p$lR))

  cat("\nCalculating standard errors via exact Delta Method...\n")
  hess_precise <- numDeriv::hessian(
    func = objective_map_3d,
    x = theta_opt,
    dy = dy,
    energy_ratio = p$energy_ratio
  )

  # Safely invert Hessian for covariance matrix
  inv_hess <- tryCatch(MASS::ginv(hess_precise), error = function(e) matrix(NA, 3, 3))

  # Jacobian for Delta Method mapping from Theta to Physical (ks, kp, lR)
  J <- matrix(0, nrow = 3, ncol = 3)
  J[1, 1] <- p$ks
  J[2, 1] <- p$ks
  J[2, 2] <- p$kp - p$ks
  J[3, 3] <- p$lR

  cov_physical <- J %*% inv_hess %*% t(J)

  # Extract explicit standard errors for the numerical parameters
  se_ks <- sqrt(max(cov_physical[1, 1], 0))
  se_kp <- sqrt(max(cov_physical[2, 2], 0))
  se_lR <- sqrt(max(cov_physical[3, 3], 0))

  # Delta Method extension for biological time constants (tau = 1/kappa)
  # d(1/k)/dk = -1/k^2 -> Var(tau) = Var(k) / k^4
  se_taus <- sqrt(max(cov_physical[1, 1] / (p$ks^4), 0))
  se_taup <- sqrt(max(cov_physical[2, 2] / (p$kp^4), 0))

  # Assembly of the Comprehensive Report Dictionary
  params_list <- list(
    list(name = "Nu_0", group = "Baseline", interp = "Intrinsic SA node pacing rate absent autonomic tone (Hz)", status = "Profiled (GLS)", est = p$nu0, se = NA),
    list(name = "Sigma_sys^2", group = "Scale", interp = "Total stochastic system variance multiplier", status = "Profiled (GLS)", est = p$sig2, se = NA),

    list(name = "Kappa_S", group = "Kinetics", interp = "Clearance rate of sympathetic neurotransmitters (Hz)", status = "Estimated (MAP)", est = p$ks, se = se_ks),
    list(name = "Kappa_P", group = "Kinetics", interp = "Clearance rate of vagal neurotransmitters (Hz)", status = "Estimated (MAP)", est = p$kp, se = se_kp),
    list(name = "Tau_S", group = "Kinetics", interp = "Relaxation half-life of sympathetic response (s)", status = "Derived (MAP)", est = 1/p$ks, se = se_taus),
    list(name = "Tau_P", group = "Kinetics", interp = "Relaxation half-life of parasympathetic response (s)", status = "Derived (MAP)", est = 1/p$kp, se = se_taup),

    list(name = "Sigma_S", group = "Volatility", interp = "Absolute amplitude of sympathetic neural drive", status = "Profiled (GLS)", est = p$sig_s, se = NA),
    list(name = "Sigma_P", group = "Volatility", interp = "Absolute amplitude of parasympathetic neural drive", status = "Derived (Spectral)", est = p$sig_p, se = NA),

    list(name = "Lambda_P", group = "Ratios", interp = "Vagal-to-Sympathetic energy balance (Parseval anchor)", status = "Profiled (Spectral)", est = p$lp, se = NA),
    list(name = "Lambda_R", group = "Ratios", interp = "Fraction of variance from SA node threshold jitter", status = "Estimated (MAP)", est = p$lR, se = se_lR)
  )

  # Format into a strictly aligned data frame
  df <- do.call(rbind, lapply(params_list, function(x) {
    if (is.na(x$se)) {
      ci_str <- "-"
      se_str <- "-"
    } else {
      ci_low <- x$est - 1.96 * x$se
      ci_upp <- x$est + 1.96 * x$se
      ci_str <- sprintf("[%.4f, %.4f]", ci_low, ci_upp)
      se_str <- sprintf("%.5f", x$se)
    }

    data.frame(
      Parameter = x$name,
      Estimate = sprintf("%.5f", x$est),
      SE = se_str,
      `CI_95%` = ci_str,
      Status = x$status,
      Interpretation = x$interp,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  # Console Print Layout
  cat("\n=========================================================================================================\n")
  cat("                              HRV PHASE-DOMAIN SDE-IPFM FITTING REPORT\n")
  cat("=========================================================================================================\n")
  print(df, row.names = FALSE, right = FALSE)
  cat("=========================================================================================================\n")
  cat("* Note: Parameters marked as 'Profiled' are solved exactly via analytical analytical integration \n")
  cat("  (GLS/Concentrated Likelihood) or deterministic continuous wavelet mapping. Because they are integrated \n")
  cat("  out of the optimization search space, they do not possess numerical MLE standard errors.\n\n")
}

# ==============================================================================
# VISUALIZATION (Phase Reconstruction & Topology)
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

  # 1. PHASE-DOMAIN RECONSTRUCTION
  # Uses the mathematical inverse of the IPFM integration: dt = (1.0 + Phase_p - Phase_s) / nu0
  implied_rr <- (1.0 + Phase_p - Phase_s) / p$nu0

  df_fit <- data.frame(Time = time_cum, Observed = dy, Implied = implied_rr)

  pA <- ggplot(df_fit, aes(x = Time)) +
    geom_line(aes(y = Observed, color = "Observed RR"), linewidth = 1, alpha = 0.5) +
    geom_line(aes(y = Implied, color = "Filtered RR"), linewidth = 1/3) +
    scale_color_manual(values = c("Observed RR" = "gray50", "Filtered RR" = "darkred")) +
    labs(title = "Phase-Domain Filtering & Predictive Fit", y = "RR Interval (s)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  Z_bal <- X_p - X_s
  Z_tot <- X_s + X_p
  df_Z <- data.frame(
    Time = rep(time_cum, 2), Value = c(Z_bal, Z_tot),
    State = factor(rep(c("Autonomic Balance (p - s)", "Total Tone (p + s)"), each = N))
  )

  pB <- ggplot(df_Z, aes(x = Time, y = Value, color = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_line(linewidth = 1/3) +
    scale_color_manual(values = c("Autonomic Balance (p - s)" = "#4DBBD5", "Total Tone (p + s)" = "#00A087")) +
    labs(title = "Geometric Projections", y = "Amplitude (Hz)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  df_X <- data.frame(
    Time = rep(time_cum, 2), Value = c(X_p, X_s),
    State = factor(rep(c("Parasympathetic Drive", "Sympathetic Drive"), each = N),
                   levels = c("Sympathetic Drive", "Parasympathetic Drive"))
  )

  pC <- ggplot(df_X, aes(x = Time, y = Value, color = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_line(linewidth = 1/3) +
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
# DIAGNOSTICS (Pre-standardized innovations)
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
  df_acf <- data.frame(lag = acf_data$lag, acf = acf_data$acf)

  pC <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0), color = "black", linewidth = 0.8) + geom_hline(yintercept = c(-1.96/sqrt(N), 1.96/sqrt(N)), color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black") + labs(title = "Autocorrelation of Innovations", subtitle = "Tests for Unmodeled Kinetics", x = "Lag (Heartbeats)", y = "ACF") + theme_classic()

  pD <- ggplot(df_qq, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)), fill = "#7E6148", alpha = 0.7, bins = 30, color = "white") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linewidth = 1.2, color = "black", linetype = "dashed") +
    labs(title = "Exact Gaussian Innovations", subtitle = "Validating Phase-Domain Linearity", x = "Standardized Residuals", y = "Density") + theme_classic()

  (pA | pB) / (pC | pD)
}
