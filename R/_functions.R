
# Prepare workspace -------------------------------------------------------

## Import libraries
library(ggplot2)
library(biwavelet)
library(patchwork)

# ==============================================================================
# 1. SPECTRAL INITIALIZATION (Amplitude & Bandwidth Mapping via CWT)
# ==============================================================================
extract_spectral_priors <- function(dy) {
  cat("Extracting Amplitude and Bandwidth priors via Morlet Wavelet...\n")
  time_cum <- cumsum(dy)
  hr_hz <- 1 / dy

  fs <- 4
  t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)
  hr_interp <- spline(time_cum, hr_hz, xout = t_grid)$y
  hr_centered <- hr_interp - mean(hr_interp)

  data_mat <- cbind(t_grid, hr_centered)
  invisible(capture.output(wt_res <- biwavelet::wt(data_mat, do.sig = FALSE)))

  global_power <- rowMeans(wt_res$power)
  freqs <- 1 / wt_res$period
  df <- abs(c(diff(freqs), 0))

  lf_idx <- which(freqs >= 0.04 & freqs < 0.15)
  hf_idx <- which(freqs >= 0.15 & freqs <= 0.60)

  # 1. AMPLITUDE MAPPING (Total Band Energy -> V)
  p_lf <- sum(global_power[lf_idx] * df[lf_idx])
  p_hf <- sum(global_power[hf_idx] * df[hf_idx])
  energy_ratio <- p_hf / p_lf

  # 2. BEHAVIOR MAPPING (Spectral Bandwidth -> Kappa)
  fc_lf <- sum(freqs[lf_idx] * global_power[lf_idx] * df[lf_idx]) / p_lf
  fc_hf <- sum(freqs[hf_idx] * global_power[hf_idx] * df[hf_idx]) / p_hf

  # Calculate spectral standard deviation (Bandwidth)
  bw_lf <- sqrt(sum((freqs[lf_idx] - fc_lf)^2 * global_power[lf_idx] * df[lf_idx]) / p_lf)
  bw_hf <- sqrt(sum((freqs[hf_idx] - fc_hf)^2 * global_power[hf_idx] * df[hf_idx]) / p_hf)

  # Convert bandwidth to angular corner frequency (Kappa MAP mode)
  ks_mode <- 2 * pi * bw_lf
  kp_mode <- 2 * pi * bw_hf

  # Safety fallbacks for ultra-clean signals
  if(is.na(ks_mode) || ks_mode < 0.01) ks_mode <- 0.10
  if(is.na(kp_mode) || kp_mode < 0.01) kp_mode <- 1.25

  cat(sprintf("Wavelet Ratio (V_p/V_s): %.4f | Priors: k_s=%.3f, k_p=%.3f\n",
              energy_ratio, ks_mode, kp_mode))

  list(
    nu0_init = mean(hr_hz),
    energy_ratio = energy_ratio,
    ks_mode = ks_mode,
    kp_mode = kp_mode
  )
}

# ==============================================================================
# 2. OU INTEGRALS HELPER (Reparameterized for Variance/Behavior Orthogonality)
# ==============================================================================
.compute_ou_block <- function(k, dt, V) {
  x <- k * dt

  # Dynamic diffusion calculation to conserve Amplitude (V)
  sigma2 <- 2 * k * V

  if (x < 1e-3) {
    phi_int <- dt - (k * dt^2)/2 + (k^2 * dt^3)/6
    q11 <- sigma2 * (dt - k * dt^2 + (2 * k^2 * dt^3)/3)
    q13 <- sigma2 * ((dt^2)/2 - (k * dt^3)/2 + (7 * k^2 * dt^4)/24)
    q33 <- sigma2 * ((dt^3)/3 - (k * dt^4)/4 + (7 * k^2 * dt^5)/60)
  } else {
    phi_int <- (1 - exp(-x)) / k
    q11 <- sigma2 * (1 - exp(-2 * x)) / (2 * k)
    q13 <- sigma2 * (1 - 2*exp(-x) + exp(-2*x)) / (2 * k^2)
    q33 <- sigma2 * (2*x - 3 + 4*exp(-x) - exp(-2*x)) / (2 * k^3)
  }

  list(phi_int = phi_int, q11 = q11, q13 = q13, q33 = q33)
}

# ==============================================================================
# 3. DUAL-FILTER GLS (Robustified with Asymptotic Null-Gain Coasting + RTS)
# ==============================================================================
run_gls_filter <- function(dy, kp, ks, lp, lR, jump_threshold, jump_power, smooth) {
  N <- length(dy)

  # ASYMPTOTIC NULL-GAIN GATE (Numerically Safe)
  dy_safe <- pmax(dy[-N], 1e-6)

  jump_back <- c(0, abs(diff(dy)) / dy_safe)
  jump_fwd  <- c(abs(diff(dy)) / dy_safe, 0)

  J_k <- pmax(jump_back, jump_fwd, na.rm = TRUE)
  gate_multiplier <- 1 + (J_k / jump_threshold)^jump_power

  # Cap the multiplier at 10 billion. This mathematically severs the Kalman
  # update but prevents R from evaluating (0 * Inf) to NaN in the covariance step.
  gate_multiplier[!is.finite(gate_multiplier)] <- 1e10
  gate_multiplier <- pmin(gate_multiplier, 1e10)

  X0 <- matrix(0, 4, 1); Xv <- matrix(0, 4, 1)
  P  <- diag(c(1, 1, 0, 0)); H  <- matrix(c(0, 0, -1, 1), 1, 4)
  S_vec <- numeric(N); v0_vec <- numeric(N); vv_vec <- numeric(N)
  S_ungated_vec <- numeric(N)

  # Storage for forward state reconstruction
  X0_store <- matrix(0, N, 4); Xv_store <- matrix(0, N, 4)

  # +++ RTS SMOOTHER CONDITIONAL STORAGE +++
  if (smooth) {
    if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
    X0_pred_store <- matrix(0, N, 4)
    Xv_pred_store <- matrix(0, N, 4)
    P_pred_store  <- array(0, dim=c(4, 4, N))
    P_post_store  <- array(0, dim=c(4, 4, N))
    Phi_store     <- array(0, dim=c(4, 4, N))
  }
  # ++++++++++++++++++++++++++++++++++++++++

  Phi <- matrix(0, 4, 4); Q <- matrix(0, 4, 4)

  for(k in 1:N) {
    dt <- dy[k]
    dynamic_lR <- min(lR * gate_multiplier[k], 1e12, na.rm = TRUE)

    blk_p <- .compute_ou_block(kp, dt, lp)
    blk_s <- .compute_ou_block(ks, dt, 1.0)
    Phi[1,1] <- exp(-kp * dt); Phi[3,1] <- blk_p$phi_int
    Phi[2,2] <- exp(-ks * dt); Phi[4,2] <- blk_s$phi_int
    Q[1,1] <- blk_p$q11; Q[1,3] <- blk_p$q13; Q[3,1] <- blk_p$q13; Q[3,3] <- blk_p$q33
    Q[2,2] <- blk_s$q11; Q[2,4] <- blk_s$q13; Q[4,2] <- blk_s$q13; Q[4,4] <- blk_s$q33

    X0_pred <- Phi %*% X0; Xv_pred <- Phi %*% Xv
    P_pred  <- Phi %*% P %*% t(Phi) + Q

    S <- as.numeric(H %*% P_pred %*% t(H) + dynamic_lR * dt)
    S <- max(S, 1e-12)
    K <- P_pred %*% t(H) / S

    S_true <- as.numeric(H %*% P_pred %*% t(H) + lR * dt)
    S_ungated_vec[k] <- max(S_true, 1e-12)

    v0 <- 1.0 - as.numeric(H %*% X0_pred)
    vv <- dt  - as.numeric(H %*% Xv_pred)

    X0 <- X0_pred + K * v0; Xv <- Xv_pred + K * vv
    I_KH <- diag(4) - K %*% H

    P <- I_KH %*% P_pred %*% t(I_KH) + K %*% (dynamic_lR * dt) %*% t(K)
    P <- (P + t(P)) / 2

    X0_store[k, ] <- as.numeric(X0)
    Xv_store[k, ] <- as.numeric(Xv)

    # +++ RTS SMOOTHER CONTINUOUS LOGGING +++
    if (smooth) {
      X0_pred_store[k, ] <- as.numeric(X0_pred)
      Xv_pred_store[k, ] <- as.numeric(Xv_pred)
      P_pred_store[,,k]  <- P_pred
      P_post_store[,,k]  <- P
      Phi_store[,,k]     <- Phi
    }
    # +++++++++++++++++++++++++++++++++++++++

    # BOUNDARY RESET
    X0[3:4,] <- 0; Xv[3:4,] <- 0
    P[3:4, ] <- 0; P[, 3:4] <- 0

    S_vec[k] <- S; v0_vec[k] <- v0; vv_vec[k] <- vv
  }

  # Analytical Marginalizations (Underflow Guarded)
  den_nu <- sum((vv_vec^2) / S_vec, na.rm = TRUE)
  nu0_hat <- if (is.finite(den_nu) && den_nu > 1e-12) sum(v0_vec * vv_vec / S_vec, na.rm = TRUE) / den_nu else 0

  v_final  <- v0_vec - nu0_hat * vv_vec
  sig2_hat <- max(mean(v_final^2 / S_vec), 1e-12)
  LL_star  <- -0.5 * sum(log(S_vec)) - (N / 2) * log(sig2_hat)

  # Reconstruct forward biological states
  states_final <- (X0_store - nu0_hat * Xv_store)
  std_innov <- v_final / sqrt(S_vec * sig2_hat)

  std_innov_diagnostic <- v_final / sqrt(S_ungated_vec * sig2_hat)
  implied_rr <- (1.0 + states_final[,3] - states_final[,4]) / nu0_hat
  time_residuals <- dy - implied_rr

  # RTS BACKWARD PASS
  if (smooth) {
    # Superimpose predictions utilizing the profiled baseline
    X_pred_combined <- X0_pred_store - nu0_hat * Xv_pred_store
    X_smooth <- matrix(0, N, 4)

    # Initialize the backward recursion with the final forward state
    X_smooth[N, ] <- states_final[N, ]

    for (k in (N-1):1) {
      # 1. Safe Inversion Fallback
      P_pred_inv <- tryCatch({
        MASS::ginv(P_pred_store[,, k+1])
      }, error = function(e) {
        # If SVD fails due to extreme scaling or an isolated NaN, return a zero matrix
        matrix(0, 4, 4)
      })

      # 2. Backward Gain
      J_k <- P_post_store[,, k] %*% t(Phi_store[,, k+1]) %*% P_pred_inv

      # 3. Guard against NA propagation from the forward states
      state_diff <- matrix(X_smooth[k+1, ] - X_pred_combined[k+1, ], 4, 1)
      if (any(!is.finite(state_diff)) || any(!is.finite(J_k))) {
        J_k <- matrix(0, 4, 4)
      }

      # 4. Smooth the state (If J_k is 0, it safely accepts the forward filtered state)
      X_smooth[k, ] <- states_final[k, ] + as.numeric(J_k %*% state_diff)
    }

    # Append to output array
    out_states <- cbind(states_final, X_smooth[, 1:2])
    colnames(out_states) <- c("X_p", "X_s", "Phase_p", "Phase_s", "X_p_smooth", "X_s_smooth")

    return(list(LL_star = LL_star, nu0 = nu0_hat, sig2 = sig2_hat,
                states = out_states, std_innov = std_innov,
                innovations_diagnostic = std_innov_diagnostic,
                time_residuals = time_residuals,
                gate_multiplier = gate_multiplier))
  }
  # ++++++++++++++++++++++++++++++++++++++++

  list(LL_star = LL_star, nu0 = nu0_hat, sig2 = sig2_hat,
       states = states_final, std_innov = std_innov,
       innovations_diagnostic = std_innov_diagnostic,
       time_residuals = time_residuals,
       gate_multiplier = gate_multiplier)
}

# ==============================================================================
# MAP OBJECTIVE FUNCTION
# ==============================================================================
objective_map_3d <- function(theta, dy, priors, jump_threshold, jump_power, smooth) {
  ks <- exp(theta[1])
  kp <- ks + exp(theta[2])
  lR <- exp(theta[3])

  # Vp is the pure amplitude ratio. NO longer entangled with (kp / ks)!
  Vp <- priors$energy_ratio

  # Pass hyperparameters to the filter (Note: signature of run_gls_filter does not change)
  res <- run_gls_filter(dy, kp, ks, Vp, lR, jump_threshold, jump_power, smooth)

  # Data-Driven Maximum A Posteriori (MAP) Priors mapped from Wavelet Bandwidth
  ks_mode <- priors$ks_mode; ks_strength <- 3.0
  pen_ks  <- (ks_strength) * log(ks) - (ks_strength/ks_mode) * ks

  kp_mode <- priors$kp_mode; kp_strength <- 5.0
  pen_kp  <- (kp_strength) * log(kp) - (kp_strength/kp_mode) * kp

  lR_mode <- 0.0005; lR_shape <- 2.0
  pen_lR  <- -(lR_shape + 1) * log(lR) - (lR_mode * (lR_shape + 1)) / lR

  return(-(res$LL_star + pen_ks + pen_kp + pen_lR))
}

# ==============================================================================
# 4. OPTIMIZATION ROUTINE
# ==============================================================================
fit_model <- function(dy, jump_threshold = 0.1, jump_power = 10, smooth = TRUE) {
  if (any(is.na(dy))) dy <- na.omit(dy)
  priors <- extract_spectral_priors(dy)

  # Initialize near the empirical bandwidths
  theta_init <- c(log(max(priors$ks_mode, 0.01)), log(max(priors$kp_mode - priors$ks_mode, 0.05)), log(0.001))

  cat("\nRunning 3D BFGS Optimization with Reparameterized MAP Penalties...\n")
  opt_res <- optim(par = theta_init, fn = objective_map_3d, dy = dy,
                   priors = priors,
                   jump_threshold = jump_threshold,
                   jump_power = jump_power,
                   smooth = smooth,
                   method = "BFGS",
                   control = list(maxit = 1000, trace = 1))

  ks_opt <- exp(opt_res$par[1])
  kp_opt <- ks_opt + exp(opt_res$par[2])
  lR_opt <- exp(opt_res$par[3])
  Vp_opt <- priors$energy_ratio

  final_res <- run_gls_filter(dy, kp_opt, ks_opt, Vp_opt, lR_opt, jump_threshold, jump_power, smooth)

  return(list(
    dy = dy, time = cumsum(dy),
    priors_obj = priors, # Required for the reporting Hessian
    params = list(
      ks = ks_opt, kp = kp_opt, lR = lR_opt, lp = Vp_opt,
      nu0 = final_res$nu0, sig2 = final_res$sig2,
      # Convert the Vp amplitude back into physical SDE diffusion (sigma) for the report
      sig_p = sqrt(2 * kp_opt * Vp_opt * final_res$sig2),
      sig_s = sqrt(2 * ks_opt * 1.0 * final_res$sig2),
      energy_ratio = priors$energy_ratio,
      jump_threshold = jump_threshold, jump_power = jump_power
    ),
    states = final_res$states,
    innovations = final_res$std_innov,
    innovations_diagnostic = final_res$innovations_diagnostic,
    time_residuals = final_res$time_residuals,
    gate_multiplier = final_res$gate_multiplier,
    opt_raw = opt_res
  ))
}

# ==============================================================================
# 5. REPORT PARAMETERS
# ==============================================================================
report_model <- function(fit_obj, smooth = TRUE) {
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
    priors = fit_obj$priors_obj,
    jump_threshold = p$jump_threshold,
    jump_power = p$jump_power,
    smooth = smooth
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
    list(name = "Lambda_R", group = "Ratios", interp = "Fraction of variance from SA node threshold jitter", status = "Estimated (MAP)", est = p$lR, se = se_lR),

    # +++ ADD HYPERPARAMETERS TO REPORT +++
    list(name = "Theta_Jump", group = "Robustness", interp = "Relative fractional jump threshold for outlier rejection", status = "Fixed (Hyper)", est = p$jump_threshold, se = NA),
    list(name = "Power_Jump", group = "Robustness", interp = "Polynomial steepness of asymptotic null-gain gate", status = "Fixed (Hyper)", est = p$jump_power, se = NA)
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
# 6. VISUALIZATION (Phase Reconstruction & Topology)
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
    geom_line(aes(y = Observed, color = "Observed RR"), linewidth = 1/3, alpha = 0.5) +
    geom_line(aes(y = Implied, color = "Filtered RR"), linewidth = 1/3) +
    scale_color_manual(values = c("Observed RR" = "gray50",
                                  "Filtered RR" = "darkred")) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "Phase-Domain Filtering & Predictive Fit",
         y = "RR Interval (s)",
         x = "") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "top")

  df_X <- data.frame(
    Time = rep(time_cum, 2), Value = c(X_p, X_s),
    State = factor(rep(c("Parasympathetic Drive", "Sympathetic Drive"), each = N),
                   levels = c("Parasympathetic Drive", "Sympathetic Drive"))
  )

  pB <- ggplot(df_X, aes(x = Time, y = Value, color = State)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "black", alpha = 0.5) +
    geom_line(linewidth = 1/3) +
    scale_color_manual(values = c("Parasympathetic Drive" = "#4DB0D0",
                                  "Sympathetic Drive" = "#DC5050")) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "Native Autonomic Drivers (Direct Tracking)",
         y = "Amplitude (Hz)",
         x = "Time (seconds)") +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "top")

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

  pC <- ggplot() +
    geom_raster(data = grid, aes(x = X_p, y = X_s, fill = Energy)) +
    geom_contour(data = grid, aes(x = X_p, y = X_s, z = Energy), color = "white", alpha = 0.2, bins = 20) +
    geom_path(data = df_traj, aes(x = X_p, y = X_s, color = Time, alpha = Time),
              arrow = ggplot2::arrow(type = "closed", length = unit(0.06, "inches")),
              show.legend = FALSE) +
    scale_fill_viridis_c(option = "mako", name = "Energy (U)") +
    scale_color_viridis_c(option = "plasma", guide = "none") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(title = "Autonomic Phase Space Topology",
         x = "Parasympathetic Drive (Hz)",
         y = "Sympathetic Drive (Hz)") +
    theme_classic() +
    theme(legend.position = "top")

  (pA / pB) | pC
}

# ==============================================================================
# 7. DIAGNOSTICS (Pre-standardized innovations)
# ==============================================================================
diagnose_model <- function(fit_obj) {
  z_raw <- fit_obj$innovations_diagnostic
  gate <- fit_obj$gate_multiplier

  # Only run diagnostics on physiologically valid (non-gated) beats
  valid_idx <- which(gate < 1.5)
  z <- z_raw[valid_idx]
  U_k <- pnorm(z)
  N <- length(U_k)

  ks_res <- ks.test(U_k, "punif")
  print(ks_res)
  ks_stat <- ks_res$statistic
  p_val <- ks_res$p.value
  bound <- 1.36 / sqrt(N)

  df_qq <- data.frame(z = z)
  pA <- ggplot(df_qq, aes(sample = z)) +
    stat_qq(color = "darkblue", alpha = 0.5) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = "Q-Q Plot of Phase Innovations",
         subtitle = "Tests Linear Observation Gaussianity",
         x = "Theoretical Normal",
         y = "Empirical") +
    theme_classic()

  empirical_cdf <- ecdf(U_k)
  df_cdf <- data.frame(U = sort(U_k), CDF = empirical_cdf(sort(U_k)))

  pB <- ggplot(df_cdf, aes(x = U, y = CDF)) +
    geom_step(color = "darkgreen", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_abline(slope = 1, intercept = bound, color = "gray", linetype = "dotted", linewidth=0.8) +
    geom_abline(slope = 1, intercept = -bound, color = "gray", linetype = "dotted", linewidth=0.8) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
    labs(title = "Time-Rescaling KS-Plot",
         subtitle = "Validates Exact IPFM Boundary Crossings",
         x = "Theoretical Uniform(0,1)",
         y = "Empirical CDF") +
    theme_classic()

  acf_data <- acf(z, plot = FALSE, lag.max = 40)
  df_acf <- data.frame(lag = acf_data$lag, acf = acf_data$acf)

  pC <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0),
                 color = "black",
                 linewidth = 0.8) +
    geom_hline(yintercept = c(-1.96/sqrt(N), 1.96/sqrt(N)),
               color = "red",
               linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black") +
    labs(title = "Autocorrelation of Innovations",
         subtitle = "Tests for Unmodeled Kinetics",
         x = "Lag (Heartbeats)",
         y = "ACF") +
    theme_classic()

  pD <- ggplot(df_qq, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = "#7E6148",
                   alpha = 0.7, bins = 30, color = "white") +
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = 1),
                  color = "black",
                  linetype = "dashed") +
    labs(title = "Exact Gaussian Innovations",
         subtitle = "Validating Phase-Domain Linearity",
         x = "Standardized Residuals",
         y = "Density") +
    theme_classic()

  (pA | pB) / (pC | pD)
}

# ==============================================================================
# 8. BATCH PROCESSING & AGGREGATION ENGINE
# ==============================================================================
batch_process_hrv <- function(rr_list, jump_threshold = 0.15, jump_power = 8, dataset_tag = "Unknown") {
  N_batch <- length(rr_list)
  subj_names <- names(rr_list)
  if (is.null(subj_names)) subj_names <- paste0("Subject_", seq_len(N_batch))

  # Storage structures
  summary_rows <- list()
  states_list <- list()
  raw_fits <- list()

  cat(sprintf("\nInitializing Batch Processing for Dataset: %s (%d recordings)\n", dataset_tag, N_batch))

  for (i in seq_len(N_batch)) {
    dy <- rr_list[[i]]
    s_name <- subj_names[i]
    cat(sprintf(" Processing [%d/%d]: %s... ", i, N_batch, s_name))

    # Safe execution wrapper
    fit_res <- tryCatch({
      fit_model(dy, jump_threshold = jump_threshold, jump_power = jump_power)
    }, error = function(e) {
      cat("FAILED\n")
      return(NULL)
    })

    if (is.null(fit_res)) next
    cat("SUCCESS\n")

    p <- fit_res$params
    raw_fits[[s_name]] <- fit_res

    # Extract Statistical Diagnostics natively matching diagnose_model()
    z <- fit_res$innovations
    N_beats <- length(z)

    # 1. Kolmogorov-Smirnov Test (Time-Rescaling Theorem)
    U_k <- pnorm(z)
    ks_res <- ks.test(U_k, "punif")

    # 2. Autocorrelation Check (Up to 40 lags)
    acf_res <- acf(z, lag.max = 40, plot = FALSE)
    acf_vals <- acf_res$acf[-1] # Remove lag 0 (always 1.0)
    acf_bound <- 1.96 / sqrt(N_beats)
    significant_lags <- sum(abs(acf_vals) > acf_bound)

    # Pack Scalar Outputs
    summary_rows[[s_name]] <- data.frame(
      Dataset = dataset_tag,
      SubjectID = s_name,
      Converged = (fit_res$opt_raw$convergence == 0),
      Iterations = fit_res$opt_raw$counts[1],
      NumBeats = N_beats,
      Nu_0 = p$nu0,
      Sigma_sys2 = p$sig2,
      Kappa_S = p$ks,
      Kappa_P = p$kp,
      Tau_S = 1 / p$ks,
      Tau_P = 1 / p$kp,
      Lambda_R = p$lR,
      KS_Stat = as.numeric(ks_res$statistic),
      KS_p_value = as.numeric(ks_res$p.value),
      ACF_Violations = significant_lags,
      stringsAsFactors = FALSE
    )

    # Pack Time-Varying States for Long-Form Pooling
    states_list[[s_name]] <- data.frame(
      Dataset = dataset_tag,
      SubjectID = s_name,
      Beat = seq_len(N_beats),
      Time = fit_res$time,
      RR_Interval = dy,
      X_p = fit_res$states[, 1],
      X_s = fit_res$states[, 2],
      Phase_p = fit_res$states[, 3],
      Phase_s = fit_res$states[, 4],
      Innovation = z,
      stringsAsFactors = FALSE
    )
  }

  # Bind into tidy metrics structures
  summary_df <- do.call(rbind, summary_rows)
  states_df  <- do.call(rbind, states_list)
  rownames(summary_df) <- NULL
  rownames(states_df)  <- NULL

  return(list(
    summary = summary_df,
    states = states_df,
    models = raw_fits
  ))
}

# ==============================================================================
# 9. GROUP METRICS REPORTING UTILITY
# ==============================================================================
report_batch_metrics <- function(batch_obj) {
  df <- batch_obj$summary
  N <- nrow(df)

  cat("\n=========================================================================\n")
  cat(sprintf("               BATCH VALIDATION REPORT (Dataset: %s)\n", df$Dataset[1]))
  cat("=========================================================================\n")
  cat(sprintf("Total Managed Series      : %d\n", N))
  cat(sprintf("Optimization Convergence  : %.1f%%\n", (sum(df$Converged) / N) * 100))
  cat(sprintf("Average Sequence Length   : %.1f beats\n", mean(df$NumBeats)))

  cat("\n--- Statistical Diagnostics (Point-Process Validity) ---\n")
  alpha_passing <- sum(df$KS_p_value > 0.05)
  cat(sprintf("Time-Rescaling KS Pass Rate (p > 0.05)   : %.1f%% (Mean Stat: %.4f)\n",
              (alpha_passing / N) * 100, mean(df$KS_Stat)))
  cat(sprintf("Mean Whiteness Violations (out of 40)    : %.2f lags\n",
              mean(df$ACF_Violations)))

  cat("\n--- Kinetic Parameter Estimations (Mean ± SD) ---\n")
  cat(sprintf("Baseline Pacing (Nu_0)         : %.3f ± %.3f Hz\n", mean(df$Nu_0), sd(df$Nu_0)))
  cat(sprintf("Sympathetic Clearance (Kappa_S): %.4f ± %.4f Hz (Tau: %.2fs)\n",
              mean(df$Kappa_S), sd(df$Kappa_S), mean(df$Tau_S)))
  cat(sprintf("Parasympathetic Clear. (Kappa_P): %.4f ± %.4f Hz (Tau: %.2fs)\n",
              mean(df$Kappa_P), sd(df$Kappa_P), mean(df$Tau_P)))
  cat(sprintf("Observation Jitter (Lambda_R)  : %.5f ± %.5f\n", mean(df$Lambda_R), sd(df$Lambda_R)))
  cat("=========================================================================\n\n")
}
