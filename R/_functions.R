#' Run Kalman Filter with Exact Van Loan Discretization and Joseph Stabilized Update
#'
#' @param theta Numeric vector of 7 transformed parameters.
#' @param log_rr_ts Numeric vector of log-transformed RR intervals.
#' @param dt Time step (e.g., 1 for beat-to-beat).
#' @param return_states Boolean; if TRUE, returns state variables and innovations.
run_kalman_filter <- function(theta, log_rr_ts, dt = 1, return_states = FALSE) {

  # ==============================================================================
  # 1. Autonomic Kinetics (7-Parameter Architecture)
  # ==============================================================================
  tau_v_min <- 0.3; tau_v_max <- 2.0
  rho_min <- 5.0; rho_max <- 15.0

  # Hyperspherical mapping replaces plogis
  tau_v <- tau_v_min + (tau_v_max - tau_v_min) * (sin(theta[1])^2)

  log_rho <- log(rho_min) + (log(rho_max) - log(rho_min)) * (sin(theta[2])^2)
  tau_s <- tau_v * exp(log_rho)

  a11 <- 1 / tau_v
  a22 <- 1 / tau_s
  max_a12 <- sqrt(a11 * a22)
  a12 <- max_a12 * (sin(theta[3])^2)

  sigma_stat_v <- exp(theta[4])
  sigma_stat_s <- exp(theta[5])

  # 2. Intrinsic SA Node Memory (Strictly Positive AR1)
  k1 <- (sin(theta[6])^2)

  # AR Variance Bounding: Parameterize stationary variance, derive driving noise conditionally
  sigma_stat_eta <- exp(theta[7])
  sig_eta <- sigma_stat_eta * sqrt(1 - k1^2)

  # Removing mu_rr from optimization
  mu_rr <- mean(log_rr_ts)
  N <- length(log_rr_ts)

  # ==============================================================================
  # 3. EXACT DISCRETIZATION VIA VAN LOAN'S METHOD
  # ==============================================================================
  A_2x2 <- matrix(c(-a11, -a12, -a12, -a22), nrow = 2, byrow = TRUE)

  # Stationary Covariance Reparameterization: Solve continuous Lyapunov algebraically
  Sigma_stat <- matrix(c(sigma_stat_v^2, 0, 0, sigma_stat_s^2), nrow = 2, ncol = 2)
  Sigma_c <- -(A_2x2 %*% Sigma_stat + Sigma_stat %*% t(A_2x2))

  if (any(!is.finite(A_2x2)) || any(!is.finite(Sigma_c))) {
    if (return_states) stop("Filter diverged: Infinite matrix entries explored.")
    return(1e9)
  }

  # PSD Check for algebraically derived continuous noise
  if (any(eigen(Sigma_c, symmetric = TRUE)$values < 0)) {
    if (return_states) stop("Filter diverged: Continuous process noise not PSD.")
    return(1e9)
  }

  Z <- rbind(
    cbind(-A_2x2, Sigma_c),
    cbind(matrix(0, 2, 2), t(A_2x2))
  )

  expZ <- expm::expm(Z * dt)
  F_2x2 <- t(expZ[3:4, 3:4])
  Q_2x2 <- F_2x2 %*% expZ[1:2, 3:4]
  Q_2x2 <- (Q_2x2 + t(Q_2x2)) / 2

  # ==============================================================================
  # 3-DIMENSIONAL STATE MATRICES SETUP
  # ==============================================================================
  F_mat <- matrix(0, nrow = 3, ncol = 3)
  F_mat[1:2, 1:2] <- F_2x2
  F_mat[3, 3] <- k1

  Q <- matrix(0, nrow = 3, ncol = 3)
  Q[1:2, 1:2] <- Q_2x2
  Q[3, 3] <- sig_eta^2

  # Observation maps directly to [Vagal, Sympathetic, SA_Memory]
  H <- matrix(c(1, -1, 1), nrow = 1)

  R <- 1e-5
  R_mat <- matrix(R, 1, 1)

  x_hat <- matrix(0, nrow = 3, ncol = N)
  P_var <- matrix(0, nrow = 3, ncol = N)
  innovations <- numeric(N)

  # ==============================================================================
  # 4. EXACT STATIONARY INITIALIZATION (Discrete Lyapunov Equation)
  # ==============================================================================
  I_9 <- diag(9)
  F_kron_F <- kronecker(F_mat, F_mat)
  vec_Q <- as.vector(Q)

  P_init <- tryCatch({
    matrix(solve(I_9 - F_kron_F, vec_Q), 3, 3)
  }, error = function(e) {
    matrix(MASS::ginv(I_9 - F_kron_F) %*% vec_Q, 3, 3)
  })

  # Enforce strict symmetry and Positive Definiteness
  P_init <- (P_init + t(P_init)) / 2
  eigs <- eigen(P_init, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs < 1e-10)) {
    P_init <- as.matrix(Matrix::nearPD(P_init, ensureSymmetry = TRUE)$mat)
  }

  x_pred <- matrix(0, nrow = 3, ncol = 1)
  P_pred <- P_init

  log_likelihood <- 0
  I_3 <- diag(3)

  # ==============================================================================
  # 5. KALMAN FILTER LOOP
  # ==============================================================================
  for (t in 1:N) {

    y_pred <- (H %*% x_pred)[1,1] + mu_rr
    y_obs <- log_rr_ts[t]

    y_err <- y_obs - y_pred
    innovations[t] <- y_err

    S <- (H %*% P_pred %*% t(H))[1,1] + R

    if (!is.finite(S) || S < 1e-9) {
      if (return_states) stop("Filter diverged.")
      return(1e9)
    }

    K <- P_pred %*% t(H) / S
    x_hat[, t] <- x_pred + K * y_err

    IKH <- I_3 - K %*% H
    P <- IKH %*% P_pred %*% t(IKH) + K %*% R_mat %*% t(K)
    P <- (P + t(P)) / 2

    P_var[, t] <- diag(P)

    log_likelihood <- log_likelihood - 0.5 * (log(2 * pi * S) + (y_err^2) / S)

    if (t < N) {
      x_pred <- F_mat %*% x_hat[, t]
      P_pred <- F_mat %*% P %*% t(F_mat) + Q
    }
  }

  nll <- -log_likelihood
  if (!is.finite(nll)) return(1e9)

  if (return_states) {
    return(list(x_hat = x_hat, P_var = P_var, innovations = innovations))
  } else {
    return(nll)
  }
}

#' L2 Penalized Negative Log-Likelihood Objective Function
#'
#' Evaluates the exact linear Gaussian state-space model and applies
#' a smooth L2 (Ridge) regularization penalty to targeted auxiliary parameters.
#'
#' @param theta Numeric vector of 7 unconstrained optimization parameters.
#' @param log_rr_ts Numeric vector of the log-transformed RR interval time series.
#' @param lambda Scalar hyperparameter dictating the severity of the L2 shrinkage.
#' @return The L2 penalized negative log-likelihood scalar.
evaluate_penalized_nll <- function(theta, log_rr_ts, lambda = 0) {
  base_nll <- run_kalman_filter(theta, log_rr_ts, dt = 1, return_states = FALSE)

  if (is.infinite(base_nll) || base_nll == 1e9) return(1e9)

  # Direct quadratic penalty on the unbounded parameter space
  l2_penalty <- lambda * (theta[3]^2)

  return(base_nll + l2_penalty)
}

#' Fit Adaptive Autonomic State-Space Model (Deterministic Initialization)
#'
#' Executes a single deep numerical minimization on the penalized log-likelihood surface
#' using physiologically anchored starting coordinates, bypassing global LHS exploration.
#'
#' @param rr_ts Numeric vector of the raw RR interval time series.
#' @param lambda Scalar hyperparameter dictating the L2 (Ridge) shrinkage severity.
#' @return A fitted model object containing the optimal parameter vector and the exact Hessian.
fit_autonomic_model <- function(rr_ts, lambda = 0) {
  log_rr_ts <- log(rr_ts)

  # 7-Parameter Deterministic Anchoring (Arcsine Square Root transformations)
  init_theta1 <- asin(sqrt((1.0 - 0.3) / (2.0 - 0.3)))
  init_theta2 <- asin(sqrt((log(10.0) - log(5.0)) / (log(15.0) - log(5.0)))) # Anchored at rho = 10
  init_theta3 <- asin(sqrt(0.10))
  init_theta4 <- log(0.05)
  init_theta5 <- log(0.05)
  init_theta6 <- asin(sqrt(0.50))  # Positive SA memory baseline
  init_theta7 <- log(0.01)

  theta_init <- c(init_theta1, init_theta2, init_theta3, init_theta4,
                  init_theta5, init_theta6, init_theta7)

  fit_try <- tryCatch({
    optim(
      par = theta_init,
      fn = evaluate_penalized_nll,
      log_rr_ts = log_rr_ts,
      lambda = lambda,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = 10000, trace = 1)
    )
  }, error = function(e) { return(NULL) })

  if (is.null(fit_try) || is.infinite(fit_try$value)) {
    return(list(par = rep(NA, 7), value = Inf, hessian = matrix(NA, 7, 7), convergence = FALSE))
  }

  return(list(par = fit_try$par, value = fit_try$value,
              hessian = fit_try$hessian, convergence = (fit_try$convergence == 0)))
}

get_model_ci <- function(fit_result) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Please install 'numDeriv' for Delta Method standard errors.")

  cov_mat_raw <- tryCatch(MASS::ginv(fit_result$hessian), error = function(e) matrix(NA, 7, 7))
  cov_mat_pd <- tryCatch(as.matrix(Matrix::nearPD(cov_mat_raw, ensureSymmetry = TRUE)$mat),
                         error = function(e) cov_mat_raw)

  # Define the transformation mapping function identical to the exact model boundaries
  transform_theta <- function(th) {
    t_v <- 0.3 + 1.7 * (sin(th[1])^2)
    rho <- exp(log(5.0) + (log(15.0) - log(5.0)) * (sin(th[2])^2))
    t_s <- t_v * rho
    a12_f <- sin(th[3])^2
    s_v <- exp(th[4])
    s_s <- exp(th[5])
    k1_val <- sin(th[6])^2
    s_e <- exp(th[7])
    return(c(t_v, t_s, a12_f, s_v, s_s, k1_val, s_e))
  }

  # Project unbounded point estimates directly into physiological space
  phys_est <- transform_theta(fit_result$par)

  # Delta Method Execution: J * Sigma * J^T
  J <- tryCatch(numDeriv::jacobian(transform_theta, fit_result$par), error = function(e) matrix(NA, 7, 7))
  phys_cov <- J %*% cov_mat_pd %*% t(J)

  # Extract valid projected standard errors
  phys_se <- sqrt(pmax(diag(phys_cov), 0))

  ci_lower <- phys_est - 1.96 * phys_se
  ci_upper <- phys_est + 1.96 * phys_se

  param_names <- c("tau_v", "tau_s", "a12_fraction", "sigma_stat_v", "sigma_stat_s", "k1", "sigma_stat_eta")

  return(data.frame(Parameter = param_names, Estimate = phys_est,
                    Std_Error = phys_se, CI_Lower = ci_lower, CI_Upper = ci_upper))
}

#' Parse and Untransform Optimized Parameters into Physical Units
report_fit <- function(fit_result) {
  if (all(is.na(fit_result$par))) return(data.frame(Message = "Optimization Failed"))

  # get_model_ci now contains Delta-mapped physiological estimates and CIs
  phys_ci <- get_model_ci(fit_result)

  interpretations <- c(
    "Vagal Recovery (beats)", "Sympathetic Recovery (beats)", "Cross-talk Coupling Fraction",
    "Vagal Stationary Amplitude", "Sympathetic Stationary Amplitude",
    "SA Node Memory (AR1)", "SA Node Memory Variance"
  )

  return(data.frame(
    Parameter_Symbol = phys_ci$Parameter,
    Interpretation = interpretations,
    Estimate = round(phys_ci$Estimate, 4),
    Std_Error = round(phys_ci$Std_Error, 4),
    CI_Lower = round(phys_ci$CI_Lower, 4),
    CI_Upper = round(phys_ci$CI_Upper, 4)
  ))
}


#' Advanced Autonomic Model Diagnostics and 3D Visualization
#'
#' @param fit_result The optimized model object returned by fit_autonomic_model()
#' @param rr_ts The numeric vector of the original observed RR intervals
#' @return A single patchwork object containing the complete 5-panel dashboard.
visualize_autonomic_model <- function(fit_result, rr_ts) {

  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install 'patchwork'")
  require(ggplot2)
  require(patchwork)

  log_rr_ts <- log(rr_ts)
  opt_theta <- fit_result$par

  # Extract System Kinetics using Log-Ratio mapping
  tau_v_min <- 0.3; tau_v_max <- 2.0
  rho_min   <- 5.0; rho_max   <- 15.0

  tau_v   <- tau_v_min + (tau_v_max - tau_v_min) * (sin(opt_theta[1])^2)
  log_rho <- log(rho_min) + (log(rho_max) - log(rho_min)) * (sin(opt_theta[2])^2)
  tau_s   <- tau_v * exp(log_rho)
  a11     <- 1 / tau_v
  a22     <- 1 / tau_s
  a12     <- sqrt(a11 * a22) * (sin(opt_theta[3])^2)
  mu_rr_log <- mean(log_rr_ts)

  # Run Exact Kalman Filter
  kf_out <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)
  states <- kf_out$x_hat
  p_var <- kf_out$P_var
  innovations <- kf_out$innovations

  time_steps <- 1:length(rr_ts)
  valid_idx <- 2:length(rr_ts)

  # =========================================================================
  # SIGNAL RECONSTRUCTION & PREDICTION
  # =========================================================================

  log_rr_pred <- log_rr_ts - innovations
  rr_pred <- exp(log_rr_pred)

  H <- matrix(c(1, -1, 1), nrow = 1)
  log_rr_rec <- as.numeric(H %*% states) + mu_rr_log
  rr_rec <- exp(log_rr_rec)

  # =========================================================================
  # GGPLOT2 COMPONENTS
  # =========================================================================

  df_ts <- data.frame(Time = time_steps, Observed = rr_ts, Predicted = rr_pred)
  p_ts <- ggplot(df_ts) +
    geom_line(aes(x = Time, y = Observed, color = "Observed"), linewidth = 0.6, alpha = 0.5) +
    geom_line(aes(x = Time, y = Predicted, color = "Kalman Predicted"), linewidth = 0.4, alpha = 0.9) +
    scale_color_manual(values = c("Observed" = "gray50", "Kalman Predicted" = "#E74C3C")) +
    theme_classic(base_size = 11) +
    labs(title = "A. Signal Prediction Tracking", x = "Beat Sequence", y = "RR Interval (ms)") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"))

  r_squared <- round(cor(rr_ts[valid_idx], rr_pred[valid_idx])^2, 3)
  df_scatter <- data.frame(Obs = rr_ts[valid_idx], Pred = rr_pred[valid_idx])
  p_scatter <- ggplot(df_scatter, aes(x = Obs, y = Pred)) +
    geom_point(color = "#34495E", alpha = 0.3, pch = 16) +
    geom_abline(slope = 1, intercept = 0, color = "#E74C3C") +
    annotate("text", x = min(df_scatter$Obs), y = max(df_scatter$Pred),
             label = paste0("R² = ", r_squared), hjust = 0, vjust = 1, fontface = "bold") +
    theme_classic(base_size = 11) +
    labs(title = "B. Prediction Accuracy", x = "Observed RR (ms)", y = "Predicted RR (ms)") +
    theme(plot.title = element_text(face = "bold"))

  df_states <- data.frame(
    Time = rep(time_steps, 2),
    Activity = c(states[1, ], states[2, ]),
    Lower = c(states[1, ] - 1.96*sqrt(p_var[1,]), states[2, ] - 1.96*sqrt(p_var[2,])),
    Upper = c(states[1, ] + 1.96*sqrt(p_var[1,]), states[2, ] + 1.96*sqrt(p_var[2,])),
    Branch = factor(rep(c("Vagal (x1)", "Sympathetic (x2)"), each = length(time_steps)),
                    levels = c("Vagal (x1)", "Sympathetic (x2)"))
  )
  p_states <- ggplot(df_states, aes(x = Time, y = Activity, color = Branch, fill = Branch)) +
    geom_hline(yintercept = 0, color = "gray80", linewidth = 1) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c("Vagal (x1)" = "#2980B9", "Sympathetic (x2)" = "#D35400")) +
    scale_fill_manual(values = c("Vagal (x1)" = "#2980B9", "Sympathetic (x2)" = "#D35400")) +
    theme_classic(base_size = 11) +
    labs(title = "C. Neural Latent States", x = "Beat Sequence", y = "Amplitude") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"))

  sp_obs <- spectrum(rr_ts - mean(rr_ts), plot = FALSE, spans = c(3, 3))
  sp_rec <- spectrum(rr_rec - mean(rr_rec), plot = FALSE, spans = c(3, 3))
  df_psd <- data.frame(
    Freq = c(sp_obs$freq, sp_rec$freq),
    Power = c(sp_obs$spec, sp_rec$spec),
    Signal = factor(rep(c("Observed", "SDE Reconstructed"), c(length(sp_obs$freq), length(sp_rec$freq))))
  )
  p_psd <- ggplot(df_psd, aes(x = Freq, y = Power, color = Signal)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    scale_y_log10() +
    scale_color_manual(values = c("Observed" = "gray70", "SDE Reconstructed" = "#8E44AD")) +
    theme_classic(base_size = 11) +
    labs(title = "D. Spectral Preservation", x = "Frequency (Cycles/Beat)", y = "Power (Log Scale)") +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(face = "bold"))

  # =========================================================================
  # NATIVE 3D PHASE-SPACE (Wrapped for Patchwork)
  # =========================================================================
  margin_x1 <- max(0.05, max(abs(states[1, ])) * 0.2)
  margin_x2 <- max(0.05, max(abs(states[2, ])) * 0.2)

  x1_seq <- seq(min(states[1, ]) - margin_x1, max(states[1, ]) + margin_x1, length.out = 45)
  x2_seq <- seq(min(states[2, ]) - margin_x2, max(states[2, ]) + margin_x2, length.out = 45)

  sigma_stat_v_vis <- exp(opt_theta[4])
  sigma_stat_s_vis <- exp(opt_theta[5])
  Sigma_stat_vis <- diag(c(sigma_stat_v_vis^2, sigma_stat_s_vis^2))
  A_2x2_vis <- matrix(c(-a11, -a12, -a12, -a22), nrow = 2, byrow = TRUE)
  Sigma_c_vis <- -(A_2x2_vis %*% Sigma_stat_vis + Sigma_stat_vis %*% t(A_2x2_vis))

  Z_vis   <- rbind(cbind(-A_2x2_vis, Sigma_c_vis), cbind(matrix(0,2,2), t(A_2x2_vis)))
  expZ_v  <- expm::expm(Z_vis)
  F_vis   <- t(expZ_v[3:4, 3:4])
  Q_vis   <- F_vis %*% expZ_v[1:2, 3:4]; Q_vis <- (Q_vis + t(Q_vis)) / 2

  P_stat_v <- tryCatch(
    matrix(solve(diag(4) - kronecker(F_vis, F_vis), as.vector(Q_vis)), 2, 2),
    error = function(e) matrix(MASS::ginv(diag(4) - kronecker(F_vis, F_vis)) %*% as.vector(Q_vis), 2, 2)
  )
  P_stat_v <- (P_stat_v + t(P_stat_v)) / 2
  Pi       <- tryCatch(solve(P_stat_v), error = function(e) MASS::ginv(P_stat_v))

  U_mat  <- outer(x1_seq, x2_seq, FUN = function(x1, x2) {
    Pi[1,1]*x1^2/2 + Pi[2,2]*x2^2/2 + Pi[1,2]*x1*x2
  })
  U_traj <- (Pi[1,1]*states[1,]^2 + Pi[2,2]*states[2,]^2 +
               2*Pi[1,2]*states[1,]*states[2,]) / 2

  # Capture Base R plotting system into a patchwork grid element
  p_3d <- wrap_elements(panel = ~{
    par(mar = c(1, 2, 2, 1))

    # Render the 3D potential energy surface
    pmat <- persp(x1_seq, x2_seq, U_mat,
                  theta = 35, phi = 25, expand = 0.6,
                  col = "#EBF5FB", border = "#34495E", lwd = 0.3,
                  shade = 0.1,
                  xlab = "\nVagal State (x1)",
                  ylab = "\nSympathetic State (x2)",
                  zlab = "\nEnergy U(x)",
                  ticktype = "detailed",
                  main = "E. 3D Phase-Space Energy Attractor",
                  font.main = 2, cex.main = 1.1)

    # Float the trajectory slightly above the surface to prevent Z-fighting overlap
    z_float <- U_traj + max(U_mat) * 0.02

    # Project 3D coordinates onto the 2D plane via the transformation matrix (pmat)
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#2C3E50", lwd = 1)
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#E74C3C", lwd = 0.5)
  })

  # =========================================================================
  # ASSEMBLE DASHBOARD
  # =========================================================================

  layout <- "
  AABB
  CCDD
  EEEE
  EEEE
  "

  dashboard <- p_ts + p_scatter + p_states + p_psd + p_3d +
    plot_layout(design = layout)

  return(dashboard)
}

#' Advanced State-Space Filter Diagnostics and Residual Analysis
#'
#' @param fit_result The optimized model object returned by fit_autonomic_model()
#' @param rr_ts The numeric vector of the original observed RR intervals
#' @return A single patchwork object containing the 5-panel diagnostic dashboard.
diagnose_autonomic_model <- function(fit_result, rr_ts) {

  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install 'patchwork'")
  require(ggplot2)
  require(patchwork)

  log_rr_ts <- log(rr_ts)
  opt_theta <- fit_result$par

  # 1. Run Exact Kalman Filter to Extract Innovations
  kf_out <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)

  # Remove the first beat (t=1) as it represents the unconditioned prior
  raw_innovations <- kf_out$innovations[-1]

  # Standardize the innovations for normalized statistical evaluation
  std_res <- as.numeric(scale(raw_innovations))
  time_steps <- 1:length(std_res)

  df_res <- data.frame(Time = time_steps, Residual = std_res)

  # =========================================================================
  # A. STANDARDIZED RESIDUALS TIME SERIES
  # =========================================================================
  # Evaluates homoscedasticity and detects isolated outliers or regime shifts
  p_time <- ggplot(df_res, aes(x = Time, y = Residual)) +
    geom_hline(yintercept = 0, color = "#2C3E50", linewidth = 0.8) +
    geom_hline(yintercept = c(-2, 2), color = "#E74C3C", linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = c(-3, 3), color = "#E74C3C", linetype = "dotted", alpha = 0.5) +
    geom_line(color = "#7F8C8D", linewidth = 0.5, alpha = 0.8) +
    theme_classic(base_size = 11) +
    labs(title = "A. Standardized Filter Innovations", x = "Beat Sequence", y = "Standardized Error") +
    theme(plot.title = element_text(face = "bold"))

  # =========================================================================
  # B. RESIDUAL DISTRIBUTION
  # =========================================================================
  # Evaluates the Gaussian noise assumption via Kernel Density vs. N(0,1)
  p_dist <- ggplot(df_res, aes(x = Residual)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#BDC3C7", color = "white", alpha = 0.7) +
    geom_density(color = "#2980B9", linewidth = 1) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
    theme_classic(base_size = 11) +
    labs(title = "B. Innovation Distribution", x = "Standardized Error", y = "Density") +
    theme(plot.title = element_text(face = "bold"))

  # =========================================================================
  # C. NORMAL Q-Q PLOT
  # =========================================================================
  # Evaluates tail behaviors and skewness
  p_qq <- ggplot(df_res, aes(sample = Residual)) +
    stat_qq(color = "#34495E", alpha = 0.4, size = 1.5) +
    stat_qq_line(color = "#E74C3C", linewidth = 1, linetype = "dashed") +
    theme_classic(base_size = 11) +
    labs(title = "C. Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme(plot.title = element_text(face = "bold"))

  # =========================================================================
  # D. AUTOCORRELATION FUNCTION (ACF)
  # =========================================================================
  # Evaluates structural leakage (unmodeled dynamics remaining in the error)
  acf_res <- acf(std_res, lag.max = 20, plot = FALSE)
  df_acf <- data.frame(Lag = acf_res$lag[-1], ACF = acf_res$acf[-1])
  ci_acf <- qnorm((1 + 0.95)/2) / sqrt(length(std_res))

  p_acf <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "#2C3E50", linewidth = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-ci_acf, ci_acf), color = "#E74C3C", linetype = "dashed") +
    theme_classic(base_size = 11) +
    labs(title = "D. Autocorrelation Function (ACF)", x = "Lag (Beats)", y = "ACF") +
    theme(plot.title = element_text(face = "bold"))

  # =========================================================================
  # E. LJUNG-BOX TEST INDEPENDENCE TRAJECTORY
  # =========================================================================
  lags <- 11:30
  p_vals <- sapply(lags, function(l) {
    tryCatch({
      Box.test(std_res, lag = l, type = "Ljung-Box", fitdf = length(opt_theta))$p.value
    }, error = function(e) NA)
  })

  df_lb <- data.frame(Lag = lags, P_Value = p_vals)

  p_lb <- ggplot(df_lb, aes(x = Lag, y = P_Value)) +
    geom_point(color = "#2980B9", size = 2) +
    geom_line(color = "#2980B9", alpha = 0.5) +
    geom_hline(yintercept = 0.05, color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme_classic(base_size = 11) +
    labs(title = "E. Ljung-Box Independence Test", x = "Lag (Beats)", y = "P-Value") +
    theme(plot.title = element_text(face = "bold"))

  # =========================================================================
  # ASSEMBLE DASHBOARD
  # =========================================================================

  # Design: Time Series spans the entire top row.
  # Distribution and Q-Q share the middle row.
  # Correlogram and Hypothesis Tests share the bottom row.
  layout <- "
  AAAAAA
  BBBCCC
  DDDEEE
  "

  dashboard <- p_time + p_dist + p_qq + p_acf + p_lb +
    plot_layout(design = layout)

  return(dashboard)
}

#' Evaluate Goodness-of-Fit Metrics (Log-Linear Fixed)
evaluate_single_fit <- function(fit_result, rr_ts) {
  # FIX: Directly pass the boolean flag from the fit_result object
  converged <- fit_result$convergence
  nll <- fit_result$value

  # AIC correctly scales dynamically with the 7 parameter length
  aic <- 2 * length(fit_result$par) + 2 * nll

  log_rr_ts <- log(rr_ts)
  kf_out <- suppressWarnings(run_kalman_filter(fit_result$par, log_rr_ts, return_states = TRUE))
  innovations <- kf_out$innovations[-1]

  # Evaluated at lag 20 to ensure valid degrees of freedom against 7 estimated parameters
  lb_test <- tryCatch(
    Box.test(innovations, lag = 20, type = "Ljung-Box", fitdf = length(fit_result$par)),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  return(data.frame(
    Converged = converged,
    NLL = round(nll, 2),
    AIC = round(aic, 2),
    LB_Stat = round(as.numeric(lb_test$statistic), 2),
    LB_p_value = as.numeric(lb_test$p.value)
  ))
}

#' Batch Validate Autonomic Model Across Multiple Recordings
#'
#' @param rr_batch_list A list of numeric vectors, where each vector is an RR time series.
#' @param subject_ids A character/numeric vector of IDs corresponding to the batch list.
#' @return A data frame summarizing the validation metrics for all subjects.
batch_validate_models <- function(rr_batch_list, subject_ids) {

  num_records <- length(rr_batch_list)
  results_list <- vector("list", num_records)

  for (i in seq_len(num_records)) {
    id <- subject_ids[i]
    rr_ts <- rr_batch_list[[i]]

    run_status <- tryCatch({

      fit <- fit_autonomic_model(rr_ts)
      metrics <- evaluate_single_fit(fit, rr_ts)

      data.frame(
        Subject_ID = id,
        Status = "Success",
        Converged = metrics$Converged,
        NLL = metrics$NLL,
        AIC = metrics$AIC,
        LB_Stat = metrics$LB_Stat,
        LB_p_value = metrics$LB_p_value,
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      data.frame(
        Subject_ID = id,
        Status = paste("Failed:", e$message),
        Converged = NA,
        NLL = NA,
        AIC = NA,
        LB_Stat = NA,
        LB_p_value = NA,
        stringsAsFactors = FALSE
      )
    })

    results_list[[i]] <- run_status

    if (i %% 10 == 0) cat(sprintf("Processed %d / %d datasets...\n", i, num_records))
  }

  final_results_df <- do.call(rbind, results_list)
  return(final_results_df)
}
