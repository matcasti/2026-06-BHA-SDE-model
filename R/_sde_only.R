# ==============================================================================
# MODULE 1: PARAMETER MAPPING & PHYSICAL TRANSLATIONS
# ==============================================================================

#' Map Unconstrained Optimization Vector to Physical Model Parameters
#'
#' @param theta Numeric vector of 7 unconstrained optimization parameters.
#' @return A named list containing bounding-protected physical constants.
map_theta_to_physical <- function(theta) {
  tau_v_min <- 0.3; tau_v_max <- 2.0
  rho_min   <- 5.0; rho_max   <- 15.0

  # Hyperspherical and exponential mappings
  tau_v   <- tau_v_min + (tau_v_max - tau_v_min) * (sin(theta[1])^2)
  log_rho <- log(rho_min) + (log(rho_max) - log(rho_min)) * (sin(theta[2])^2)
  tau_s   <- tau_v * exp(log_rho)

  a11     <- 1 / tau_v
  a22     <- 1 / tau_s
  max_a12 <- sqrt(a11 * a22)
  a12     <- max_a12 * 0.999 * (sin(theta[3])^2)

  sigma_proc_v   <- exp(theta[4])
  sigma_proc_s   <- exp(theta[5])

  # Sinoatrial node memory parameters
  k1             <- 0.999 * (sin(theta[6])^2)
  sigma_stat_eta <- exp(theta[7])
  sig_eta        <- sigma_stat_eta * sqrt(1 - k1^2)

  list(
    tau_v = tau_v, tau_s = tau_s, a11 = a11, a22 = a22, a12 = a12,
    sigma_proc_v = sigma_proc_v, sigma_proc_s = sigma_proc_s,
    k1 = k1, sigma_stat_eta = sigma_stat_eta, sig_eta = sig_eta
  )
}

# ==============================================================================
# MODULE 2: SYSTEM DISCRETIZATION VIA VAN LOAN'S METHOD
# ==============================================================================

#' Discretize Continuous-Time 2x2 Autonomic SDE Kinetics
#'
#' @param p Named list of physical parameters derived from map_theta_to_physical().
#' @param dt Time step interval scalar.
#' @return A list containing the discrete transitions F_2x2 and process noise Q_2x2.
discretize_kinetics <- function(p, dt = 1) {
  A_2x2 <- matrix(c(-p$a11, -p$a12, -p$a12, -p$a22), nrow = 2, byrow = TRUE)

  # Explicit continuous process noise mapping
  sig1 <- p$sigma_proc_v * sqrt(2 * p$a11)
  sig2 <- p$sigma_proc_s * sqrt(2 * p$a22)
  Sigma_c <- matrix(c(sig1^2, 0, 0, sig2^2), nrow = 2, ncol = 2)

  if (any(!is.finite(A_2x2)) || any(!is.finite(Sigma_c))) return(NULL)

  # Build block matrix for automatic matrix exponential integration
  Z <- rbind(
    cbind(-A_2x2, Sigma_c),
    cbind(matrix(0, 2, 2), t(A_2x2))
  )

  expZ <- expm::expm(Z * dt)
  F_2x2 <- t(expZ[3:4, 3:4])
  Q_2x2 <- F_2x2 %*% expZ[1:2, 3:4]
  Q_2x2 <- (Q_2x2 + t(Q_2x2)) / 2  # Enforce numerical symmetry

  list(F_2x2 = F_2x2, Q_2x2 = Q_2x2)
}

# ==============================================================================
# MODULE 3: STATE-SPACE MATRIX COUPLING & INITIALIZATION
# ==============================================================================

#' Assemble Full 3-Dimensional System Transition and Noise Covariance Matrices
#'
#' @param disc List containing the 2x2 discretized SDE matrices.
#' @param p Named list of physical parameters.
#' @return A list containing the full 3x3 F_mat and Q matrices.
assemble_3d_system <- function(disc, p) {
  F_mat <- matrix(0, nrow = 3, ncol = 3)
  F_mat[1:2, 1:2] <- disc$F_2x2
  F_mat[3, 3]     <- p$k1

  Q_mat <- matrix(0, nrow = 3, ncol = 3)
  Q_mat[1:2, 1:2] <- disc$Q_2x2
  Q_mat[3, 3]     <- p$sig_eta^2

  list(F_mat = F_mat, Q_mat = Q_mat)
}

#' Compute Exact Stationary Initial State Covariance Matrix
#'
#' @param F_mat Full 3x3 system transition matrix.
#' @param Q_mat Full 3x3 process noise matrix.
#' @return A stable positive semi-definite 3x3 covariance matrix.
initialize_stationary_covariance <- function(F_mat, Q_mat) {
  I_9 <- diag(9)
  F_kron_F <- kronecker(F_mat, F_mat)
  vec_Q <- as.vector(Q_mat)

  P_init <- tryCatch({
    matrix(solve(I_9 - F_kron_F, vec_Q), 3, 3)
  }, error = function(e) {
    matrix(MASS::ginv(I_9 - F_kron_F) %*% vec_Q, 3, 3)
  })

  P_init <- (P_init + t(P_init)) / 2
  eigs <- eigen(P_init, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs < 1e-10)) {
    P_init <- as.matrix(Matrix::nearPD(P_init, ensureSymmetry = TRUE)$mat)
  }
  return(P_init)
}

# ==============================================================================
# MODULE 4: INTERCHANGEABLE MEASUREMENT DOCKING LAYER
# ==============================================================================

#' Execute Single-Step Linear Gaussian Measurement Update
#'
#' @note DOCKING NODE FOR POINT PROCESSES: Replace this specific block with your
#'   Inverse Gaussian Point Process local EKF gradient/Hessian calculation loop
#'   when transitioning models later.
#'
#' @return List containing log-likelihood contribution and updated state properties.
measurement_step_gaussian <- function(y_obs, x_pred, P_pred, mu_rr, H, R, R_mat, I_3) {
  y_pred <- (H %*% x_pred)[1,1] + mu_rr
  y_err  <- y_obs - y_pred
  S      <- (H %*% P_pred %*% t(H))[1,1] + R

  if (!is.finite(S) || S < 1e-9) return(NULL)

  # Joseph-Stabilized Filter Calculations
  K       <- P_pred %*% t(H) / S
  x_hat_t <- x_pred + K * y_err

  IKH <- I_3 - K %*% H
  P_t <- IKH %*% P_pred %*% t(IKH) + K %*% R_mat %*% t(K)
  P_t <- (P_t + t(P_t)) / 2

  log_lik_contrib <- -0.5 * (log(2 * pi * S) + (y_err^2) / S)

  list(
    x_hat_t = x_hat_t, P_t = P_t, log_lik_contrib = log_lik_contrib,
    y_err = y_err, innov_sd = sqrt(S)
  )
}

# ==============================================================================
# MODULE 5: CORE ESTIMATION PASS ROUTINE
# ==============================================================================

#' Run Decentralized State-Space Estimation Pass
run_kalman_filter <- function(theta, log_rr_ts, dt = 1, return_states = FALSE) {
  # Pipe through isolated preparation subroutines
  p    <- map_theta_to_physical(theta)
  disc <- discretize_kinetics(p, dt)
  if (is.null(disc)) {
    if (return_states) stop("Filter diverged during integration.")
    return(1e9)
  }

  sys    <- assemble_3d_system(disc, p)
  P_init <- initialize_stationary_covariance(sys$F_mat, sys$Q_mat)

  N      <- length(log_rr_ts)
  mu_rr  <- mean(log_rr_ts)

  # Observation parameters
  H <- matrix(c(1, -1, 1), nrow = 1)
  R <- 1e-5; R_mat <- matrix(R, 1, 1); I_3 <- diag(3)

  # Arrays allocation
  x_hat       <- matrix(0, nrow = 3, ncol = N)
  P_var       <- matrix(0, nrow = 3, ncol = N)
  innovations <- numeric(N)
  innov_sd    <- numeric(N)

  x_pred         <- matrix(0, nrow = 3, ncol = 1)
  P_pred         <- P_init
  log_likelihood <- 0

  # Temporal filter propagation sequence
  for (t in 1:N) {
    upd <- measurement_step_gaussian(log_rr_ts[t], x_pred, P_pred, mu_rr, H, R, R_mat, I_3)
    if (is.null(upd)) {
      if (return_states) stop("Filter collapsed during measurement step.")
      return(1e9)
    }

    x_hat[, t]     <- upd$x_hat_t
    P_var[, t]     <- diag(upd$P_t)
    innovations[t] <- upd$y_err
    innov_sd[t]    <- upd$innov_sd
    log_likelihood <- log_likelihood + upd$log_lik_contrib

    if (t < N) {
      x_pred <- sys$F_mat %*% x_hat[, t]
      P_pred <- sys$F_mat %*% upd$P_t %*% t(sys$F_mat) + sys$Q_mat
    }
  }

  nll <- -log_likelihood
  if (!is.finite(nll)) return(1e9)

  if (return_states) {
    return(list(x_hat = x_hat, P_var = P_var, innovations = innovations, innov_sd = innov_sd))
  } else {
    return(nll)
  }
}

# ==============================================================================
# MODULE 6: CALIBRATION & INFERENCE ENGINE
# ==============================================================================

evaluate_penalized_nll <- function(theta, log_rr_ts, lambda = 0) {
  base_nll <- run_kalman_filter(theta, log_rr_ts, dt = 1, return_states = FALSE)
  if (is.infinite(base_nll) || base_nll == 1e9) return(1e9)
  return(base_nll + lambda * (theta[3]^2))
}

fit_model <- function(rr_ts, lambda = 0) {
  log_rr_ts <- log(rr_ts)

  init_theta1 <- asin(sqrt((1.0 - 0.3) / (2.0 - 0.3)))
  init_theta2 <- asin(sqrt((log(10.0) - log(5.0)) / (log(15.0) - log(5.0))))
  init_theta3 <- asin(sqrt(0.10))
  init_theta4 <- log(0.05)
  init_theta5 <- log(0.05)
  init_theta6 <- asin(sqrt(0.50))
  init_theta7 <- log(0.01)

  theta_init <- c(init_theta1, init_theta2, init_theta3, init_theta4,
                  init_theta5, init_theta6, init_theta7)

  fit_try <- tryCatch({
    optim(
      par = theta_init, fn = evaluate_penalized_nll, log_rr_ts = log_rr_ts,
      lambda = lambda, method = "BFGS", hessian = TRUE, control = list(maxit = 10000, trace = 1)
    )
  }, error = function(e) { return(NULL) })

  if (is.null(fit_try) || is.infinite(fit_try$value) || is.na(fit_try$value)) {
    return(list(par = rep(NA, 7), value = Inf, hessian = matrix(NA, 7, 7), convergence = FALSE))
  }

  return(list(par = fit_try$par, value = fit_try$value,
              hessian = fit_try$hessian, convergence = (fit_try$convergence == 0)))
}

get_model_ci <- function(fit_result) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Install 'numDeriv' library.")

  cov_mat_raw <- tryCatch(MASS::ginv(fit_result$hessian), error = function(e) matrix(NA, 7, 7))
  cov_mat_pd  <- tryCatch(as.matrix(Matrix::nearPD(cov_mat_raw, ensureSymmetry = TRUE)$mat),
                          error = function(e) cov_mat_raw)

  # Centralized Delta Method transformation function utilizing Module 1 mapping
  transform_theta <- function(th) {
    p <- map_theta_to_physical(th)
    max_a12 <- sqrt(p$a11 * p$a22)
    a12_f   <- p$a12 / max_a12
    return(c(p$tau_v, p$tau_s, a12_f, p$sigma_proc_v, p$sigma_proc_s, p$k1, p$sigma_stat_eta))
  }

  phys_est <- transform_theta(fit_result$par)
  J        <- tryCatch(numDeriv::jacobian(transform_theta, fit_result$par), error = function(e) matrix(NA, 7, 7))
  phys_cov <- J %*% cov_mat_pd %*% t(J)
  phys_se  <- sqrt(pmax(diag(phys_cov), 0))

  ci_lower <- phys_est - 1.96 * phys_se
  ci_upper <- phys_est + 1.96 * phys_se
  param_names <- c("tau_v", "tau_s", "a12_fraction", "sigma_proc_v", "sigma_proc_s", "k1", "sigma_stat_eta")

  return(data.frame(Parameter = param_names, Estimate = phys_est,
                    Std_Error = phys_se, CI_Lower = ci_lower, CI_Upper = ci_upper))
}

report_model <- function(fit_result) {
  if (all(is.na(fit_result$par))) return(data.frame(Message = "Optimization Failed"))
  phys_ci <- get_model_ci(fit_result)
  interpretations <- c("Vagal Recovery (beats)", "Sympathetic Recovery (beats)", "Cross-talk Coupling Fraction",
                       "Vagal Diffusion Scale", "Sympathetic Diffusion Scale", "SA Node Memory (AR1)", "SA Node Memory Variance")

  return(data.frame(Parameter_Symbol = phys_ci$Parameter, Interpretation = interpretations,
                    Estimate = round(phys_ci$Estimate, 4), Std_Error = round(phys_ci$Std_Error, 4),
                    CI_Lower = round(phys_ci$CI_Lower, 4), CI_Upper = round(phys_ci$CI_Upper, 4)))
}

# ==============================================================================
# MODULE 7: DIAGNOSTIC GRAPHICS & POST-PASS ANALYTICS
# ==============================================================================

visualize_model <- function(fit_result, rr_ts) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install 'patchwork'")
  require(ggplot2); require(patchwork)

  log_rr_ts <- log(rr_ts)
  opt_theta <- fit_result$par
  mu_rr_log <- mean(log_rr_ts)

  # Access shared modules
  p    <- map_theta_to_physical(opt_theta)
  disc <- discretize_kinetics(p, dt = 1)

  kf_out      <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)
  states      <- kf_out$x_hat
  p_var       <- kf_out$P_var
  innovations <- kf_out$innovations

  time_steps <- 1:length(rr_ts)
  valid_idx  <- 2:length(rr_ts)

  log_rr_pred <- log_rr_ts - innovations
  rr_pred     <- exp(log_rr_pred)
  H           <- matrix(c(1, -1, 1), nrow = 1)
  log_rr_rec  <- as.numeric(H %*% states) + mu_rr_log
  rr_rec      <- exp(log_rr_rec)

  df_ts <- data.frame(Time = time_steps, Observed = rr_ts, Predicted = rr_pred)
  p_ts <- ggplot(df_ts) +
    geom_line(aes(x = Time, y = Observed, color = "Observed"), linewidth = 0.6, alpha = 0.5) +
    geom_line(aes(x = Time, y = Predicted, color = "Kalman Predicted"), linewidth = 0.4, alpha = 0.9) +
    scale_color_manual(values = c("Observed" = "gray50", "Kalman Predicted" = "#E74C3C")) +
    theme_classic(base_size = 11) + labs(title = "A. Signal Prediction Tracking", x = "Beat Sequence", y = "RR Interval (ms)") +
    theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold"))

  r_squared  <- round(cor(rr_ts[valid_idx], rr_pred[valid_idx])^2, 3)
  df_scatter <- data.frame(Obs = rr_ts[valid_idx], Pred = rr_pred[valid_idx])
  p_scatter  <- ggplot(df_scatter, aes(x = Obs, y = Pred)) +
    geom_point(color = "#34495E", alpha = 0.3, pch = 16) + geom_abline(slope = 1, intercept = 0, color = "#E74C3C") +
    annotate("text", x = min(df_scatter$Obs), y = max(df_scatter$Pred), label = paste0("R² = ", r_squared), hjust = 0, vjust = 1, fontface = "bold") +
    theme_classic(base_size = 11) + labs(title = "B. Prediction Accuracy", x = "Observed RR (ms)", y = "Predicted RR (ms)") +
    theme(plot.title = element_text(face = "bold"))

  df_states <- data.frame(
    Time = rep(time_steps, 2), Activity = c(states[1, ], states[2, ]),
    Lower = c(states[1, ] - 1.96*sqrt(p_var[1,]), states[2, ] - 1.96*sqrt(p_var[2,])),
    Upper = c(states[1, ] + 1.96*sqrt(p_var[1,]), states[2, ] + 1.96*sqrt(p_var[2,])),
    Branch = factor(rep(c("Vagal (x1)", "Sympathetic (x2)"), each = length(time_steps)), levels = c("Vagal (x1)", "Sympathetic (x2)"))
  )
  p_states <- ggplot(df_states, aes(x = Time, y = Activity, color = Branch, fill = Branch)) +
    geom_hline(yintercept = 0, color = "gray80", linewidth = 1) + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.6) + scale_color_manual(values = c("Vagal (x1)" = "#2980B9", "Sympathetic (x2)" = "#D35400")) +
    scale_fill_manual(values = c("Vagal (x1)" = "#2980B9", "Sympathetic (x2)" = "#D35400")) +
    theme_classic(base_size = 11) + labs(title = "C. Neural Latent States", x = "Beat Sequence", y = "Amplitude") +
    theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold"))

  sp_obs <- spectrum(rr_ts - mean(rr_ts), plot = FALSE, spans = c(3, 3))
  sp_rec <- spectrum(rr_rec - mean(rr_rec), plot = FALSE, spans = c(3, 3))
  df_psd <- data.frame(Freq = c(sp_obs$freq, sp_rec$freq), Power = c(sp_obs$spec, sp_rec$spec),
                       Signal = factor(rep(c("Observed", "SDE Reconstructed"), c(length(sp_obs$freq), length(sp_rec$freq)))))
  p_psd <- ggplot(df_psd, aes(x = Freq, y = Power, color = Signal)) +
    geom_line(linewidth = 0.8, alpha = 0.8) + scale_y_log10() +
    scale_color_manual(values = c("Observed" = "gray70", "SDE Reconstructed" = "#8E44AD")) +
    theme_classic(base_size = 11) + labs(title = "D. Spectral Preservation", x = "Frequency (Cycles/Beat)", y = "Power (Log Scale)") +
    theme(legend.position = "bottom", legend.title = element_blank(), plot.title = element_text(face = "bold"))

  # Phase space parameters using synchronized continuous process noise parameters
  margin_x1 <- max(0.05, max(abs(states[1, ])) * 0.2)
  margin_x2 <- max(0.05, max(abs(states[2, ])) * 0.2)
  x1_seq    <- seq(min(states[1, ]) - margin_x1, max(states[1, ]) + margin_x1, length.out = 45)
  x2_seq    <- seq(min(states[2, ]) - margin_x2, max(states[2, ]) + margin_x2, length.out = 45)

  P_stat_v <- tryCatch(
    matrix(solve(diag(4) - kronecker(disc$F_2x2, disc$F_2x2), as.vector(disc$Q_2x2)), 2, 2),
    error = function(e) matrix(MASS::ginv(diag(4) - kronecker(disc$F_2x2, disc$F_2x2)) %*% as.vector(disc$Q_2x2), 2, 2)
  )
  P_stat_v <- (P_stat_v + t(P_stat_v)) / 2
  Pi       <- tryCatch(solve(P_stat_v), error = function(e) MASS::ginv(P_stat_v))

  U_mat  <- outer(x1_seq, x2_seq, FUN = function(x1, x2) Pi[1,1]*x1^2/2 + Pi[2,2]*x2^2/2 + Pi[1,2]*x1*x2)
  U_traj <- (Pi[1,1]*states[1,]^2 + Pi[2,2]*states[2,]^2 + 2*Pi[1,2]*states[1,]*states[2,]) / 2

  p_3d <- wrap_elements(panel = ~{
    par(mar = c(1, 2, 2, 1))
    pmat <- persp(x1_seq, x2_seq, U_mat, theta = 35, phi = 25, expand = 0.6,
                  col = "#EBF5FB", border = "#34495E", lwd = 0.3, shade = 0.1,
                  xlab = "\nVagal State (x1)", ylab = "\nSympathetic State (x2)", zlab = "\nEnergy U(x)",
                  ticktype = "detailed", main = "E. 3D Phase-Space Energy Attractor", font.main = 2, cex.main = 1.1)
    z_float <- U_traj + max(U_mat) * 0.02
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#2C3E50", lwd = 1)
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#E74C3C", lwd = 0.5)
  })

  layout    <- "AABB\nCCDD\nEEEE\nEEEE"
  dashboard <- p_ts + p_scatter + p_states + p_psd + p_3d + plot_layout(design = layout)
  return(dashboard)
}

diagnose_model <- function(fit_result, rr_ts) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Install 'patchwork'")
  require(ggplot2); require(patchwork)

  kf_out  <- run_kalman_filter(fit_result$par, log(rr_ts), return_states = TRUE)
  std_res <- kf_out$innovations[-1] / pmax(kf_out$innov_sd[-1], .Machine$double.eps)
  df_res  <- data.frame(Time = 1:length(std_res), Residual = std_res)

  p_time <- ggplot(df_res, aes(x = Time, y = Residual)) +
    geom_hline(yintercept = 0, color = "#2C3E50", linewidth = 0.8) +
    geom_hline(yintercept = c(-2, 2), color = "#E74C3C", linetype = "dashed", alpha = 0.7) +
    geom_line(color = "#7F8C8D", linewidth = 0.5, alpha = 0.8) +
    theme_classic(base_size = 11) + labs(title = "A. Standardized Filter Innovations", x = "Beat Sequence", y = "Standardized Error") +
    theme(plot.title = element_text(face = "bold"))

  p_dist <- ggplot(df_res, aes(x = Residual)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#BDC3C7", color = "white", alpha = 0.7) +
    geom_density(color = "#2980B9", linewidth = 1) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
    theme_classic(base_size = 11) + labs(title = "B. Innovation Distribution", x = "Standardized Error", y = "Density") +
    theme(plot.title = element_text(face = "bold"))

  p_qq <- ggplot(df_res, aes(sample = Residual)) +
    stat_qq(color = "#34495E", alpha = 0.4, size = 1.5) + stat_qq_line(color = "#E74C3C", linewidth = 1, linetype = "dashed") +
    theme_classic(base_size = 11) + labs(title = "C. Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme(plot.title = element_text(face = "bold"))

  acf_res <- acf(std_res, lag.max = 20, plot = FALSE)
  df_acf  <- data.frame(Lag = acf_res$lag[-1], ACF = acf_res$acf[-1])
  p_acf   <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "#2C3E50", linewidth = 1.5, alpha = 0.8) +
    geom_hline(yintercept = c(-1, 1) * (qnorm(0.975)/sqrt(length(std_res))), color = "#E74C3C", linetype = "dashed") +
    theme_classic(base_size = 11) + labs(title = "D. Autocorrelation Function (ACF)", x = "Lag (Beats)", y = "ACF") +
    theme(plot.title = element_text(face = "bold"))

  lags   <- 11:30
  p_vals <- sapply(lags, function(l) {
    tryCatch(Box.test(std_res, lag = l, type = "Ljung-Box", fitdf = length(fit_result$par))$p.value, error = function(e) NA)
  })
  p_lb <- ggplot(data.frame(Lag = lags, P_Value = p_vals), aes(x = Lag, y = P_Value)) +
    geom_point(color = "#2980B9", size = 2) + geom_line(color = "#2980B9", alpha = 0.5) +
    geom_hline(yintercept = 0.05, color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
    scale_y_continuous(limits = c(0, 1)) + theme_classic(base_size = 11) + labs(title = "E. Ljung-Box Independence Test", x = "Lag (Beats)", y = "P-Value") +
    theme(plot.title = element_text(face = "bold"))

  layout    <- "AAAAAA\nBBBCCC\nDDDEEE"
  dashboard <- p_time + p_dist + p_qq + p_acf + p_lb + plot_layout(design = layout)
  return(dashboard)
}

evaluate_single_fit <- function(fit_result, rr_ts) {
  converged <- fit_result$convergence
  nll       <- run_kalman_filter(fit_result$par, log(rr_ts), dt = 1, return_states = FALSE)
  aic       <- 2 * length(fit_result$par) + 2 * nll

  kf_out    <- suppressWarnings(run_kalman_filter(fit_result$par, log(rr_ts), return_states = TRUE))
  std_innov <- kf_out$innovations[-1] / pmax(kf_out$innov_sd[-1], .Machine$double.eps)
  lb_test   <- tryCatch(Box.test(std_innov, lag = 20, type = "Ljung-Box", fitdf = length(fit_result$par)), error = function(e) list(statistic = NA, p.value = NA))

  return(data.frame(Converged = converged, NLL = round(nll, 2), AIC = round(aic, 2),
                    LB_Stat = round(as.numeric(lb_test$statistic), 2), LB_p_value = as.numeric(lb_test$p.value)))
}

batch_validate_models <- function(rr_batch_list, subject_ids) {
  num_records  <- length(rr_batch_list)
  results_list <- vector("list", num_records)

  for (i in seq_len(num_records)) {
    id    <- subject_ids[i]
    rr_ts <- rr_batch_list[[i]]

    results_list[[i]] <- tryCatch({
      fit     <- fit_autonomic_model(rr_ts)
      metrics <- evaluate_single_fit(fit, rr_ts)
      data.frame(Subject_ID = id, Status = "Success", Converged = metrics$Converged, NLL = metrics$NLL,
                 AIC = metrics$AIC, LB_Stat = metrics$LB_Stat, LB_p_value = metrics$LB_p_value, stringsAsFactors = FALSE)
    }, error = function(e) {
      data.frame(Subject_ID = id, Status = paste("Failed:", e$message), Converged = NA, NLL = NA, AIC = NA, LB_Stat = NA, LB_p_value = NA, stringsAsFactors = FALSE)
    })
    if (i %% 10 == 0) cat(sprintf("Processed %d / %d records...\n", i, num_records))
  }
  return(do.call(rbind, results_list))
}
