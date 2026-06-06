#' Run Kalman Filter with Exact Van Loan Discretization and Joseph Stabilized Update
#'
#' @param theta Numeric vector of 10 transformed parameters.
#' @param log_rr_ts Numeric vector of log-transformed RR intervals.
#' @param dt Time step (e.g., 1 for beat-to-beat).
#' @param return_states Boolean; if TRUE, returns state variables and innovations.
run_kalman_filter <- function(theta, log_rr_ts, dt = 1, return_states = FALSE) {

  # 1. Autonomic Kinetics
  tau_v_min <- 0.3; tau_v_max <- 2.0
  tau_s_min <- 10.0; tau_s_max <- 30.0

  tau_v <- tau_v_min + (tau_v_max - tau_v_min) * plogis(theta[1])
  tau_s <- tau_s_min + (tau_s_max - tau_s_min) * plogis(theta[2])

  a11 <- 1 / tau_v
  a22 <- 1 / tau_s
  max_a12 <- sqrt(a11 * a22)
  a12 <- max_a12 * plogis(theta[3])

  sig1 <- exp(theta[4])
  sig2 <- exp(theta[5])

  # 2. Intrinsic SA Node Memory (AR3 Persistent State) via Durbin-Levinson
  k1 <- 2 * plogis(theta[6]) - 1
  k2 <- 2 * plogis(theta[7]) - 1
  k3 <- 2 * plogis(theta[8]) - 1

  phi1_1 <- k1
  phi2_2 <- k2
  phi2_1 <- phi1_1 - k2 * phi1_1
  phi3_3 <- k3
  phi3_2 <- phi2_2 - k3 * phi2_1
  phi3_1 <- phi2_1 - k3 * phi2_2

  phi1 <- phi3_1
  phi2 <- phi3_2
  phi3 <- phi3_3

  sig_eta <- exp(theta[9])
  mu_rr <- theta[10]

  N <- length(log_rr_ts)

  # ==============================================================================
  # 3. EXACT DISCRETIZATION VIA VAN LOAN'S METHOD (Solves the fast tau_v problem)
  # ==============================================================================
  A_2x2 <- matrix(c(-a11, -a12, -a12, -a22), nrow = 2, byrow = TRUE)
  Sigma_c <- matrix(c(sig1^2, 0, 0, sig2^2), nrow = 2, ncol = 2)

  # Build the 4x4 augmented block matrix
  Z <- rbind(
    cbind(-A_2x2, Sigma_c),
    cbind(matrix(0, 2, 2), t(A_2x2))
  )

  # A single matrix exponential yields both exact F and exact Q
  expZ <- expm::expm(Z * dt)
  F_2x2 <- t(expZ[3:4, 3:4])
  Q_2x2 <- F_2x2 %*% expZ[1:2, 3:4]
  # ==============================================================================

  # Augment into a 5x5 Block Matrix (Autonomic Dynamics + AR3 Memory)
  F_mat <- matrix(0, nrow = 5, ncol = 5)
  F_mat[1:2, 1:2] <- F_2x2
  F_mat[3, 3:5] <- c(phi1, phi2, phi3)
  F_mat[4, 3] <- 1
  F_mat[5, 4] <- 1

  # 5x5 Exact Process Noise Covariance
  Q <- matrix(0, nrow = 5, ncol = 5)
  Q[1:2, 1:2] <- Q_2x2
  Q[3, 3] <- sig_eta^2

  H <- matrix(c(1, -1, 1, 0, 0), nrow = 1)
  R <- 1e-8
  R_mat <- matrix(R, 1, 1)

  x_hat <- matrix(0, nrow = 5, ncol = N)
  P_var <- matrix(0, nrow = 5, ncol = N)
  innovations <- numeric(N)

  # ==============================================================================
  # 4. EXACT STATIONARY INITIALIZATION (Discrete Lyapunov Equation)
  # ==============================================================================
  # Solves vec(P) = (I - F x F)^-1 * vec(Q) to instantly eliminate burn-in
  P_init <- tryCatch({
    I_25 <- diag(25)
    F_kron_F <- kronecker(F_mat, F_mat)
    vec_Q <- as.vector(Q)
    matrix(solve(I_25 - F_kron_F, vec_Q), 5, 5)
  }, error = function(e) diag(5) * 1.0) # Fallback to Identity if optimizer proposes singular bounds

  P <- P_init
  P_var[, 1] <- diag(P)
  # ==============================================================================

  log_likelihood <- 0
  I_5 <- diag(5)

  for (t in 2:N) {
    x_pred <- F_mat %*% x_hat[, t-1]
    P_pred <- F_mat %*% P %*% t(F_mat) + Q

    y_pred <- (H %*% x_pred)[1,1] + mu_rr
    y_obs <- log_rr_ts[t]

    y_err <- y_obs - y_pred
    innovations[t] <- y_err

    S <- (H %*% P_pred %*% t(H))[1,1] + R

    if (is.na(S) || S <= 1e-10) {
      if (return_states) stop("Filter diverged.")
      return(1e9)
    }

    K <- P_pred %*% t(H) / S
    x_hat[, t] <- x_pred + K * y_err

    # ==============================================================================
    # 5. JOSEPH STABILIZED COVARIANCE UPDATE (Guarantees Positive Definiteness)
    # ==============================================================================
    IKH <- I_5 - K %*% H
    P <- IKH %*% P_pred %*% t(IKH) + K %*% R_mat %*% t(K)
    P <- (P + t(P)) / 2
    # ==============================================================================

    P_var[, t] <- diag(P)

    log_likelihood <- log_likelihood - 0.5 * (log(2 * pi * S) + (y_err^2) / S)
  }

  nll <- -log_likelihood
  if (!is.finite(nll)) return(1e9)

  if (return_states) {
    return(list(x_hat = x_hat, P_var = P_var, innovations = innovations))
  } else {
    return(nll)
  }
}

#' Fit the 10-Parameter Augmented Autonomic SDE Model
fit_autonomic_model <- function(rr_ts) {

  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Package 'lhs' is required.")
  }

  log_rr_ts <- log(rr_ts)
  mu_log_rr <- mean(log_rr_ts)

  # LHS already implements Claude's recommendation for global multi-start robustness
  n_starts <- 20
  lhc <- lhs::randomLHS(n_starts, 10)
  theta_starts <- matrix(0, nrow = n_starts, ncol = 10)

  theta_starts[, 1] <- qlogis(0.1) + lhc[, 1] * (qlogis(0.9) - qlogis(0.1))
  theta_starts[, 2] <- qlogis(0.1) + lhc[, 2] * (qlogis(0.9) - qlogis(0.1))
  theta_starts[, 3] <- qlogis(0.01) + lhc[, 3] * (qlogis(0.5) - qlogis(0.01))
  theta_starts[, 4] <- log(0.01) + lhc[, 4] * (log(0.1) - log(0.01))
  theta_starts[, 5] <- log(0.01) + lhc[, 5] * (log(0.1) - log(0.01))
  theta_starts[, 6] <- qlogis(0.5) + lhc[, 6] * (qlogis(0.95) - qlogis(0.5))
  theta_starts[, 7] <- qlogis(0.1) + lhc[, 7] * (qlogis(0.9) - qlogis(0.1))
  theta_starts[, 8] <- qlogis(0.1) + lhc[, 8] * (qlogis(0.9) - qlogis(0.1))
  theta_starts[, 9] <- log(0.001) + lhc[, 9] * (log(0.05) - log(0.001))
  theta_starts[, 10] <- (mu_log_rr - 0.15) + lhc[, 10] * (0.30)

  best_nll <- Inf
  best_theta <- NULL

  for (i in 1:n_starts) {
    fit_try <- tryCatch({
      optim(par = theta_starts[i, ],
            fn = run_kalman_filter,
            log_rr_ts = log_rr_ts,
            dt = 1,
            method = "BFGS",
            hessian = FALSE,
            control = list(maxit = 300))
    }, error = function(e) list(value = Inf))

    if (fit_try$value < best_nll) {
      best_nll <- fit_try$value
      best_theta <- fit_try$par
    }
  }

  if (is.infinite(best_nll)) stop("Filter divergence: LHS failed.")

  final_fit <- optim(par = best_theta,
                     fn = run_kalman_filter,
                     log_rr_ts = log_rr_ts,
                     dt = 1,
                     method = "BFGS",
                     hessian = TRUE,
                     control = list(maxit = 2500, trace = 1))

  return(final_fit)
}


#' Extract Raw Model Parameters
get_model_ci <- function(fit_result) {
  cov_mat_raw <- tryCatch(MASS::ginv(fit_result$hessian), error = function(e) matrix(NA, 10, 10))
  cov_mat_pd <- tryCatch(as.matrix(Matrix::nearPD(cov_mat_raw, ensureSymmetry = TRUE)$mat),
                         error = function(e) cov_mat_raw)

  std_errs <- sqrt(pmax(diag(cov_mat_pd), 0))
  param_names <- c("logit_tau_v", "logit_tau_s", "logit_a12_frac", "log_sig1", "log_sig2", "logit_k1", "logit_k2", "logit_k3", "log_sig_eta", "mu_rr")

  ci_lower <- fit_result$par - 1.96 * std_errs
  ci_upper <- fit_result$par + 1.96 * std_errs

  results_df <- data.frame(Parameter = param_names, Estimate = fit_result$par,
                           Std_Error = std_errs, CI_Lower = ci_lower, CI_Upper = ci_upper)
  return(results_df)
}

#' Generate a Bounded Physiological Interpretation Report
report_fit <- function(fit_result) {

  raw_ci <- get_model_ci(fit_result)
  get_val <- function(param) raw_ci$Estimate[raw_ci$Parameter == param]

  tau_v_est <- 0.5 + 3.5 * plogis(get_val("logit_tau_v"))
  tau_s_est <- 8.0 + 27.0 * plogis(get_val("logit_tau_s"))

  mu_0_est_ms <- exp(get_val("mu_rr"))
  sig_v_frac <- exp(get_val("log_sig1"))
  sig_s_frac <- exp(get_val("log_sig2"))

  k1_est <- 2 * plogis(get_val("logit_k1")) - 1
  k2_est <- 2 * plogis(get_val("logit_k2")) - 1
  k3_est <- 2 * plogis(get_val("logit_k3")) - 1
  sig_eta_frac <- exp(get_val("log_sig_eta"))

  report <- data.frame(
    Parameter_Symbol = c("tau_v", "tau_s", "sigma_v", "sigma_s", "k1", "k2", "k3", "sigma_eta", "mu_0"),
    Interpretation = c(
      "Vagal Recovery Time Constant", "Sympathetic Recovery Time Constant",
      "Vagal Process Noise (Fractional)", "Sympathetic Process Noise (Fractional)",
      "SA Memory PACF 1 (Lag 1)", "SA Memory PACF 2 (Lag 2)", "SA Memory PACF 3 (Lag 3)",
      "SA Memory Driving Noise (Broadband)", "Intrinsic Pacemaker Baseline"
    ),
    Estimate = round(c(tau_v_est, tau_s_est, sig_v_frac, sig_s_frac, k1_est, k2_est, k3_est, sig_eta_frac, mu_0_est_ms), 4),
    Units = c("Beats", "Beats", "Unitless (%)", "Unitless (%)", "Correlation", "Correlation", "Correlation", "Unitless (%)", "ms")
  )
  return(report)
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

  # Extract System Kinetics
  tau_v <- 0.5 + 3.5 * plogis(opt_theta[1])
  tau_s <- 8.0 + 27.0 * plogis(opt_theta[2])
  a11 <- 1 / tau_v
  a22 <- 1 / tau_s
  a12 <- sqrt(a11 * a22) * plogis(opt_theta[3])
  mu_rr_log <- opt_theta[10]

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

  H <- matrix(c(1, -1, 1, 0, 0), nrow = 1)
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
  p_psd <- ggplot(df_psd, aes(x = Freq, y = Power, color = Signal, linetype = Signal)) +
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

  U_mat <- outer(x1_seq, x2_seq, FUN = function(x1, x2) {
    0.5 * a11 * x1^2 + 0.5 * a22 * x2^2 + a12 * x1 * x2
  })

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

    # Extract patient trajectory
    U_traj <- 0.5 * a11 * states[1,]^2 + 0.5 * a22 * states[2,]^2 + a12 * states[1,] * states[2,]

    # Float the trajectory slightly above the surface to prevent Z-fighting overlap
    z_float <- U_traj + max(U_mat) * 0.02

    # Project 3D coordinates onto the 2D plane via the transformation matrix (pmat)
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#2C3E50", lwd = 2)
    lines(trans3d(states[1,], states[2,], z_float, pmat = pmat), col = "#E74C3C", lwd = 1)
  })

  # =========================================================================
  # ASSEMBLE DASHBOARD
  # =========================================================================

  # Design forces the Time Series and States to take 3 columns,
  # while Scatter and PSD take 1 column.
  # The 3D plot takes all 4 columns at the bottom.
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

#' Elegant ggplot2 Diagnostic Dashboard (Log-Linear Fixed)
diagnose_autonomic_model <- function(fit_result, rr_ts) {
  require(ggplot2)
  require(gridExtra)

  log_rr_ts <- log(rr_ts)

  opt_theta <- fit_result$par
  kf_out <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)

  innovations <- kf_out$innovations[-1]

  df_res <- data.frame(Time = 1:length(innovations), Innovation = innovations)

  p1 <- ggplot(df_res, aes(x = Time, y = Innovation)) +
    geom_line(color = "gray40", linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "#D85A30", linetype = "dashed", linewidth = 0.8) +
    theme_classic() + labs(title = "A. Innovation Time Series", x = "Beat Number", y = "Error (Fractional)")

  p2 <- ggplot(df_res, aes(sample = Innovation)) +
    stat_qq(color = "#378ADD", alpha = 0.6, size = 1.5) +
    stat_qq_line(color = "#D85A30", linetype = "dashed", linewidth = 0.8) +
    theme_classic() + labs(title = "B. Q-Q Plot", x = "Theoretical", y = "Sample")

  acf_res <- acf(innovations, plot = FALSE)
  df_acf <- data.frame(Lag = acf_res$lag[-1], ACF = acf_res$acf[-1])
  ci_acf <- qnorm((1 + 0.95)/2) / sqrt(length(innovations))

  p3 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "gray50", linewidth = 1) +
    geom_point(color = "#378ADD", size = 2) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-ci_acf, ci_acf), color = "#D85A30", linetype = "dashed") +
    theme_classic() + labs(title = "C. Autocorrelation Function", x = "Lag", y = "ACF")

  gridExtra::grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 1), c(2, 3)))
}

#' Evaluate Goodness-of-Fit Metrics (Log-Linear Fixed)
evaluate_single_fit <- function(fit_result, rr_ts) {
  converged <- ifelse(fit_result$convergence == 0, TRUE, FALSE)
  nll <- fit_result$value

  # AIC correctly scales dynamically with the 10 parameter length
  aic <- 2 * length(fit_result$par) + 2 * nll

  log_rr_ts <- log(rr_ts)
  kf_out <- suppressWarnings(run_kalman_filter(fit_result$par, log_rr_ts, return_states = TRUE))
  innovations <- kf_out$innovations[-1]

  lb_test <- tryCatch(
    Box.test(innovations, lag = 5, type = "Ljung-Box"),
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
