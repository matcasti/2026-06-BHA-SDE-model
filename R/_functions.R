#' Run Kalman Filter with 3D State-Augmented AR(1) SA Node Memory
#'
#' @param theta Numeric vector of 8 transformed parameters.
#' @param log_rr_ts Numeric vector of log-transformed RR intervals.
#' @param dt Time step (e.g., 1 for beat-to-beat).
#' @param return_states Boolean; if TRUE, returns state variables and innovations.
run_kalman_filter <- function(theta, log_rr_ts, dt = 1, return_states = FALSE) {

  # 1. Autonomic Kinetics (Bounds: ACh = Fast, NE = Slow)
  tau_v_min <- 0.5; tau_v_max <- 4.0
  tau_s_min <- 8.0; tau_s_max <- 35.0

  tau_v <- tau_v_min + (tau_v_max - tau_v_min) * plogis(theta[1])
  tau_s <- tau_s_min + (tau_s_max - tau_s_min) * plogis(theta[2])

  a11 <- 1 / tau_v
  a22 <- 1 / tau_s
  max_a12 <- sqrt(a11 * a22)
  a12 <- max_a12 * plogis(theta[3])

  sig1 <- exp(theta[4])
  sig2 <- exp(theta[5])

  # 2. Intrinsic SA Node Memory (AR1 State)
  # Map theta[6] from real line to (-1, 1) for strictly stationary AR(1) memory
  phi <- 2 * plogis(theta[6]) - 1
  sig_eta <- exp(theta[7]) # Driving noise of the AR(1) process

  mu_rr <- theta[8] # Intrinsic HR (in log-space)

  N <- length(log_rr_ts)

  # 3. 3D State Augmentation
  # Build the 2D Autonomic Transition Matrix
  A_2x2 <- matrix(c(-a11, -a12, -a12, -a22), nrow = 2, byrow = TRUE)
  F_2x2 <- as.matrix(expm::expm(A_2x2 * dt))

  # Augment into a 3x3 Block Matrix (Autonomic Dynamics + SA Memory)
  F_mat <- matrix(0, nrow = 3, ncol = 3)
  F_mat[1:2, 1:2] <- F_2x2
  F_mat[3, 3] <- phi

  M_vec <- c(0, 0, 0)

  # 3x3 Process Noise Covariance
  Q <- diag(c(sig1^2 * dt, sig2^2 * dt, sig_eta^2))

  # Observation matrix now sums: Vagal - Sympathetic + SA_Memory
  H <- matrix(c(1, -1, 1), nrow = 1)

  # Measurement noise is effectively 0 because the noise is explicitly modeled in state 3.
  # We use 1e-8 strictly to prevent numerical singularity in matrix inversion.
  R <- 1e-8

  x_hat <- matrix(0, nrow = 3, ncol = N)
  P_var <- matrix(0, nrow = 3, ncol = N)
  innovations <- numeric(N)

  x_hat[, 1] <- M_vec
  P <- diag(3) * 1.0

  log_likelihood <- 0

  for (t in 2:N) {
    x_pred <- F_mat %*% (x_hat[, t-1] - M_vec) + M_vec
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
    P <- (diag(3) - K %*% H) %*% P_pred
    P <- (P + t(P)) / 2

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

#' Fit the 8-Parameter Augmented Autonomic SDE Model
fit_autonomic_model <- function(rr_ts) {

  log_rr_ts <- log(rr_ts)

  init_tau_v <- 1.5
  init_tau_s <- 15.0
  logit_init_v <- qlogis((init_tau_v - 0.5) / (4.0 - 0.5))
  logit_init_s <- qlogis((init_tau_s - 8.0) / (35.0 - 8.0))

  # Initialize SA memory phi at ~0.5 (plogis(0.5+1)/2 = 0.75 in logit space)
  logit_init_phi <- qlogis(0.75)

  # 8 Parameters
  init_theta <- c(logit_init_v,
                  logit_init_s,
                  qlogis(0.1),
                  log(0.05), log(0.05),
                  logit_init_phi,
                  log(0.01),
                  mean(log_rr_ts))

  fit <- optim(par = init_theta,
               fn = run_kalman_filter,
               log_rr_ts = log_rr_ts,
               dt = 1,
               method = "BFGS",
               hessian = TRUE,
               control = list(maxit = 2000, trace = 1)) # Increased maxit for 8 parameters

  return(fit)
}


#' Extract Raw Model Parameters
get_model_ci <- function(fit_result) {
  cov_mat_raw <- tryCatch(MASS::ginv(fit_result$hessian), error = function(e) matrix(NA, 8, 8))
  cov_mat_pd <- tryCatch(as.matrix(Matrix::nearPD(cov_mat_raw, ensureSymmetry = TRUE)$mat),
                         error = function(e) cov_mat_raw)

  std_errs <- sqrt(pmax(diag(cov_mat_pd), 0))
  param_names <- c("logit_tau_v", "logit_tau_s", "logit_a12_frac", "log_sig1", "log_sig2", "logit_phi", "log_sig_eta", "mu_rr")

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

  # Extract the new Intrinsic SA Node Memory Parameters
  phi_est <- 2 * plogis(get_val("logit_phi")) - 1
  sig_eta_frac <- exp(get_val("log_sig_eta"))

  # Calculate the functional half-life of the SA memory (how long before the mechanical stretch washes out)
  memory_halflife <- ifelse(phi_est > 0, log(0.5) / log(phi_est), 0)

  report <- data.frame(
    Parameter_Symbol = c("tau_v", "tau_s", "sigma_v", "sigma_s", "phi_SA", "tau_SA_memory", "sigma_eta", "mu_0"),
    Interpretation = c(
      "Vagal Recovery Time Constant", "Sympathetic Recovery Time Constant",
      "Vagal Process Noise (Fractional)", "Sympathetic Process Noise (Fractional)",
      "Intrinsic SA Node Memory (AR1)", "SA Node Memory Half-Life",
      "Memory Driving Noise (Fractional)", "Intrinsic Pacemaker Baseline"
    ),
    Estimate = round(c(tau_v_est, tau_s_est, sig_v_frac, sig_s_frac, phi_est, memory_halflife, sig_eta_frac, mu_0_est_ms), 4),
    Units = c("Beats", "Beats", "Unitless (%)", "Unitless (%)", "Correlation", "Beats", "Unitless (%)", "ms")
  )
  return(report)
}

#' Elegant ggplot2 Visualization (Log-Linear Fixed)
visualize_autonomic_model <- function(fit_result, rr_ts) {
  require(ggplot2)
  require(gridExtra)

  log_rr_ts <- log(rr_ts)
  opt_theta <- fit_result$par

  tau_v <- 0.5 + 3.5 * plogis(opt_theta[1])
  tau_s <- 8.0 + 27.0 * plogis(opt_theta[2])
  a11 <- 1 / tau_v
  a22 <- 1 / tau_s
  a12 <- sqrt(a11 * a22) * plogis(opt_theta[3])

  mu_rr_ms <- exp(opt_theta[8])

  kf_out <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)
  states <- kf_out$x_hat
  p_var <- kf_out$P_var
  time_steps <- 1:length(rr_ts)

  df_rr <- data.frame(Time = time_steps, RR = rr_ts, Baseline = mu_rr_ms)
  p1 <- ggplot(df_rr, aes(x = Time)) +
    geom_line(aes(y = RR), color = "black", linewidth = 0.4, alpha = 0.6) +
    geom_hline(yintercept = mu_rr_ms, color = "#D85A30", linewidth = 1.2, linetype = "dashed") +
    theme_minimal() + labs(title = "Observed RR and Intrinsic Baseline", x = "Beat", y = "RR (ms)")

  # Plot B: We only extract rows 1 and 2 (the neural states), filtering out the 3rd AR(1) SA memory state for visual clarity
  df_states <- data.frame(
    Time = rep(time_steps, 2),
    Activity = c(states[1, ], states[2, ]),
    Lower = c(states[1, ] - 1.96*sqrt(p_var[1,]), states[2, ] - 1.96*sqrt(p_var[2,])),
    Upper = c(states[1, ] + 1.96*sqrt(p_var[1,]), states[2, ] + 1.96*sqrt(p_var[2,])),
    Branch = factor(rep(c("Vagal (x1)", "Sympathetic (x2)"), each = length(time_steps)),
                    levels = c("Vagal (x1)", "Sympathetic (x2)"))
  )

  p2 <- ggplot(df_states, aes(x = Time, y = Activity, color = Branch, fill = Branch)) +
    geom_hline(yintercept = 0, color = "gray80", linewidth = 1) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    scale_color_manual(values = c("Vagal (x1)" = "#378ADD", "Sympathetic (x2)" = "#D85A30")) +
    scale_fill_manual(values = c("Vagal (x1)" = "#378ADD", "Sympathetic (x2)" = "#D85A30")) +
    theme_minimal() + labs(title = "Inferred Neural States (SA Node Memory Filtered)", x = "Beat", y = "Amplitude (Fractional)") +
    theme(legend.position = "bottom")

  margin <- max(0.05, max(abs(states[1:2, ])) * 0.5)
  grid_df <- expand.grid(x1 = seq(min(states[1, ])-margin, max(states[1, ])+margin, length.out = 100),
                         x2 = seq(min(states[2, ])-margin, max(states[2, ])+margin, length.out = 100))
  grid_df$U <- 0.5 * a11 * grid_df$x1^2 + 0.5 * a22 * grid_df$x2^2 + a12 * grid_df$x1 * grid_df$x2
  df_traj <- data.frame(x1 = states[1, ], x2 = states[2, ])

  p3 <- ggplot() +
    geom_contour_filled(data = grid_df, aes(x = x1, y = x2, z = U), bins = 12, alpha = 0.9) +
    scale_fill_viridis_d(option = "magma", direction = -1, guide = "none") +
    geom_path(data = df_traj, aes(x = x1, y = x2), color = "white", alpha = 0.6, linewidth = 0.4) +
    theme_minimal() + coord_fixed(ratio = 1) +
    labs(title = "Phase-Space Energy Map U(x1, x2)", x = "Vagal State (x1)", y = "Sympathetic State (x2)")

  gridExtra::grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 3, 3), c(2, 3, 3)), widths = c(1.2, 1, 1))
}

#' Elegant ggplot2 Diagnostic Dashboard (Log-Linear Fixed)
diagnose_autonomic_model <- function(fit_result, rr_ts) {
  require(ggplot2)
  require(gridExtra)

  # CRITICAL FIX: Transform the data to log-space to match the fitted model
  log_rr_ts <- log(rr_ts)

  opt_theta <- fit_result$par
  kf_out <- run_kalman_filter(opt_theta, log_rr_ts, return_states = TRUE)

  innovations <- kf_out$innovations[-1]

  df_res <- data.frame(Time = 1:length(innovations), Innovation = innovations)

  p1 <- ggplot(df_res, aes(x = Time, y = Innovation)) +
    geom_line(color = "gray40", linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "#D85A30", linetype = "dashed", linewidth = 0.8) +
    theme_minimal() + labs(title = "A. Innovation Time Series", x = "Beat Number", y = "Error (Fractional)")

  p2 <- ggplot(df_res, aes(sample = Innovation)) +
    stat_qq(color = "#378ADD", alpha = 0.6, size = 1.5) +
    stat_qq_line(color = "#D85A30", linetype = "dashed", linewidth = 0.8) +
    theme_minimal() + labs(title = "B. Q-Q Plot", x = "Theoretical", y = "Sample")

  acf_res <- acf(innovations, plot = FALSE)
  df_acf <- data.frame(Lag = acf_res$lag[-1], ACF = acf_res$acf[-1])
  ci_acf <- qnorm((1 + 0.95)/2) / sqrt(length(innovations))

  p3 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
    geom_segment(aes(xend = Lag, yend = 0), color = "gray50", linewidth = 1) +
    geom_point(color = "#378ADD", size = 2) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-ci_acf, ci_acf), color = "#D85A30", linetype = "dashed") +
    theme_minimal() + labs(title = "C. Autocorrelation Function", x = "Lag", y = "ACF")

  gridExtra::grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 1), c(2, 3)))
}

#' Evaluate Goodness-of-Fit Metrics (Log-Linear Fixed)
evaluate_single_fit <- function(fit_result, rr_ts) {
  converged <- ifelse(fit_result$convergence == 0, TRUE, FALSE)
  nll <- fit_result$value

  # AIC correctly scales with the new 8 parameter length
  aic <- 2 * length(fit_result$par) + 2 * nll

  log_rr_ts <- log(rr_ts)
  kf_out <- suppressWarnings(run_kalman_filter(fit_result$par, log_rr_ts, return_states = TRUE))
  innovations <- kf_out$innovations[-1]

  lb_test <- tryCatch(
    Box.test(innovations, lag = 10, type = "Ljung-Box"),
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

  # Pre-allocate a list to store results for speed
  num_records <- length(rr_batch_list)
  results_list <- vector("list", num_records)

  for (i in seq_len(num_records)) {
    id <- subject_ids[i]
    rr_ts <- rr_batch_list[[i]]

    # Use tryCatch to prevent a single bad dataset from stopping the loop
    run_status <- tryCatch({

      # 1. Fit the model (using the function defined previously)
      fit <- fit_autonomic_model(rr_ts)

      # 2. Evaluate the fit
      metrics <- evaluate_single_fit(fit, rr_ts)

      # 3. Store success data
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
      # Fallback for datasets that cause the optimizer or filter to crash
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

    # Optional: Print progress for large batches
    if (i %% 10 == 0) cat(sprintf("Processed %d / %d datasets...\n", i, num_records))
  }

  # Combine the list of single-row data frames into one large table
  final_results_df <- do.call(rbind, results_list)
  return(final_results_df)
}
