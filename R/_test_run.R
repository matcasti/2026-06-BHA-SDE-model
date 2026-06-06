
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions.R")

## Load test data
rr_data <- readLines("data-raw/rri-jabf.txt") |>
  as.numeric() |>
  CardioCurveR::clean_outlier()

fit <- fit_autonomic_model(rr_ts = rr_data)

get_model_ci(fit)

visualize_autonomic_model(fit_result = fit, rr_ts = rr_data)

diagnose_autonomic_model(fit, rr_ts = rr_data)

evaluate_single_fit(fit, rr_ts = rr_data)

report_fit(fit_result = fit)
