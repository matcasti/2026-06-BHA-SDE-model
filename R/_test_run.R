
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions.R")

## Load test data
rr_data <- (readLines("data/rri-jabf.txt") |>
  as.numeric() |>
  CardioCurveR::clean_outlier(threshold = 3))/1000

plot(rr_data, type = "o", pch = 16, cex = 0.5)

fit <- fit_model(rr_data)

visualize_model(fit)

diagnose_model(fit)

report_model(fit)
