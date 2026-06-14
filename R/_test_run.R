
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions.R")

## Load test data
rr_data <- readLines("data/rri-jabf.txt") |>
  as.numeric()

plot(rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data)

visualize_model(fit)

diagnose_model(fit)

report_model(fit)
