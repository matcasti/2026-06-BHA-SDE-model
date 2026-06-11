
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions_v3.R")

## Load test data
rr_data <- readLines("data-raw/rri-prcp.txt")[2:1650] |>
  as.numeric() |>
  CardioCurveR::clean_outlier(threshold = 3)

rr_data <- (readLines("data-raw/rri-jabf.txt") |>
  as.numeric() |>
  CardioCurveR::clean_outlier(threshold = 3))/1000

rr_data <- (readLines("data-raw/rri-maca-hiit.txt") |>
  as.numeric() |>
  CardioCurveR::clean_outlier(threshold = 3))/1000

rr_data <- (read.csv("data-raw/rri-rcmt.csv", comment.char = "#")[,2] |>
  as.numeric() |>
  CardioCurveR::clean_outlier(threshold = 3))/1000

plot(rr_data, type = "o", pch = 16, cex = 0.5)

fit <- fit_model(rr_data)

visualize_model(fit)

diagnose_model(fit)

report_model(fit)
