
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions.R")

## Load test data
prcp_data <-
  lapply(
    X = list.files("data/prcp", pattern = "\\.txt", full.names = TRUE),
    FUN = function(x) as.numeric(x = readLines(x))
  )

prcp_batch <- batch_process_hrv(prcp_data, dataset_tag = "PRCP dataset")

report_batch_metrics(prcp_batch)

plot(rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data, jump_power = 0.1)

visualize_model(fit)

diagnose_model(fit)

report_model(fit)
