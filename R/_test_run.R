
# Prepare workspace -------------------------------------------------------

## Load function
source("R/01_model_engine.R")

# -------------------------------------------------------------------------

## Load test data
rr_data <- (readLines("data/prcp/12734.txt") |> as.numeric())

plot(rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data)

fig <- (visualize_model(fit) | diagnose_model(fit)) +
  patchwork::plot_layout(guides = "keep", widths = c(1,1,0.75))

ggsave(filename = "manuscript/_example_fig.png", plot = fig,
       device = "png", scale = 6, height = 500, width = 500, units = "px")

report_model(fit)

# -------------------------------------------------------------------------

## Load data
rr_data <- (readLines("data/rri-jabf.txt") |> as.numeric())

rr_data <- (readLines("data/tmst/007.txt") |> as.numeric())
rr_data <- (readLines("data/tmst/008.txt") |> as.numeric())

plot(cumsum(rr_data)/60, rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data)

fig <- (visualize_model(fit) | diagnose_model(fit)) +
  patchwork::plot_layout(guides = "keep", widths = c(1,1,0.75))

ggsave(filename = "manuscript/_example_fig2.png", plot = fig,
       device = "png", scale = 6, height = 500, width = 500, units = "px")

report_model(fit)
