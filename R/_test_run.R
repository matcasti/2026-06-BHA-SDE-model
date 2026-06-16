
# Prepare workspace -------------------------------------------------------

## Load function
source("R/_functions.R")

# -------------------------------------------------------------------------

## Load test data
rr_data <- (readLines("data/prcp/12726.txt") |> as.numeric())

plot(rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data, jump_threshold = 0.1, jump_power = 12)

fig <- visualize_model(fit) / diagnose_model(fit) +
  patchwork::plot_layout(heights = c(1,0.75))

ggsave(filename = "manuscript/_example_fig.png", plot = fig,
       device = "png", scale = 6, height = 500, width = 500, units = "px")

report_model(fit)

# -------------------------------------------------------------------------

## Load data
rr_data <- (readLines("data/rri-jabf.txt") |> as.numeric())/1000

plot(rr_data, type = "o", pch = "•")

fit <- fit_model(rr_data, jump_threshold = 0.1, jump_power = 12)

fig <- visualize_model(fit) / diagnose_model(fit) +
  patchwork::plot_layout(heights = c(1,0.75))

ggsave(filename = "manuscript/_example_fig2.png", plot = fig,
       device = "png", scale = 6, height = 500, width = 500, units = "px")

report_model(fit)
