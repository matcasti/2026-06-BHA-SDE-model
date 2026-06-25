# ==============================================================================
# This script executes the high-performance batch optimization across the three
# clinical cohorts: Fantasia (Baseline), Orthostatic Tilt (Transient), and
# Exercise Recovery (Metabolic). It extracts the fit statistics, parameters,
# and PACF matrices, exporting them for visualization.
# ==============================================================================

# Prepare Workspace ------------------------------------------------------------

# Source the core mathematical engine (contains zero execution code, only functions)
source("R/01_model_engine.R")

# 2. Data Loading --------------------------------------------------------------

# --- COHORT 1: Fantasia (Resting Baseline) ---
files <- list.files("data/fantasia")
rr_fantasia <- lapply(files, function(x) {
  readLines(paste0("data/fantasia/",x)) |> as.numeric()
})
names(rr_fantasia) <- gsub("\\.txt", "", files)


# --- COHORT 2: Postural Response (Orthostatic Tilt) ---
files <- list.files("data/prcp")
rr_tilt <- lapply(files, function(x) {
  readLines(paste0("data/prcp/",x)) |> as.numeric()
})
names(rr_tilt) <- gsub("\\.txt", "", files)

# --- COHORT 3: Exercise Recovery (Metabolic Stress) ---
files <- list.files("data/tmst")
rr_exercise <- lapply(files, function(x) {
  readLines(paste0("data/tmst/",x)) |> as.numeric()
})
names(rr_exercise) <- gsub("\\.txt", "", files)

# 3. Execute Batch Processing --------------------------------------------------

# We use the updated batch_process_hrv() which calculates RMSE, MAPE, and PACF.

# A. Process Fantasia
res_fantasia <- batch_process_hrv(
  rr_list = rr_fantasia,
  dataset_tag = "Fantasia_Baseline"
)
report_batch_metrics(res_fantasia)

# B. Process Orthostatic Tilt
res_tilt <- batch_process_hrv(
  rr_list = rr_tilt,
  dataset_tag = "Orthostatic_Tilt"
)
report_batch_metrics(res_tilt)

# C. Process Exercise Recovery
res_exercise <- batch_process_hrv(
  rr_list = rr_exercise,
  dataset_tag = "Exercise_Recovery"
)
report_batch_metrics(res_exercise)


# 4. Global Aggregation --------------------------------------------------------

# A. Population Summary (Parameters, Fit Statistics, KS, RMSE, MAPE)
df_population_summary <- rbind(
  res_fantasia$summary,
  res_tilt$summary,
  res_exercise$summary
)

# B. Population PACF Matrix (For Figure 1B)
mat_population_pacf <- rbind(
  res_fantasia$pacf_matrix,
  res_tilt$pacf_matrix,
  res_exercise$pacf_matrix
)

# C. Population States (Optional: for extracting specific examples later)
df_population_states <- rbind(
  res_fantasia$states,
  res_tilt$states,
  res_exercise$states
)

# 5. Export Results for Visualization Scripts ----------------------------------

## Exporting highly compressed datasets to 'results/'...
saveRDS(df_population_summary, file = "results/df_population_summary.rds")
saveRDS(mat_population_pacf, file = "results/mat_population_pacf.rds")
saveRDS(df_population_states, file = "results/df_population_states.rds")

# Save the raw models (useful if you need to run diagnose_model() on a specific subject)
raw_models_all <- c(res_fantasia$models, res_tilt$models, res_exercise$models)
saveRDS(raw_models_all, file = "results/list_raw_models.rds")
