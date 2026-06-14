
# Prepare workspace -------------------------------------------------------

## Import libraries
library(data.table)
library(RHRV)
library(reticulate)
library(CardioCurveR)


# Importing Fantasia dataset ----------------------------------------------

subjects <- list.files("data-raw/fantasia/") |>
  gsub(pattern = "\\.hea|\\.dat|\\.ecg|-", replacement = "") |>
  unique()

for (id in subjects) {
  rr_intervals <- RHRV::LoadBeatWFDB(
    HRVData = RHRV::CreateHRVData(),
    RecordName = paste0("data-raw/fantasia/",id),
    annotator = "ecg"
  )$Beat$Time |>
    diff() |> CardioCurveR::clean_outlier(threshold = 5)

cat(rr_intervals, sep = "\n", file = paste0("data/fantasia/",id,".txt"))
}

plot(rr_intervals, type = "o", pch = "•")


# Importing PRCP dataset ---------------------------------------------------

subjects <- list.files("data-raw/prcp/") |>
  gsub(pattern = "\\.hea|\\.dat|\\.anI|\\.wqrs|\\.wabp|SHA256SUMS\\.txt|-", replacement = "") |>
  unique()

subjects <- subjects[nchar(subjects) > 0]

py_require("wfdb")
wfdb <- import("wfdb")

for (id in subjects) {
  rr_intervals <- wfdb$processing$ann2rr(
    record_name = paste0("data-raw/prcp/",id),
    extension = "wqrs",
    format = "s",
    as_array = TRUE
  )

  cat(rr_intervals, sep = "\n", file = paste0("data/prcp/",id,".txt"))
}

plot(rr_intervals, type = "o", pch = "•")


# Importing TMST dataset --------------------------------------------------

subjects <- list.files("data-raw/tmst/hrv_data/") |>
  unique()

for (id in subjects) {
  rr_intervals <- import_RRi_txt(file = paste0("data-raw/tmst/hrv_data/",id), F, F, warn = FALSE)$RRi
  rr_intervals <- clean_outlier(rr_intervals, threshold = 3)
  rr_intervals <- filter_signal(rr_intervals, n = 1, W = 0.8)
  rr_intervals <- rr_intervals[!is.na(rr_intervals)]/1000
  cat(rr_intervals, sep = "\n", file = paste0("data/tmst/",id))
}

plot(rr_intervals, type = "o", pch = "•")
