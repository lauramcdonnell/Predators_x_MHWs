# User Details ====================================================================
# 01_detect_marine_heatwaves.R

# Purpose
#   Detect daily marine heatwaves (MHWs) from gridded daily SST and produce
#   one binary MHW occurrence file per day.
#
# Method
#   1. Calculate a grid-cell-specific, seasonally varying 90th-percentile SST
#      threshold from an 11-day centered window during a 2002-2022 baseline.
#   2. Smooth the 365 daily threshold fields with a 31-day circular running mean.
#   3. Define threshold exceedance as SST > the smoothed day-of-year threshold.
#   4. At each grid cell, join exceedance days into candidate events when no more
#      than two consecutive sub-threshold days occur between exceedance days.
#   5. Retain candidate events containing at least five ACTUAL threshold-
#      exceedance days. Bridged gap days are part of the retained event, but do
#      not count toward the five exceedance days.
#   6. Write one binary MHW occurrence file per day.
#
#   The binary files are the only MHW product passed to script 02. Script 02
#   calculates regional MHW extent directly from these files.
#
# Input assumptions
#   - One NetCDF SST file per day.
#   - All files are already cropped/masked to the same study region and grid.
#   - File names contain YYYY-MM-DD or YYYY_MM_DD.
#   - The SST variable is named "analysed_sst" by default.
#
# Packages
#   install.packages(c("terra", "dplyr", "stringr"))
#
# Notes
#   - This script uses raw SST, not detrended SST.
#   - The baseline defaults below reproduce the regional analysis script
#     (June 2002 through June 2022).
#   - February 29 is excluded when estimating the 365-day climatology. Its
#     threshold during detection is interpolated as the mean of Feb 28 and Mar 1.
#   - For very large 1-km domains, terra will use disk-backed processing, but the
#     custom time-series event step can still be computationally intensive.

library(terra)
library(dplyr)
library(stringr)

# 1. USER SETTINGS =============================================================

# Folder containing daily, region-masked SST NetCDF files.
sst_dir <- "PATH/TO/DAILY_SST"

# Folder for climatology, intermediate files, and daily binary MHW files.
out_dir <- "PATH/TO/OUTPUT"

# SST variable inside each NetCDF.
sst_var <- "analysed_sst"

# Date range used to estimate the climatological threshold.
# Exact climatology interval used in the regional analysis.
baseline_start <- as.Date("2002-06-01")
baseline_end   <- as.Date("2022-06-01")

# Study period for MHW detection.
analysis_start <- as.Date("2002-06-01")
analysis_end   <- as.Date("2024-12-31")

# Threshold/event parameters.
threshold_prob        <- 0.90
climatology_window    <- 11L  # centered 11-day window: +/- 5 days
smoothing_window      <- 31L  # centered 31-day circular running mean
min_exceedance_days   <- 5L   # actual SST-threshold exceedance days required
max_gap_days          <- 2L   # sub-threshold days allowed BETWEEN exceedances

# Number of CPU cores passed to terra::app().
# Increase cautiously on shared systems.
n_cores <- 1L

# Optional terra temporary directory.
tmp_dir <- file.path(tempdir(), "mhw_terra")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir, progress = 1, memfrac = 0.8)

# Output folders.
clim_raw_dir    <- file.path(out_dir, "climatology_p90_raw")
clim_smooth_dir <- file.path(out_dir, "climatology_p90_smoothed")
mhw_binary_dir  <- file.path(out_dir, "mhw_binary")

dir.create(out_dir,         recursive = TRUE, showWarnings = FALSE)
dir.create(clim_raw_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(clim_smooth_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mhw_binary_dir,  recursive = TRUE, showWarnings = FALSE)

# 2. HELPER FUNCTIONS ===============================================================


extract_date <- function(x) {
  z <- stringr::str_extract(basename(x), "\\d{4}[-_]\\d{2}[-_]\\d{2}")
  as.Date(gsub("_", "-", z))
}

# Map dates to a 365-day calendar without shifting dates after Feb 29.
# Feb 29 itself is NA here and is handled separately during detection.
calendar_day_365 <- function(x) {
  md <- format(as.Date(x), "%m-%d")
  out <- rep(NA_integer_, length(md))
  keep <- md != "02-29"
  ref_dates <- as.Date(paste0("2001-", md[keep]))  # 2001 is not a leap year
  out[keep] <- as.integer(format(ref_dates, "%j"))
  out
}

circular_window <- function(day, half_width, n_days = 365L) {
  z <- (day - half_width):(day + half_width)
  ((z - 1L) %% n_days) + 1L
}

open_sst <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  terra::rast(paste0("NETCDF:", path, ":", sst_var))
}

open_sst_stack <- function(paths) {
  src <- paste0(
    "NETCDF:",
    normalizePath(paths, winslash = "/", mustWork = TRUE),
    ":", sst_var
  )
  r <- terra::rast(src)
  if (terra::nlyr(r) != length(paths)) {
    stop("Expected one SST layer per input file. Check sst_var and NetCDF structure.")
  }
  r
}

# Detect retained MHW events in one grid-cell time series.
# Candidate events are defined by groups of exceedance dates separated by no more
# than max_gap_days intervening sub-threshold days. A candidate is retained only
# if it contains >= min_exceedance_days TRUE exceedances.
detect_mhw_series <- function(x,
                              min_exceedance_days = 5L,
                              max_gap_days = 2L) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))

  exceed <- which(!is.na(x) & x > 0)
  out <- integer(length(x))

  if (length(exceed) == 0L) return(out)

  # A difference of max_gap_days + 1 still belongs to the same event.
  # Example: max_gap_days = 2 allows exceedances on days 1 and 4.
  grp <- cumsum(c(TRUE, diff(exceed) > (max_gap_days + 1L)))
  groups <- split(exceed, grp)

  for (g in groups) {
    if (length(g) >= min_exceedance_days) {
      out[min(g):max(g)] <- 1L
    }
  }

  out
}

# Matrix-aware wrapper for terra::app(). Rows are cells, columns are dates.
detect_mhw_app <- function(x) {
  if (is.null(dim(x))) {
    return(detect_mhw_series(
      x,
      min_exceedance_days = min_exceedance_days,
      max_gap_days = max_gap_days
    ))
  }

  ans <- t(apply(
    x,
    1,
    detect_mhw_series,
    min_exceedance_days = min_exceedance_days,
    max_gap_days = max_gap_days
  ))

  storage.mode(ans) <- "integer"
  ans
}

# 3. INDEX DAILY SST FILES ==============================================================

sst_files <- list.files(
  sst_dir,
  pattern = "\\.nc$",
  full.names = TRUE
)

sst_index <- tibble(
  file = sst_files,
  date = extract_date(sst_files)
) %>%
  filter(!is.na(date)) %>%
  arrange(date)

if (nrow(sst_index) == 0L) {
  stop("No dated NetCDF files found in sst_dir.")
}

if (anyDuplicated(sst_index$date)) {
  stop("More than one SST file was found for at least one date.")
}

analysis_index <- sst_index %>%
  filter(date >= analysis_start, date <= analysis_end)

expected_dates <- seq.Date(analysis_start, analysis_end, by = "day")
missing_dates <- setdiff(expected_dates, analysis_index$date)

if (length(missing_dates) > 0L) {
  stop(
    "Daily SST series is incomplete. Missing ", length(missing_dates),
    " dates in the requested analysis period."
  )
}

baseline_index <- sst_index %>%
  filter(date >= baseline_start, date <= baseline_end) %>%
  mutate(cal_day = calendar_day_365(date)) %>%
  filter(!is.na(cal_day))  # exclude Feb 29 from the 365-day climatology

if (nrow(baseline_index) == 0L) {
  stop("No SST files fall within the requested climatology baseline.")
}

# Confirm common geometry using the first and last SST files.
r_first <- open_sst(analysis_index$file[1])
r_last  <- open_sst(analysis_index$file[nrow(analysis_index)])
if (!terra::compareGeom(r_first, r_last, stopOnError = FALSE)) {
  stop("SST files do not share a common raster geometry.")
}
rm(r_first, r_last)
gc()

# 4. GRID-CELL 90TH-PERCENTILE CLIMATOLOGY ======================================================

message("Calculating raw 90th-percentile climatology...")

half_clim <- (climatology_window - 1L) / 2L
if (climatology_window %% 2L != 1L) stop("climatology_window must be odd.")

p_fun <- function(x) {
  stats::quantile(
    x,
    probs = threshold_prob,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
}

raw_clim_files <- character(365)

for (d in seq_len(365)) {
  out_file <- file.path(clim_raw_dir, sprintf("p90_raw_day%03d.tif", d))
  raw_clim_files[d] <- out_file

  if (file.exists(out_file)) next

  days_use <- circular_window(d, half_width = half_clim)
  use <- baseline_index %>% filter(cal_day %in% days_use)

  if (nrow(use) == 0L) stop("No baseline SST observations for climatology day ", d)

  s <- open_sst_stack(use$file)

  thr <- terra::app(
    s,
    fun = p_fun,
    filename = out_file,
    overwrite = TRUE
  )

  rm(s, thr)
  gc()
}

# 5. 31-DAY CIRCULAR SMOOTHING OF THE CLIMATOLOGICAL THRESHOLD =========================================================

message("Smoothing climatological threshold with a 31-day circular mean...")

half_smooth <- (smoothing_window - 1L) / 2L
if (smoothing_window %% 2L != 1L) stop("smoothing_window must be odd.")

smooth_clim_files <- character(365)

for (d in seq_len(365)) {
  out_file <- file.path(clim_smooth_dir, sprintf("p90_smoothed_day%03d.tif", d))
  smooth_clim_files[d] <- out_file

  if (file.exists(out_file)) next

  days_use <- circular_window(d, half_width = half_smooth)
  s <- terra::rast(raw_clim_files[days_use])

  sm <- terra::mean(
    s,
    na.rm = TRUE,
    filename = out_file,
    overwrite = TRUE
  )

  rm(s, sm)
  gc()
}

clim365 <- terra::rast(smooth_clim_files)
names(clim365) <- sprintf("day_%03d", seq_len(365))

# 6. BUILD SST AND MATCHED THRESHOLD TIME STACKS ======================================================

message("Building SST and day-matched threshold stacks...")

dates <- analysis_index$date
sst <- open_sst_stack(analysis_index$file)
names(sst) <- format(dates, "%Y-%m-%d")

# Day-matched threshold. Feb 29 uses the mean of Feb 28 (day 59) and Mar 1 (day 60).
threshold_layers <- lapply(dates, function(d) {
  if (format(d, "%m-%d") == "02-29") {
    terra::mean(clim365[[c(59, 60)]], na.rm = TRUE)
  } else {
    clim365[[calendar_day_365(d)]]
  }
})

threshold <- terra::rast(threshold_layers)
names(threshold) <- names(sst)

if (!terra::compareGeom(sst[[1]], threshold[[1]], stopOnError = FALSE)) {
  stop("SST and climatological threshold grids do not align.")
}

# 7. THRESHOLD EXCEEDANCE, MHW EVENT DETECTION, AND BINARY OUTPUT =====================================================


message("Identifying daily threshold exceedances...")
exceed <- sst > threshold

# terra::app() evaluates each grid-cell time series across the full analysis
# period. The resulting stack contains one binary layer per day.
#
# Values:
#   1  = grid cell belongs to a qualifying MHW event
#   0  = grid cell does not belong to a qualifying MHW event
#   NA = land / cells without an SST time series
#
# Allowed sub-threshold gap days are retained within a qualifying event span.
# They allow the event to continue but do not count toward the five required
# threshold-exceedance days.

mhw_flag_stack_file <- file.path(tmp_dir, "mhw_flag_stack.tif")

message("Applying the >=5 exceedance-day / <=2 gap-day MHW definition...")
mhw_flag <- terra::app(
  exceed,
  fun = detect_mhw_app,
  cores = n_cores,
  filename = mhw_flag_stack_file,
  overwrite = TRUE,
  datatype = "INT1U"
)

names(mhw_flag) <- format(dates, "%Y-%m-%d")

message("Writing daily binary MHW files...")

for (i in seq_along(dates)) {
  d <- dates[i]
  message(i, "/", length(dates), "  ", d)

  flag_i <- mhw_flag[[i]]
  names(flag_i) <- "mhw_flag"

  out_file <- file.path(
    mhw_binary_dir,
    sprintf("MHW_binary_%s.nc", format(d, "%Y_%m_%d"))
  )

  terra::writeCDF(
    flag_i,
    filename = out_file,
    overwrite = TRUE
  )
}

message("Done.")
message("Daily binary MHW files: ", mhw_binary_dir)
