# User Details =====================================================
# 02_build_event_control_windows.R

# Purpose
#   Calculate daily regional marine heatwave (MHW) extent from the binary files
#   produced by 01_detect_marine_heatwaves.R and build the event-control catalog
#   used for matched habitat comparisons.
#
# Method implemented here
#   1. Read daily binary MHW occurrence files from script 01.
#   2. Calculate regional MHW extent (frac_mhw) as the proportion of ocean pixels
#      classified as MHW on each day.
#   3. Define qualifying event days as days at or above the regional 70th
#      percentile of frac_mhw.
#   4. Group qualifying event days into an event window when no more than ONE
#      intervening calendar day separates successive qualifying event days.
#   5. Retain event windows spanning at least seven calendar days. This is a
#      separate filter from the >=5 exceedance-day criterion used to DETECT MHWs
#      at individual grid cells in script 01.
#   6. Define control candidates as days at or below the regional 30th percentile
#      of frac_mhw.
#   7. For each retained event window, select control observations from the same
#      day-of-year range in other years that satisfy the control criterion.
#   8. Exclude event windows without corresponding control observations.
#   9. Retain only days of year represented in at least three years under BOTH
#      event and control conditions.
#
# Input
#   Daily files written by script 01, named:
#     MHW_binary_YYYY_MM_DD.nc
#   Each file must contain a binary layer named "mhw_flag".
#
# Output
#   - regional_mhw_extent.csv: daily total pixels, MHW pixels, and frac_mhw
#   - event_windows.csv: one row per retained regional event window
#   - event_control_lookup.csv: one row per event/control date and window ID
#   - event_control_thresholds.csv: q70 and q30 regional extent thresholds
#   - event_control_doy_year_counts.csv: DOY representation diagnostics
#
# Packages
#   install.packages(c("terra", "dplyr", "tidyr", "lubridate", "readr", "stringr"))

library(terra)
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)

# 1. USER SETTINGS ==========================================================


# Folder containing daily binary MHW files produced by script 01.
mhw_binary_dir <- "PATH/TO/OUTPUT/mhw_binary"

# Output folder.
out_dir <- "PATH/TO/EVENT_CONTROL_OUTPUT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Name of the binary MHW variable inside each NetCDF.
mhw_var <- "mhw_flag"

# Event-control thresholds.
event_quantile   <- 0.70
control_quantile <- 0.30

# Regional event-window rules.
max_intervening_days <- 1L  # one non-event day may occur between event days
min_window_span_days <- 7L  # retained windows must span >= 7 calendar days

# Minimum number of years representing a calendar DOY in BOTH groups.
min_years_per_doy <- 3L

#  2. INDEX DAILY BINARY MHW FILES ========================================================

extract_date <- function(x) {
  z <- stringr::str_extract(
    basename(x),
    "\\d{4}[-_]\\d{2}[-_]\\d{2}"
  )
  as.Date(gsub("_", "-", z))
}

mhw_files <- list.files(
  mhw_binary_dir,
  pattern = "^MHW_binary_\\d{4}_\\d{2}_\\d{2}\\.nc$",
  full.names = TRUE
)

mhw_index <- tibble(
  file = mhw_files,
  date = extract_date(mhw_files)
) %>%
  filter(!is.na(date)) %>%
  arrange(date)

if (nrow(mhw_index) == 0L) {
  stop("No daily binary MHW files were found in mhw_binary_dir.")
}

if (anyDuplicated(mhw_index$date)) {
  stop("More than one binary MHW file was found for at least one date.")
}

# 3. CALCULATE DAILY REGIONAL MHW EXTENT ================================================

# Regional extent is the proportion of valid ocean pixels classified as MHW.
# Script 01 writes land/non-ocean cells as NA, so total_pixels is simply the
# number of non-NA cells in each daily binary file.

summarize_binary_file <- function(path, date) {
  r <- try(terra::rast(paste0("NETCDF:", normalizePath(path, winslash = "/", mustWork = TRUE), ":", mhw_var)), silent = TRUE)

  # Some NetCDF readers expose a single-layer file without requiring the subdataset
  # syntax. Fall back to a normal rast() call if needed.
  if (inherits(r, "try-error")) {
    r <- terra::rast(path)
    if (terra::nlyr(r) > 1L && mhw_var %in% names(r)) r <- r[[mhw_var]]
    if (terra::nlyr(r) > 1L) r <- r[[1]]
  }

  vals <- terra::values(r, mat = FALSE)

  total_pixels <- sum(!is.na(vals))
  n_mhw_any    <- sum(vals == 1, na.rm = TRUE)

  if (total_pixels == 0L) {
    stop("No valid ocean pixels found in: ", basename(path))
  }

  tibble(
    date = as.Date(date),
    total_pixels = total_pixels,
    n_mhw_any = n_mhw_any,
    frac_mhw = n_mhw_any / total_pixels
  )
}

message("Calculating daily regional MHW extent...")

regional_list <- vector("list", nrow(mhw_index))
for (i in seq_len(nrow(mhw_index))) {
  message(i, "/", nrow(mhw_index), "  ", mhw_index$date[i])
  regional_list[[i]] <- summarize_binary_file(
    mhw_index$file[i],
    mhw_index$date[i]
  )
}

mhw_daily <- bind_rows(regional_list) %>%
  arrange(date) %>%
  mutate(
    year = year(date),
    doy  = yday(date),
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                           ~ "Fall"
    )
  )

write_csv(
  mhw_daily,
  file.path(out_dir, "regional_mhw_extent.csv")
)
# 4. REGIONAL 70TH / 30TH PERCENTILE THRESHOLDS ==============================================================================


q70 <- unname(quantile(mhw_daily$frac_mhw, event_quantile, na.rm = TRUE))
q30 <- unname(quantile(mhw_daily$frac_mhw, control_quantile, na.rm = TRUE))

write_csv(
  tibble(
    metric = "frac_mhw",
    event_quantile = event_quantile,
    event_threshold = q70,
    control_quantile = control_quantile,
    control_threshold = q30
  ),
  file.path(out_dir, "event_control_thresholds.csv")
)

# 5. BUILD REGIONAL EVENT WINDOWS===============================================

qualifying_event_days <- mhw_daily %>%
  filter(frac_mhw >= q70) %>%
  arrange(date)

if (nrow(qualifying_event_days) == 0L) {
  stop("No days satisfy the regional event criterion.")
}

# With max_intervening_days = 1, dates separated by 1 or 2 days remain in the
# same event window. A new event begins when the gap between qualifying event
# dates is > 2 days.
event_windows <- qualifying_event_days %>%
  mutate(
    date_gap = as.numeric(date - lag(date)),
    new_window = is.na(date_gap) |
      date_gap > (max_intervening_days + 1L),
    block = cumsum(new_window)
  ) %>%
  group_by(block) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    n_qualifying_days = n_distinct(date),
    window_span_days = as.integer(end_date - start_date) + 1L,
    .groups = "drop"
  ) %>%
  filter(window_span_days >= min_window_span_days) %>%
  arrange(start_date) %>%
  mutate(window_id = row_number()) %>%
  select(
    window_id,
    start_date,
    end_date,
    window_span_days,
    n_qualifying_days
  )

if (nrow(event_windows) == 0L) {
  stop("No regional event windows satisfy the minimum seven-day span.")
}

# Expand each event window to every calendar day between its first and last
# qualifying event day. This includes the allowed single intervening day.
event_lookup <- event_windows %>%
  rowwise() %>%
  mutate(date = list(seq.Date(start_date, end_date, by = "day"))) %>%
  ungroup() %>%
  unnest(date) %>%
  mutate(window_type = "Event") %>%
  select(date, window_id, window_type)

# 6. MATCH CONTROLS BY DAY-OF-YEAR==============================================

control_candidates <- mhw_daily %>%
  filter(frac_mhw <= q30) %>%
  mutate(
    year = year(date),
    doy = yday(date)
  )

# For each event window, draw low-MHW observations from the corresponding DOY
# range in OTHER years. Cross-year event windows are handled by the OR condition
# for wrapped DOY ranges.
control_lookup <- event_windows %>%
  mutate(
    start_doy = yday(start_date),
    end_doy   = yday(end_date)
  ) %>%
  rowwise() %>%
  mutate(
    event_years = list(unique(year(seq.Date(start_date, end_date, by = "day")))),
    control_dates = list({
      sy <- start_doy
      ey <- end_doy
      yrs_exclude <- unlist(event_years)

      z <- control_candidates %>%
        filter(!year %in% yrs_exclude)

      if (ey >= sy) {
        z <- z %>% filter(doy >= sy, doy <= ey)
      } else {
        z <- z %>% filter(doy >= sy | doy <= ey)
      }

      z %>% pull(date)
    })
  ) %>%
  ungroup() %>%
  select(window_id, control_dates) %>%
  unnest(control_dates) %>%
  transmute(
    date = as.Date(control_dates),
    window_id = window_id,
    window_type = "Control"
  )

# Remove event windows that have no corresponding controls.
valid_event_ids <- intersect(
  unique(event_lookup$window_id),
  unique(control_lookup$window_id)
)

event_lookup <- event_lookup %>%
  filter(window_id %in% valid_event_ids)

control_lookup <- control_lookup %>%
  filter(window_id %in% valid_event_ids)

event_windows <- event_windows %>%
  filter(window_id %in% valid_event_ids)

# 7. REQUIRE >=3 YEARS OF BOTH EVENT AND CONTROL REPRESENTATION PER DOY =========

contrast_lookup <- bind_rows(event_lookup, control_lookup) %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    doy = yday(date),
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Winter",
      month(date) %in% c(3, 4, 5)  ~ "Spring",
      month(date) %in% c(6, 7, 8)  ~ "Summer",
      TRUE                           ~ "Fall"
    )
  ) %>%
  arrange(window_id, window_type, date)

doy_year_counts <- contrast_lookup %>%
  distinct(window_type, doy, year) %>%
  count(window_type, doy, name = "n_years") %>%
  pivot_wider(
    names_from = window_type,
    values_from = n_years,
    values_fill = 0
  )

# Ensure both expected columns exist even in unusual datasets.
if (!"Event" %in% names(doy_year_counts)) doy_year_counts$Event <- 0L
if (!"Control" %in% names(doy_year_counts)) doy_year_counts$Control <- 0L

valid_doys <- doy_year_counts %>%
  filter(
    Event >= min_years_per_doy,
    Control >= min_years_per_doy
  ) %>%
  pull(doy)

contrast_lookup_filtered <- contrast_lookup %>%
  filter(doy %in% valid_doys) %>%
  select(date, window_id, window_type, doy, season)

# Retain only window IDs still represented after the >=3-year DOY filter.
final_window_ids <- unique(contrast_lookup_filtered$window_id)
event_windows_final <- event_windows %>%
  filter(window_id %in% final_window_ids)

# 8. WRITE OUTPUTS==============================================================


write_csv(
  event_windows_final,
  file.path(out_dir, "event_windows.csv")
)

write_csv(
  contrast_lookup_filtered,
  file.path(out_dir, "event_control_lookup.csv")
)

write_csv(
  doy_year_counts,
  file.path(out_dir, "event_control_doy_year_counts.csv")
)

message("Done.")
message("Regional MHW extent: ", file.path(out_dir, "regional_mhw_extent.csv"))
message("Event windows: ", file.path(out_dir, "event_windows.csv"))
message("Event-control lookup: ", file.path(out_dir, "event_control_lookup.csv"))
