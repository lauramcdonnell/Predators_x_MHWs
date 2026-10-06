# 03_calculate_core_habitat_area.R

# Purpose
#   Convert daily species distribution model (SDM) predictions to binary core
#   habitat, calculate daily core-habitat area, and link habitat estimates to the
#   MHW event-control windows produced by 02_build_event_control_windows.R.
#
# Inputs
#   - Daily SDM rasters, organized by species.
#   - event_control_lookup.csv from Script 02.
#   - Species-specific habitat thresholds.
#   - Retained occupancy months for each species.
#   - Daily binary MHW files from Script 01, used only to define regional extent.
#
# Outputs
#   - daily_core_habitat_area.csv
#   - event_control_habitat_area.csv
#
# Packages
#   install.packages(c("terra", "dplyr", "tidyr", "readr", "stringr", "lubridate", "purrr"))

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(purrr)

# 1. USER SETTINGS =============================================================

# Root folder containing daily SDMs.
# Expected structure by default:
#   SDM_ROOT/species/glorys/species_YYYYMMDD.grd
sdm_root <- "PATH/TO/DAILY_SDM_FILES"
sdm_subdir <- "glorys"
sdm_extension <- "grd"

# Daily binary MHW files from Script 01.
mhw_binary_dir <- "PATH/TO/OUTPUT/mhw_binary"

# Event-control lookup from Script 02.
event_control_file <- "PATH/TO/EVENT_CONTROL_OUTPUT/event_control_lookup.csv"

# CSV with one threshold per species.
# Required columns: species, threshold
threshold_file <- "PATH/TO/sdm_thresholds.csv"

# CSV listing retained occupancy months.
# Required columns: species, month
# Include one row for every species-month combination retained for analysis.
occupancy_file <- "PATH/TO/species_occupancy_months.csv"

# Output folder.
out_dir <- "PATH/TO/HABITAT_OUTPUT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Species to process. Set to NULL to use all species in threshold_file.
species_vec <- NULL

# Optional terra temporary directory.
tmp_dir <- file.path(tempdir(), "sdm_terra")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = tmp_dir, progress = 1, memfrac = 0.8)

# 2. READ LOOKUP TABLES =========================================================

event_control <- read_csv(event_control_file, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    month = month(date)
  )

thresholds <- read_csv(threshold_file, show_col_types = FALSE)
occupancy <- read_csv(occupancy_file, show_col_types = FALSE)

stopifnot(all(c("date", "window_id", "window_type", "doy", "season") %in%
                names(event_control)))
stopifnot(all(c("species", "threshold") %in% names(thresholds)))
stopifnot(all(c("species", "month") %in% names(occupancy)))

if (is.null(species_vec)) {
  species_vec <- thresholds %>% distinct(species) %>% pull(species)
}

# 3. INDEX DAILY MHW FILES ======================================================

extract_mhw_date <- function(x) {
  z <- str_extract(basename(x), "\\d{4}[-_]\\d{2}[-_]\\d{2}")
  as.Date(gsub("_", "-", z))
}

mhw_files <- list.files(
  mhw_binary_dir,
  pattern = "^MHW_binary_\\d{4}_\\d{2}_\\d{2}\\.nc$",
  full.names = TRUE
)

mhw_index <- tibble(
  file = mhw_files,
  date = extract_mhw_date(mhw_files)
) %>%
  filter(!is.na(date)) %>%
  distinct(date, .keep_all = TRUE)

if (nrow(mhw_index) == 0L) {
  stop("No daily binary MHW files were found in mhw_binary_dir.")
}

# 4. HELPER FUNCTIONS ===========================================================

extract_sdm_date <- function(x) {
  z <- str_extract(basename(x), "\\d{8}")
  as.Date(z, format = "%Y%m%d")
}

# Calculate daily core-habitat area from one SDM raster.
# A day with no cells above threshold is retained as 0 km2.
core_area_km2 <- function(sdm_file, region_file, threshold) {
  sdm_r <- terra::rast(sdm_file)
  region_r <- terra::rast(region_file)[[1]]

  r <- terra::crop(sdm_r, terra::ext(region_r))
  if (terra::ncell(r) == 0L) return(NA_real_)

  vals <- terra::values(r, mat = FALSE)
  if (all(is.na(vals))) return(NA_real_)

  core_cells <- which(!is.na(vals) & vals >= threshold)
  if (length(core_cells) == 0L) return(0)

  cell_area <- terra::values(
    terra::cellSize(r, unit = "km"),
    mat = FALSE
  )

  sum(cell_area[core_cells], na.rm = TRUE)
}

process_species <- function(species) {
  message("Processing species: ", species)

  threshold <- thresholds %>%
    filter(.data$species == species) %>%
    pull(threshold)

  if (length(threshold) != 1L || is.na(threshold)) {
    stop("Missing or ambiguous threshold for species: ", species)
  }

  retained_months <- occupancy %>%
    filter(.data$species == species) %>%
    pull(month) %>%
    unique()

  if (length(retained_months) == 0L) {
    warning("No retained occupancy months for species: ", species)
    return(NULL)
  }

  lookup <- event_control %>%
    filter(month %in% retained_months)

  sdm_dir <- file.path(sdm_root, species, sdm_subdir)
  pattern <- paste0("^", species, "_\\d{8}\\.", sdm_extension, "$")

  sdm_files <- list.files(
    sdm_dir,
    pattern = pattern,
    full.names = TRUE
  )

  sdm_index <- tibble(
    sdm_file = sdm_files,
    date = extract_sdm_date(sdm_files)
  ) %>%
    filter(!is.na(date)) %>%
    distinct(date, .keep_all = TRUE)

  lookup %>%
    left_join(sdm_index, by = "date") %>%
    left_join(mhw_index %>% rename(region_file = file), by = "date") %>%
    filter(!is.na(sdm_file), !is.na(region_file)) %>%
    mutate(
      species = species,
      area_km2 = map2_dbl(
        sdm_file,
        region_file,
        ~ core_area_km2(.x, .y, threshold)
      )
    ) %>%
    select(
      species, date, doy, season,
      window_id, window_type, area_km2
    )
}

# 5. CALCULATE DAILY CORE-HABITAT AREA =========================================

daily_core <- map_dfr(species_vec, process_species) %>%
  arrange(species, window_id, window_type, date)

write_csv(
  daily_core,
  file.path(out_dir, "daily_core_habitat_area.csv")
)

# 6. EVENT-CONTROL HABITAT-AREA CHANGE =========================================

window_area <- daily_core %>%
  group_by(species, season, window_id, window_type) %>%
  summarise(
    area_km2 = median(area_km2, na.rm = TRUE),
    n_days = sum(!is.na(area_km2)),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = window_type,
    values_from = c(area_km2, n_days),
    names_sep = "_"
  ) %>%
  filter(!is.na(area_km2_Event), !is.na(area_km2_Control)) %>%
  mutate(
    delta_area_km2 = area_km2_Event - area_km2_Control,
    pct_change_area = if_else(
      area_km2_Control > 0,
      100 * delta_area_km2 / area_km2_Control,
      NA_real_
    )
  ) %>%
  arrange(species, season, window_id)

write_csv(
  window_area,
  file.path(out_dir, "event_control_habitat_area.csv")
)

message("Done.")
message("Daily habitat area: ", file.path(out_dir, "daily_core_habitat_area.csv"))
message("Event-control habitat area: ", file.path(out_dir, "event_control_habitat_area.csv"))
