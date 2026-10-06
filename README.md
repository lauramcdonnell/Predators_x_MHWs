# Predators × Marine Heatwaves

Code supporting the manuscript:

**McDonnell et al. _Marine heatwaves expose multiple categories of predator habitat disruption._**

This repository contains code used to identify marine heatwaves (MHWs), construct regional event-control windows, and quantify changes in predicted predator core-habitat area during MHW events.

## Repository contents

### `01_detect_marine_heatwaves.R`

Detects daily MHWs from gridded daily sea surface temperature (SST) and writes one binary MHW occurrence file per day.

The script:

- calculates a grid-cell-specific 90th-percentile SST threshold using a centered 11-day window over the 1 June 2002 to 1 June 2022 climatological baseline;
- smooths the daily threshold climatology using a centered 31-day circular running mean;
- allows up to two consecutive sub-threshold days within a candidate MHW event;
- requires at least five actual threshold-exceedance days for an event to qualify; and
- writes daily binary MHW fields (`1 = MHW`, `0 = non-MHW`).

### `02_build_event_control_windows.R`

Uses the daily binary MHW files from Script 01 to calculate regional MHW extent and construct matched event-control windows.

The script:

- calculates daily regional MHW extent as the proportion of ocean pixels classified as MHW;
- defines qualifying event days as days at or above the regional 70th percentile of MHW extent;
- groups qualifying event days while allowing one intervening non-event day;
- retains event windows spanning at least seven calendar days;
- defines control observations as days at or below the regional 30th percentile of MHW extent;
- matches controls to the same day-of-year range as each event window; and
- retains days of year represented in at least three years under both event and control conditions.

## Workflow

Run the scripts in numerical order:

```text
01_detect_marine_heatwaves.R
        |
        v
daily binary MHW files
        |
        v
02_build_event_control_windows.R
        |
        v
regional MHW extent and event-control lookup tables
```

Users should edit the paths in the **USER SETTINGS** section at the top of each script before running.

## `03_calculate_core_habitat_area.R`
Calculates daily species core-habitat area from species distribution model predictions and links those estimates to the matched MHW event-control windows.

The script:
- applies a species-specific habitat-suitability threshold to each daily SDM prediction to identify core habitat;
- restricts calculations to species-specific retained occupancy months;
- calculates total daily core-habitat area in km² using grid-cell area;
- links daily habitat estimates to the event and control dates produced by Script 02;
- summarizes habitat area within each event and matched control period using the median; and
- calculates percent habitat-area change as \((A_{event}-A_{control})/A_{control}\times100\).



## Data requirements

### Sea surface temperature

Daily SST fields are from the NASA/JPL **GHRSST Level 4 MUR Global Foundation Sea Surface Temperature Analysis, version 4.1**.

DOI: https://doi.org/10.5067/GHGMR-4FJ04

The scripts assume one daily NetCDF SST file per day on a common grid and study-region mask.

### Species distribution models

The species distribution models used in the manuscript were adapted from:

Braun et al. (2023), _Building use-inspired species distribution models: Using multiple data types to examine and improve model performance_, **Ecological Applications**.

Code associated with development of those models is archived at:

https://doi.org/10.5281/zenodo.7971532

## R packages

Script 01 requires:

- `terra`
- `dplyr`
- `stringr`

Script 02 requires:

- `terra`
- `dplyr`
- `tidyr`
- `lubridate`
- `readr`
- `stringr`

## Reproducibility scope

The scripts currently provided here reproduce the MHW detection and event-control construction portions of the analysis. 

## Citation

Citation information will be added upon publication.
