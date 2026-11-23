# 01 Setup and inputs
# Packages, model settings, study sites, mangrove extent and wind data.
# Outputs : objects used by scripts 02-07
# Inputs  : Vector_data/lamu_coastline_polygon.shp, lamu_wind_2024_2025.csv

# Packages 

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)
SITE_ORDER <- c("Mkunumbi", "Matondoni", "Manda", "Mtangawanda", "Faza")

# Settings 

# Working directory
# Run all scripts from the repository root (e.g. open the RStudio project)
BASE_DIR <- getwd()

# Output label: results files are named with this label and figures are saved
# in plots/<RUN_TAG>/. Change it to keep outputs from different runs separate.
RUN_TAG  <- "main"
PLOT_DIR <- file.path("plots", RUN_TAG)

dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create("checkpoints", showWarnings = FALSE)

# Physical constants
G <- 9.81
WAVE_BREAK_KD <- 0.78  # Breaker index (McCowan, 1894)

# Wind-stress adjustment (UA) parameters
UA_A <- 0.71  # Wind-stress factor UA = 0.71 * U10^1.23 (USACE, 2002)
UA_B <- 1.23

# Ray-casting parameters (fetch calculation)
MAX_FETCH_M <- 50000  # Maximum fetch (m)
RAY_SPACING <- 10  # Ray spacing (degrees); 36 rays per site

DIRECTIONS <- seq(0, 350, by = RAY_SPACING)

MIN_DIST_EXCL <- 50  # Ignore shoreline hits within 50 m of the site

# Mangrove structural type
MANGROVE_TYPE <- "aerial"

# Tidal inundation percentages (TESSA M4): plots flooded at high tide and
# exposed at low tide, placing all sites in zone D (between MHWN and MLWN)
PCT_INUNDATION <- c(
  Manda       = 50,
  Mkunumbi    = 50,
  Matondoni   = 50,
  Faza        = 50,
  Mtangawanda = 50
)

# Upper limit on realistic belt width
MAX_REALISTIC_WIDTH <- 1500  # Upper limit on belt width (m)

# Spatial data 

# Study sites
site_coords <- data.frame(
  site_name = c("Manda", "Mkunumbi", "Matondoni", "Faza", "Mtangawanda"),
  lon       = c(40.935987, 40.732707, 40.842305, 41.109528, 40.971728),
  lat       = c(-2.291773, -2.302789, -2.261917, -2.048144, -2.122702)
)

lamu_sites <- st_as_sf(site_coords, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = 32737)

# Mangrove polygon
lamu_mangrove <- st_read("Vector_data/lamu_coastline_polygon.shp") %>%
  st_transform(crs = 32737) %>%
  st_union() %>%
  st_as_sf()

mangrove_boundary <- st_boundary(lamu_mangrove)

mangrove_plot <- st_simplify(lamu_mangrove, preserveTopology = TRUE, dTolerance = 10)

# Community (village) locations
comm_coords <- data.frame(
  site_name = c("Manda","Mkunumbi","Matondoni","Faza","Mtangawanda"),
  lon       = c(40.936131, 40.728842, 40.844122, 41.106481, 40.974952),
  lat       = c(-2.293624, -2.308761, -2.263894, -2.050565, -2.120821)
)

communities_sf <- st_as_sf(comm_coords, coords = c("lon","lat"), crs = 4326) %>%
  st_transform(crs = 32737)

# Wind data 

wind_raw <- read_csv("lamu_wind_2024_2025.csv", show_col_types = FALSE) %>%
  mutate(
    datetime       = dmy_hm(datetime),
    wind_speed     = wind_speed_mps,
    wind_direction = wind_dir_deg
  ) %>%
  dplyr::select(datetime, wind_speed, wind_direction) %>%
  drop_na() %>%
  filter(wind_speed > 0) %>%
  mutate(wind_direction = if_else(wind_direction == 0, 360, wind_direction))

# Save session for the next script
save.image("checkpoints/after_01.RData")
