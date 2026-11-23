# 03 Wave attenuation and WASP index
# Seaward edge, mangrove belt width, water depth, wave attenuation
# and the WASP index for each site (values reported in Table 1).
# Requires: 01-02
# Outputs : lamu_wemo_results.csv and supporting tables

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("wave_exp")) load("checkpoints/after_02.RData")

# Seaward edge of the mangroves ----

comm_with_dir <- communities_sf %>%
  left_join(dplyr::select(wave_summary, site_name, dominant_dir), by = "site_name") %>%
  mutate(bearing_landward = (dominant_dir + 180) %% 360)

true_seaward_list <- vector("list", nrow(lamu_sites))

for (i in seq_len(nrow(lamu_sites))) {
  nm <- lamu_sites$site_name[i]
  interior_pt <- lamu_sites[i, ]
  dom_dir <- wave_summary$dominant_dir[wave_summary$site_name == nm]

  theta <- (90 - dom_dir) * pi / 180

  ox <- st_coordinates(interior_pt)[1]
  oy <- st_coordinates(interior_pt)[2]

  short_ray <- st_sfc(
    st_linestring(matrix(c(ox, oy,
                           ox + 300 * cos(theta),
                           oy + 300 * sin(theta)),
                         ncol = 2, byrow = TRUE)),
    crs = st_crs(interior_pt)
  )

  crossings <- st_intersection(short_ray, mangrove_boundary)

  if (length(crossings) > 0 && !all(st_is_empty(crossings))) {
    pts <- st_coordinates(crossings)
    dists <- sqrt((pts[, "X"] - ox)^2 + (pts[, "Y"] - oy)^2)
    valid <- dists > 1
    if (any(valid)) {
      first <- which.min(dists[valid])
      true_x <- pts[valid, ][first, "X"]
      true_y <- pts[valid, ][first, "Y"]
    } else {
      true_x <- ox; true_y <- oy
    }
  } else {
    true_x <- ox; true_y <- oy
  }

  true_seaward_list[[i]] <- st_sf(
    site_name = nm,
    geometry = st_sfc(st_point(c(true_x, true_y)), crs = st_crs(interior_pt))
  )

  message(sprintf("  %-12s : seaward edge moved %0.0f m", nm,
                  sqrt((true_x-ox)^2 + (true_y-oy)^2)))
}

true_seaward_sf <- do.call(rbind, true_seaward_list)

# Mangrove belt width ----

measure_belt_width <- function(seaward_pt, bearing_landward, village_pt,
                               mangrove_poly, max_dist = 5000) {
  theta <- (90 - bearing_landward) * pi / 180

  ox <- st_coordinates(seaward_pt)[1]
  oy <- st_coordinates(seaward_pt)[2]

  ray <- st_sfc(
    st_linestring(matrix(c(ox, oy,
                           ox + max_dist * cos(theta),
                           oy + max_dist * sin(theta)),
                         ncol = 2, byrow = TRUE)),
    crs = st_crs(seaward_pt)
  )

  inside_segs <- tryCatch(
    st_intersection(ray, mangrove_poly),
    error = function(e) NULL
  )

  if (is.null(inside_segs) || length(inside_segs) == 0) {
    return(list(width_m = 0, note = "No mangroves along transect"))
  }

  if (grepl("MULTILINESTRING|GEOMETRYCOLLECTION",
            as.character(st_geometry_type(inside_segs)))) {
    segs <- st_collection_extract(inside_segs, "LINESTRING")
  } else {
    segs <- inside_segs
  }

  seg_info <- lapply(seq_len(length(segs)), function(k) {
    seg    <- segs[k]
    coords <- st_coordinates(seg)
    if (nrow(coords) < 2) return(NULL)
    start_dist <- sqrt((coords[1, "X"] - ox)^2 + (coords[1, "Y"] - oy)^2)
    list(start = start_dist, length = as.numeric(st_length(seg)))
  })

  seg_info <- Filter(Negate(is.null), seg_info)
  if (length(seg_info) == 0) return(list(width_m = 0, note = "No valid segments"))

  first <- seg_info[[which.min(sapply(seg_info, `[[`, "start"))]]

  width_mangrove <- first$length

  vx <- st_coordinates(village_pt)[1]
  vy <- st_coordinates(village_pt)[2]
  dist_to_village <- sqrt((vx - ox)^2 + (vy - oy)^2)

  effective <- min(width_mangrove, dist_to_village)

  if (effective > MAX_REALISTIC_WIDTH) {
    warning(sprintf("%s: width %.0f m exceeds realistic limit. Capping.",
                    seaward_pt$site_name, effective))
    effective <- MAX_REALISTIC_WIDTH
  }

  return(list(width_m = round(effective, 1),
              note = sprintf("mangrove belt=%.0f m, village=%.0f m",
                             width_mangrove, dist_to_village)))
}

width_results <- vector("list", nrow(comm_with_dir))
for (i in seq_len(nrow(comm_with_dir))) {
  nm         <- comm_with_dir$site_name[i]
  bearing_lw <- comm_with_dir$bearing_landward[i]
  seaward_pt <- true_seaward_sf[true_seaward_sf$site_name == nm, ]
  village_pt <- communities_sf[communities_sf$site_name == nm, ]
  res        <- measure_belt_width(seaward_pt, bearing_lw, village_pt, lamu_mangrove)
  message(sprintf("  %-12s : width = %.0f m", nm, res$width_m))
  width_results[[i]] <- tibble(site_name = nm, mangrove_width_m = res$width_m)
}

width_summary <- bind_rows(width_results)

# Water depth at the seaward edge

site_depths <- c(
  Manda       = 1.70,  # Deepest site - channels at Manda Bay
  Mkunumbi    = 0.51,  # Shallow intertidal flat
  Matondoni   = 0.92,
  Faza        = 1.33,
  Mtangawanda = 1.69
)

# Wave attenuation and WASP 

# Depth-dependent attenuation coefficient (after Mazda et al., 2006)
r_mangrove <- function(d) {
  r_raw <- 0.0086 * d^2 - 0.0128 * d + 0.0076
  return(max(r_raw, 0))
}

wave_summary <- wave_summary %>%
  left_join(width_summary, by = "site_name") %>%
  mutate(
    depth_m = site_depths[site_name],
    H0_eff_m = pmin(Hs_p95_max, WAVE_BREAK_KD * depth_m),
    r_value = vapply(depth_m, r_mangrove, numeric(1)),
    H_landward_m = H0_eff_m * exp(-r_value * mangrove_width_m),  # Mazda et al. (2006)
    attenuation_pct = (1 - H_landward_m / H0_eff_m) * 100,
    # Onshore-wave frequency: % of wind records from within +/-90 deg of the
    # dominant wave direction. Winds exactly 90 deg from it (alongshore) are
    # excluded by the strict upper bound (< hi).
    pct_onshore_waves = map_dbl(site_name, function(nm) {
      dom <- wave_summary$dominant_dir[wave_summary$site_name == nm]
      lo  <- (dom - 90) %% 360
      hi  <- (dom + 90) %% 360
      n_on <- if (lo < hi) {
        sum(wind_raw$wind_direction >= lo & wind_raw$wind_direction < hi)
      } else {
        sum(wind_raw$wind_direction >= lo | wind_raw$wind_direction < hi)
      }
      100 * n_on / nrow(wind_raw)
    }),
    pct_inundation = PCT_INUNDATION[site_name],
    # WASP = % inundation x % onshore waves x % attenuation / 10^6 (TESSA)
    WASP = pct_inundation * pct_onshore_waves * attenuation_pct / 1e6
  )

# Save results 

write_csv(wave_summary, paste0("lamu_wemo_results_", RUN_TAG, ".csv"))

write_csv(width_summary, paste0("lamu_mangrove_widths_", RUN_TAG, ".csv"))

write_csv(wave_exp, paste0("lamu_directional_exposure_", RUN_TAG, ".csv"))

write_csv(fetch_all, paste0("lamu_fetch_distances_", RUN_TAG, ".csv"))

# Save session for the next script
save.image("checkpoints/after_03.RData")
