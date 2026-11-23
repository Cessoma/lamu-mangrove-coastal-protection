# 02 Fetch and wave heights
# Fetch by ray-casting, SMB wave heights and directional wave exposure.
# Requires: 01

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("wind_raw")) load("checkpoints/after_01.RData")

# Fetch calculation 

cast_ray <- function(origin, bearing_deg, coast_geom, max_dist = MAX_FETCH_M) {
  theta <- (90 - bearing_deg) * pi / 180

  ox <- st_coordinates(origin)[1]
  oy <- st_coordinates(origin)[2]

  ray <- st_sfc(
    st_linestring(matrix(c(ox, oy,
                           ox + max_dist * cos(theta),
                           oy + max_dist * sin(theta)),
                         ncol = 2, byrow = TRUE)),
    crs = st_crs(origin)
  )

  ints <- tryCatch(
    st_intersection(ray, coast_geom),
    error = function(e) st_sfc(crs = st_crs(origin))
  )

  if (length(ints) == 0 || all(st_is_empty(ints))) return(max_dist)

  pts <- st_coordinates(ints)
  if (nrow(pts) == 0) return(max_dist)

  dists <- sqrt((pts[, "X"] - ox)^2 + (pts[, "Y"] - oy)^2)

  valid <- dists[dists > MIN_DIST_EXCL]

  if (length(valid) == 0) return(max_dist)

  return(min(valid, max_dist))
}

message("Calculating fetch (open-water distances) - may take a minute...")

fetch_all <- map_dfr(seq_len(nrow(lamu_sites)), function(i) {
  site <- lamu_sites[i, ]

  fetches <- vapply(DIRECTIONS,
                    function(d) cast_ray(site, d, mangrove_boundary),
                    numeric(1))

  tibble(site_name = site$site_name,
         direction = DIRECTIONS,
         fetch_m   = fetches)
})

# Wave height model (WEMo, SMB equations) 

wind_stress_ua <- function(U10) UA_A * U10^UA_B

smb_Hs <- function(U10, F_m, d = Inf) {
  if (U10 <= 0 || F_m <= 0) return(0)

  UA <- wind_stress_ua(U10)

  gF <- G * F_m / UA^2

  # SMB deep-water significant wave height (USACE, 2002)
  Hs <- 0.283 * (UA^2 / G) * tanh(0.0125 * gF^0.42)

  return(max(Hs, 0))
}

smb_Hs_vec <- function(U10_v, F_m_v) {
  mapply(smb_Hs, U10_v, F_m_v, MoreArgs = list(d = Inf))
}

# Directional wave exposure 

compute_exposure <- function(site_nm, fetch_df, wind_df) {
  half_bin <- RAY_SPACING / 2  # 5 deg; each ray represents a +/-5 deg sector
  n_total  <- nrow(wind_df)

  map_dfr(seq_len(nrow(fetch_df)), function(i) {
    d   <- fetch_df$direction[i]
    F_m <- fetch_df$fetch_m[i]

    lo <- (d - half_bin) %% 360
    hi <- (d + half_bin) %% 360

    wsub <- if (lo < hi) {
      filter(wind_df, wind_direction >= lo, wind_direction < hi)
    } else {
      filter(wind_df, wind_direction >= lo | wind_direction < hi)
    }

    n_obs <- nrow(wsub)
    freq  <- n_obs / n_total

    if (n_obs == 0) {
      return(tibble(site_name = site_nm, direction = d, fetch_m = F_m,
                    n_obs = 0L, freq = 0, U_mean = 0, U_p95 = 0,
                    Hs_mean = 0, Hs_p95 = 0, RWE = 0))
    }

    U_mean <- mean(wsub$wind_speed)

    U_p95  <- as.numeric(quantile(wsub$wind_speed, 0.95))

    Hs_obs  <- smb_Hs_vec(wsub$wind_speed, rep(F_m, n_obs))

    Hs_mean <- mean(Hs_obs)

    Hs_p95  <- smb_Hs(U_p95, F_m)

    RWE <- freq * mean(Hs_obs^2)

    tibble(site_name = site_nm, direction = d, fetch_m = F_m,
           n_obs = as.integer(n_obs), freq, U_mean, U_p95,
           Hs_mean, Hs_p95, RWE)
  })
}

message("Computing wave exposure for all directions...")

wave_exp <- map_dfr(unique(fetch_all$site_name), function(nm) {
  fd <- filter(fetch_all, site_name == nm)
  compute_exposure(nm, fd, wind_raw)
})

# Site-level wave summary 

wave_summary <- wave_exp %>%
  group_by(site_name) %>%
  summarise(
    max_fetch_km  = max(fetch_m) / 1000,
    mean_fetch_km = mean(fetch_m) / 1000,
    dominant_dir  = direction[which.max(Hs_p95)],
    # Mean wave height across all wind records: each direction's mean Hs is
    # weighted by how often the wind blew from that direction (freq)
    Hs_mean_all   = sum(freq * Hs_mean),
    Hs_p95_max    = max(Hs_p95),
    RWE_total     = sum(RWE),
    .groups = "drop"
  )

# Save session for the next script
save.image("checkpoints/after_02.RData")
