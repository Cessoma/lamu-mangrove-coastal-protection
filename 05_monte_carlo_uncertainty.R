# 05 Monte Carlo uncertainty analysis
# Error distributions (Supplementary Table 4), Monte Carlo simulation
# and 90% uncertainty intervals.
# Requires: 01-03
# Outputs : lamu_mc_uncertainty_summary.csv, lamu_mc_raw_iterations.csv

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("width_summary")) load("checkpoints/after_03.RData")

# Uncertainty parameters 
set.seed(42)

N_MC <- 5000              # Monte Carlo iterations

WIND_SPEED_SD <- 0.5      # Wind speed SD (m/s; WMO, 2023)

WIND_DIR_SD <- 5.0        # Wind direction SD (degrees; WMO, 2023)

SMB_CV <- 0.15            # SMB structural error (CV; USACE, 2002)

WIDTH_SD <- 15.0          # Belt width SD (m)

R_COEF_CV <- 0.20         # Attenuation coefficient CV

cat("\n=== Uncertainty parameters used in the Monte Carlo simulation ===\n")
cat(sprintf("  Wind speed SD  : +/-%.1f m/s  [WMO2018]\n",  WIND_SPEED_SD))
cat(sprintf("  Wind dir SD    : +/-%.0f deg  [WMO2018]\n",  WIND_DIR_SD))
cat(sprintf("  SMB model CV   : %.0f%%        [USACE-CEM]\n", SMB_CV * 100))
cat(sprintf("  Belt width SD  : +/-%.0f m    [McIvor2012]\n", WIDTH_SD))
cat(sprintf("  Mazda r CV     : %.0f%%        [Mazda2006]\n",  R_COEF_CV * 100))

# Helper functions for the monte carlo loop 
ua_fn <- function(U10) UA_A * U10^UA_B

hs_fn <- function(U10, F_m) {
  if (U10 <= 0 || F_m <= 0) return(0)
  UA <- ua_fn(U10)
  gF <- G * F_m / UA^2
  max(0.283 * (UA^2 / G) * tanh(0.0125 * gF^0.42), 0)
}

r_fn <- function(d) {
  max(0.0086 * d^2 - 0.0128 * d + 0.0076, 0)
}

# Monte carlo simulation 

message("Running Monte Carlo simulation (", N_MC,
        " iterations x ", nrow(wave_summary), " sites) ...")
message("  Expected run time: 2-5 minutes.")

mc_results <- map_dfr(wave_summary$site_name, function(nm) {
  # Best-estimate inputs for this site (from scripts 02-03)
  dom_dir   <- wave_summary$dominant_dir[wave_summary$site_name == nm]
  F_dom     <- fetch_all$fetch_m[
    fetch_all$site_name == nm & fetch_all$direction == dom_dir]
  W_base    <- wave_summary$mangrove_width_m[wave_summary$site_name == nm]
  d_site    <- site_depths[nm]
  r_base    <- r_fn(d_site)
  pct_inund <- PCT_INUNDATION[nm]
  n_total   <- nrow(wind_raw)

  # Wind observations in the dominant direction sector
  lo_dir <- (dom_dir - RAY_SPACING / 2) %% 360
  hi_dir <- (dom_dir + RAY_SPACING / 2) %% 360
  wsub <- if (lo_dir < hi_dir) {
    filter(wind_raw, wind_direction >= lo_dir, wind_direction < hi_dir)
  } else {
    filter(wind_raw, wind_direction >= lo_dir | wind_direction < hi_dir)
  }
  base_speeds <- wsub$wind_speed
  n_obs       <- length(base_speeds)

  # Pre-generate all random draws for this site

  spd_err <- matrix(rnorm(N_MC * n_obs, 0, WIND_SPEED_SD),
                    nrow = N_MC, ncol = n_obs)

  smb_fac <- pmax(rnorm(N_MC, 1, SMB_CV), 0.05)

  wid_err <- rnorm(N_MC, 0, WIDTH_SD)

  sigma_ln <- sqrt(log(1 + R_COEF_CV^2))
  r_fac    <- rlnorm(N_MC, meanlog = -sigma_ln^2 / 2, sdlog = sigma_ln)

  dir_err <- rnorm(N_MC, 0, WIND_DIR_SD)

  # Run N_MC iterations
  rows <- vector("list", N_MC)

  for (k in seq_len(N_MC)) {
    sp_k <- pmax(base_speeds + spd_err[k, ], 0)

    Hs_obs_k <- vapply(sp_k, function(u) hs_fn(u, F_dom), numeric(1))

    Hs_p95_k <- hs_fn(as.numeric(quantile(sp_k, 0.95)), F_dom) * smb_fac[k]

    H0_eff_k <- min(Hs_p95_k, WAVE_BREAK_KD * d_site)

    W_k <- max(W_base + wid_err[k], 0)

    r_k <- max(r_base * r_fac[k], 0)

    H_land_k <- H0_eff_k * exp(-r_k * W_k)

    atten_k <- (1 - H_land_k / H0_eff_k) * 100

    # Onshore-wave frequency with wind-direction instrument error. The error
    # shifts the recorded wind directions only; the onshore sector (+/-90 deg of
    # the dominant wave direction) is fixed by the coastline, as in script 03.
    lo_k   <- (dom_dir - 90) %% 360
    hi_k   <- (dom_dir + 90) %% 360
    dirs_k <- (wind_raw$wind_direction + dir_err[k]) %% 360
    n_on_k <- if (lo_k < hi_k) {
      sum(dirs_k >= lo_k & dirs_k < hi_k)
    } else {
      sum(dirs_k >= lo_k | dirs_k < hi_k)
    }
    pct_on_k <- 100 * n_on_k / n_total

    WASP_k <- pct_inund * pct_on_k * atten_k / 1e6

    rows[[k]] <- list(
      site_name = nm,
      iteration = k,
      Hs_p95    = Hs_p95_k,
      H0_eff    = H0_eff_k,
      H_land    = H_land_k,
      atten_pct = atten_k,
      WASP      = WASP_k
    )
  }

  bind_rows(rows)
})

message("  Monte Carlo complete: ", nrow(mc_results), " rows  |  ",
        "Columns: ", paste(names(mc_results), collapse = ", "))

stopifnot(all(c("Hs_p95","H0_eff","H_land","atten_pct","WASP") %in%
                names(mc_results)))
stopifnot(!anyNA(mc_results$WASP))
message("  Sanity checks passed.")

# Summarise uncertainty intervals 

q <- function(x, p) unname(quantile(x, p))

mc_summary <- lapply(
  split(mc_results, mc_results$site_name),
  function(df) {
    data.frame(
      site_name = df$site_name[1],
      Hs_lo  = q(df$Hs_p95,    0.05),
      Hs_med = q(df$Hs_p95,    0.50),
      Hs_hi  = q(df$Hs_p95,    0.95),
      H0_lo  = q(df$H0_eff,    0.05),
      H0_med = q(df$H0_eff,    0.50),
      H0_hi  = q(df$H0_eff,    0.95),
      Hl_lo  = q(df$H_land,    0.05),
      Hl_med = q(df$H_land,    0.50),
      Hl_hi  = q(df$H_land,    0.95),
      at_lo  = q(df$atten_pct, 0.05),
      at_med = q(df$atten_pct, 0.50),
      at_hi  = q(df$atten_pct, 0.95),
      wp_lo  = q(df$WASP,      0.05),
      wp_med = q(df$WASP,      0.50),
      wp_hi  = q(df$WASP,      0.95),
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }
) |>
  bind_rows() |>
  left_join(
    data.frame(
      site_name = wave_summary$site_name,
      Hs_det    = wave_summary$Hs_p95_max,
      H0_det    = wave_summary$H0_eff_m,
      Hl_det    = wave_summary$H_landward_m,
      at_det    = wave_summary$attenuation_pct,
      wp_det    = wave_summary$WASP
    ),
    by = "site_name"
  )

cat("\n=== Uncertainty summary - 90% credible intervals ===\n")
mc_summary |>
  dplyr::select(site_name, Hs_lo, Hs_det, Hs_hi,
                at_lo, at_det, at_hi, wp_lo, wp_det, wp_hi) |>
  as.data.frame() |>
  print()

write.csv(mc_summary, paste0("lamu_mc_uncertainty_summary_", RUN_TAG, ".csv"), row.names = FALSE)
write.csv(mc_results, paste0("lamu_mc_raw_iterations_", RUN_TAG, ".csv"), row.names = FALSE)
message("Uncertainty tables saved.")

saveRDS(list(mc_results = mc_results, mc_summary = mc_summary,
             wave_summary = wave_summary), paste0("mc_cache_", RUN_TAG, ".rds"))
message("Cache saved for faster reloading.")

file.exists(paste0("mc_cache_", RUN_TAG, ".rds"))

mc_results |>
  filter(site_name == "Mkunumbi") |>
  summarise(mean_WASP = mean(WASP),
            q10 = quantile(WASP, 0.10),
            q50 = quantile(WASP, 0.50),
            q89 = quantile(WASP, 0.89),
            q90 = quantile(WASP, 0.90))

# Save session for the next script
save.image("checkpoints/after_05.RData")
