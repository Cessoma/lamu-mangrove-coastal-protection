# 06 Sensitivity analysis
# Contribution of each uncertain input to WASP variance.
# Requires: 01-05

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("mc_results")) load("checkpoints/after_05.RData")

# Fixed seed for reproducible results
set.seed(43)

# Sensitivity analysis 

run_one_source <- function(nm, src) {
  dom_dir   <- wave_summary$dominant_dir[wave_summary$site_name == nm]
  F_dom     <- fetch_all$fetch_m[
    fetch_all$site_name == nm & fetch_all$direction == dom_dir]
  W_base    <- wave_summary$mangrove_width_m[wave_summary$site_name == nm]
  d_site    <- site_depths[nm]
  r_base    <- r_fn(d_site)
  pct_inund <- PCT_INUNDATION[nm]
  n_total   <- nrow(wind_raw)

  lo_dir <- (dom_dir - RAY_SPACING / 2) %% 360
  hi_dir <- (dom_dir + RAY_SPACING / 2) %% 360
  wsub <- if (lo_dir < hi_dir) {
    filter(wind_raw, wind_direction >= lo_dir, wind_direction < hi_dir)
  } else {
    filter(wind_raw, wind_direction >= lo_dir | wind_direction < hi_dir)
  }
  base_speeds <- wsub$wind_speed
  n_obs       <- length(base_speeds)
  sigma_ln    <- sqrt(log(1 + R_COEF_CV^2))

  vapply(seq_len(N_MC), function(k) {
    sp_k <- if (src == "wind_speed") {
      pmax(base_speeds + rnorm(n_obs, 0, WIND_SPEED_SD), 0)
    } else base_speeds

    smb_k <- if (src == "smb") pmax(rnorm(1, 1, SMB_CV), 0.05) else 1.0
    w_k   <- if (src == "width") max(W_base + rnorm(1, 0, WIDTH_SD), 0) else W_base
    r_k   <- if (src == "r_coef") {
      max(r_base * rlnorm(1, -sigma_ln^2/2, sigma_ln), 0)
    } else r_base
    d_err <- if (src == "wind_dir") rnorm(1, 0, WIND_DIR_SD) else 0

    Hs_k <- hs_fn(as.numeric(quantile(sp_k, 0.95)), F_dom) * smb_k
    H0_k <- min(Hs_k, WAVE_BREAK_KD * d_site)
    Hl_k <- H0_k * exp(-r_k * w_k)
    at_k <- (1 - Hl_k / H0_k) * 100

    # Direction error shifts the wind records only; onshore sector stays fixed
    lo_k   <- (dom_dir - 90) %% 360
    hi_k   <- (dom_dir + 90) %% 360
    dirs_k <- (wind_raw$wind_direction + d_err) %% 360
    n_on_k <- if (lo_k < hi_k) sum(dirs_k >= lo_k & dirs_k < hi_k) else
      sum(dirs_k >= lo_k | dirs_k < hi_k)

    pct_inund * (100 * n_on_k / n_total) * at_k / 1e6
  }, numeric(1))
}

message("Running sensitivity analysis ...")

src_ids    <- c("wind_speed", "wind_dir", "smb", "width", "r_coef")
# Labels match the uncertain inputs in Supplementary Table 4
src_labels <- c("Wind speed (SD = 0.5 m s\u207B\u00B9)",
                "Wind direction (SD = 5\u00B0)",
                "Wave growth model, SMB (CV = 15%)",
                "Mangrove belt width (SD = 15 m)",
                "Wave-decay coefficient, r (CV = 20%)")

sens_results <- map_dfr(wave_summary$site_name, function(nm) {
  map_dfr(seq_along(src_ids), function(j) {
    vals <- run_one_source(nm, src_ids[j])
    data.frame(
      site_name = nm,
      source    = src_labels[j],
      src_id    = src_ids[j],
      variance  = var(vals),
      sd_WASP   = sd(vals),
      range_90  = unname(diff(quantile(vals, c(0.05, 0.95)))),
      stringsAsFactors = FALSE
    )
  })
}) |>
  group_by(site_name) |>
  mutate(pct_var = 100 * variance / sum(variance)) |>
  ungroup()

message("  Sensitivity analysis complete.")

# Save sensitivity results (variance share of each uncertainty source)
write.csv(sens_results, paste0("lamu_sensitivity_", RUN_TAG, ".csv"), row.names = FALSE)

# Save session for the next script
save.image("checkpoints/after_06.RData")
