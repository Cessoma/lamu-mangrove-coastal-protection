# Runs the wave model and uncertainty analysis (scripts 01-07) in order.
# Script 08 (household survey) is independent and can be run on its own.

scripts <- c(
  "01_setup_and_inputs.R",
  "02_fetch_and_wave_heights.R",
  "03_attenuation_and_WASP.R",
  "04_Figure1_wind_and_fetch_roses.R",
  "05_monte_carlo_uncertainty.R",
  "06_sensitivity_analysis.R",
  "07_uncertainty_figures.R"
)

for (s in scripts) {
  message("Running ", s)
  source(s)
}
