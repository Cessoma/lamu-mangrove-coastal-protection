# 07 Figures 2 and 3, Supplementary Figures 1 and 2
# Requires: 01-06
# Outputs : plots/Figure2_*, plots/Figure3_*, plots/SupplementaryFigure1_*, plots/SupplementaryFigure2_*

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("sens_results")) load("checkpoints/after_06.RData")

# Shared settings 

# Publication theme: Arial, black text, no minor gridlines
theme_pub <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      text             = element_text(colour = "black", family = "Arial"),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size, colour = "black"),
      legend.title     = element_text(size = base_size, colour = "black"),
      legend.text      = element_text(size = base_size, colour = "black"),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(t = 4, r = 4, b = 4, l = 4)
    )
}

# Bold panel letters (a, b) for combined figures
tag_theme <- theme(plot.tag = element_text(face = "bold", size = 16,
                                           colour = "black", family = "Arial"))

site_colors <- c("Faza" = "#1D2F6F", "Manda" = "#7B4FA6",
                 "Matondoni" = "#A64D79", "Mkunumbi" = "#E8934A",
                 "Mtangawanda" = "#E8D84A")

# Save a figure in one or more formats
save_fig <- function(plot, name, width, height, formats = c("png", "pdf")) {
  for (fmt in formats) {
    args <- list(filename = file.path(PLOT_DIR, paste0(name, ".", fmt)),
                 plot = plot, width = width, height = height, bg = "white")
    if (fmt != "pdf") args$dpi <- 400
    if (fmt == "tiff") args$compression <- "lzw"
    do.call(ggsave, args)
  }
  print(plot)
}

# Interquartile range of a Monte Carlo output for each site
iqr_by_site <- function(var) {
  mc_results |>
    group_by(site_name) |>
    summarise(q25 = unname(quantile(.data[[var]], 0.25)),
              q75 = unname(quantile(.data[[var]], 0.75)),
              .groups = "drop")
}

# Interval plot used for Figures 2b and 3a: 90% interval (light bar),
# interquartile range (dark bar), median (white tick) and best estimate (dot)
interval_plot <- function(df, label_fun, x_title, x_labels = waiver()) {
  df |>
    mutate(site_name = factor(site_name, levels = SITE_ORDER)) |>
    ggplot(aes(y = site_name)) +
    geom_linerange(aes(xmin = lo, xmax = hi),
                   linewidth = 4, colour = "#5E4FA2", alpha = 0.30) +
    geom_linerange(aes(xmin = q25, xmax = q75),
                   linewidth = 4, colour = "#5E4FA2", alpha = 0.70) +
    geom_point(aes(x = med), shape = 124, size = 4,
               colour = "white", stroke = 1) +
    geom_point(aes(x = det), size = 4, shape = 21,
               fill = "#F46D43", colour = "white", stroke = 1.5) +
    geom_text(aes(x = hi + max(hi) * 0.04, label = label_fun(det)),
              hjust = 0, size = 4.6, fontface = "bold", family = "Arial",
              colour = "black") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.20)),
                       labels = x_labels) +
    labs(x = x_title, y = NULL) +
    theme_pub()
}

# Figure 2: attenuation capacity

# Panel a: seaward and landward wave heights with 90% intervals
heights_df <- mc_summary |>
  dplyr::select(site_name,
                `Seaward edge (H0_eff)_lo`  = H0_lo,
                `Seaward edge (H0_eff)_det` = H0_det,
                `Seaward edge (H0_eff)_hi`  = H0_hi,
                `Landward edge (H1)_lo`     = Hl_lo,
                `Landward edge (H1)_det`    = Hl_det,
                `Landward edge (H1)_hi`     = Hl_hi) |>
  pivot_longer(-site_name,
               names_to  = c("position", "stat"),
               names_sep = "_(?=[^_]+$)",
               values_to = "value") |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(position  = factor(position, levels = c("Seaward edge (H0_eff)",
                                                 "Landward edge (H1)")),
         site_name = factor(site_name, levels = SITE_ORDER))

p_fig2a <- heights_df |>
  ggplot(aes(x = site_name, y = det, fill = position)) +
  geom_col(position = position_dodge(0.7), width = 0.62,
           colour = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7),
                width = 0.22, linewidth = 0.7, colour = "grey25") +
  geom_text(aes(label = sprintf("%.2f", det), y = hi + 0.03),
            position = position_dodge(0.7), vjust = 0,
            size = 4.6, fontface = "bold", family = "Arial", colour = "black") +
  scale_fill_manual(
    values = c("Seaward edge (H0_eff)" = "#AED6F2",
               "Landward edge (H1)"    = "#B2D8B2"),
    labels = c("Seaward edge (H0_eff)" = expression("Seaward edge (" * H[0 * "," * eff] * ")"),
               "Landward edge (H1)"    = expression("Landward edge (" * H[landward] * ")")),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Wave height (m)") +
  theme_pub() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# Panel b: wave height attenuation (%) with uncertainty intervals
p_fig2b <- mc_summary |>
  transmute(site_name, lo = at_lo, med = at_med, hi = at_hi, det = at_det) |>
  left_join(iqr_by_site("atten_pct"), by = "site_name") |>
  interval_plot(label_fun = function(x) paste0(round(x, 1), "%"),
                x_title   = "Wave height attenuation (%)",
                x_labels  = label_number(suffix = "%"))

figure2 <- (p_fig2a + p_fig2b) + plot_annotation(tag_levels = "a") & tag_theme
save_fig(figure2, "Figure2_attenuation_capacity", 13, 6, c("tiff", "pdf"))

# Figure 3: realised service provision 

# Panel a: WASP with uncertainty intervals
p_fig3a <- mc_summary |>
  transmute(site_name, lo = wp_lo, med = wp_med, hi = wp_hi, det = wp_det) |>
  left_join(iqr_by_site("WASP"), by = "site_name") |>
  interval_plot(label_fun = function(x) sprintf("%.2f", x),
                x_title   = "WASP (0-1)")

# Panel b: attenuation versus WASP, sized by onshore-wave frequency
p_fig3b <- wave_summary |>
  dplyr::select(site_name, attenuation_pct, WASP, pct_onshore_waves) |>
  mutate(site_name = factor(site_name, levels = SITE_ORDER)) |>
  ggplot(aes(x = attenuation_pct, y = WASP, size = pct_onshore_waves,
             fill = site_name)) +
  geom_point(shape = 21, colour = "white", alpha = 0.9, stroke = 1) +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_size_continuous(name = "Onshore wave\nfrequency (%)", range = c(6, 16),
                        breaks = c(40, 60, 80, 100), limits = c(40, 100)) +
  scale_x_continuous(labels = label_number(suffix = "%"),
                     expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 5)),
         size = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
  labs(x = "Wave height attenuation (%)", y = "WASP (0-1)") +
  theme_pub()

figure3 <- (p_fig3a + p_fig3b) + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_fig(figure3, "Figure3_realised_service", 13, 6)

# Supplementary Figure 1: distributions of attenuation and WASP 

# Best estimates from the deterministic model, shown as red diamonds. For
# bimodal sites the simulation mean falls between the modes, so it is not shown.
best_est <- wave_summary |>
  transmute(site_name = factor(site_name, levels = SITE_ORDER),
            atten_pct = attenuation_pct, WASP)

violin_plot <- function(var, y_title, y_labels = waiver()) {
  mc_results |>
    mutate(site_name = factor(site_name, levels = SITE_ORDER)) |>
    ggplot(aes(x = site_name, y = .data[[var]], fill = site_name)) +
    geom_violin(colour = "grey40", alpha = 0.5, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, colour = "grey20") +
    geom_point(data = best_est, aes(x = site_name, y = .data[[var]]),
               inherit.aes = FALSE, shape = 18, size = 3.5, colour = "#D62728") +
    scale_fill_manual(values = site_colors, guide = "none") +
    scale_y_continuous(labels = y_labels) +
    labs(x = NULL, y = y_title) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

figure_s1 <- (violin_plot("atten_pct", "Wave height attenuation (%)",
                          label_number(suffix = "%")) +
                violin_plot("WASP", "WASP (0–1)")) +
  plot_annotation(tag_levels = "a") & tag_theme
save_fig(figure_s1, "SupplementaryFigure1_distributions", 12, 6)

# Supplementary Figure 2: sensitivity analysis ----

# Colours in the same order as src_labels (script 06)
figure_s2 <- sens_results |>
  mutate(source    = factor(source, levels = src_labels),
         site_name = factor(site_name, levels = SITE_ORDER)) |>
  ggplot(aes(x = site_name, y = pct_var, fill = source)) +
  geom_col(colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = setNames(c("#1565C0", "#AED6F2", "#F4511E", "#B2D8B2", "#E8D84A"),
                      src_labels),
    name = "Uncertainty source"
  ) +
  scale_y_continuous(labels = function(x) paste0(round(x, 0), "%")) +
  labs(x = NULL, y = "% of total WASP variance") +
  theme_pub() +
  theme(legend.text          = element_text(size = 12),
        legend.title         = element_text(size = 13),
        legend.key.spacing.y = unit(4, "pt"),
        panel.grid.major.x   = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

save_fig(figure_s2, "SupplementaryFigure2_sensitivity", 11, 6, c("png", "tiff"))

# Values reported in the text

# WASP mode ranges at the three east-facing sites (bimodal distributions).
# Iterations are split at the midpoint of each site's WASP range.
mc_results |>
  filter(site_name %in% c("Mkunumbi", "Manda", "Faza")) |>
  group_by(site_name) |>
  mutate(mode = ifelse(WASP > (min(WASP) + max(WASP)) / 2, "upper", "lower")) |>
  group_by(site_name, mode) |>
  summarise(share = n() / N_MC,
            lo = quantile(WASP, 0.05), hi = quantile(WASP, 0.95),
            .groups = "drop") |>
  print()

# Sensitivity: percentage of WASP variance from each uncertainty source
sens_results |>
  mutate(pct_var = round(pct_var, 1)) |>
  dplyr::select(site_name, src_id, pct_var) |>
  pivot_wider(names_from = src_id, values_from = pct_var) |>
  print()
