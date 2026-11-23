# 04 Figure 1b-g: wind rose and fetch roses
# Panel (a), the study area map, was produced separately in ArcGIS.
# Panel letters and the final layout were added when combining the panels.
# Requires: 01-03
# Outputs : plots/Figure1_wind_rose.png, plots/Figure1_fetch_<site>.png

library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

# Load results of earlier scripts if they are not already in the session
if (!exists("width_summary")) load("checkpoints/after_03.RData")

# Figure 1b: wind rose 

# Directions are plotted on a continuous 0-360 degree axis so that each
# 10-degree sector is centred on its compass bearing (N at the top).
p_wind_rose <- wind_raw %>%
  mutate(
    dir_centre = (round(wind_direction / 10) * 10) %% 360,
    speed_bin = cut(wind_speed,
                    breaks = c(2, 4, 6, 8, Inf),
                    labels = c("2-4","4-6","6-8",">8"))
  ) %>%
  drop_na(dir_centre, speed_bin) %>%
  count(dir_centre, speed_bin) %>%
  ggplot(aes(x = dir_centre, y = n, fill = speed_bin)) +
  geom_col(width = 10, colour = "white", linewidth = 0.2) +
  coord_polar(start = -pi / 36) +
  scale_x_continuous(limits = c(-5, 355), breaks = c(0, 90, 180, 270),
                     labels = c("N", "E", "S", "W")) +
  scale_fill_viridis_d(option = "plasma", name = "Speed (m s\u207b\u00b9)", direction = -1,
                       drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank())

plot_fetch_rose_arrow <- function(df, nm) {
  dom <- wave_summary$dominant_dir[wave_summary$site_name == nm]

  dir_levels <- sort(unique(df$direction))
  dom_bin <- dir_levels[which.min(abs(((dir_levels - dom + 180) %% 360) - 180))]

  rmax <- max(df$fetch_m, na.rm = TRUE) / 1000

  ggplot(df, aes(x = direction %% 360, y = fetch_m / 1000)) +
    geom_col(fill = "#AED6F2", colour = "white", linewidth = 0.15, width = 10) +
    geom_segment(data = data.frame(dom_bin = dom_bin, rmax = rmax),
                 aes(x = dom_bin %% 360, xend = dom_bin %% 360,
                     y = 0, yend = rmax * 0.95),
                 colour = "#C9536B", linewidth = 1.1,
                 arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
                 inherit.aes = FALSE) +
    coord_polar(start = -pi / 36) +
    scale_x_continuous(limits = c(-5, 355), breaks = c(0, 90, 180, 270),
                       labels = c("N", "E", "S", "W")) +
    scale_y_continuous(labels = label_number(suffix = " km")) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      text = element_text(colour = "black", family = "Arial"),
      axis.text = element_text(size = 14, colour = "black", family = "Arial"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5,
                                 colour = "black", family = "Arial"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
    )
}

fetch_rose_arrow_list <- wave_exp %>%
  mutate(site_name = factor(site_name, levels = SITE_ORDER)) %>%
  group_by(site_name) %>%
  group_map(~ plot_fetch_rose_arrow(.x, .y$site_name), .keep = TRUE)

# Wind rose panel formatting
wind_panel <- p_wind_rose +
  labs(title = NULL, subtitle = NULL, x = NULL) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.text = element_text(size = 14, colour = "black", family = "Arial"),
    axis.title = element_text(size = 14, colour = "black", family = "Arial"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 14, colour = "black", family = "Arial"),
    legend.text = element_text(size = 14, colour = "black", family = "Arial"),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4)
  )

site_order <- wave_exp %>%
  mutate(site_name = factor(site_name, levels = SITE_ORDER)) %>%
  group_by(site_name) %>% group_keys() %>% pull(site_name)

# Save panels 
p_wind <- wind_panel
ggsave(file.path(PLOT_DIR, "Figure1_wind_rose.png"), p_wind, width = 5, height = 5, dpi = 400, bg = "white")
print(p_wind)

for (i in seq_along(fetch_rose_arrow_list)) {
  nm <- site_order[i]
  p_i <- fetch_rose_arrow_list[[i]] + labs(title = nm) +
    theme(plot.title = element_text(size = 16, colour = "black", family = "Arial"))
  fname <- paste0(file.path(PLOT_DIR, "Figure1_fetch_"), tolower(nm), ".png")
  ggsave(fname, p_i, width = 5, height = 5, dpi = 300, bg = "white")
  print(p_i)
}
