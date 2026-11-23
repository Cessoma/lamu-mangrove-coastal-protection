# 08 Household survey: Figure 4 and Supplementary Table 2
# Perceived coastal protection and mangrove effectiveness
# Descriptive statistics by village (panels a-c) and regression forest plots
# (panels d-f). Runs independently of scripts 01-07.
# Inputs  : Survey_New.xlsx (anonymised household survey)
# Outputs : Final_outputs/ (Figure 4 panels, Supplementary_regression_table.csv)

suppressMessages({
  library(readxl)
  library(dplyr)
  library(forcats)
  library(MASS)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

options(scipen = 999)

if (!dir.exists("Final_outputs")) dir.create("Final_outputs")

# Load and prepare data ----
coastal_data <- read_excel("Survey_New.xlsx")

level_village <- c("Manda", "Mkunumbi", "Matondoni", "Faza", "Mtangawanda")
level_education <- c("No formal education", "Primary", "Secondary", "College and University")
level_land_owner <- c("Yes", "No")
level_experience_flooding <- c("No", "Yes")
level_years_lived <- c("0_10", "11_20", "21_30", "31_40", "41_50", ">50")
level_effectiveness_mangrove <- c("Less effective", "Equally effective", "More effective")
level_coast_protected <- c("No", "Yes")

coastal_data <- coastal_data %>%
  mutate(village = factor(village, levels = level_village),
         education = factor(education, levels = level_education),
         land_owner = factor(land_owner, levels = level_land_owner),
         experience_flooding = factor(experience_flooding, levels = level_experience_flooding),
         years_lived = factor(years_lived, levels = level_years_lived),
         effectiveness_mangrove = factor(effectiveness_mangrove, levels = level_effectiveness_mangrove),
         coast_protected = factor(coast_protected, levels = level_coast_protected))

coastal_data$effectiveness_mangrove <- factor(coastal_data$effectiveness_mangrove,
                                              levels = c("Less effective", "Equally effective", "More effective"), ordered = TRUE)
coastal_data$experience_flooding <- factor(coastal_data$experience_flooding)
coastal_data$village <- factor(coastal_data$village)
coastal_data$land_owner <- factor(coastal_data$land_owner)
coastal_data$land_owner <- relevel(coastal_data$land_owner, ref = "No")
coastal_data$years_lived <- factor(coastal_data$years_lived)
coastal_data$education <- factor(coastal_data$education, ordered = FALSE)
coastal_data$coast_protected <- factor(coastal_data$coast_protected, levels = c("No", "Yes"))

# Order of villages in the figures (southwest to northeast). Manda remains the
# reference category in the regression models (level_village).
site_order <- c("Mkunumbi", "Matondoni", "Manda", "Mtangawanda", "Faza")

# Colours: panels a-c (response categories), panels d-f (significance)
col_positive <- "#AED6F1"
col_negative <- "#F4A9AE"
col_middle   <- "#B2D8B2"
col_sig      <- "#C9536B"
col_ns       <- "#999999"

# Shared theme (180 mm print width, 8 pt text)
base_theme <- theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
    plot.title = element_text(face = "bold", size = 10, color = "black", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    legend.position = "bottom",
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(t = 0, b = 0, l = 0, r = 0),
    legend.box.margin = margin(t = -5)
  )

# Panels a-c: descriptive statistics by village 
panel_a_data <- coastal_data %>%
  mutate(village = factor(village, levels = site_order),
         flood_cat = factor(if_else(experience_flooding == "Yes", "Yes", "No"),
                            levels = c("Yes", "No"))) %>%
  count(village, flood_cat) %>%
  group_by(village) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

panel_a <- ggplot(panel_a_data, aes(x = village, y = pct, fill = flood_cat)) +
  geom_col(position = "stack", width = 0.75) +
  geom_text(aes(label = paste0(round(pct), "%")), position = position_stack(vjust = 0.5),
            size = 2.2, color = "black") +
  scale_fill_manual(values = c("Yes" = col_positive, "No" = col_negative), name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "a", x = NULL, y = "Percentage of households (%)") +
  base_theme

panel_b_data <- coastal_data %>%
  mutate(village = factor(village, levels = site_order),
         protection_cat = factor(if_else(coast_protected == "No", "Not well protected", "Well protected"),
                                 levels = c("Well protected", "Not well protected"))) %>%
  count(village, protection_cat) %>%
  group_by(village) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

panel_b <- ggplot(panel_b_data, aes(x = village, y = pct, fill = protection_cat)) +
  geom_col(position = "stack", width = 0.75) +
  geom_text(aes(label = paste0(round(pct), "%")), position = position_stack(vjust = 0.5),
            size = 2.2, color = "black") +
  scale_fill_manual(values = c("Not well protected" = col_negative, "Well protected" = col_positive), name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "b", x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank())

panel_c_data <- coastal_data %>%
  mutate(village = factor(village, levels = site_order),
         effectiveness_mangrove = factor(effectiveness_mangrove,
                                         levels = c("More effective", "Equally effective", "Less effective"))) %>%
  count(village, effectiveness_mangrove) %>%
  group_by(village) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

panel_c <- ggplot(panel_c_data, aes(x = village, y = pct, fill = effectiveness_mangrove)) +
  geom_col(position = "stack", width = 0.75) +
  geom_text(aes(label = paste0(round(pct), "%")), position = position_stack(vjust = 0.5),
            size = 2.2, color = "black") +
  scale_fill_manual(values = c("Less effective" = col_negative, "Equally effective" = col_middle,
                               "More effective" = col_positive), name = NULL,
                    breaks = c("Less effective", "Equally effective", "More effective")) +
  guides(fill = guide_legend(nrow = 3)) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "c", x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank())

# Panels d-f: regression forest plots 
adj_flood <- glm(experience_flooding ~ village + land_owner + years_lived + education,
                 data = coastal_data, family = "binomial")

adj_protected <- glm(coast_protected ~ village + land_owner + years_lived + education + experience_flooding,
                     data = coastal_data, family = "binomial")

model_effectiveness <- polr(effectiveness_mangrove ~ experience_flooding + village + land_owner +
                              years_lived + education,
                            data = coastal_data, Hess = TRUE)

tidy_glm <- function(model) {
  ct <- coef(summary(model))
  ci <- exp(confint(model))
  data.frame(term = names(coef(model)), aOR = exp(coef(model)),
             LB = ci[, 1], UB = ci[, 2], p_value = ct[, 4], row.names = NULL)
}

tidy_polr <- function(model) {
  ct <- coef(summary(model))
  coef_rows <- rownames(ct)[!grepl("\\|", rownames(ct))]
  ci <- exp(confint(model))
  p_values <- pnorm(abs(ct[coef_rows, "t value"]), lower.tail = FALSE) * 2
  data.frame(term = coef_rows, aOR = exp(coef(model))[coef_rows],
             LB = ci[coef_rows, 1], UB = ci[coef_rows, 2], p_value = p_values, row.names = NULL)
}

flood_tidy <- tidy_glm(adj_flood) %>% filter(term != "(Intercept)")
protected_tidy <- tidy_glm(adj_protected) %>% filter(term != "(Intercept)")
effectiveness_tidy <- tidy_polr(model_effectiveness)

flood_tidy$outcome <- "d"
protected_tidy$outcome <- "e"
effectiveness_tidy$outcome <- "f"

forest_data <- bind_rows(flood_tidy, protected_tidy, effectiveness_tidy)

label_map <- c(
  "villageMkunumbi" = "Mkunumbi", "villageMatondoni" = "Matondoni",
  "villageFaza" = "Faza", "villageMtangawanda" = "Mtangawanda",
  "land_ownerYes" = "Land owner",
  "years_lived11_20" = "Lived 11-20 yrs", "years_lived21_30" = "Lived 21-30 yrs",
  "years_lived31_40" = "Lived 31-40 yrs", "years_lived41_50" = "Lived 41-50 yrs",
  "years_lived>50" = "Lived >50 yrs",
  "educationPrimary" = "Primary", "educationSecondary" = "Secondary",
  "educationCollege and University" = "College/University",
  "experience_floodingYes" = "Experienced flooding"
)

group_map <- c(
  "villageMkunumbi" = "Village", "villageMatondoni" = "Village",
  "villageFaza" = "Village", "villageMtangawanda" = "Village",
  "land_ownerYes" = "Land ownership",
  "years_lived11_20" = "Years lived", "years_lived21_30" = "Years lived",
  "years_lived31_40" = "Years lived", "years_lived41_50" = "Years lived",
  "years_lived>50" = "Years lived",
  "educationPrimary" = "Education", "educationSecondary" = "Education",
  "educationCollege and University" = "Education",
  "experience_floodingYes" = "Flood experience"
)

row_order <- unname(label_map[c(
  "villageMkunumbi", "villageMatondoni", "villageMtangawanda", "villageFaza", "land_ownerYes",
  "years_lived11_20", "years_lived21_30", "years_lived31_40", "years_lived41_50", "years_lived>50",
  "educationPrimary", "educationSecondary", "educationCollege and University",
  "experience_floodingYes"
)])

n <- length(row_order)

# Y-axis labels with variable prefix (e.g. Village: Mkunumbi)
row_order_display <- c(
  "Village: Mkunumbi", "Village: Matondoni", "Village: Mtangawanda", "Village: Faza",
  "Land ownership: Yes",
  "Years lived: 11-20 yrs", "Years lived: 21-30 yrs", "Years lived: 31-40 yrs",
  "Years lived: 41-50 yrs", "Years lived: >50 yrs",
  "Education: Primary", "Education: Secondary", "Education: College/University",
  "Flood experience: Yes"
)

forest_data <- forest_data %>%
  mutate(
    label = label_map[term],
    group = group_map[term],
    significant = ifelse(LB > 1 | UB < 1, "Significant", "Not significant"),
    y_pos = n - match(label, row_order) + 1,
    outcome = factor(outcome, levels = c("d", "e", "f"))
  )

term_for_label <- setNames(names(label_map), label_map)
group_seq <- group_map[term_for_label[row_order]]
change_at <- which(group_seq[-1] != group_seq[-n])
group_bounds <- (n - change_at) + 0.5

forest_plot <- ggplot(forest_data, aes(x = aOR, y = y_pos, color = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = group_bounds, linetype = "dotted", color = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = LB, xmax = UB), height = 0, linewidth = 1.1) +
  geom_point(size = 1.6) +
  scale_color_manual(values = c("Not significant" = col_ns, "Significant" = col_sig)) +
  scale_x_log10(breaks = c(0.1, 1, 10), labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(breaks = n:1, labels = row_order_display, limits = c(0.5, n + 0.5)) +
  facet_wrap(~outcome, nrow = 1, scales = "fixed") +
  labs(x = "Adjusted odds ratio (log scale)", y = "Predictor", color = NULL) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    strip.text = element_text(face = "bold", size = 10, color = "black")
  )

# Save panels 
descriptive_panels <- panel_a + panel_b + panel_c + plot_layout(nrow = 1)

ggsave("Final_outputs/descriptive_panels_abc.tiff", descriptive_panels,
       device = "tiff", dpi = 300, width = 7.09, height = 2.6,
       units = "in", compression = "lzw", bg = "white")

ggsave("Final_outputs/forest_plot_regression.tiff", forest_plot,
       device = "tiff", dpi = 300, width = 7.09, height = 3.3,
       units = "in", compression = "lzw", bg = "white")

outcome_full_names <- c("d" = "Flood exposure", "e" = "Perceived as not well protected",
                        "f" = "Mangrove effectiveness")

supp_table <- forest_data %>%
  transmute(
    Outcome = factor(outcome_full_names[as.character(outcome)],
                     levels = outcome_full_names[c("d", "e", "f")]),
    Predictor = factor(label, levels = row_order),
    Group = group,
    aOR = round(aOR, 2), `95% CI` = paste0(round(LB, 2), "-", round(UB, 2)),
    p_value = round(p_value, 3)
  ) %>%
  arrange(Outcome, Predictor)

write.csv(supp_table, "Final_outputs/Supplementary_regression_table.csv", row.names = FALSE)

# Combine into one six-panel figure 
combined_figure <- descriptive_panels / forest_plot +
  plot_layout(heights = c(1, 1.3))

print(combined_figure)

ggsave("Final_outputs/combined_figure_abcdef.tiff", combined_figure,
       device = "tiff", dpi = 300, width = 7.09, height = 6.5,
       units = "in", compression = "lzw", bg = "white")
