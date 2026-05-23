# =============================================================================
# fibre_csa_analysis.R
# Skeletal Muscle Fibre Cross-Sectional Area (CSA) Analysis
#
# Author:  B. Shamim
# Dataset: fibre-cross-sectional-area-dataset.csv
#
# Overview:
#   Loads immunohistochemistry-derived fibre CSA data, applies a circularity
#   filter to exclude non-circular (i.e. obliquely sectioned) fibres, computes
#   summary statistics, runs a linear mixed-effects model, and produces a
#   suite of visualisations.
#
# Conditions: Condition A, Condition B, Condition C
# Legs:       L (Control), R (Intervention)
#
# Circularity filter: > 0.6 (standard threshold for IHC fibre CSA analysis)
#
# Output files saved to: outputs/
#   - figures/  → PNG plots
#   - tables/   → CSV summary files
# =============================================================================


# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------

# Install any missing packages before loading
# install.packages(c("tidyverse", "nlme", "lsmeans", "patchwork", "janitor"))

library(tidyverse)   # dplyr, ggplot2, tidyr, readr
library(janitor)     # clean_names()
library(nlme)        # lme() — linear mixed-effects models
library(lsmeans)     # pairwise contrasts from mixed models
library(patchwork)   # combining ggplot panels

# Working directory — set to project root
setwd("~/Documents/fibre-cross-sectional-area/data")

# Create output folders if they don't exist
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

# Consistent colour palette
COL_L    <- "#1B6B7B"   # Left leg
COL_R    <- "#C47A1E"   # Right leg
COL_A    <- "#1B6B7B"   # Condition A
COL_B    <- "#C47A1E"   # Condition B
COL_C    <- "#8B1A2E"   # Condition C

THEME_CSA <- theme_minimal(base_size = 12) +
  theme(
    panel.background  = element_blank(),
    axis.line         = element_line(colour = "black"),
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(colour = "grey40", size = 10),
    strip.text        = element_text(face = "bold")
  )


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------

csa_raw <- read_csv(
  "fibre-cross-sectional-area-dataset.csv",
  show_col_types = FALSE
) %>%
  clean_names()   # standardises column names to snake_case

# Confirm structure
cat("=== Raw data ===\n")
cat("Rows:", nrow(csa_raw), "\n")
cat("Columns:", paste(names(csa_raw), collapse = ", "), "\n\n")

# Confirm factor levels
cat("Conditions:", paste(unique(csa_raw$condition), collapse = ", "), "\n")
cat("Legs:      ", paste(unique(csa_raw$leg), collapse = ", "), "\n")
cat("Sex:       ", paste(unique(csa_raw$sex), collapse = ", "), "\n")
cat("Subjects:  ", n_distinct(csa_raw$subject), "unique\n\n")


# -----------------------------------------------------------------------------
# 2. FILTER & PREPARE
# -----------------------------------------------------------------------------

# Apply circularity filter (> 0.6) to exclude obliquely sectioned fibres
# This is a standard threshold in IHC fibre morphology analysis
csa_filtered <- csa_raw %>%
  filter(circularity > 0.6) %>%
  # Set condition order for all plots
  mutate(
    condition = factor(condition,
                       levels = c("Condition A", "Condition B", "Condition C")),
    leg       = factor(leg, levels = c("L", "R"),
                       labels = c("Left", "Right")),
    sex       = factor(sex),
    subject   = factor(subject)
  ) %>%
  arrange(condition)

cat("=== After circularity filter (> 0.6) ===\n")
cat("Rows retained:", nrow(csa_filtered),
    paste0("(", round(nrow(csa_filtered) / nrow(csa_raw) * 100, 1), "% of raw)\n"))
cat("Rows excluded:", nrow(csa_raw) - nrow(csa_filtered), "\n\n")

# Split by condition for condition-specific plots
csa_A <- filter(csa_filtered, condition == "Condition A")
csa_B <- filter(csa_filtered, condition == "Condition B")
csa_C <- filter(csa_filtered, condition == "Condition C")


# -----------------------------------------------------------------------------
# 3. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

# Per-subject median and IQR (fibre-level → subject-level)
summary_subject <- csa_filtered %>%
  group_by(condition, sex, subject, leg) %>%
  summarise(
    n_fibres = n(),
    median   = median(fibre_area),
    iqr      = IQR(fibre_area),
    mean     = mean(fibre_area),
    sd       = sd(fibre_area),
    .groups  = "drop"
  )

# Group-level summary (condition × leg)
summary_group <- csa_filtered %>%
  group_by(condition, leg) %>%
  summarise(
    n_fibres   = n(),
    mean_area  = mean(fibre_area),
    sd_area    = sd(fibre_area),
    median_area = median(fibre_area),
    iqr_area   = IQR(fibre_area),
    .groups    = "drop"
  )

# Median per condition × leg for density plot vertical lines
median_by_group <- csa_filtered %>%
  group_by(condition, leg) %>%
  summarise(median_area = median(fibre_area), .groups = "drop")

cat("=== Group-level summary ===\n")
print(summary_group)
cat("\n")

# Save tables
write_csv(summary_subject,
          "~/Documents/fibre-cross-sectional-area/outputs/tables/summary_by_subject.csv")
write_csv(summary_group,
          "~/Documents/fibre-cross-sectional-area/outputs/tables/summary_by_group.csv")
write_csv(csa_filtered,
          "~/Documents/fibre-cross-sectional-area/outputs/tables/csa_filtered_0.6.csv")

cat("Tables saved to outputs/tables/\n\n")


# -----------------------------------------------------------------------------
# 4. VISUALISATIONS
# -----------------------------------------------------------------------------

# ── 4a. Density plots — one panel per condition ────────────────────────────
# Shows full distribution of fibre areas by leg within each condition.
# Dashed vertical lines mark the median per leg.

make_density_plot <- function(data, condition_name, median_data, colours) {
  med <- median_data %>% filter(condition == condition_name)
  ggplot(data, aes(x = fibre_area, fill = leg)) +
    geom_density(alpha = 0.5, colour = NA) +
    geom_vline(
      data     = med,
      aes(xintercept = median_area, colour = leg),
      linetype = "dashed", linewidth = 0.8
    ) +
    scale_fill_manual(values   = colours, name = "Leg") +
    scale_colour_manual(values = colours, name = "Leg") +
    labs(
      title    = condition_name,
      subtitle = "Dashed lines indicate median fibre area per leg",
      x        = expression("Fibre area" ~ (mu * m^{2})),
      y        = "Density"
    ) +
    THEME_CSA
}

leg_colours <- c("Left" = COL_L, "Right" = COL_R)

p_density_A <- make_density_plot(csa_A, "Condition A", median_by_group, leg_colours)
p_density_B <- make_density_plot(csa_B, "Condition B", median_by_group, leg_colours)
p_density_C <- make_density_plot(csa_C, "Condition C", median_by_group, leg_colours)

# Combined density panel
p_density_combined <- p_density_A / p_density_B / p_density_C +
  plot_annotation(
    title   = "Fibre CSA Distributions by Condition and Leg",
    caption = "Circularity filter > 0.6 applied"
  )

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/01_density_by_condition.png",
       p_density_combined, width = 9, height = 12, dpi = 150)
cat("Saved: 01_density_by_condition.png\n")


# ── 4b. Bar chart with error bars — mean ± SD by condition and leg ─────────
# Replicates the facet_wrap bar chart from the original script, updated for
# new column names and condition labels.

p_bar <- ggplot(summary_group,
                aes(x = leg, y = mean_area, fill = leg)) +
  facet_wrap(~ condition, strip.position = "bottom") +
  geom_col(position = position_dodge(), alpha = 0.85, width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_area, ymax = mean_area + sd_area),
    width    = 0.2,
    position = position_dodge(0.9),
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = leg_colours,
    name   = "Leg"
  ) +
  labs(
    title = "Mean Fibre CSA by Condition and Leg",
    x     = NULL,
    y     = expression("Fibre CSA" ~ (mu * m^{2}))
  ) +
  THEME_CSA +
  theme(strip.placement = "outside")

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/02_bar_mean_sd.png",
       p_bar, width = 10, height = 5, dpi = 150)
cat("Saved: 02_bar_mean_sd.png\n")


# ── 4c. Box plots — fibre-level distribution by condition and leg ──────────
# More informative than bar charts for skewed fibre area distributions.

p_box <- ggplot(csa_filtered,
                aes(x = leg, y = fibre_area, fill = leg)) +
  facet_wrap(~ condition) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21,
               outlier.size = 1, outlier.alpha = 0.4) +
  scale_fill_manual(values = leg_colours, name = "Leg") +
  labs(
    title    = "Fibre CSA Distribution by Condition and Leg",
    subtitle = "Box = IQR, line = median, whiskers = 1.5× IQR",
    x        = NULL,
    y        = expression("Fibre area" ~ (mu * m^{2}))
  ) +
  THEME_CSA

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/03_boxplot_by_condition.png",
       p_box, width = 10, height = 5, dpi = 150)
cat("Saved: 03_boxplot_by_condition.png\n")


# ── 4d. Violin + jitter — subject-level medians overlaid ──────────────────
# Shows fibre distribution shape AND individual subject data points.

p_violin <- ggplot(csa_filtered,
                   aes(x = leg, y = fibre_area, fill = leg)) +
  facet_wrap(~ condition) +
  geom_violin(alpha = 0.5, colour = NA, trim = TRUE) +
  geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA,
               fill = "white", colour = "grey30") +
  # Overlay subject-level medians
  geom_point(
    data     = summary_subject,
    aes(x = leg, y = median),
    shape    = 21, size = 3,
    fill     = "white", colour = "grey20",
    position = position_jitter(width = 0.05, seed = 42)
  ) +
  scale_fill_manual(values = leg_colours, name = "Leg") +
  labs(
    title    = "Fibre CSA — Violin Plot with Subject Medians",
    subtitle = "White points = per-subject median fibre area",
    x        = NULL,
    y        = expression("Fibre area" ~ (mu * m^{2}))
  ) +
  THEME_CSA

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/04_violin_with_subjects.png",
       p_violin, width = 10, height = 5, dpi = 150)
cat("Saved: 04_violin_with_subjects.png\n")


# ── 4e. Histogram grid — fibre count distribution ─────────────────────────
# Useful for reporting n fibres per bin and identifying multimodal distributions.

p_hist <- ggplot(csa_filtered,
                 aes(x = fibre_area, fill = leg)) +
  facet_grid(leg ~ condition) +
  geom_histogram(bins = 40, alpha = 0.8, colour = "white") +
  geom_vline(
    data     = median_by_group,
    aes(xintercept = median_area),
    linetype = "dashed", colour = "grey30", linewidth = 0.7
  ) +
  scale_fill_manual(values = leg_colours, name = "Leg") +
  labs(
    title    = "Fibre CSA Frequency Distribution",
    subtitle = "Dashed line = median per group",
    x        = expression("Fibre area" ~ (mu * m^{2})),
    y        = "Fibre count"
  ) +
  THEME_CSA +
  theme(legend.position = "none")

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/05_histogram_grid.png",
       p_hist, width = 11, height = 6, dpi = 150)
cat("Saved: 05_histogram_grid.png\n")


# ── 4f. Cumulative distribution — condition comparison ────────────────────
# Shows the proportion of fibres below a given area threshold.
# Useful for comparing fibre size distributions across conditions.

p_ecdf <- ggplot(csa_filtered,
                 aes(x = fibre_area, colour = condition, linetype = leg)) +
  stat_ecdf(linewidth = 0.9) +
  scale_colour_manual(
    values = c("Condition A" = COL_A,
               "Condition B" = COL_B,
               "Condition C" = COL_C),
    name = "Condition"
  ) +
  scale_linetype_manual(
    values = c("Left" = "solid", "Right" = "dashed"),
    name   = "Leg"
  ) +
  labs(
    title    = "Cumulative Distribution of Fibre CSA",
    subtitle = "Solid = Left leg, Dashed = Right leg",
    x        = expression("Fibre area" ~ (mu * m^{2})),
    y        = "Cumulative proportion"
  ) +
  THEME_CSA

ggsave("~/Documents/fibre-cross-sectional-area/outputs/figures/06_cumulative_distribution.png",
       p_ecdf, width = 9, height = 5, dpi = 150)
cat("Saved: 06_cumulative_distribution.png\n\n")


# -----------------------------------------------------------------------------
# 5. STATISTICAL ANALYSIS
# -----------------------------------------------------------------------------

cat("=== Statistical Analysis ===\n\n")

# ── 5a. Shapiro-Wilk normality test on fibre area ──────────────────────────
# NOTE: With n = 5000+ fibres, the Shapiro-Wilk test will almost always
# reject normality due to sample size sensitivity — this is expected for
# fibre CSA data. The mixed model below is robust to this.
# We test on a random subsample (n = 500) as a practical check.

set.seed(42)
sw_sample <- csa_filtered %>% slice_sample(n = 500)
sw_result <- shapiro.test(sw_sample$fibre_area)
cat("Shapiro-Wilk (n=500 subsample):\n")
cat("  W =", round(sw_result$statistic, 4),
    "  p =", format.pval(sw_result$p.value, digits = 3), "\n\n")


# ── 5b. Linear mixed-effects model ─────────────────────────────────────────
# Fixed effects:  Condition, Leg, Condition × Leg interaction
# Random effect:  Subject (accounts for repeated measures per animal)
# Method:         REML

cat("Fitting linear mixed-effects model...\n")

model_csa <- lme(
  fibre_area ~ condition * leg,
  data   = csa_filtered,
  random = ~ 1 | subject,
  method = "REML"
)

cat("\n--- Model Summary ---\n")
print(summary(model_csa))


# ── 5c. Pairwise contrasts ─────────────────────────────────────────────────

cat("\n--- Pairwise contrasts: Leg ---\n")
print(lsmeans(model_csa, pairwise ~ leg, adjust = "none"))

cat("\n--- Pairwise contrasts: Condition ---\n")
print(lsmeans(model_csa, pairwise ~ condition, adjust = "none"))

cat("\n--- Pairwise contrasts: Condition × Leg (LSD) ---\n")
print(lsmeans(model_csa, pairwise ~ condition * leg, adjust = "none"))


# -----------------------------------------------------------------------------
# 6. DONE
# -----------------------------------------------------------------------------

cat("\n=================================================\n")
cat("  Analysis complete.\n")
cat("  Figures saved to:  outputs/figures/\n")
cat("  Tables saved to:   outputs/tables/\n")
cat("=================================================\n")
