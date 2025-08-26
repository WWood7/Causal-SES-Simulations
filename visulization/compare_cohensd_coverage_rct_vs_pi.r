# Compare CI Coverage: RCT (pi=0.5) vs RCT (pi!=0.5) for Cohen's d and One-step

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# Output directory for comparison
save_dir <- 
"/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results"

# Source CSVs produced by the two visualization scripts
rct_summary_csv <- 
  "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct_visualization/simulation_summary_statistics_rct.csv"
rct_pi_summary_csv <- 
  "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct_pi_visualization/simulation_summary_statistics_rct_pi.csv"

cat("Loading summary CSVs for coverage comparison...\n")

rct_stats <- read.csv(rct_summary_csv, stringsAsFactors = FALSE)
rct_pi_stats <- read.csv(rct_pi_summary_csv, stringsAsFactors = FALSE)

# Keep only the columns we need and tag the setting
rct_stats <- rct_stats %>%
  select(n_factor, vtype_factor, cohens_d_coverage, es_one_step_coverage) %>%
  mutate(setting = "RCT (pi=0.5)")

rct_pi_stats <- rct_pi_stats %>%
  select(n_factor, vtype_factor, cohens_d_coverage, es_one_step_coverage) %>%
  mutate(setting = "RCT (pi!=0.5)")

# Combine
coverage_compare <- bind_rows(rct_stats, rct_pi_stats) %>%
  mutate(
    n_factor = factor(n_factor, levels = unique(n_factor)),
    vtype_factor = factor(vtype_factor, levels = c("Small ES", "Medium ES", "Large ES", "Super Large ES")),
    setting = factor(setting, levels = c("RCT (pi=0.5)", "RCT (pi!=0.5)"))
  )

# Dynamic y-axis limits
min_cov <- min(coverage_compare$cohens_d_coverage, na.rm = TRUE)
max_cov <- max(coverage_compare$cohens_d_coverage, na.rm = TRUE)
y_lower <- max(0, min_cov - 0.05)
y_upper <- min(1, max_cov + 0.02)

# Theme
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )

p <- ggplot(coverage_compare, aes(x = n_factor, y = cohens_d_coverage, color = setting, group = setting)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_wrap(~ vtype_factor, ncol = 3) +
  scale_color_viridis_d(name = "Setting") +
  scale_y_continuous(limits = c(y_lower, y_upper), labels = scales::percent) +
  labs(title = "Cohen's d 95% CI Coverage: RCT vs RCT (pi!=0.5)",
       subtitle = "Red dashed line denotes nominal 95% coverage",
       x = "Sample Size",
       y = "Coverage Probability") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(save_dir, "cohensd_coverage_comparison_rct_vs_rct_pi.png")
ggsave(outfile, p, width = 12, height = 8, dpi = 300, bg = "white")

cat("Saved coverage comparison to:", outfile, "\n")

# One-step coverage plot
min_cov_os <- min(coverage_compare$es_one_step_coverage, na.rm = TRUE)
max_cov_os <- max(coverage_compare$es_one_step_coverage, na.rm = TRUE)
y_lower_os <- max(0, min_cov_os - 0.05)
y_upper_os <- min(1, max_cov_os + 0.02)

p_os <- ggplot(coverage_compare, aes(x = n_factor, y = es_one_step_coverage, color = setting, group = setting)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_wrap(~ vtype_factor, ncol = 3) +
  scale_color_viridis_d(name = "Setting") +
  scale_y_continuous(limits = c(y_lower_os, y_upper_os), labels = scales::percent) +
  labs(title = "One-step 95% CI Coverage: RCT vs RCT (pi!=0.5)",
       subtitle = "Red dashed line denotes nominal 95% coverage",
       x = "Sample Size",
       y = "Coverage Probability") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

outfile_os <- file.path(save_dir, "one_step_coverage_comparison_rct_vs_rct_pi.png")
ggsave(outfile_os, p_os, width = 12, height = 8, dpi = 300, bg = "white")

cat("Saved one-step coverage comparison to:", outfile_os, "\n")


