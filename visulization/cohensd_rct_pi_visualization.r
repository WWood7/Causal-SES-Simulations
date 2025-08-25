# Causal Effect Size Simulation Results Analysis and Visualization - RCT PI VERSION
# Comparing Cohen's d, Robust Cohen's d, and Causal Effect Size Estimators (pi!=0.5)

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

# set the saving directory
save_dir <- 
"/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct_pi_visualization"

# Set results directory for RCT pi (change if you stage files elsewhere)
results_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct_pi"

# =============================================================================
# 1. LOAD AND COMBINE ALL RESULTS
# =============================================================================

cat("Loading RCT pi simulation results...\n")

# Get all RDS files
rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "result files\n")

# Read and combine all results
all_results <- data.frame()

for (file in rds_files) {
    tryCatch({
        result <- readRDS(file)
        all_results <- rbind(all_results, result)
    }, error = function(e) {
        cat("Error reading file:", file, "\n")
    })
}

cat("Combined", nrow(all_results), "simulation results\n")

# Sort and clean data
all_results <- all_results[order(all_results$seed, all_results$n, all_results$effect_size_type), ]

# Add factor labels for better plotting
n_levels <- sort(unique(all_results$n))
all_results$n_factor <- factor(all_results$n, levels = n_levels,
                              labels = paste0("n=", n_levels))

all_results$vtype_factor <- factor(all_results$effect_size_type, 
                                  levels = c(1, 2, 3, 4),
                                  labels = c("Small ES", "Medium ES", "Large ES", "Super Large ES"))

# =============================================================================
# 2. SUMMARY STATISTICS
# =============================================================================

cat("\n=== RCT PI SUMMARY STATISTICS ===\n")

# Calculate bias, MSE, and coverage for each estimator
summary_stats <- all_results %>%
  group_by(n, effect_size_type, vtype_factor, n_factor) %>%
  summarise(
    # Sample info
    n_sims = n(),
    true_es_mean = mean(true_es),
    
    # Cohen's d statistics
    cohens_d_mean = mean(cohens_d, na.rm = TRUE),
    cohens_d_bias = mean(cohens_d - true_es, na.rm = TRUE),
    cohens_d_mse = mean((cohens_d - true_es)^2, na.rm = TRUE),
    cohens_d_coverage = mean(cohens_d_lb <= true_es & true_es <= cohens_d_ub, na.rm = TRUE),
    cohens_d_ci_width = mean(cohens_d_ub - cohens_d_lb, na.rm = TRUE),

    # Robust Cohen's d statistics
    robust_cohens_d_mean = mean(robust_cohens_d, na.rm = TRUE),
    robust_cohens_d_bias = mean(robust_cohens_d - true_es, na.rm = TRUE),
    robust_cohens_d_mse = mean((robust_cohens_d - true_es)^2, na.rm = TRUE),
    robust_cohens_d_coverage = mean(robust_cohens_d_lb <= true_es & true_es <= robust_cohens_d_ub, na.rm = TRUE),
    robust_cohens_d_ci_width = mean(robust_cohens_d_ub - robust_cohens_d_lb, na.rm = TRUE),
    
    # One-step causal ES statistics
    es_one_step_mean = mean(es_one_step, na.rm = TRUE),
    es_one_step_bias = mean(es_one_step - true_es, na.rm = TRUE),
    es_one_step_mse = mean((es_one_step - true_es)^2, na.rm = TRUE),
    es_one_step_coverage = mean(es_one_step_lb <= true_es & true_es <= es_one_step_ub, na.rm = TRUE),
    es_one_step_ci_width = mean(es_one_step_ub - es_one_step_lb, na.rm = TRUE),
    
    .groups = "drop"
  )

# Print summary table
print(summary_stats)

# Save summary statistics
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_stats, paste0(save_dir, "/simulation_summary_statistics_rct_pi.csv"), row.names = FALSE)
cat("Summary statistics saved to: simulation_summary_statistics_rct_pi.csv\n")

# =============================================================================
# 3. DATA PREPARATION FOR VISUALIZATION
# =============================================================================

# Create long format for bias comparison
bias_data <- all_results %>%
  mutate(
    cohens_d_bias = cohens_d - true_es,
    robust_cohens_d_bias = robust_cohens_d - true_es,
    es_one_step_bias = es_one_step - true_es
  ) %>%
  select(seed, n, effect_size_type, n_factor, vtype_factor, true_es,
         cohens_d_bias, robust_cohens_d_bias, es_one_step_bias) %>%
  pivot_longer(cols = c(cohens_d_bias, robust_cohens_d_bias, es_one_step_bias),
               names_to = "estimator", values_to = "bias") %>%
  mutate(estimator = factor(estimator, 
                           levels = c("cohens_d_bias", "robust_cohens_d_bias", "es_one_step_bias"),
                           labels = c("Cohen's d", "Robust Cohen's d", "One-step Causal ES")))

# Create long format for estimates comparison
estimates_data <- all_results %>%
  select(seed, n, effect_size_type, n_factor, vtype_factor, true_es,
         cohens_d, robust_cohens_d, es_one_step) %>%
  pivot_longer(cols = c(cohens_d, robust_cohens_d, es_one_step),
               names_to = "estimator", values_to = "estimate") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d", "robust_cohens_d", "es_one_step"),
                           labels = c("Cohen's d", "Robust Cohen's d", "One-step Causal ES")))

# =============================================================================
# 4. VISUALIZATION FUNCTIONS
# =============================================================================

# Theme for consistent plotting
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold")
  )

# =============================================================================
# 5. DASHBOARD VISUALIZATION
# =============================================================================

# Plot 1: Bias by Sample Size and Variability Type
p1 <- ggplot(bias_data, aes(x = n_factor, y = bias, fill = estimator)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_wrap(~ vtype_factor, ncol = 3) +
  scale_fill_viridis_d(name = "Estimator") +
  labs(title = "RCT (pi!=0.5): Bias Comparison Across Conditions",
       subtitle = "Bias = Estimate - True Effect Size",
       x = "Sample Size", 
       y = "Bias") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2: MSE Comparison (from summary stats)
mse_data <- summary_stats %>%
  select(n_factor, vtype_factor, cohens_d_mse, robust_cohens_d_mse, es_one_step_mse) %>%
  pivot_longer(cols = c(cohens_d_mse, robust_cohens_d_mse, es_one_step_mse),
               names_to = "estimator", values_to = "mse") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_mse", "robust_cohens_d_mse", "es_one_step_mse"),
                           labels = c("Cohen's d", "Robust Cohen's d", "One-step Causal ES")))

p2 <- ggplot(mse_data, aes(x = n_factor, y = mse, color = estimator, group = estimator)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ vtype_factor, ncol = 3, scales = "free_y") +
  scale_color_viridis_d(name = "Estimator") +
  labs(title = "RCT (pi!=0.5): Mean Squared Error by Sample Size",
       subtitle = "Lower is better",
       x = "Sample Size", 
       y = "MSE") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Coverage Probability
coverage_data <- summary_stats %>%
  select(n_factor, vtype_factor, cohens_d_coverage, robust_cohens_d_coverage, es_one_step_coverage) %>%
  pivot_longer(cols = c(cohens_d_coverage, robust_cohens_d_coverage, es_one_step_coverage),
               names_to = "estimator", values_to = "coverage") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_coverage", "robust_cohens_d_coverage", "es_one_step_coverage"),
                           labels = c("Cohen's d", "Robust Cohen's d", "One-step Causal ES")))

# Dynamic y-axis limits
min_coverage <- min(coverage_data$coverage, na.rm = TRUE)
max_coverage <- max(coverage_data$coverage, na.rm = TRUE)
y_lower <- max(0, min_coverage - 0.05)
y_upper <- min(1, max_coverage + 0.02)

p3 <- ggplot(coverage_data, aes(x = n_factor, y = coverage, color = estimator, group = estimator)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_wrap(~ vtype_factor, ncol = 3) +
  scale_color_viridis_d(name = "Estimator") +
  scale_y_continuous(limits = c(y_lower, y_upper), labels = scales::percent) +
  labs(title = "RCT (pi!=0.5): 95% Confidence Interval Coverage",
       subtitle = "Should be close to 95% (red dashed line)",
       x = "Sample Size", 
       y = "Coverage Probability") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 4: Confidence Interval Width Comparison
ci_width_data <- summary_stats %>%
  select(n_factor, vtype_factor, cohens_d_ci_width, robust_cohens_d_ci_width, es_one_step_ci_width) %>%
  pivot_longer(cols = c(cohens_d_ci_width, robust_cohens_d_ci_width, es_one_step_ci_width),
               names_to = "estimator", values_to = "ci_width") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_ci_width", "robust_cohens_d_ci_width", "es_one_step_ci_width"),
                           labels = c("Cohen's d", "Robust Cohen's d", "One-step Causal ES")))

p4 <- ggplot(ci_width_data, aes(x = n_factor, y = ci_width, color = estimator, group = estimator)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ vtype_factor, ncol = 3, scales = "free_y") +
  scale_color_viridis_d(name = "Estimator") +
  labs(title = "RCT (pi!=0.5): 95% Confidence Interval Width",
       subtitle = "Narrower intervals indicate higher precision",
       x = "Sample Size", 
       y = "CI Width") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =============================================================================
# 6. SAVE DASHBOARD
# =============================================================================

# Create and save summary dashboard
dashboard <- grid.arrange(p1, p2, p3, p4, ncol = 1)
ggsave(paste0(save_dir, "/simulation_dashboard_rct_pi.png"), dashboard, width = 12, height = 24, dpi = 300)

# =============================================================================
# 7. NUMERICAL SUMMARY TABLE
# =============================================================================

# Create a nice summary table for reporting
final_summary <- summary_stats %>%
  select(n, vtype_factor, true_es_mean, 
         cohens_d_bias, cohens_d_mse, cohens_d_coverage, cohens_d_ci_width,
         robust_cohens_d_bias, robust_cohens_d_mse, robust_cohens_d_coverage, robust_cohens_d_ci_width,
         es_one_step_bias, es_one_step_mse, es_one_step_coverage, es_one_step_ci_width) %>%
  mutate(across(c(true_es_mean, cohens_d_bias, cohens_d_mse, robust_cohens_d_bias, robust_cohens_d_mse,
                  es_one_step_bias, es_one_step_mse, cohens_d_ci_width, robust_cohens_d_ci_width, es_one_step_ci_width), 
                ~ round(.x, 4)),
         across(c(cohens_d_coverage, robust_cohens_d_coverage, es_one_step_coverage), 
                ~ round(.x, 3))) %>%
  arrange(vtype_factor, n)

print("=== RCT PI FINAL SUMMARY TABLE ===")
print(final_summary)

write.csv(final_summary, paste0(save_dir, "/final_summary_table_rct_pi.csv"), row.names = FALSE)

# =============================================================================
# 8. KEY FINDINGS SUMMARY
# =============================================================================

cat("\n=== RCT PI KEY FINDINGS ===\n")

# Overall bias comparison
overall_bias <- all_results %>%
  summarise(
    cohens_d_abs_bias = mean(abs(cohens_d - true_es), na.rm = TRUE),
    robust_cohens_d_abs_bias = mean(abs(robust_cohens_d - true_es), na.rm = TRUE),
    es_one_step_abs_bias = mean(abs(es_one_step - true_es), na.rm = TRUE),
    cohens_d_coverage = mean(cohens_d_lb <= true_es & true_es <= cohens_d_ub, na.rm = TRUE),
    robust_cohens_d_coverage = mean(robust_cohens_d_lb <= true_es & true_es <= robust_cohens_d_ub, na.rm = TRUE),
    es_one_step_coverage = mean(es_one_step_lb <= true_es & true_es <= es_one_step_ub, na.rm = TRUE)
  )

cat("Cohen's d - Average absolute bias:", round(overall_bias$cohens_d_abs_bias, 4), "\n")
cat("Robust Cohen's d - Average absolute bias:", round(overall_bias$robust_cohens_d_abs_bias, 4), "\n")
cat("Causal ES - Average absolute bias:", round(overall_bias$es_one_step_abs_bias, 4), "\n")
cat("Cohen's d - Coverage probability:", round(overall_bias$cohens_d_coverage, 3), "\n")
cat("Robust Cohen's d - Coverage probability:", round(overall_bias$robust_cohens_d_coverage, 3), "\n")
cat("Causal ES - Coverage probability:", round(overall_bias$es_one_step_coverage, 3), "\n")

cat("\nRCT PI Files created:\n")
cat("- simulation_dashboard_rct_pi.png\n")
cat("- simulation_summary_statistics_rct_pi.csv\n")
cat("- final_summary_table_rct_pi.csv\n")

cat("\nRCT PI Analysis complete! 🎉\n")


