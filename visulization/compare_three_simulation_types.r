# Comprehensive Comparison of Three Simulation Types
# Creates 4 separate images: Bias, MSE, Coverage, and CI Width
# Each image contains three panels comparing Regular, RCT, and RCT PI results

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

# =============================================================================
# SETUP DIRECTORIES AND PATHS
# =============================================================================

# Output directory for comparison plots
save_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/three_way_comparison"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# Input directories for the three simulation types
regular_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d"
rct_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct"
rct_pi_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/cohens_d_rct_pi"

cat("Loading and processing data from three simulation types...\n")

# =============================================================================
# HELPER FUNCTION TO LOAD AND PROCESS DATA
# =============================================================================

load_simulation_data <- function(results_dir, simulation_type) {
  cat(paste0("Loading ", simulation_type, " simulation results...\n"))
  
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
                                    labels = c("Small", "Medium", "Large", "Huge"))
  
  # Add simulation type identifier
  all_results$simulation_type <- simulation_type
  
  return(all_results)
}

# =============================================================================
# LOAD ALL THREE DATASETS
# =============================================================================

regular_data <- load_simulation_data(regular_dir, "Confounded")
rct_data <- load_simulation_data(rct_dir, "RCT (balanced)")
rct_pi_data <- load_simulation_data(rct_pi_dir, "RCT (imbalanced)")

# Combine all data
combined_data <- bind_rows(regular_data, rct_data, rct_pi_data)

# Make simulation_type a factor with proper ordering
# Order panels for clarity in the article figures
# Requested order (top to bottom): RCT (balanced), RCT (imbalanced), Confounded
combined_data$simulation_type <- factor(
  combined_data$simulation_type,
  levels = c("RCT (balanced)", "RCT (imbalanced)", "Confounded")
)

cat("Combined data from all three simulation types\n")
cat("Total rows:", nrow(combined_data), "\n")

# =============================================================================
# COMPUTE SUMMARY STATISTICS FOR ALL THREE TYPES
# =============================================================================

summary_stats <- combined_data %>%
  group_by(simulation_type, n, effect_size_type, vtype_factor, n_factor) %>%
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
    
    # Plugin causal ES statistics  
    es_plugin_mean = mean(es_plugin, na.rm = TRUE),
    es_plugin_bias = mean(es_plugin - true_es, na.rm = TRUE),
    es_plugin_mse = mean((es_plugin - true_es)^2, na.rm = TRUE),
    
    # One-step causal ES statistics (available in all three types)
    es_one_step_mean = mean(es_one_step, na.rm = TRUE),
    es_one_step_bias = mean(es_one_step - true_es, na.rm = TRUE),
    es_one_step_mse = mean((es_one_step - true_es)^2, na.rm = TRUE),
    es_one_step_coverage = mean(es_one_step_lb <= true_es & true_es <= es_one_step_ub, na.rm = TRUE),
    es_one_step_ci_width = mean(es_one_step_ub - es_one_step_lb, na.rm = TRUE),
    
    .groups = "drop"
  )

# Save combined summary statistics
write.csv(summary_stats, paste0(save_dir, "/combined_summary_statistics.csv"), row.names = FALSE)

# =============================================================================
# PLOTTING THEME
# =============================================================================

my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(linewidth = 0.6),
    axis.ticks.length = unit(3, "pt"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  )

# =============================================================================
# 1. BIAS COMPARISON PLOT
# =============================================================================

cat("Creating bias comparison plot...\n")

# Prepare bias data for plotting
bias_data <- combined_data %>%
  mutate(
    cohens_d_bias = cohens_d - true_es,
    robust_cohens_d_bias = robust_cohens_d - true_es,
    es_plugin_bias = es_plugin - true_es,
    es_one_step_bias = es_one_step - true_es
  ) %>%
  select(simulation_type, n_factor, vtype_factor, 
         cohens_d_bias, robust_cohens_d_bias, es_plugin_bias, es_one_step_bias) %>%
  pivot_longer(cols = c(cohens_d_bias, robust_cohens_d_bias, es_plugin_bias, es_one_step_bias),
               names_to = "estimator", values_to = "bias") %>%
  filter(!is.na(bias)) %>%  # Handle any missing plugin data
  filter(estimator != "es_plugin_bias") %>%
  mutate(estimator = factor(estimator, 
                           levels = c("cohens_d_bias", "robust_cohens_d_bias", "es_one_step_bias"),
                           labels = c("Cohen's d", "rescaled RESI", "One-step causal Cohen's d")))

# Create bias plot
p_bias <- ggplot(bias_data, aes(x = n_factor, y = bias, fill = estimator)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_grid(simulation_type ~ vtype_factor) +
  scale_fill_viridis_d(name = "Estimator") +
  labs(title = "Bias Comparison Across Three Simulation Types",
       subtitle = "Bias = Estimate - True Effect Size",
       x = "Sample Size", 
       y = "Bias") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save bias plot
ggsave(paste0(save_dir, "/bias_comparison_three_types.png"), p_bias, 
       width = 16, height = 12, dpi = 300)

# =============================================================================
# 2. MSE COMPARISON PLOT
# =============================================================================

cat("Creating MSE comparison plot...\n")

# Prepare MSE data
mse_data <- summary_stats %>%
  select(simulation_type, n_factor, vtype_factor, 
         cohens_d_mse, robust_cohens_d_mse, es_plugin_mse, es_one_step_mse) %>%
  pivot_longer(cols = c(cohens_d_mse, robust_cohens_d_mse, es_plugin_mse, es_one_step_mse),
               names_to = "estimator", values_to = "mse") %>%
  filter(!is.na(mse)) %>%  # Handle any missing plugin data
  filter(estimator != "es_plugin_mse") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_mse", "robust_cohens_d_mse", "es_one_step_mse"),
                           labels = c("Cohen's d", "rescaled RESI", "One-step causal Cohen's d")))

# Create MSE plot
p_mse <- ggplot(
  mse_data,
  aes(x = n_factor, y = mse, color = estimator, linetype = estimator, shape = estimator, group = estimator)
) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black", alpha = 0.9) +
  facet_grid(simulation_type ~ vtype_factor, scales = "free_y", switch = "y") +
  scale_color_viridis_d(name = "Estimator") +
  scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotdash"), guide = "none") +
  scale_shape_manual(values = c(16, 17, 15, 18), guide = "none") +
  labs(
    title = "Mean Squared Error Comparison Across Three Simulation Types",
    x = "Sample Size",
    y = "MSE"
  ) +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Save MSE plot
ggsave(paste0(save_dir, "/mse_comparison_three_types.png"), p_mse, width = 12, height = 8, dpi = 300)
ggsave(paste0(save_dir, "/mse_comparison_three_types.pdf"), p_mse, width = 12, height = 8)

# =============================================================================
# 3. COVERAGE COMPARISON PLOT
# =============================================================================

cat("Creating coverage comparison plot...\n")

# Prepare coverage data
coverage_data <- summary_stats %>%
  select(simulation_type, n_factor, vtype_factor, 
         cohens_d_coverage, robust_cohens_d_coverage, es_one_step_coverage) %>%
  pivot_longer(cols = c(cohens_d_coverage, robust_cohens_d_coverage, es_one_step_coverage),
               names_to = "estimator", values_to = "coverage") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_coverage", "robust_cohens_d_coverage", "es_one_step_coverage"),
                           labels = c("Cohen's d", "rescaled RESI", "One-step causal Cohen's d")))

# Calculate dynamic y-axis limits
min_coverage <- min(coverage_data$coverage, na.rm = TRUE)
max_coverage <- max(coverage_data$coverage, na.rm = TRUE)
y_lower <- max(0, min_coverage - 0.05)
y_upper <- min(1, max_coverage + 0.02)

# Create coverage plot
p_coverage <- ggplot(coverage_data, aes(x = n_factor, y = coverage, color = estimator, group = estimator)) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_point(size = 2.5, stroke = 0.4) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black", alpha = 0.9) +
  facet_grid(simulation_type ~ vtype_factor, switch = "y") +
  scale_color_viridis_d(name = "Estimator") +
  scale_y_continuous(limits = c(0, y_upper), labels = scales::percent) +
  labs(title = "95% Confidence Interval Coverage Comparison",
       subtitle = "Should be close to 95% (red dashed line)",
       x = "Sample Size", 
       y = "Coverage Probability") +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = NA, color = NA),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# Save coverage plot
ggsave(paste0(save_dir, "/coverage_comparison_three_types.png"), p_coverage, 
       width = 16, height = 12, dpi = 300)

# =============================================================================
# 4. CI WIDTH COMPARISON PLOT
# =============================================================================

cat("Creating CI width comparison plot...\n")

# Prepare CI width data
ci_width_data <- summary_stats %>%
  select(simulation_type, n_factor, vtype_factor, 
         cohens_d_ci_width, robust_cohens_d_ci_width, es_one_step_ci_width) %>%
  pivot_longer(cols = c(cohens_d_ci_width, robust_cohens_d_ci_width, es_one_step_ci_width),
               names_to = "estimator", values_to = "ci_width") %>%
  mutate(estimator = factor(estimator,
                           levels = c("cohens_d_ci_width", "robust_cohens_d_ci_width", "es_one_step_ci_width"),
                           labels = c("Cohen's d", "rescaled RESI", "One-step Causal ES")))

# Create CI width plot
p_ci_width <- ggplot(ci_width_data, aes(x = n_factor, y = ci_width, color = estimator, group = estimator)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  facet_grid(simulation_type ~ vtype_factor, scales = "free_y") +
  scale_color_viridis_d(name = "Estimator") +
  labs(title = "95% Confidence Interval Width Comparison",
       subtitle = "Narrower intervals indicate higher precision",
       x = "Sample Size", 
       y = "CI Width") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save CI width plot
ggsave(paste0(save_dir, "/ci_width_comparison_three_types.png"), p_ci_width, 
       width = 16, height = 12, dpi = 300)

# =============================================================================
# SUMMARY AND COMPLETION
# =============================================================================

cat("\n=== THREE-WAY COMPARISON COMPLETE ===\n")
cat("Created 4 comparison plots:\n")
cat("1. bias_comparison_three_types.png\n")
cat("2. mse_comparison_three_types.png\n")
cat("3. coverage_comparison_three_types.png\n")
cat("4. ci_width_comparison_three_types.png\n")
cat("\nAll plots saved in:", save_dir, "\n")
cat("Combined summary statistics saved as: combined_summary_statistics.csv\n")

# Print summary of what we compared
cat("\nSimulation types compared:\n")
cat("- RCT (balanced): Randomized controlled trial with balanced treatment assignment\n")
cat("- RCT (imbalanced): Randomized controlled trial with unbalanced treatment assignment\n")
cat("- Confounded: Observational data simulation with confounding\n")

cat("\nEach plot shows three panels (by simulation type) with four effect size conditions.\n")
cat("Analysis complete! 🎉\n")
