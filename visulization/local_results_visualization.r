# Visualize results from the local causal effect-size simulation.

required_packages <- c("ggplot2", "dplyr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1])))
} else {
  getwd()
}
repository_directory <- dirname(script_directory)

# Optional command-line arguments:
# 1. input RDS file
# 2. output directory
arguments <- commandArgs(trailingOnly = TRUE)
default_input <- file.path(
  repository_directory,
  "local",
  "results",
  "simulation_results.rds"
)
input_file <- if (length(arguments) >= 1) arguments[1] else default_input
output_directory <- if (length(arguments) >= 2) {
  arguments[2]
} else {
  file.path(repository_directory, "local", "results", "figures")
}

if (!file.exists(input_file)) {
  stop("Simulation results file does not exist: ", input_file)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

results <- readRDS(input_file)
required_columns <- c(
  "trial", "n", "effect_size_type", "propensity_score_type",
  "true_es", "true_resi", "cohens_d", "cohens_d_lb", "cohens_d_ub",
  "causal_es", "causal_es_lb", "causal_es_ub",
  "causal_resi", "causal_resi_lb", "causal_resi_ub", "error"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop(
    "Simulation results are missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

effect_labels <- c(
  "1" = "Small effect",
  "2" = "Medium effect",
  "3" = "Large effect",
  "4" = "Very large effect"
)
propensity_labels <- c(
  "1" = "Randomized: P(A=1) = 0.5",
  "2" = "Randomized: P(A=1) = 0.3",
  "3" = "Confounded treatment"
)
sample_size_levels <- sort(unique(results$n))

results <- results %>%
  mutate(
    effect_type = factor(
      effect_size_type,
      levels = as.integer(names(effect_labels)),
      labels = unname(effect_labels)
    ),
    propensity_type = factor(
      propensity_score_type,
      levels = as.integer(names(propensity_labels)),
      labels = unname(propensity_labels)
    ),
    sample_size = factor(
      n,
      levels = sample_size_levels,
      labels = paste0("n = ", sample_size_levels)
    )
  )

successful_results <- results %>%
  filter(is.na(error))

if (nrow(successful_results) == 0) {
  stop("No successful simulation trials are available for visualization")
}

coverage_rate <- function(lower, truth, upper) {
  valid <- is.finite(lower) & is.finite(truth) & is.finite(upper)
  if (!any(valid)) {
    return(NA_real_)
  }
  mean(lower[valid] <= truth[valid] & truth[valid] <= upper[valid])
}

coverage_summary <- results %>%
  group_by(
    n, sample_size, effect_size_type, effect_type,
    propensity_score_type, propensity_type
  ) %>%
  summarise(
    n_attempted = n(),
    n_fit_failed = sum(!is.na(error)),
    n_success = sum(is.na(error)),
    cohens_d_valid_ci = sum(
      is.na(error) & is.finite(cohens_d_lb) & is.finite(cohens_d_ub)
    ),
    causal_es_valid_ci = sum(
      is.na(error) & is.finite(causal_es_lb) & is.finite(causal_es_ub)
    ),
    causal_resi_valid_ci = sum(
      is.na(error) &
        is.finite(causal_resi_lb) & is.finite(causal_resi_ub)
    ),
    cohens_d_coverage = coverage_rate(
      cohens_d_lb[is.na(error)],
      true_es[is.na(error)],
      cohens_d_ub[is.na(error)]
    ),
    causal_es_coverage = coverage_rate(
      causal_es_lb[is.na(error)],
      true_es[is.na(error)],
      causal_es_ub[is.na(error)]
    ),
    causal_resi_coverage = coverage_rate(
      causal_resi_lb[is.na(error)],
      true_resi[is.na(error)],
      causal_resi_ub[is.na(error)]
    ),
    .groups = "drop"
  )

coverage_output <- coverage_summary %>%
  select(
    n, effect_size_type, propensity_score_type,
    n_attempted, n_fit_failed, n_success,
    cohens_d_valid_ci, causal_es_valid_ci, causal_resi_valid_ci,
    cohens_d_coverage, causal_es_coverage, causal_resi_coverage
  )
write.csv(
  coverage_output,
  file.path(output_directory, "coverage_summary.csv"),
  row.names = FALSE
)

coverage_plot_data <- coverage_summary %>%
  select(
    n, sample_size, effect_type, propensity_type,
    cohens_d_coverage, causal_es_coverage, causal_resi_coverage
  ) %>%
  pivot_longer(
    cols = c(
      cohens_d_coverage,
      causal_es_coverage,
      causal_resi_coverage
    ),
    names_to = "estimator",
    values_to = "coverage"
  ) %>%
  mutate(
    estimator = factor(
      estimator,
      levels = c(
        "cohens_d_coverage",
        "causal_es_coverage",
        "causal_resi_coverage"
      ),
      labels = c("Cohen's d", "Causal effect size", "Causal RESI")
    )
  )

plot_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.7, "lines"),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

estimator_colors <- c(
  "Cohen's d" = "#0072B2",
  "Causal effect size" = "#D55E00",
  "Causal RESI" = "#009E73"
)

coverage_plot <- ggplot(
  coverage_plot_data,
  aes(x = sample_size, y = coverage, color = estimator, group = estimator)
) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_line(linewidth = 0.75, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  facet_grid(effect_type ~ propensity_type) +
  scale_color_manual(values = estimator_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Confidence-interval coverage across simulation scenarios",
    subtitle = "Undefined intervals and failed fits are excluded; dashed line is 95% nominal coverage",
    x = "Sample size",
    y = "Coverage"
  ) +
  plot_theme

estimate_plot_data <- successful_results %>%
  select(
    sample_size, effect_type, propensity_type, true_es,
    cohens_d, causal_es
  ) %>%
  pivot_longer(
    cols = c(cohens_d, causal_es),
    names_to = "estimator",
    values_to = "estimate"
  ) %>%
  filter(is.finite(estimate)) %>%
  mutate(
    estimator = factor(
      estimator,
      levels = c("cohens_d", "causal_es"),
      labels = c("Cohen's d", "Causal effect size")
    )
  ) %>%
  group_by(
    sample_size, effect_type, propensity_type, true_es, estimator
  ) %>%
  summarise(
    estimate = mean(estimate),
    .groups = "drop"
  )

true_es_lines <- successful_results %>%
  distinct(effect_type, propensity_type, true_es)

estimate_plot <- ggplot(
  estimate_plot_data,
  aes(x = sample_size, y = estimate, color = estimator, group = estimator)
) +
  geom_hline(
    data = true_es_lines,
    aes(yintercept = true_es),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey25"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(effect_type ~ propensity_type, scales = "free_y") +
  scale_color_manual(
    values = estimator_colors[c("Cohen's d", "Causal effect size")]
  ) +
  labs(
    title = "Mean Cohen's d and causal effect-size estimates",
    subtitle = "Points are trial means; dashed line is the true causal effect size",
    x = "Sample size",
    y = "Mean estimate"
  ) +
  plot_theme

allocation_plot_data <- successful_results %>%
  filter(is.finite(causal_es), is.finite(causal_resi))

allocation_means <- allocation_plot_data %>%
  group_by(
    sample_size, propensity_type, effect_size_type, effect_type,
    true_es, true_resi
  ) %>%
  summarise(
    causal_es = mean(causal_es),
    causal_resi = mean(causal_resi),
    .groups = "drop"
  )

allocation_x_limits <- quantile(
  allocation_plot_data$causal_es,
  probs = c(0.005, 0.995),
  na.rm = TRUE,
  names = FALSE
)
allocation_x_limits <- range(
  allocation_x_limits,
  allocation_means$true_es
)

propensity_colors <- c(
  "Randomized: P(A=1) = 0.5" = "#0072B2",
  "Randomized: P(A=1) = 0.3" = "#D55E00",
  "Confounded treatment" = "#009E73"
)

allocation_plot <- ggplot(
  allocation_plot_data,
  aes(
    x = causal_es,
    y = causal_resi,
    color = propensity_type,
    shape = propensity_type
  )
) +
  geom_point(alpha = 0.18, size = 0.9) +
  geom_point(
    data = allocation_means,
    size = 2.7,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  geom_point(
    data = allocation_means,
    aes(x = true_es, y = true_resi),
    shape = 4,
    size = 2.8,
    stroke = 0.9
  ) +
  facet_grid(sample_size ~ effect_type) +
  scale_color_manual(values = propensity_colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  coord_cartesian(xlim = allocation_x_limits) +
  labs(
    title = "Treatment-allocation sensitivity of causal RESI",
    subtitle = paste(
      "Transparent points are trials, circles are scenario means,",
      "and crosses are population targets; x-axis zooms to the",
      "0.5th–99.5th percentile"
    ),
    x = "Estimated causal effect size",
    y = "Estimated causal RESI",
    color = "Treatment assignment",
    shape = "Treatment assignment"
  ) +
  plot_theme

standardized_rmse_summary <- successful_results %>%
  group_by(
    n, sample_size, effect_size_type, effect_type,
    propensity_score_type, propensity_type, true_es, true_resi
  ) %>%
  summarise(
    causal_es_n = sum(is.finite(causal_es)),
    causal_resi_n = sum(is.finite(causal_resi)),
    causal_es_standardized_rmse = sqrt(
      mean((causal_es - true_es) ^ 2, na.rm = TRUE)
    ) / abs(first(true_es)),
    causal_resi_standardized_rmse = sqrt(
      mean((causal_resi - true_resi) ^ 2, na.rm = TRUE)
    ) / abs(first(true_resi)),
    .groups = "drop"
  )

write.csv(
  standardized_rmse_summary %>%
    select(
      n, effect_size_type, propensity_score_type,
      true_es, true_resi, causal_es_n, causal_resi_n,
      causal_es_standardized_rmse, causal_resi_standardized_rmse
    ),
  file.path(output_directory, "standardized_rmse_summary.csv"),
  row.names = FALSE
)

standardized_rmse_plot_data <- standardized_rmse_summary %>%
  select(
    sample_size, effect_type, propensity_type,
    causal_es_standardized_rmse,
    causal_resi_standardized_rmse
  ) %>%
  pivot_longer(
    cols = c(
      causal_es_standardized_rmse,
      causal_resi_standardized_rmse
    ),
    names_to = "estimator",
    values_to = "standardized_rmse"
  ) %>%
  filter(is.finite(standardized_rmse), standardized_rmse > 0) %>%
  mutate(
    estimator = factor(
      estimator,
      levels = c(
        "causal_es_standardized_rmse",
        "causal_resi_standardized_rmse"
      ),
      labels = c("Causal effect size", "Causal RESI")
    )
  )

standardized_rmse_plot <- ggplot(
  standardized_rmse_plot_data,
  aes(
    x = sample_size,
    y = standardized_rmse,
    color = estimator,
    group = estimator
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(effect_type ~ propensity_type) +
  scale_color_manual(
    values = estimator_colors[c("Causal effect size", "Causal RESI")]
  ) +
  scale_y_log10(
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(
    title = "Standardized RMSE of the causal estimators",
    subtitle = paste(
      "Standardized RMSE = RMSE / |true value|;",
      "logarithmic y-axis and dashed reference at 1"
    ),
    x = "Sample size",
    y = "Standardized root mean squared error"
  ) +
  plot_theme

save_plot <- function(plot, filename, width, height) {
  ggsave(
    file.path(output_directory, paste0(filename, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(output_directory, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
}

save_plot(coverage_plot, "coverage_comparison", width = 15, height = 12)
save_plot(estimate_plot, "cohens_d_vs_causal_es", width = 15, height = 12)
save_plot(allocation_plot, "resi_allocation_sensitivity", width = 15, height = 12)
save_plot(
  standardized_rmse_plot,
  "standardized_rmse_comparison",
  width = 15,
  height = 12
)

cat(
  "Created four visualizations and summary CSV files in:\n",
  normalizePath(output_directory),
  "\n"
)
