# Illustrate why Glass's delta ignores changes in Var(Y(1)).

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2 before running this script")
}

library(ggplot2)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1])))
} else {
  getwd()
}
repository_directory <- dirname(script_directory)

arguments <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(arguments) >= 1) {
  arguments[1]
} else {
  file.path(repository_directory, "local", "results", "figures")
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

mu_0 <- 0
mu_1 <- 2
ate <- mu_1 - mu_0
sigma_0 <- 1

scenario_parameters <- data.frame(
  scenario = factor(
    c("Scenario 1: smaller Var(Y(1))", "Scenario 2: larger Var(Y(1))"),
    levels = c(
      "Scenario 1: smaller Var(Y(1))",
      "Scenario 2: larger Var(Y(1))"
    )
  ),
  sigma_1 = c(1, 3)
)
scenario_parameters$glass_delta <- ate / sigma_0
scenario_parameters$causal_es <- ate / sqrt(
  0.5 * sigma_0 ^ 2 + 0.5 * scenario_parameters$sigma_1 ^ 2
)
scenario_parameters$ate_start <- mu_0
scenario_parameters$ate_end <- mu_1
scenario_parameters$ate_midpoint <- (mu_0 + mu_1) / 2
scenario_parameters$ate_label <- sprintf("ATE = %.0f", ate)

x_values <- seq(-7, 11, length.out = 1200)
density_data <- do.call(
  rbind,
  lapply(seq_len(nrow(scenario_parameters)), function(index) {
    scenario <- scenario_parameters[index, ]
    rbind(
      data.frame(
        scenario = scenario$scenario,
        potential_outcome = "Y(0)",
        x = x_values,
        density = dnorm(x_values, mean = mu_0, sd = sigma_0)
      ),
      data.frame(
        scenario = scenario$scenario,
        potential_outcome = "Y(1)",
        x = x_values,
        density = dnorm(x_values, mean = mu_1, sd = scenario$sigma_1)
      )
    )
  })
)

density_data$potential_outcome <- factor(
  density_data$potential_outcome,
  levels = c("Y(0)", "Y(1)")
)

scenario_annotations <- transform(
  scenario_parameters,
  label = sprintf(
    "SD[Y(0)] = %.0f    SD[Y(1)] = %.0f\nCausal Effect Size (beta=0.5) = %.2f\nCausal Effect Size (beta=0.5) = %.2f",
    sigma_0,
    sigma_1,
    glass_delta,
    causal_es
  )
)

mean_annotations <- do.call(
  rbind,
  lapply(scenario_parameters$scenario, function(scenario) {
    data.frame(
      scenario = scenario,
      x = c(mu_0, mu_1),
      y = c(0.43, 0.43),
      label = c(
        "mu[0] == E[Y(0)]",
        "mu[1] == E[Y(1)]"
      )
    )
  })
)

plot_colors <- c("Y(0)" = "#7FBF3F", "Y(1)" = "#F87621")

illustration <- ggplot(
  density_data,
  aes(x = x, y = density, color = potential_outcome)
) +
  geom_area(
    aes(fill = potential_outcome),
    position = "identity",
    alpha = 0.12,
    color = NA
  ) +
  geom_line(linewidth = 1.25) +
  geom_vline(
    xintercept = c(mu_0, mu_1),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_segment(
    data = scenario_parameters,
    aes(
      x = ate_start,
      xend = ate_end,
      y = -0.035,
      yend = -0.035
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    arrow = arrow(
      ends = "both",
      type = "closed",
      length = grid::unit(0.08, "inches")
    )
  ) +
  geom_text(
    data = scenario_parameters,
    aes(
      x = ate_midpoint,
      y = -0.063,
      label = ate_label
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.6
  ) +
  geom_text(
    data = mean_annotations,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    size = 3.5
  ) +
  geom_text(
    data = scenario_annotations,
    aes(x = -6.7, y = 0.35, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    lineheight = 1.1,
    size = 3.6
  ) +
  facet_grid(scenario ~ ., scales = "fixed") +
  scale_color_manual(values = plot_colors) +
  scale_fill_manual(values = plot_colors) +
  coord_cartesian(xlim = c(-7, 11), ylim = c(-0.08, 0.47), clip = "off") +
  labs(
    x = "Potential outcome value",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11),
    plot.caption = element_text(hjust = 0, size = 10),
    strip.text.y = element_text(face = "bold", size = 11),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.spacing = grid::unit(0.8, "lines"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

png_file <- file.path(
  output_directory,
  "glass_delta_variance_illustration.png"
)
pdf_file <- file.path(
  output_directory,
  "glass_delta_variance_illustration.pdf"
)

ggsave(
  png_file,
  illustration,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  pdf_file,
  illustration,
  width = 12,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

cat(
  "Created Glass's delta variance illustration:\n",
  normalizePath(png_file),
  "\n",
  normalizePath(pdf_file),
  "\n"
)
