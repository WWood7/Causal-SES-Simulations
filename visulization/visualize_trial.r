# Visualization for RCT simulation results


library(ggplot2)
library(dplyr)
library(tidyr)


results_dir <- 
    "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/trial_results"
output_dir <- 
    "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/trial_visualization"

# Helper to parse seed, n, effect_size_type from filename
parse_meta_from_filename <- function(filepath) {
    fname <- basename(filepath)
    # expected: result_seed<seed>_n<n>_vtype<effect>.rds
    pattern <- "^result_seed([0-9]+)_n([0-9]+)_vtype([0-9]+)\\.rds$"
    m <- regexec(pattern, fname)
    mm <- regmatches(fname, m)
    if (length(mm) == 0 || length(mm[[1]]) != 4) return(NULL)
    vals <- mm[[1]]
    return(list(
        seed = as.numeric(vals[2]),
        n = as.numeric(vals[3]),
        effect_size_type = as.numeric(vals[4])
    ))
}

# Load all result files
files <- list.files(results_dir, pattern = "^result_seed[0-9]+_n[0-9]+_vtype[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) {
    stop(paste0("No result files found in ", results_dir))
}

results_list <- lapply(files, function(f) {
    meta <- parse_meta_from_filename(f)
    if (is.null(meta)) return(NULL)
    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res) || !is.data.frame(res)) return(NULL)
    res$seed <- meta$seed
    res$n <- meta$n
    res$effect_size_type <- meta$effect_size_type
    # true ES by effect_size_type: c(0.209, 0.557, 0.876)
    res$true_es <- c(0.209, 0.557, 0.876)[res$effect_size_type][1]
    res
})

results <- bind_rows(results_list)
if (nrow(results) == 0) {
    stop("Parsed zero result rows from files.")
}

results <- results %>%
    mutate(
        plugin_bias = es_plugin - true_es,
        plugin_rf_bias = es_plugin_rf - true_es,
        plugin_sq_err = (es_plugin - true_es)^2,
        plugin_rf_sq_err = (es_plugin_rf - true_es)^2,
        n_factor = factor(n),
        effect_factor = factor(effect_size_type)
    )

# Save summary table by n and effect size type
summary_by_n_vtype <- results %>%
    group_by(n, effect_size_type) %>%
    summarise(
        mean_es_plugin = mean(es_plugin),
        mean_es_plugin_rf = mean(es_plugin_rf),
        sd_es_plugin = sd(es_plugin),
        sd_es_plugin_rf = sd(es_plugin_rf),
        bias_plugin = mean(plugin_bias),
        bias_plugin_rf = mean(plugin_rf_bias),
        rmse_plugin = sqrt(mean(plugin_sq_err)),
        rmse_plugin_rf = sqrt(mean(plugin_rf_sq_err)),
        .groups = "drop"
    )

write.csv(summary_by_n_vtype, file.path(output_dir, "summary_by_n_vtype.csv"), row.names = FALSE)

# Prepare long format for ES plots
plot_df <- bind_rows(
    results %>% transmute(n_factor, effect_factor, method = "SL", es = es_plugin, true_es),
    results %>% transmute(n_factor, effect_factor, method = "RF", es = es_plugin_rf, true_es)
)

true_lines <- plot_df %>% distinct(effect_factor, true_es)

# Boxplot of estimated ES by n and method, faceted by effect size type
p_box <- ggplot(plot_df, aes(x = n_factor, y = es, fill = method)) +
    geom_boxplot(outlier.alpha = 0.2) +
    geom_hline(data = true_lines, aes(yintercept = true_es), linetype = "dashed", color = "red", inherit.aes = FALSE) +
    facet_wrap(~ effect_factor, labeller = labeller(effect_factor = label_both)) +
    labs(x = "Sample size (n)", y = "Estimated effect size", fill = "Method",
         title = "Estimated effect sizes by n and method",
         subtitle = "Dashed red line = true effect size per effect_size_type") +
    theme_bw()
ggsave(filename = file.path(output_dir, "boxplot_es_by_n_method.png"), plot = p_box, width = 10, height = 6, dpi = 300)

# Bias boxplot (ES - true)
plot_bias <- plot_df %>% mutate(bias = es - true_es)
p_bias <- ggplot(plot_bias, aes(x = n_factor, y = bias, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~ effect_factor, labeller = labeller(effect_factor = label_both)) +
    labs(x = "Sample size (n)", y = "Bias (estimate - true)", fill = "Method",
         title = "Bias of estimated effect sizes by n and method") +
    theme_bw()
ggsave(filename = file.path(output_dir, "boxplot_bias_by_n_method.png"), plot = p_bias, width = 10, height = 6, dpi = 300)

# Scatter: SL vs RF
p_scatter <- ggplot(results, aes(x = es_plugin, y = es_plugin_rf, color = n_factor)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.5) +
    facet_wrap(~ effect_factor, labeller = labeller(effect_factor = label_both)) +
    labs(x = "SL estimate (es_plugin)", y = "RF estimate (es_plugin_rf)", color = "n",
         title = "SL vs RF plug-in estimates") +
    theme_bw()
ggsave(filename = file.path(output_dir, "scatter_sl_vs_rf.png"), plot = p_scatter, width = 10, height = 6, dpi = 300)

# Learner coefficient summaries
coef_outcome <- results %>%
    select(effect_factor, starts_with("coef_reg_")) %>%
    group_by(effect_factor) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(cols = starts_with("coef_reg_"), names_to = "learner", values_to = "mean_coef") %>%
    mutate(task = "Outcome (y)")

coef_outcome_sq <- results %>%
    select(effect_factor, starts_with("coef_regsq_")) %>%
    group_by(effect_factor) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(cols = starts_with("coef_regsq_"), names_to = "learner", values_to = "mean_coef") %>%
    mutate(task = "Outcome squared (y_sq)")

coef_long <- bind_rows(coef_outcome, coef_outcome_sq) %>%
    mutate(learner = gsub("^coef_regsq_|^coef_reg_", "", learner))

p_coef <- ggplot(coef_long, aes(x = learner, y = mean_coef, fill = effect_factor)) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_wrap(~ task, ncol = 1, scales = "free_y") +
    labs(x = "Learner", y = "Mean SL coefficient", fill = "effect_size_type",
         title = "Average SuperLearner coefficients across simulations") +
    theme_bw()
ggsave(filename = file.path(output_dir, "mean_sl_coefficients.png"), plot = p_coef, width = 10, height = 8, dpi = 300)

# Also write the aggregated results to CSV for ad-hoc analyses
write.csv(results, file.path(output_dir, "all_results_flat.csv"), row.names = FALSE)

message(paste0("Wrote plots and summaries to ", output_dir))


