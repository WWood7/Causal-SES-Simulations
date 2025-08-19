# Visualization for RF hyperparameter tuning results


library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Args: results_dir [output_dir]
args <- commandArgs(trailingOnly = TRUE)
results_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/trial_results"

output_dir <- "/Users/winnwu/Emory_Projects/causal_effect_size/CausalEffectSize/simulation_results/trial_visualization"

if (!dir.exists(results_dir)) {
    stop(paste0("Results directory not found: ", results_dir))
}
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Parse seed, n, effect_size_type from filename
parse_meta_from_filename <- function(filepath) {
    fname <- basename(filepath)
    pattern <- "^result_seed([0-9]+)_n([0-9]+)_vtype([0-9]+)\\.rds$"
    m <- regexec(pattern, fname)
    mm <- regmatches(fname, m)
    if (length(mm) == 0 || length(mm[[1]]) != 4) return(NULL)
    vals <- mm[[1]]
    list(
        seed = as.numeric(vals[2]),
        n = as.numeric(vals[3]),
        effect_size_type = as.numeric(vals[4])
    )
}

files <- list.files(results_dir, pattern = "^result_seed[0-9]+_n[0-9]+_vtype[0-9]+\\.rds$", full.names = TRUE)
if (length(files) == 0) stop(paste0("No result files found in ", results_dir))

results_list <- lapply(files, function(f) {
    meta <- parse_meta_from_filename(f)
    if (is.null(meta)) return(NULL)
    df <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(df) || !is.data.frame(df)) return(NULL)
    # Expect columns: mtry, depth, bias
    if (!all(c("mtry", "depth", "bias") %in% names(df))) return(NULL)
    df$seed <- meta$seed
    df$n <- meta$n
    df$effect_size_type <- meta$effect_size_type
    df
})

rf <- bind_rows(results_list)
if (nrow(rf) == 0) stop("Parsed zero rows from RF result files.")

rf <- rf %>%
    mutate(
        n_factor = factor(n),
        depth_factor = factor(depth),
        mtry_factor = factor(mtry),
        effect_factor = factor(effect_size_type)
    )

# Summary by n, mtry, depth
summary_ndm <- rf %>%
    group_by(n, mtry, depth) %>%
    summarise(
        mean_bias = mean(bias, na.rm = TRUE),
        sd_bias = sd(bias, na.rm = TRUE),
        rmse = sqrt(mean(bias^2, na.rm = TRUE)),
        .groups = "drop"
    )
write.csv(summary_ndm, file.path(output_dir, "summary_bias_by_n_mtry_depth.csv"), row.names = FALSE)

# Boxplot of bias by depth, colored by mtry, faceted by n
p_box_bias <- ggplot(rf, aes(x = depth_factor, y = bias, fill = mtry_factor)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(outlier.alpha = 0.2) +
    facet_wrap(~ n_factor) +
    labs(x = "max_depth", y = "Bias (estimate - true)", fill = "mtry",
         title = "RF hyperparameters: bias by depth and mtry across seeds",
         subtitle = "Faceted by sample size n") +
    theme_bw()
ggsave(file.path(output_dir, "bias_boxplot_by_depth_mtry_facet_n.png"), p_box_bias, width = 12, height = 6, dpi = 300)

# Heatmap of RMSE by (mtry, depth) for each n
p_heat <- ggplot(summary_ndm, aes(x = factor(depth), y = factor(mtry), fill = rmse)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +
    facet_wrap(~ factor(n), labeller = label_both) +
    labs(x = "max_depth", y = "mtry", fill = "RMSE",
         title = "RMSE of bias by RF hyperparameters") +
    theme_bw()
ggsave(file.path(output_dir, "rmse_heatmap_by_mtry_depth_facet_n.png"), p_heat, width = 10, height = 6, dpi = 300)

# Best hyperparameters per n by RMSE
best_by_n <- summary_ndm %>%
    group_by(n) %>%
    slice_min(order_by = rmse, n = 1, with_ties = FALSE) %>%
    ungroup()
write.csv(best_by_n, file.path(output_dir, "best_hyperparams_by_n.csv"), row.names = FALSE)

message(paste0("Wrote RF plots and summaries to ", output_dir))


