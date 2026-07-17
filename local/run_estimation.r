library(parallel)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_argument[1])))
} else {
  getwd()
}
source(file.path(script_directory, "data_generation.r"))
source(file.path(script_directory, "estimation_functions.r"))

# Command-line arguments:
# 1. number of trials per scenario (default: 100)
# 2. number of parallel workers (default: detected cores minus one)
# 3. output RDS path (default: local/results/simulation_results.rds)
arguments <- commandArgs(trailingOnly = TRUE)
n_trials <- if (length(arguments) >= 1) as.integer(arguments[1]) else 100L
detected_cores <- parallel::detectCores()
default_cores <- if (is.na(detected_cores)) 1L else max(1L, detected_cores - 1L)
n_cores <- if (length(arguments) >= 2) as.integer(arguments[2]) else default_cores
default_output <- file.path(script_directory, "results", "simulation_results.rds")
output_file <- if (length(arguments) >= 3) arguments[3] else default_output

if (is.na(n_trials) || n_trials < 1) {
  stop("n_trials must be a positive integer")
}
if (is.na(n_cores) || n_cores < 1) {
  stop("n_cores must be a positive integer")
}

sample_sizes <- c(100L, 500L, 1000L, 2000L)
effect_size_types <- 1:4
propensity_score_types <- 1:3
true_effect_sizes <- c(0.209, 0.557, 0.876, 1.868)
# Rows are effect-size types and columns are propensity-score types.
# Values use the known DGP and the population variance of the ATE EIF.
true_resi <- matrix(
  c(
    0.104807527, 0.081307998, 0.089359379,
    0.285567719, 0.239021586, 0.239141267,
    0.471015216, 0.370474379, 0.409765360,
    1.517646588, 1.418305385, 1.358407610
  ),
  nrow = 4,
  byrow = TRUE
)

simulation_grid <- expand.grid(
  trial = seq_len(n_trials),
  n = sample_sizes,
  effect_size_type = effect_size_types,
  propensity_score_type = propensity_score_types,
  KEEP.OUT.ATTRS = FALSE
)

run_single_trial <- function(row_index) {
  scenario <- simulation_grid[row_index, ]
  set.seed(scenario$trial)

  base_result <- list(
    trial = scenario$trial,
    n = scenario$n,
    effect_size_type = scenario$effect_size_type,
    propensity_score_type = scenario$propensity_score_type,
    true_es = true_effect_sizes[scenario$effect_size_type],
    true_resi = true_resi[
      scenario$effect_size_type,
      scenario$propensity_score_type
    ]
  )

  tryCatch({
    simulated_data <- generate_data(
      n = scenario$n,
      effect_size_type = scenario$effect_size_type,
      propensity_score_type = scenario$propensity_score_type
    )
    cohens_result <- estimate_cohens_d(simulated_data)
    causal_result <- estimate_ces_and_resi(simulated_data)

    data.frame(
      base_result,
      cohens_d = cohens_result$cohens_d,
      cohens_d_lb = cohens_result$cohens_d_lb,
      cohens_d_ub = cohens_result$cohens_d_ub,
      causal_es = causal_result$est,
      causal_es_se = causal_result$est_se,
      causal_es_lb = causal_result$est_lb,
      causal_es_ub = causal_result$est_ub,
      causal_resi = causal_result$resi,
      causal_resi_lb = causal_result$resi_lb,
      causal_resi_ub = causal_result$resi_ub,
      error = NA_character_
    )
  }, error = function(error) {
    data.frame(
      base_result,
      cohens_d = NA_real_,
      cohens_d_lb = NA_real_,
      cohens_d_ub = NA_real_,
      causal_es = NA_real_,
      causal_es_se = NA_real_,
      causal_es_lb = NA_real_,
      causal_es_ub = NA_real_,
      causal_resi = NA_real_,
      causal_resi_lb = NA_real_,
      causal_resi_ub = NA_real_,
      error = conditionMessage(error)
    )
  })
}

n_cores <- min(n_cores, nrow(simulation_grid))
cat(
  "Running", nrow(simulation_grid), "trials across",
  nrow(unique(simulation_grid[c("n", "effect_size_type", "propensity_score_type")])),
  "scenarios with", n_cores, "workers...\n"
)

results_list <- parallel::mclapply(
  seq_len(nrow(simulation_grid)),
  run_single_trial,
  mc.cores = n_cores,
  mc.set.seed = FALSE
)
results <- do.call(rbind, results_list)

output_directory <- dirname(output_file)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, output_file)

successful_results <- results[is.na(results$error), , drop = FALSE]
if (nrow(successful_results) > 0) {
  scenario_groups <- split(
    successful_results,
    list(
      successful_results$n,
      successful_results$effect_size_type,
      successful_results$propensity_score_type
    ),
    drop = TRUE
  )
  summary_results <- do.call(rbind, lapply(scenario_groups, function(group) {
    data.frame(
      n = group$n[1],
      effect_size_type = group$effect_size_type[1],
      propensity_score_type = group$propensity_score_type[1],
      true_es = group$true_es[1],
      true_resi = group$true_resi[1],
      n_success = nrow(group),
      causal_es_mean = mean(group$causal_es),
      causal_es_bias = mean(group$causal_es - group$true_es),
      causal_es_mse = mean((group$causal_es - group$true_es) ^ 2),
      causal_es_coverage = mean(
        group$causal_es_lb <= group$true_es &
          group$true_es <= group$causal_es_ub
      ),
      causal_resi_mean = mean(group$causal_resi),
      causal_resi_bias = mean(group$causal_resi - group$true_resi),
      causal_resi_mse = mean((group$causal_resi - group$true_resi) ^ 2),
      causal_resi_coverage = mean(
        group$causal_resi_lb <= group$true_resi &
          group$true_resi <= group$causal_resi_ub
      )
    )
  }))
  rownames(summary_results) <- NULL
  summary_file <- sub("\\.rds$", "_summary.csv", output_file, ignore.case = TRUE)
  if (identical(summary_file, output_file)) {
    summary_file <- paste0(output_file, "_summary.csv")
  }
  write.csv(summary_results, summary_file, row.names = FALSE)
  print(summary_results)
}

cat(
  "Completed:", nrow(successful_results), "successful and",
  sum(!is.na(results$error)), "failed trials.\n",
  "Results:", normalizePath(output_file), "\n"
)