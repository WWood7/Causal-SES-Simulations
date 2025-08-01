library(parallel)
source("data_generation.r")
source("estimation_functions.r")

# Simple parallel simulation function
n_trials <- 100
run_trials_parallel <- function(n_trials = 100, n = 200, variability_type = 1) {
  
  # Single trial function
  single_trial <- function(i) {
    data <- generate_data(n = n, variability_type = variability_type)
    cohens_result <- estimate_cohens_d(data)
    causal_result <- estimate_causal_es(data)
    
    # Calculate true value
    if (variability_type == 1) {
      true_es <- 0.528
    } else if (variability_type == 2) {
      true_es <- 1.626
    } else {
      true_es <- 1.752
    }
    
    return(data.frame(
      trial = i,
      true_es = true_es,
      cohens_d = cohens_result$cohens_d,
      cohens_d_lb = cohens_result$cohens_d_lb,
      cohens_d_ub = cohens_result$cohens_d_ub,
      es_plugin = causal_result$es_plugin,
      es_one_step = causal_result$es_one_step,
      es_one_step_lb = causal_result$es_one_step_lb,
      es_one_step_ub = causal_result$es_one_step_ub
    ))
  }
  
  # Run in parallel (Mac/Linux)
  n_cores <- min(detectCores() - 1, n_trials)
  results_list <- mclapply(1:n_trials, single_trial, mc.cores = n_cores)
  
  # Combine results
  results <- do.call(rbind, results_list)
  return(results)
}

# Run the simulation
results <- run_trials_parallel(n_trials = 100, n = 200, variability_type = 1)
print(head(results))

# Add progress tracking (optional)
cat("Running", n_trials, "trials with", n_cores, "cores...\n")

# Add simple analysis at the end
cat("\n=== RESULTS SUMMARY ===\n")
cat("Bias - Cohen's d:", round(mean(results$cohens_d - results$true_es), 4), "\n")
cat("Bias - Plugin:", round(mean(results$es_plugin - results$true_es), 4), "\n") 
cat("Bias - One-step:", round(mean(results$es_one_step - results$true_es), 4), "\n")

# Coverage rates
cohens_coverage <- mean(results$true_es >= results$cohens_d_lb & results$true_es <= results$cohens_d_ub)
onestep_coverage <- mean(results$true_es >= results$es_one_step_lb & results$true_es <= results$es_one_step_ub)
cat("Coverage - Cohen's d:", round(cohens_coverage, 3), "\n")
cat("Coverage - One-step:", round(onestep_coverage, 3), "\n")