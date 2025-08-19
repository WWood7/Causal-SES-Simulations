# Source the required functions
source("data_generation.r")
source("estimation_functions.r")

# Define a grid of parameters for the simulation
params <- expand.grid(seed = 1:500,
                     n = c(100, 500, 1000, 2000, 5000),
                     effect_size_type = c(1, 2, 3))

# Import environment parameters
# Get iter from SLURM array
iter <- Sys.getenv("SLURM_ARRAY_TASK_ID")
iter <- as.numeric(iter)

# Get nloop from command line arguments
args <- commandArgs(trailingOnly = TRUE)
nloop <- as.numeric(args[1])

# Calculate task id
max_jobs <- 500  # Must match max_jobs in run_simulation.sh
task_id <- max_jobs * (nloop-1) + iter

# Extract parameters for this job
seed <- as.numeric(params[task_id, "seed"])
n <- as.numeric(params[task_id, "n"])
effect_size_type <- as.numeric(params[task_id, "effect_size_type"])

# Set seed for reproducibility
set.seed(seed)

# Generate data
data <- generate_data(n, effect_size_type)

# get the true effect size
true_es <- c(0.191, 0.527, 0.804)[effect_size_type]

# Estimate Cohen's d
cohens_d_results <- estimate_cohens_d(data)

# Estimate causal effect size
causal_es_results <- estimate_causal_es(data)

# Combine all results
results <- data.frame(
    seed = seed,
    n = n,
    effect_size_type = effect_size_type,
    true_es = true_es,
    cohens_d = cohens_d_results$cohens_d,
    cohens_d_lb = cohens_d_results$cohens_d_lb,
    cohens_d_ub = cohens_d_results$cohens_d_ub,
    es_plugin = causal_es_results$es_plugin,
    es_one_step = causal_es_results$es_one_step,
    es_one_step_lb = causal_es_results$es_one_step_lb,
    es_one_step_ub = causal_es_results$es_one_step_ub
)

# Create results directory if it doesn't exist
results_dir <- "/home/wwu227/CSES_results/results"
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
}

# Store the results
filename <- paste0(results_dir, "/result_seed", seed, "_n", n, "_vtype", effect_size_type, ".rds")
saveRDS(results, file = filename)
