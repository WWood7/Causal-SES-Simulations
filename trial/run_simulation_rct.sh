#!/bin/bash
# Causal Effect Size Simulation Script
# Usage: chmod u+x run_simulation_rct.sh
# ./run_simulation_rct.sh trial_rct.r es_rct

analysis=$2      # change for every analysis you run (2nd arg)
maildom='@emory.edu'   # your email domain (for receiving error messages)
myscratch="/home/wwu227/CSES_results/scratch_rct"  # location of your persistent scratch dir
resultdir="/home/wwu227/CSES_results/scratch_rct/out"  # This is a folder in permanent storage
script=$1      # your code as (R or Python) script (1st arg)
workdir=$(dirname "$script")
script_base=$(basename "$script")
max_jobs=500  # max number of jobs per loop
loops=5   # total loops so 5 * 500 = 2500 tasks (500 seeds * 5 sample sizes * 1 effect type)


username=$(id -nu)

# if scratch directory doesn't exist, make it
[ ! -d ${myscratch} ] && mkdir ${myscratch}
[ ! -d ${myscratch}/out ] && mkdir ${myscratch}/out
[ ! -d ${myscratch}/err ] && mkdir ${myscratch}/err

# submit simulation jobs and collect job IDs
job_ids=()
echo "Submitting ${loops} job arrays with ${max_jobs} tasks each..."

for i in $(seq 1 ${loops}); do
	echo "#!/bin/bash" > cses_cohensd$i.sh
	echo "#SBATCH --array=1-$max_jobs" >> cses_cohensd$i.sh
	echo "#SBATCH --partition=day-long-cpu" >> cses_cohensd$i.sh
	echo "#SBATCH --error=${myscratch}/err/${analysis}$i.err" >> cses_cohensd$i.sh
	echo "#SBATCH --output=${myscratch}/out/${analysis}$i.out" >> cses_cohensd$i.sh
	echo "#SBATCH --job-name=${analysis}_batch$i" >> cses_cohensd$i.sh

    echo "source ~/miniconda3/etc/profile.d/conda.sh" >> cses_cohensd$i.sh
	echo "conda activate r_fresh" >> cses_cohensd$i.sh
	echo "Rscript ${script} $i" >> cses_cohensd$i.sh

    # Submit job and capture job ID
    job_output=$(sbatch cses_cohensd$i.sh)
    job_id=$(echo $job_output | awk '{print $4}')
    job_ids+=($job_id)
    echo "Submitted batch $i: Job ID $job_id"
done

# Create dependency string for all jobs
dependency_string=$(IFS=:; echo "${job_ids[*]}")
echo "All simulation jobs submitted. Job IDs: ${job_ids[*]}"

# Create and submit final notification job
echo "Creating final notification job..."
echo "#!/bin/bash" > ${analysis}_complete.sh
echo "#SBATCH --partition=day-long-cpu" >> ${analysis}_complete.sh
echo "#SBATCH --dependency=afterany:$dependency_string" >> ${analysis}_complete.sh
echo "#SBATCH --mail-user=$username$maildom" >> ${analysis}_complete.sh
echo "#SBATCH --mail-type=END,FAIL" >> ${analysis}_complete.sh
echo "#SBATCH --job-name=${analysis}_notification" >> ${analysis}_complete.sh
echo "#SBATCH --time=00:10:00" >> ${analysis}_complete.sh
echo "" >> ${analysis}_complete.sh
echo "echo 'Simulation ${analysis} completed at \$(date)'" >> ${analysis}_complete.sh
echo "echo 'Total jobs submitted: ${loops} batches x ${max_jobs} tasks = $((loops * max_jobs)) tasks'" >> ${analysis}_complete.sh
echo "echo 'Check results in: /home/wwu227/trial_results/results_rct'" >> ${analysis}_complete.sh

# Submit the notification job
notification_output=$(sbatch ${analysis}_complete.sh)
notification_id=$(echo $notification_output | awk '{print $4}')
echo "Submitted notification job: $notification_id"
echo "You will receive an email when all simulation jobs complete."
