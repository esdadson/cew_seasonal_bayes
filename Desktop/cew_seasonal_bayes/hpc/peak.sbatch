#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Job wrapper for the step-1 peak model; submitted by hpc/submit_peak.sh,
# which sets the resources and the PEAK_* settings. $1 = fit | loyo
cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # conda (default) or Hazel's R module
echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1} task=${SLURM_ARRAY_TASK_ID:-} model=${PEAK_MODEL:-site}: Rscript cew_peak_hpc.R $1"
Rscript cew_peak_hpc.R "$1"
