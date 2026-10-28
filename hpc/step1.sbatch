#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Job wrapper for step 1; submitted by hpc/submit_step1.sh, which sets the
# resources and the STEP1_* settings. $1 = fit | compare
cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # conda (default) or Hazel's R module
echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1} model=${STEP1_MODEL:-}: Rscript step1/step1_hpc.R $1"
Rscript step1/step1_hpc.R "$1"
