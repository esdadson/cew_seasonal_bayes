#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Job wrapper for the full stage-1 fit; submitted by hpc/submit_stage1.sh,
# which sets the resources and the STAGE1_* settings.
cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # Hazel R module (default) or conda
echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1}: Rscript cew_stage1_hpc.R fit"
Rscript cew_stage1_hpc.R fit
