#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
# Job wrapper for cew_covariates_hpc.R; submitted by hpc/submit_covariates.sh,
# which sets the resources. Usage: sbatch [options] hpc/covariates.sbatch <mode> [args]
# In array jobs the argument YEAR becomes YEAR0 + SLURM_ARRAY_TASK_ID
# (array indices stay small; Slurm caps them, so years can't be used directly).
source ~/.bashrc                      # required before conda activate (Hazel docs)
conda activate "$CONDA_ENV"
cd "$SLURM_SUBMIT_DIR"

args=("$@")
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  args=("${args[@]/#YEAR/$(( YEAR0 + SLURM_ARRAY_TASK_ID ))}")
fi
echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1}: Rscript cew_covariates_hpc.R ${args[*]}"
Rscript cew_covariates_hpc.R "${args[@]}"
