#!/bin/bash
#SBATCH --job-name=cew_dryrun
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/dryrun.%j.out
# Dry run of the covariate extraction, on a compute node. Submit from the
# project directory (cew_seasonal_bayes):
#   sbatch hpc/dryrun.sbatch
# Limits are generous so the measurement can't be cut short. When the job
# ends: tail -30 logs/dryrun.<jobid>.out (suggested settings) and seff <jobid>.
export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$SLURM_SUBMIT_DIR/covariates}

cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # Hazel R module (default) or conda
echo "$(date '+%F %T') $(hostname): dry run"
Rscript cew_covariates_hpc.R dryrun 20
