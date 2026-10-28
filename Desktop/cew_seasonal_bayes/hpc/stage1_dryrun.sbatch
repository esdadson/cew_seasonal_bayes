#!/bin/bash
#SBATCH --job-name=cew_s1_dry
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/stage1_dryrun.%j.out
# Dry run of the stage-1 fit, on a compute node: compiles the model and runs
# one short chain on the full training data. Submit from the project directory
# (cew_seasonal_bayes), with the full-run settings
# you intend to use so the suggestion matches (defaults: fit_stage1()):
#   STAGE1_CHAINS=4 STAGE1_ITER=20000 STAGE1_WARMUP=19000 sbatch hpc/stage1_dryrun.sbatch
# Limits are generous so the measurement can't be cut short. When the job
# ends: tail -25 logs/stage1_dryrun.<jobid>.out (suggested settings) and seff <jobid>.
export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$SLURM_SUBMIT_DIR/covariates}

cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # Hazel R module (default) or conda
echo "$(date '+%F %T') $(hostname): stage-1 dry run"
Rscript cew_stage1_hpc.R dryrun
