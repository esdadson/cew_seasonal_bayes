#!/bin/bash
#SBATCH --job-name=cew_s1_dry
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/step1_dryrun.%j.out
# Dry run of one step-1 model, on a compute node: compiles it and runs one
# short chain. Submit from the project directory
# (cew_seasonal_bayes) with the full-run settings
# you intend to use; --cpus-per-task must be >= STEP1_THREADS:
#   STEP1_MODEL=fourier STEP1_THREADS=4 sbatch hpc/step1_dryrun.sbatch
# (and again with STEP1_MODEL=spatial, fixed, simple). When it ends:
#   tail -25 logs/step1_dryrun.<jobid>.out (suggested settings) and seff <jobid>.
export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$SLURM_SUBMIT_DIR/covariates}
export STEP1_MODEL=${STEP1_MODEL:-fourier}
export STEP1_THREADS=${STEP1_THREADS:-4}

cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # conda (default) or Hazel's R module
echo "$(date '+%F %T') $(hostname): step-1 dry run of '$STEP1_MODEL', $STEP1_THREADS thread(s)"
Rscript step1/step1_hpc.R dryrun
