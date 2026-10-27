#!/bin/bash
#SBATCH --job-name=cew_pk_dry
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/peak_dryrun.%j.out
# Dry run of the step-1 peak model, on a compute node: compiles the model and
# runs one short chain on all years. Submit from the project directory
# (cew_seasonal_bayes/Desktop/cew_seasonal_bayes) with the full-run settings
# you intend to use; --cpus-per-task must be >= PEAK_THREADS:
#   PEAK_MODEL=zone_year PEAK_CHAINS=4 PEAK_ITER=4000 PEAK_WARMUP=2000 PEAK_THREADS=4 sbatch hpc/peak_dryrun.sbatch
# When it ends: tail -25 logs/peak_dryrun.<jobid>.out (suggested settings) and seff <jobid>.
export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$SLURM_SUBMIT_DIR/covariates}
export PEAK_THREADS=${PEAK_THREADS:-4}

cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                   # conda (default) or Hazel's R module
echo "$(date '+%F %T') $(hostname): step-1 peak model dry run, $PEAK_THREADS thread(s)"
Rscript cew_peak_hpc.R dryrun
