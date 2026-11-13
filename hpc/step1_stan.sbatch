#!/bin/bash
#SBATCH --job-name=step1_stan
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # = STAN_CHAINS (one core per chain)
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/step1_stan.%j.out
#SBATCH --error=logs/step1_stan.%j.out
# Step 1 Stan models: fixed, fourier, spatiotemporal (step1_stan/step1_stan.R).
# Submit from the project directory. Dry run first (one short chain per model:
# time per gradient and the projected time of 4 chains x 15,000 iterations),
# then the fits, then the comparison:
#   STAN_MODE=dryrun sbatch --cpus-per-task=1 --time=0-04:00:00 hpc/step1_stan.sbatch
#   sbatch --time=<1.5 x projected> hpc/step1_stan.sbatch                     # STAN_MODE=fit, all models
#   STAN_MODE=summary sbatch --cpus-per-task=1 --time=0-01:00:00 hpc/step1_stan.sbatch
# One model per job: STAN_MODEL=spatiotemporal sbatch ... (fixed and fourier take minutes).
# Settings pass through, e.g. STAN_ITER=15000 STAN_WARMUP=7500 STAN_THIN=30 (the defaults).
# Laptop (Nov 2026): fixed ~5 min, fourier ~20 min, spatiotemporal ~18 h per chain at 15,000 iterations.
# Longer than 4 days: add --qos=long.
# Watch: squeue -u $USER   Progress: tail -f logs/step1_stan.<jobid>.out   After: seff <jobid>
cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                                   # conda (default) or Hazel's R module

export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
export STAN_CHAINS=${STAN_CHAINS:-${SLURM_CPUS_PER_TASK:-4}}
MODE=${STAN_MODE:-fit}

echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1}: step1_stan.R $MODE (${STAN_MODEL:-all})"
Rscript step1_stan/step1_stan.R "$MODE"
