#!/bin/bash
#SBATCH --job-name=step1_compare
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4          # = STEP1_CHAINS x STEP1_THREADS
#SBATCH --mem=24G
#SBATCH --time=0-06:00:00
#SBATCH --output=logs/step1_compare.%j.out
#SBATCH --error=logs/step1_compare.%j.out
# Step 1 model comparison (step1/step1_compare.R): fits each model in
# STEP1_MODELS to 2008-2019, one after another, scores them with PSIS-LOO and
# on the held-out year, and writes results/step1/.
#   sbatch hpc/step1_compare.sbatch                          # all models, defaults
#   STEP1_MODELS=fixed sbatch hpc/step1_compare.sbatch       # one model
#   STEP1_ITER=4000 STEP1_WARMUP=2000 STEP1_THIN=8 sbatch hpc/step1_compare.sbatch
# Resubmitting reuses finished models (results/step1/check_<model>.rds).
# Watch: squeue -u $USER   Progress: tail -f logs/step1_compare.<jobid>.out   After: seff <jobid>
cd "$SLURM_SUBMIT_DIR"
source hpc/r_env.sh                                   # conda (default) or Hazel's R module

export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
export STEP1_HOLDOUT=${STEP1_HOLDOUT:-2020}
export STEP1_MODELS=${STEP1_MODELS:-"fixed seasonal spatial"}
export STEP1_CHAINS=${STEP1_CHAINS:-4} STEP1_ITER=${STEP1_ITER:-2000}
export STEP1_WARMUP=${STEP1_WARMUP:-1000} STEP1_THIN=${STEP1_THIN:-4} STEP1_THREADS=${STEP1_THREADS:-1}

echo "$(date '+%F %T') $(hostname) cpus=${SLURM_CPUS_PER_TASK:-1} models='$STEP1_MODELS'" \
     "chains=$STEP1_CHAINS iter=$STEP1_ITER warmup=$STEP1_WARMUP thin=$STEP1_THIN"
Rscript step1/step1_compare.R
echo "$(date '+%F %T') done"
