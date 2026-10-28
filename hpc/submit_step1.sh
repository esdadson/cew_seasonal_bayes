#!/bin/bash
# Submit step 1 (step1/step1_hpc.R) to Slurm on Hazel.
# This script only calls sbatch, so it is fine on a login node; all R code
# runs on compute nodes through hpc/step1.sbatch.
# Dry run first (its own job script):  STEP1_MODEL=fourier sbatch hpc/step1_dryrun.sbatch
#   ./hpc/submit_step1.sh fit       # one job per model in STEP1_MODELS
#                                   #   -> results/step1/fit_<model>.rds, loo_<model>.rds
#   ./hpc/submit_step1.sh compare   # loo_compare + 52-week predictions and peaks of the
#                                   #   selected model (after the fits have finished)
#   ./hpc/submit_step1.sh all       # fit, then compare once every fit has succeeded
#   STEP1_MODELS=fixed ./hpc/submit_step1.sh fit        # one model only
# Watch: squeue -u $USER     Progress: tail logs/step1_*.out     Usage after: seff <jobid>
set -euo pipefail
cd "$(dirname "$0")/.."                     # project dir
what=${1:-fit}

export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
PARTITION=${PARTITION:-compute}
MAIL_USER=${MAIL_USER:-}                    # e.g. MAIL_USER=<unityID>@ncsu.edu: email at end/failure
STEP1_MODELS=${STEP1_MODELS:-"fixed spatial fourier"}

# --- Sampler settings (the dry run must use the same ones) --------------------
export STEP1_CHAINS=${STEP1_CHAINS:-4}      # chains run in parallel
export STEP1_ITER=${STEP1_ITER:-20000}      # per chain, warmup included
export STEP1_WARMUP=${STEP1_WARMUP:-10000}
export STEP1_THIN=${STEP1_THIN:-40}         # 4 x (20000 - 10000) / 40 = 1000 draws
export STEP1_THREADS=${STEP1_THREADS:-4}    # threads per chain
CPUS=$((STEP1_CHAINS * STEP1_THREADS))

# --- Resources: replace with the values printed by the dry run / seff ----------
# One setting for every model unless overridden per model below; the GP
# models are the slow ones.
STEP1_TIME=${STEP1_TIME:-4-00:00:00}   STEP1_QOS=${STEP1_QOS:-normal}   STEP1_MEM=${STEP1_MEM:-32G}
TIME_fixed=${TIME_fixed:-1-00:00:00}        # e.g. TIME_fourier=8-00:00:00 with STEP1_QOS=long
# compare: load one fit, loo tables, 52-week predictions (~1-3 h)
COMPARE_TIME=${COMPARE_TIME:-08:00:00}   COMPARE_MEM=${COMPARE_MEM:-48G}
# ------------------------------------------------------------------------------

mkdir -p logs results/step1
fit_jobs=()
if [ "$what" = fit ] || [ "$what" = all ]; then
  for m in $STEP1_MODELS; do
    tvar="TIME_$m"; t=${!tvar:-$STEP1_TIME}
    j=$(sbatch --parsable --partition="$PARTITION" --qos="$STEP1_QOS" --job-name="s1_$m" \
          --cpus-per-task="$CPUS" --time="$t" --mem="$STEP1_MEM" --output="logs/step1_fit_$m.%j.out" \
          --export=ALL,STEP1_MODEL="$m" \
          ${MAIL_USER:+--mail-type=END,FAIL --mail-user="$MAIL_USER"} hpc/step1.sbatch fit)
    fit_jobs+=("$j")
    echo "step 1 fit ($m): $j ($STEP1_CHAINS chains x $STEP1_ITER iterations, thin $STEP1_THIN, $STEP1_THREADS threads/chain, $t, $STEP1_MEM, qos=$STEP1_QOS)"
  done
fi
if [ "$what" = compare ] || [ "$what" = all ]; then
  dep=""
  if [ ${#fit_jobs[@]} -gt 0 ]; then dep="--dependency=afterok:$(IFS=:; echo "${fit_jobs[*]}")"; fi
  j=$(sbatch --parsable --partition="$PARTITION" --qos=normal --job-name="s1_compare" \
        --cpus-per-task=1 --time="$COMPARE_TIME" --mem="$COMPARE_MEM" --output="logs/step1_compare.%j.out" \
        $dep ${MAIL_USER:+--mail-type=END,FAIL --mail-user="$MAIL_USER"} hpc/step1.sbatch compare)
  echo "step 1 compare + predict: $j ${dep:+(starts after ${fit_jobs[*]})}"
fi
