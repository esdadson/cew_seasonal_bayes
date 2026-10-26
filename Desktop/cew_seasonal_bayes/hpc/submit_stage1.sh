#!/bin/bash
# Submit the full stage-1 fit (cew_stage1_hpc.R fit) to Slurm on Hazel.
# This script only calls sbatch, so it is fine on a login node; all R code
# runs on a compute node through hpc/stage1.sbatch.
# Dry run first (its own job script):  sbatch hpc/stage1_dryrun.sbatch
#   ./hpc/submit_stage1.sh
# Watch: squeue -u $USER     Progress: tail logs/stage1.<jobid>.out     Usage after: seff <jobid>
set -euo pipefail
cd "$(dirname "$0")/.."                     # project dir (where the R scripts live)

export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
PARTITION=${PARTITION:-compute}
MAIL_USER=${MAIL_USER:-}                    # e.g. MAIL_USER=<unityID>@ncsu.edu: email at end/failure

# --- Sampler settings (the dry run must use the same ones) --------------------
export STAGE1_CHAINS=${STAGE1_CHAINS:-4}    # one core per chain, run in parallel
export STAGE1_ITER=${STAGE1_ITER:-20000}    # per chain, warmup included
export STAGE1_WARMUP=${STAGE1_WARMUP:-19000}
export STAGE1_THIN=${STAGE1_THIN:-4}

# --- Resources: replace with the values printed by the dry run / seff ----------
STAGE1_TIME=4-00:00:00; STAGE1_MEM=16G; STAGE1_QOS=normal   # qos long: up to 10 days
# ------------------------------------------------------------------------------

mkdir -p logs
jid=$(sbatch --parsable --partition="$PARTITION" --qos="$STAGE1_QOS" --job-name=cew_stage1 \
        --cpus-per-task="$STAGE1_CHAINS" --time="$STAGE1_TIME" --mem="$STAGE1_MEM" \
        --output="logs/stage1.%j.out" ${MAIL_USER:+--mail-type=END,FAIL --mail-user="$MAIL_USER"} \
        hpc/stage1.sbatch)
echo "stage 1: $jid ($STAGE1_CHAINS chains x $STAGE1_ITER iterations, $STAGE1_TIME, $STAGE1_MEM, qos=$STAGE1_QOS)"
