#!/bin/bash
# Submit the step-1 peak model (cew_peak_hpc.R) to Slurm on Hazel, for each
# candidate model in PEAK_MODELS (cew_peak_model.R: site, zone_year).
# This script only calls sbatch, so it is fine on a login node; all R code
# runs on compute nodes through hpc/peak.sbatch.
# Dry run first (its own job script):  sbatch hpc/peak_dryrun.sbatch
#   ./hpc/submit_peak.sh fit     # all years -> fit_peak_<model>.rds, loo_peak_<model>.rds,
#                                #   peak_draws_<model>.rds
#   ./hpc/submit_peak.sh loyo    # 13 leave-one-year-out fits per model -> loyo/<model>/
#   ./hpc/submit_peak.sh all     # both (they run side by side)
#   PEAK_MODELS=site ./hpc/submit_peak.sh fit     # one candidate only
# Watch: squeue -u $USER     Progress: tail logs/peak_*.out     Usage after: seff <jobid>
set -euo pipefail
cd "$(dirname "$0")/.."                     # project dir (where the R scripts live)
what=${1:-fit}

export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
PARTITION=${PARTITION:-compute}
MAIL_USER=${MAIL_USER:-}                    # e.g. MAIL_USER=<unityID>@ncsu.edu: email at end/failure
PEAK_MODELS=${PEAK_MODELS:-"site zone_year"}

# --- Sampler settings (the dry run must use the same ones) --------------------
export PEAK_CHAINS=${PEAK_CHAINS:-4}        # chains run in parallel
export PEAK_ITER=${PEAK_ITER:-2000}         # per chain, warmup included
export PEAK_WARMUP=${PEAK_WARMUP:-1000}
export PEAK_THREADS=${PEAK_THREADS:-4}      # threads per chain
CPUS=$((PEAK_CHAINS * PEAK_THREADS))

# --- Resources: replace with the values printed by the dry run / seff ----------
PEAK_TIME=1-00:00:00   PEAK_QOS=normal   PEAK_MEM=32G   # qos normal: up to 4 days
# ------------------------------------------------------------------------------

YEARS=13                                    # 2008-2020: array tasks 0-12
mkdir -p logs
sub() {   # sub <mode> <job name> <log> [extra sbatch args]; PEAK_MODEL is exported to the job
  sbatch --parsable --partition="$PARTITION" --qos="$PEAK_QOS" --job-name="$2" \
    --cpus-per-task="$CPUS" --time="$PEAK_TIME" --mem="$PEAK_MEM" --output="$3" \
    ${MAIL_USER:+--mail-type=END,FAIL --mail-user="$MAIL_USER"} "${@:4}" hpc/peak.sbatch "$1"
}
settings="$PEAK_CHAINS chains x $PEAK_ITER iterations, $PEAK_THREADS threads/chain, $PEAK_TIME, $PEAK_MEM, qos=$PEAK_QOS"
for m in $PEAK_MODELS; do
  export PEAK_MODEL=$m
  if [ "$what" = fit ] || [ "$what" = all ]; then
    echo "peak fit ($m): $(sub fit "pk_$m" "logs/peak_fit_$m.%j.out") ($settings)"
  fi
  if [ "$what" = loyo ] || [ "$what" = all ]; then
    echo "peak loyo ($m): $(sub loyo "pkloyo_$m" "logs/peak_loyo_$m.%A_%a.out" --array=0-$((YEARS - 1))) ($YEARS tasks; $settings each)"
  fi
done
