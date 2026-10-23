#!/bin/bash
# Submit the covariate extraction (cew_covariates_hpc.R) to Slurm on Hazel.
# This script only calls sbatch, so it is fine on a login node; all R code
# runs on compute nodes through hpc/covariates.sbatch.
#   ./hpc/submit_covariates.sh dry    # one job: checks inputs, times one Daymet unit + 20 CDL sites
#   ./hpc/submit_covariates.sh run    # sites -> Daymet arrays (3 x 14 jobs) + CDL array (13 jobs)
#   ./hpc/submit_covariates.sh daymet|cdl   # re-submit one array (e.g. after failures)
# Watch: squeue -u $USER     Usage after a job: seff <jobid>
set -euo pipefail
cd "$(dirname "$0")/.."                     # project dir (where the R scripts live)

# --- Paths ---------------------------------------------------------------------
# Inputs under $CEW_DATA: Original/pest-occurrence/corn_earworm/, Original/Daymet/,
# Raster/USA/landcover/
export CEW_DATA=${CEW_DATA:-/rs1/researchers/c/cmjone25/Data}
export CEW_OUT=${CEW_OUT:-$PWD/covariates}
export CONDA_ENV=${CONDA_ENV:-/usr/local/usrapps/cmjone25/$USER/cew_env}
PARTITION=${PARTITION:-compute}

# --- Resources: replace with the values printed by the dry run / seff ----------
DAYMET_TIME=01:00:00; DAYMET_MEM=8G
CDL_CPUS=4;           CDL_TIME=02:00:00;   CDL_MEM=16G   # --mem is per job (whole node share)
DAYMET_YEAR0=2007; DAYMET_N=14              # 2007-2020: min(yrs) - 1 .. max(yrs)
CDL_YEAR0=2008;    CDL_N=13                 # 2008-2020
# ------------------------------------------------------------------------------

mkdir -p logs "$CEW_OUT"
sb() { sbatch --parsable --partition="$PARTITION" "$@"; }

submit_daymet() {   # $1 = optional --dependency=...
  for v in tmax tmin prcp; do
    sb --job-name="dm_$v" ${1:+"$1"} --array=0-$((DAYMET_N - 1)) --export=ALL,YEAR0=$DAYMET_YEAR0 \
       --time="$DAYMET_TIME" --mem="$DAYMET_MEM" --output="logs/dm_$v.%A_%a.out" \
       hpc/covariates.sbatch daymet YEAR "$v"
  done
}
submit_cdl() {
  sb --job-name=cdl ${1:+"$1"} --array=0-$((CDL_N - 1)) --export=ALL,YEAR0=$CDL_YEAR0 \
     --cpus-per-task="$CDL_CPUS" --time="$CDL_TIME" --mem="$CDL_MEM" --output="logs/cdl.%A_%a.out" \
     hpc/covariates.sbatch cdl YEAR
}

case "${1:-}" in
  dry)
    # Generous limits so the measurement itself can't be cut short
    jid=$(sb --job-name=cew_dryrun --time=04:00:00 --mem=32G --output="logs/dryrun.%j.out" \
             hpc/covariates.sbatch dryrun 20)
    echo "Submitted dry run $jid. When it finishes:"
    echo "  tail -30 logs/dryrun.$jid.out   # suggested settings"
    echo "  seff $jid                       # measured time and peak memory" ;;
  run)
    jid=$(sb --job-name=cew_sites --time=00:15:00 --mem=4G --output="logs/sites.%j.out" \
             hpc/covariates.sbatch sites)
    echo "sites: $jid"
    submit_daymet "--dependency=afterok:$jid"
    submit_cdl    "--dependency=afterok:$jid" ;;
  daymet) submit_daymet ;;
  cdl)    submit_cdl ;;
  *) echo "usage: $0 dry|run|daymet|cdl"; exit 1 ;;
esac
