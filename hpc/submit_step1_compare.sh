#!/bin/bash
# Step 1 comparison as separate Slurm jobs: one per model (they run at the
# same time, each with its own time limit), then a short job that waits for
# them and builds comparison.csv, time_estimates.csv and the figures from
# their cached results (STEP1_REUSE_ONLY=1). A model that fails or runs out of
# time is left out of the comparison; resubmit just that one with
#   STEP1_MODELS=spatial bash hpc/submit_step1_compare.sh
#
#   bash hpc/submit_step1_compare.sh                      # fixed seasonal spatial
#   STEP1_ITER=4000 STEP1_WARMUP=2000 STEP1_THIN=8 bash hpc/submit_step1_compare.sh
#
# Time limits (edit below): the exact GP costs ~20 ms per gradient over 422
# sites, so at the treedepth cap (4095 steps) one iteration is ~80 s and
# 2000 iterations ~45 h; it should need fewer steps, but leave room. Check the
# partition's maximum with: sinfo -p compute -o "%P %l"
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs

declare -A TIME=( [fixed]=0-04:00:00 [seasonal]=1-00:00:00 [spatial]=3-00:00:00 [spatial_k]=1-00:00:00 )
models=${STEP1_MODELS:-"fixed seasonal spatial"}

ids=()
for m in $models; do
  t=${TIME[$m]:-1-00:00:00}
  id=$(STEP1_MODELS=$m sbatch --parsable --job-name=s1_$m --time=$t \
         --output=logs/step1_$m.%j.out --error=logs/step1_$m.%j.out hpc/step1_compare.sbatch)
  echo "$m: job $id (time limit $t)  log: logs/step1_$m.$id.out"
  ids+=("$id")
done

dep=$(IFS=:; echo "${ids[*]}")
# compare every default model that has results, plus any extra requested here
compare="fixed seasonal spatial"
for m in $models; do [[ " $compare " == *" $m "* ]] || compare="$compare $m"; done
id=$(STEP1_MODELS="$compare" STEP1_REUSE_ONLY=1 sbatch --parsable --job-name=s1_compare \
       --time=0-02:00:00 --cpus-per-task=1 --mem=16G --dependency=afterany:$dep \
       --output=logs/step1_compare.%j.out --error=logs/step1_compare.%j.out hpc/step1_compare.sbatch)
echo "compare: job $id (starts after $dep)  log: logs/step1_compare.$id.out"
