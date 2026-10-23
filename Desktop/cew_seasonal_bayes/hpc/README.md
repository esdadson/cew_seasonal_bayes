# Running the covariate extraction on Hazel

This extracts Daymet weather and CDL landscape metrics for every trap site (sections 2–3 of `cew_source_sink_pipeline.R`) as Slurm jobs on Hazel. The output files are then read by the pipeline.

> **Rule:** never run R or Python on a login node. Every step below that runs code is submitted with `sbatch`. Login nodes are only for editing, `git`, conda installs, and `sbatch`/`squeue`/`seff`.

Where each step runs:
- **[Mac]**: your laptop.
- **[Hazel]**: a login node, via `ssh`.

## 0. Push the code [Mac]

```bash
cd ~/Desktop/cew_seasonal_bayes
git add Desktop/cew_seasonal_bayes
git commit -m "HPC covariate extraction"
git push
```

## 1. Log in [Hazel]

```bash
ssh <unityID>@login.hpc.ncsu.edu
```

Approve the two-factor prompt. Idle sessions are disconnected automatically.

## 2. One-time setup [Hazel]

1. Check your group name. The scripts assume `cmjone25`:
   ```bash
   groups
   ls /usr/local/usrapps/cmjone25      # must exist; if not, ask your PI to request usrapps space
   ```
2. Initialize conda, then log out and back in:
   ```bash
   module load conda
   conda init bash                     # ignore the warnings
   exit
   ```
3. Point the conda package cache at scratch, so it doesn't fill your home quota:
   ```bash
   mkdir -p /share/cmjone25/$USER/conda/pkgs
   cat > ~/.condarc <<'EOF'
   pkgs_dirs:
     - /share/cmjone25/$USER/conda/pkgs
   channels:
     - conda-forge
   EOF
   sed -i "s/\$USER/$USER/" ~/.condarc
   ```
4. Get the code:
   ```bash
   git clone https://github.com/esdadson/cew_seasonal_bayes.git
   cd cew_seasonal_bayes/Desktop/cew_seasonal_bayes
   ```
   If the repo is private, use a GitHub personal access token as the password, or add an SSH key from Hazel to GitHub.
5. Create the R environment. Installs must run on a login node, because compute nodes can't download packages. This step takes 10–30 minutes:
   ```bash
   module load conda
   conda env create --prefix /usr/local/usrapps/cmjone25/$USER/cew_env -f hpc/cew_env.yml
   ```
   If it fails, don't patch it with `conda install`. Delete the prefix directory and create it again.

## 3. Choose where outputs go [Hazel]

By default the outputs go to `covariates/` inside the repo. Alternatively, write them to research storage so your Mac can read them through the `/Volumes/cmjone25` mount, with no copying:

```bash
export CEW_OUT=/rs1/researchers/c/cmjone25/Data/Original/pest-occurrence/corn_earworm/covariates
```

`export` only lasts for the current login. Set it again before every `submit_covariates.sh` call, or add it to `~/.bashrc`.

Don't use `/share` for outputs you need to keep, because scratch is purged automatically.

## 4. Dry run [compute node, submitted from Hazel]

The dry run runs as a Slurm job on a compute node. The only thing you type on the login node is `sbatch`. Run it from `cew_seasonal_bayes/Desktop/cew_seasonal_bayes`:

```bash
sbatch hpc/dryrun.sbatch            # prints "Submitted batch job <jobid>"
squeue -u $USER                     # PD = waiting, R = running; empty = done
```

As a safeguard, `cew_covariates_hpc.R` refuses to start on Hazel unless it is inside a Slurm job. If it is ever launched on a login node by mistake, it exits immediately without doing any work.

If you set `CEW_OUT` in step 3, `sbatch` passes it to the job.

When the job is done:

```bash
tail -30 logs/dryrun.<jobid>.out
seff <jobid>
```

The log shows the following:
- **`N/56 input files found`**: if any files are missing, the job stops and lists them. Fix the paths, then submit the dry run again.
- **Site count:** the number of sites and how many returned CDL metrics. Sites outside the US return none; that is expected.
- **Suggested settings:** the `DAYMET_*` and `CDL_*` values to use.

`seff` shows the job's **Memory Utilized** and **Job Wall-clock time**. Memory measured inside R can miss some usage, so trust `seff` if the two differ.

If the job fails, read the error at the bottom of the log. The most likely problems are a wrong `CONDA_ENV` path or a missing R package in the environment.

## 5. Set resources [Hazel]

```bash
nano hpc/submit_covariates.sh
```

Copy the suggested values into the `Resources` block (`DAYMET_TIME`, `DAYMET_MEM`, `CDL_CPUS`, `CDL_TIME`, `CDL_MEM`).

For `CDL_MEM`, use about 1.5 × the dry run's **Memory Utilized** × `CDL_CPUS`, because each core runs its own R process.

## 6. Full run [Hazel]

```bash
./hpc/submit_covariates.sh run
squeue -u $USER
```

This submits the following:
- **Site table:** 1 job.
- **Daymet:** 3 arrays of 14 jobs, one per variable (`tmax`, `tmin`, `prcp`) per year from 2007 to 2020.
- **CDL:** 1 array of 13 jobs, one per year from 2008 to 2020.

The arrays wait until the site table has finished successfully.

You can log out while the jobs run. Log in again later to check on them.

## 7. Check the results [Hazel]

```bash
ls $CEW_OUT/daymet_*.rds | wc -l    # expect 42
ls $CEW_OUT/cdl_*.rds    | wc -l    # expect 13
sacct -u $USER -S today --format=JobName%12,JobID%16,State,Elapsed,MaxRSS | grep -v COMPLETED
```

If no `export` is active, `CEW_OUT` is empty here, so use `ls covariates/...` instead.

Any line from `sacct` is a job that didn't complete. Read its log:
- Daymet: `logs/dm_<var>.<arrayid>_<index>.out`
- CDL: `logs/cdl.<arrayid>_<index>.out`

The index counts years from 0, starting at 2007 for Daymet and 2008 for CDL.

To resubmit failures:
- **One array:** `./hpc/submit_covariates.sh daymet` or `./hpc/submit_covariates.sh cdl`. This re-runs the whole array.
- **Out of time or memory** (`TIMEOUT` or `OUT_OF_MEMORY`): raise that value in step 5 first.

## 8. Run the pipeline [Mac]

If you skipped step 3, copy the outputs back:

```bash
rsync -av <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/Desktop/cew_seasonal_bayes/covariates/ \
  ~/Desktop/cew_seasonal_bayes/Desktop/cew_seasonal_bayes/covariates/
```

Then run `cew_source_sink_pipeline.R` from the project. If the outputs are in research storage, first run `Sys.setenv(CEW_OUT = "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/covariates")` in R.

The pipeline stops if its site table doesn't match the one the HPC jobs used. If that happens, re-run step 6.

## Quick reference

| Task | Command |
| --- | --- |
| My jobs | `squeue -u $USER` |
| Cancel one / all | `scancel <jobid>` / `scancel -u $USER` |
| Usage of a finished job | `seff <jobid>` |
| Today's job history | `sacct -u $USER -S today` |

Policies: [Hazel AUP](https://hpc.ncsu.edu/About/AUP.php), [LSF→Slurm guide](https://hpc.ncsu.edu/RunningJobs/Migration.php). Help: help@ncsu.edu
