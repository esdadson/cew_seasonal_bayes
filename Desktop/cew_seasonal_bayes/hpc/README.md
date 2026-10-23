# Running the covariate extraction on Hazel

This extracts Daymet weather and CDL landscape metrics for every trap site (sections 2–3 of `cew_source_sink_pipeline.R`) as Slurm jobs on Hazel. The output files are then read by the pipeline.

> **Rule:** never run analysis code on a login node. Login nodes are for editing, `git`, `module`, installing packages, and `sbatch`/`squeue`/`seff`. Every step that runs the analysis is submitted with `sbatch`. As a safeguard, `cew_covariates_hpc.R` refuses to start on Hazel unless it is inside a Slurm job.

Where each step runs:
- **[Mac]**: your laptop.
- **[Hazel]**: a login node, via `ssh`.

R comes from Hazel's R module, with packages in a personal library, so no conda environment is needed. Conda remains a fallback (appendix).

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

## 2. Get the code [Hazel]

```bash
git clone https://github.com/esdadson/cew_seasonal_bayes.git
cd cew_seasonal_bayes/Desktop/cew_seasonal_bayes
```

If the repo is private, use a GitHub personal access token as the password.

## 3. Set up R (one time) [Hazel]

1. **Check your group and package-library space.** Hazel's docs put R packages in `/usr/local/usrapps/<group>/<user>/libs/R`.
   ```bash
   groups                               # should include cmjone25
   ls /usr/local/usrapps/cmjone25       # must exist; if not, ask your PI to request usrapps space
   ```

2. **Modules are set in `hpc/r_env.sh`.** You shouldn't need to edit it.
   - `module load R` (R 4.4.0, Hazel's documented default) is always loaded.
   - `gdal`, `proj` and `geos` are loaded if they exist. `terra` needs them to compile and run.
   - The installer and every job print which modules were loaded. If GDAL is missing, they also print `no gdal-config on PATH`.

   If that warning appears, check what's available:
   ```bash
   module avail gdal ; module avail proj ; module avail geos
   ```
   - **Versioned names** (e.g. `gdal/3.8.4`): put those in `R_OPTIONAL_MODULES` in `hpc/r_env.sh`.
   - **Nothing:** GDAL on Hazel comes through the `gis` software group. Email HPC support to request access, or use the conda fallback in the appendix.

3. **Install the packages.** This runs on the login node, which Hazel's docs require: compute nodes have no internet. It only installs packages and runs no analysis.
   ```bash
   ./hpc/install_r_packages.sh
   ```
   - It installs only what is missing from `terra`, `landscapemetrics`, `dplyr`, `tidyr` and `geosphere`.
   - Compiling `terra` can take 20–40 minutes.
   - It finishes with `All packages available` or names the packages that failed.

   **If `terra` fails to compile** (errors mentioning `gdal-config`, `proj` or `geos`), see step 3.2.

## 4. Choose where outputs go (optional) [Hazel]

To save results to research storage, where your Mac can read them directly:

```bash
export CEW_OUT=/rs1/researchers/c/cmjone25/Data/Original/pest-occurrence/corn_earworm/covariates
```

Otherwise they go to `covariates/` in the repo.

`export` only lasts for the current login, so set it again before each `sbatch` or `submit_covariates.sh`. Don't use `/share` for outputs: it is purged.

## 5. Dry run [submitted from Hazel, runs on a compute node]

Run from `cew_seasonal_bayes/Desktop/cew_seasonal_bayes`:

```bash
sbatch hpc/dryrun.sbatch            # "Submitted batch job <jobid>"
squeue -u $USER                     # PD = waiting, R = running; gone = finished
```

When it is finished:

```bash
tail -30 logs/dryrun.<jobid>.out
seff <jobid>
```

What the log shows:
- **`N/56 input files found`:** if any files are missing, the job stops and lists them. Fix the paths, then resubmit.
- **Sites returned:** sites outside the US get no CDL metrics. That is expected.
- **Suggested settings:** the `DAYMET_*` and `CDL_*` values to use.
- **`seff`:** *Memory Utilized* and *Wall-clock time*. Trust these if they differ from the suggestions.

If the job fails:
- **`there is no package called ...`:** re-run step 3.3.
- **`libgdal... cannot open shared object`:** the job didn't load the GDAL module the packages were built with. See step 3.2.

## 6. Set resources [Hazel]

```bash
nano hpc/submit_covariates.sh
```

Paste the suggested values into the `Resources` block.

For `CDL_MEM`, use about 1.5 × *Memory Utilized* × `CDL_CPUS`, because each core runs its own R process.

## 7. Full run [Hazel]

```bash
./hpc/submit_covariates.sh run
squeue -u $USER
```

This submits the following:
- **Site table:** 1 job.
- **Daymet:** 42 jobs, one per variable (`tmax`, `tmin`, `prcp`) per year from 2007 to 2020.
- **CDL:** 13 jobs, one per year from 2008 to 2020.

The Daymet and CDL jobs start only after the site table finishes successfully. You can log out while they run.

## 8. Check the results [Hazel]

```bash
ls covariates/daymet_*.rds | wc -l    # expect 42  (or $CEW_OUT/... if you set it)
ls covariates/cdl_*.rds    | wc -l    # expect 13
sacct -u $USER -S today --format=JobName%12,JobID%16,State,Elapsed,MaxRSS | grep -v COMPLETED
```

Any line from `sacct` is a job that didn't complete. Read its log:
- Daymet: `logs/dm_<var>.<arrayid>_<index>.out`
- CDL: `logs/cdl.<arrayid>_<index>.out`

The index counts years from 0, starting at 2007 for Daymet and 2008 for CDL.

To fix and resubmit:
- **Out of time or memory** (`TIMEOUT` or `OUT_OF_MEMORY`): raise that value in step 6.
- **Resubmit:** `./hpc/submit_covariates.sh daymet` or `./hpc/submit_covariates.sh cdl`. This re-runs the whole array.

## 9. Run the pipeline [Mac]

If you skipped step 4, copy the outputs back:

```bash
rsync -av <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/Desktop/cew_seasonal_bayes/covariates/ \
  ~/Desktop/cew_seasonal_bayes/Desktop/cew_seasonal_bayes/covariates/
```

If you used research storage, first run this in R:

```r
Sys.setenv(CEW_OUT = "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/covariates")
```

Then run `cew_source_sink_pipeline.R`.

The pipeline stops if its site table doesn't match the one the HPC jobs used. If that happens, re-run step 7.

## Quick reference

| Task | Command |
| --- | --- |
| My jobs | `squeue -u $USER` |
| Cancel one / all | `scancel <jobid>` / `scancel -u $USER` |
| Finished job's usage | `seff <jobid>` |
| Today's history | `sacct -u $USER -S today` |

If `sbatch` rejects `--partition=compute`, change it in `hpc/dryrun.sbatch` and in `PARTITION` in `hpc/submit_covariates.sh`.

## Appendix: conda fallback

Use this only if `terra` can't be built against Hazel's modules. Do these steps instead of step 3.3. All commands run on the login node.

```bash
module load conda && conda init bash        # then log out and back in
mkdir -p /share/cmjone25/$USER/conda/pkgs
printf "pkgs_dirs:\n  - /share/cmjone25/$USER/conda/pkgs\nchannels:\n  - conda-forge\n" > ~/.condarc
conda env create --prefix /share/cmjone25/$USER/cew_env -f hpc/cew_env.yml
```

Then set `R_SETUP=${R_SETUP:-conda}` in `hpc/r_env.sh`.

`/share` is purged, so if a job fails at `conda activate`, re-run the `conda env create` line.

Policies: [Hazel AUP](https://hpc.ncsu.edu/About/AUP.php), [R on Hazel](https://hpc.ncsu.edu/Software/Apps.php?app=R), [LSF→Slurm guide](https://hpc.ncsu.edu/RunningJobs/Migration.php). Help: help@ncsu.edu
