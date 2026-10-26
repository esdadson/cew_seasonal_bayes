# Running the covariate extraction and the stage-1 fit on Hazel

Steps 1–9 extract Daymet weather and CDL landscape metrics for every trap site (sections 2–3 of `cew_source_sink_pipeline.R`) as Slurm jobs on Hazel. The output files are then read by the pipeline. Step 10 fits the stage-1 seasonal model (section 5).

> **Rule:** never run analysis code on a login node. Login nodes are for editing, `git`, `module`, installing packages, and `sbatch`/`squeue`/`seff`. Every step that runs the analysis is submitted with `sbatch`. As a safeguard, `cew_covariates_hpc.R` and `cew_stage1_hpc.R` refuse to start on Hazel unless they are inside a Slurm job.

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
   ls -ld /usr/local/usrapps/*/$USER    # your own folder, if one exists
   ls -ld /usr/local/usrapps/cmjone25   # must exist; if not, ask your PI to request usrapps space
   ```
   You need a folder you can write to there. If `mkdir /usr/local/usrapps/cmjone25/$USER` gives *Permission denied*, ask your PI to create it and give you ownership (`chown`), or email help@ncsu.edu.
   If your folder is under a different group, put it in `~/.bashrc` and log in again:
   ```bash
   echo 'export CEW_APPS=/usr/local/usrapps/<group>/$USER' >> ~/.bashrc
   ```
   `hpc/r_env.sh` puts the R packages in `$CEW_APPS/libs/R` and CmdStan in `$CEW_APPS/cmdstan`.

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

### Sites per year

By default, each year is extracted only for the sites trapped that year. The site step writes `cew_site_years.csv` next to the outputs, listing every site-year with at least one count.

- **CDL:** year Y covers the sites trapped in Y.
- **Daymet:** the file for calendar year Y covers the sites trapped in Y−1, Y or Y+1. That file holds the December before the next season, and the ISO weeks that cross New Year.

To extract every site in every year instead (the old behaviour), run `export CEW_SITES=all` before submitting.

The choice changes weather anomalies. With per-year sites, each site's climatology averages only its own extracted years. With `CEW_SITES=all`, it averages 2008–2020.

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
- **Sites per extraction year:** a table of how many sites each Daymet file and CDL year covers. Timings use the busiest year.
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

## 10. Fit the stage-1 model [Hazel]

`cew_stage1_hpc.R` fits the stage-1 seasonal model (your fit_gp2.1 model, in `cew_stage1_data.R`) on 2008–2019. The pipeline then only reads the saved fit. The exact GP makes each iteration slow: on the laptop, 3,100 iterations took about 13 hours.

Run everything from `cew_seasonal_bayes/Desktop/cew_seasonal_bayes`.

1. **Inputs.** The job rebuilds the modelling data, so it needs the covariates from steps 7–8 in `$CEW_OUT` (or `covariates/`). It also reads these files under `$CEW_DATA/Original/pest-occurrence/corn_earworm/`:
   - `Lawton_etal_data.csv` (in the dryad folder)
   - `MCD43_EVI_points.csv`
   - `cew_landscape_metrics_full_timeseries.csv`

2. **Install Stan (one time, login node).**
   ```bash
   ./hpc/install_stan.sh
   ```
   - It installs `brms`, `cmdstanr`, `posterior`, `zoo` and `lubridate`, then CmdStan into `$CEW_APPS/cmdstan` (see step 3.1).
   - `brms` compiles `rstan`, and CmdStan compiles itself, so expect an hour or more.
   - It finishes with `All available` and the CmdStan version.
   - If `hpc/install_r_packages.sh` hasn't been run yet, run it first: `terra` is needed too.

3. **Dry run.** It compiles the model and runs one short chain (150 warmup + 50 sampling iterations) on the full training data. Give it the settings you plan to use for the full run:
   ```bash
   STAGE1_CHAINS=4 STAGE1_ITER=20000 STAGE1_WARMUP=19000 STAGE1_THIN=4 sbatch hpc/stage1_dryrun.sbatch
   ```
   When it's finished:
   ```bash
   tail -25 logs/stage1_dryrun.<jobid>.out
   seff <jobid>
   ```
   What the log shows:
   - **Seconds per gradient and leapfrog steps per iteration.** Time per iteration is the number of steps × seconds per gradient. An iteration takes at most 1,023 steps (tree depth 10).
   - **Share of iterations at max tree depth.** If it's near 100%, every iteration costs the maximum, and the worst-case estimate is the realistic one.
   - **Hours per chain**, at the dry run's tree depth and in the worst case. Chains run in parallel, so this is the job's wall time.
   - **Suggested `STAGE1_TIME`, `STAGE1_QOS` and `STAGE1_MEM`.** Time is the worst case × 1.15 + 1 h. Hazel's `normal` QOS allows up to 4 days and `long` up to 10. Above 10 days, the log gives the largest `STAGE1_ITER` that fits.

4. **Set resources.** Paste the suggested values into the `Resources` block of `hpc/submit_stage1.sh`:
   ```bash
   nano hpc/submit_stage1.sh
   ```
   Keep the sampler settings there the same as in the dry run. The script requests one core per chain.

5. **Full run.**
   ```bash
   MAIL_USER=<unityID>@ncsu.edu ./hpc/submit_stage1.sh    # MAIL_USER is optional: email when it ends
   tail -f logs/stage1.<jobid>.out                         # progress every 100 iterations; Ctrl-C stops watching only
   ```
   - It writes `fit1_gp21_train.rds` and prints the summary, divergences, the tree-depth share and the largest R-hat.
   - CmdStan's draws are also kept in `stan_csv/`. If R fails after sampling, read them with `cmdstanr::as_cmdstan_fit(list.files("stan_csv", full.names = TRUE))` instead of refitting.
   - Stan can't resume a chain, so a job that runs out of time loses everything. That's why the requested time is padded.

6. **Copy the fit back [Mac].**
   ```bash
   rsync -av <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/Desktop/cew_seasonal_bayes/fit1_gp21_train.rds \
     ~/Desktop/cew_seasonal_bayes/Desktop/cew_seasonal_bayes/
   ```
   The pipeline reads it in section 5 and never refits it (`file_refit = "never"`). Delete the file to refit, for example after changing the data or the model.

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
