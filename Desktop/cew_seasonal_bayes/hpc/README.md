# Running the covariate extraction and the step-1 models on Hazel

Steps 1–9 extract Daymet weather and CDL landscape metrics for every trap site (sections 2–3 of `cew_source_sink_pipeline.R`) as Slurm jobs on Hazel. The output files are then read by the pipeline. Step 10 fits the stage-1 seasonal model (section 5).

> **Rule:** never run analysis code on a login node. Login nodes are for editing, `git`, `module`, installing packages, and `sbatch`/`squeue`/`seff`. Every step that runs the analysis is submitted with `sbatch`. As a safeguard, `cew_covariates_hpc.R` and `cew_stage1_hpc.R` refuse to start on Hazel unless they are inside a Slurm job.

Where each step runs:
- **[Mac]**: your laptop.
- **[Hazel]**: a login node, via `ssh`.

R and all packages come from a conda environment (step 3). Hazel's R module is a fallback (appendix).

## 0. Push the code [Mac]

```bash
cd ~/Desktop/cew_seasonal_bayes
git add -A
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
cd cew_seasonal_bayes
```

If the repo is private, use a GitHub personal access token as the password.

## 3. Set up R (one time) [Hazel]

R and every package come from one conda environment built from `hpc/cew_env.yml`. conda-forge ships `terra` already built with GDAL, PROJ and GEOS, plus `brms`, `cmdstanr` and CmdStan with its C++ compiler. Nothing compiles against Hazel's modules, and you don't need the `gis` software group. Hazel's R page recommends conda for packages with external dependencies.

1. **Check your software folder.** It goes in `/usr/local/usrapps/<group>/<user>`.
   ```bash
   groups                               # should include rkmeente
   ls -ld /usr/local/usrapps/*/$USER    # your own folder, if one exists
   mkdir -p /usr/local/usrapps/rkmeente/$USER
   ```
   - If `mkdir` gives *Permission denied*, ask your PI to create the folder and give you ownership (`chown`), or email help@ncsu.edu.
   - If your folder is under a different group, put it in `~/.bashrc` and log in again:
     ```bash
     echo 'export CEW_APPS=/usr/local/usrapps/<group>/$USER' >> ~/.bashrc
     ```

2. **Remove any `~/.Renviron`.** Hazel's docs say conda R must not see a personal library. If an earlier `R_LIBS` line is there, delete it:
   ```bash
   cat ~/.Renviron 2>/dev/null        # if it sets R_LIBS: rm ~/.Renviron
   ```

3. **Create the environment.** This runs on the login node, because compute nodes have no internet. It only installs software and runs no analysis. Expect 20–60 minutes.
   ```bash
   module load conda
   mkdir -p /share/rkmeente/$USER/conda/pkgs            # download cache only; /share is purged, which is fine
   export CONDA_PKGS_DIRS=/share/rkmeente/$USER/conda/pkgs
   conda env create -f hpc/cew_env.yml -p ${CEW_APPS:-/usr/local/usrapps/rkmeente/$USER}/cew_env
   ```
   `CONDA_PKGS_DIRS` keeps the downloads (several GB) out of your home directory, which is too small for them. If you see `Disk quota exceeded` under `/home/$USER/.conda/pkgs`, run `rm -rf ~/.conda/pkgs` and the folder `-p` points to, then run these lines again.
   To add packages later, edit `hpc/cew_env.yml` and run `conda env update -f hpc/cew_env.yml -p <same path>`.

4. **Check it.**
   ```bash
   source hpc/r_env.sh
   Rscript -e 'library(terra); library(brms); cat("CmdStan", cmdstanr::cmdstan_version(), "\n")'
   ```
   `r_env.sh` prints the environment, R and CmdStan paths, and the R command prints the CmdStan version. Running these two `library()` calls on the login node is allowed: they're a check, not an analysis.

   Every job script sources `hpc/r_env.sh`, so the jobs use this environment automatically.

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

Run from `cew_seasonal_bayes`:

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
- **`there is no package called ...`** or **`no conda env at ...`:** the environment from step 3 is missing or incomplete. Re-run step 3.3, or step 3.4 to check it.

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
rsync -av <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/covariates/ \
  ~/Desktop/cew_seasonal_bayes/covariates/
```

If you used research storage, first run this in R:

```r
Sys.setenv(CEW_OUT = "/Volumes/cmjone25/Data/Original/pest-occurrence/corn_earworm/covariates")
```

Then run `cew_source_sink_pipeline.R`.

The pipeline stops if its site table doesn't match the one the HPC jobs used. If that happens, re-run step 7.

## 10. Fit the stage-1 model [Hazel]

`cew_stage1_hpc.R` fits the stage-1 seasonal model (your fit_gp2.1 model, in `cew_stage1_data.R`) on 2008–2019. The pipeline then only reads the saved fit. The exact GP makes each iteration slow: on the laptop, 3,100 iterations took about 13 hours.

Run everything from `cew_seasonal_bayes`.

1. **Inputs.** The job rebuilds the modelling data, so it needs the covariates from steps 7–8 in `$CEW_OUT` (or `covariates/`). It also reads these files under `$CEW_DATA/Original/pest-occurrence/corn_earworm/`:
   - `Lawton_etal_data.csv` (in the dryad folder)
   - `MCD43_EVI_points.csv`
   - `cew_landscape_metrics_full_timeseries.csv`

2. **Software.** The conda environment from step 3 already includes `brms`, `cmdstanr` and CmdStan, so there is nothing more to install. Step 3.4 checks it. (`hpc/install_stan.sh` is only for the module setup in the appendix.)

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
   rsync -av <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/fit1_gp21_train.rds \
     ~/Desktop/cew_seasonal_bayes/
   ```
   The pipeline reads it in section 5 and never refits it (`file_refit = "never"`). Delete the file to refit, for example after changing the data or the model.

## 11. Step 1: fit and compare the weekly models [Hazel]

Step 1 lives in `step1/` and step 2 in `step2/`. The original scripts at the top level are unchanged.

| File | What it does |
| --- | --- |
| `step1/step1_data.R` | Builds the modelling data from `covariates/`: weekly counts plus a prediction grid of weeks 1–52 for every monitored site-year. |
| `step1/step1_models.R` | Defines the three candidate models, their priors and sampler settings, and the functions for the 52-week predictions and peak weeks. |
| `step1/step1_hpc.R` | Slurm driver with three modes: `dryrun`, `fit` and `compare`. |
| `step2/step2_granger.R`, `step2/step2_functions.R` | Spatial clustering and Granger causality. Run on your laptop from `results/step1/`. |

The three models are nested. Each one adds a term to the one before. `compare` ranks them with `loo_compare`.
- `fixed`: zone, year and trap type as fixed effects, plus weekly EVI, GDD and crop shares.
- `spatial`: `fixed` plus a spatial error term. This is an approximate GP with an exponential kernel on km coordinates.
- `fourier`: `spatial` plus Fourier seasonal curves. Each zone has one curve shared over space, and each site gets its own offset from it.

Sampling defaults: 4 chains × 20 000 iterations, with 10 000 warmup and thin = 40. That gives 1000 posterior draws. The sampler uses `adapt_delta = 0.95` and `max_treedepth = 12`. Change these with `STEP1_ITER`, `STEP1_WARMUP` and `STEP1_THIN` in `hpc/submit_step1.sh`. Thinning makes the saved draws smaller, but every iteration still has to run, so it doesn't shorten the job.

Run everything from `cew_seasonal_bayes`. The inputs are the same as in step 10.1: `covariates/` and the data root.

1. **Dry run each model.** Each dry run compiles the model and runs one short chain (150 warmup + 50 sampling iterations) with the full run's threads per chain:
   ```bash
   STEP1_MODEL=fourier STEP1_THREADS=4 sbatch hpc/step1_dryrun.sbatch
   STEP1_MODEL=spatial STEP1_THREADS=4 sbatch hpc/step1_dryrun.sbatch
   STEP1_MODEL=fixed   STEP1_THREADS=4 sbatch hpc/step1_dryrun.sbatch
   ```
   `--cpus-per-task` in `hpc/step1_dryrun.sbatch` must be at least `STEP1_THREADS`. When a dry run finishes:
   ```bash
   tail -25 logs/step1_dryrun.<jobid>.out
   seff <jobid>
   ```
   The log shows:
   - the data step's time
   - seconds per gradient and leapfrog steps per iteration
   - the share of iterations at max tree depth
   - hours per chain
   - suggested `STEP1_TIME`, `STEP1_QOS` and `STEP1_MEM`

2. **Set resources.** Paste the suggested values into the `Resources` block of `hpc/submit_step1.sh`. `STEP1_TIME` applies to every model, and `TIME_<model>` overrides it for one model. For example, use `TIME_fourier=8-00:00:00` with `STEP1_QOS=long`. Each fit job requests `STEP1_CHAINS × STEP1_THREADS` cores.

3. **Fit, compare, predict.**
   ```bash
   chmod +x hpc/submit_step1.sh                            # once, if git didn't keep the execute bit
   MAIL_USER=<unityID>@ncsu.edu ./hpc/submit_step1.sh all  # or: fit | compare
   squeue -u $USER
   ```
   `all` submits one fit job per model. It also submits a `compare` job, which starts only after every fit has succeeded.
   - **`fit`** writes `results/step1/fit_<model>.rds` (the brms fit) and `loo_<model>.rds` (its PSIS-LOO). It prints the summary, divergences, the share of iterations at max tree depth, the largest R-hat, and the smallest bulk and tail ESS. CmdStan's draws are kept in `results/step1/stan_csv/<model>/`. Log: `logs/step1_fit_<model>.<jobid>.out`.
   - **`compare`** prints Pareto-k tables and `loo_compare`, and writes `loo_compare.csv`. It selects the model with the best ELPD. If the difference is within about 2 `se_diff`, the simpler model is just as good: set `STEP1_SELECT=<model> ./hpc/submit_step1.sh compare` to use it instead. It then predicts the selected model over weeks 1–52 of every monitored site-year and writes:
     - `predictions.rds`: posterior mean, sd and 95% interval of E[log(CEW + 1)] for every row, plus `STEP1_PRED_DRAWS` (default 100) draws for step 2
     - `peak_draws.rds`: the peak week of every site-year in all 1000 draws
     - `peak_site_year.csv`, `peak_zone_year.csv`, `peak_zone_delay.csv` and `peak_site.csv`
     - `pred_zone_week.csv`: zone mean predicted curve per year
     - `model_data.rds`: the modelling data, for step 2

     Log: `logs/step1_compare.<jobid>.out`.
   - Stan can't resume a chain, so a job that runs out of time loses everything. That's why the requested time is padded.

4. **Copy the results back [Mac].** Step 2 and the figures need everything in `results/step1/` except the fits and `stan_csv/`:
   ```bash
   rsync -av --exclude 'fit_*.rds' --exclude 'stan_csv' \
     <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/results/step1 \
     ~/Desktop/cew_seasonal_bayes/results/
   ```

5. **Step 2 [Mac].** From the project directory, run `step2/step2_granger.R`. It does three things:
   - clusters the step-1 weekly populations (observed weeks plus predictions for unsampled weeks) into regions within each zone
   - tests Granger causality from Southern to Transitional and Northern regions at lags of 1–8 weeks, scored in-sample, leave-one-year-out and against a year-permutation null
   - compares the most predictive lags with the step-1 peak delays

   It writes `results/step2/`. Then run the "Step 1 and step 2 figures" section at the end of `CEW_JIPM_Figures.R`.

## Quick reference

| Task | Command |
| --- | --- |
| My jobs | `squeue -u $USER` |
| Cancel one / all | `scancel <jobid>` / `scancel -u $USER` |
| Finished job's usage | `seff <jobid>` |
| Today's history | `sacct -u $USER -S today` |

If `sbatch` rejects `--partition=compute`, change it in `hpc/dryrun.sbatch` and in `PARTITION` in `hpc/submit_covariates.sh`.

## Appendix: Hazel's R module instead of conda

This builds the packages from source against `module load R`, which needs GDAL, PROJ and GEOS modules (the `gis` software group). Use it only if conda isn't an option. All commands run on the login node.

```bash
export R_SETUP=module                 # or set R_SETUP=${R_SETUP:-module} in hpc/r_env.sh
./hpc/install_r_packages.sh           # terra, landscapemetrics, dplyr, tidyr, geosphere
./hpc/install_stan.sh                 # brms, cmdstanr, CmdStan (an hour or more)
```

Packages go to `$CEW_APPS/libs/R` and CmdStan to `$CEW_APPS/cmdstan`. To use the versioned GDAL/PROJ/GEOS modules that `module avail gdal` shows, put their names in `R_OPTIONAL_MODULES` in `hpc/r_env.sh`.

Policies: [Hazel AUP](https://hpc.ncsu.edu/About/AUP.php), [R on Hazel](https://hpc.ncsu.edu/Software/Apps.php?app=R), [LSF→Slurm guide](https://hpc.ncsu.edu/RunningJobs/Migration.php). Help: help@ncsu.edu
