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

Step 1 lives in `step1/step1_data.R` (data) and `step1_stan/` (models); step 2 in `step2/`. The original scripts at the top level are unchanged. The earlier brms comparison is in `archive/step1_v2/`.

| File | What it does |
| --- | --- |
| `step1_stan/fixed.stan` | Fixed-effects regression: intercept, zone intercept (Transitional and Northern against Southern; it shifts the level, not the timing), trap type, weekly EVI and GDD, six crop shares. |
| `step1_stan/fourier.stan` | `fixed` + a seasonal curve shared by all sites (harmonics of week, at most 2: `STAN_HARM`) + a year random intercept. |
| `step1_stan/spatiotemporal.stan` | `fourier` + a spatial site intercept (Matern 3/2 on the raw site locations) + spatially correlated site slopes on cos/sin of week (each site's peak shifts) + an AR(1) residual over weeks within each site-year. |
| `step1_stan/stan_common.R` | Data, design matrices, priors, held-out scoring, peak weeks, plots. |
| `step1_stan/step1_stan.R` | Driver with three modes: `dryrun`, `fit`, `summary`. |

Every model is fit to **2008–2019**. **2020 is held out** and predicted out of sample (`STEP1_HOLDOUT`). Zone enters every model as a fixed intercept only; the seasonal curve is shared, so timing differences come from the spatial site slopes and the covariates. Zones also summarise the predicted peak weeks.

Sampling defaults: 4 chains × 15 000 iterations (7 500 warmup), thin = 30, so 4 × 7 500 / 30 = **1000 posterior draws**; `adapt_delta = 0.9`, `max_treedepth = 10`. Change them with `STAN_ITER`, `STAN_WARMUP`, `STAN_THIN`, `STAN_ADAPT_DELTA`, `STAN_MAX_TREEDEPTH`. Thinning makes the saved draws smaller, but every iteration still has to run.

1. **Dry run.** One short chain per model (150 warmup + 50 sampling). It prints seconds per gradient, leapfrog steps per iteration, divergences, the hyperparameters, and the projected hours per chain for the full run, and writes `results/step1_stan/timing.csv`.
   ```bash
   git pull
   STAN_MODE=dryrun sbatch --cpus-per-task=1 --time=0-04:00:00 hpc/step1_stan.sbatch
   tail -40 logs/step1_stan.<jobid>.out
   ```
   Laptop timing (Nov 2026, full data, 422 sites, 15 000 iterations per chain, before zone was added and v dropped): fixed ~5 min, fourier ~20 min, spatiotemporal ~18 h (17 ms per gradient, ~255 leapfrog steps, no divergences; 73 h if every iteration hit max tree depth). Hazel cores may be slower, so check the dry run there before setting `--time`.

2. **Fit.** All three models in one job (fixed and fourier first, minutes; then spatiotemporal), 4 chains in parallel. Request about 1.5 × the projected time; above 4 days add `--qos=long`.
   ```bash
   sbatch --time=<1.5 x projected> hpc/step1_stan.sbatch
   STAN_MODEL=spatiotemporal sbatch --time=... hpc/step1_stan.sbatch   # one model per job instead
   squeue -u $USER
   ```
   Each model writes `results/step1_stan/<model>/draws.rds` (hyperparameters, coefficients, site effects, peak weeks, held-out means), `loo.rds` and its CmdStan CSVs in `stan_csv/`. The log prints the summary, divergences, max-tree-depth hits, E-BFMI, the largest R-hat and the smallest ESS. Stan can't resume a chain, so a job that runs out of time loses that model.

3. **Compare.**
   ```bash
   STAN_MODE=summary sbatch --cpus-per-task=1 --time=0-01:00:00 hpc/step1_stan.sbatch
   ```
   It writes to `results/step1_stan/`:
   - `comparison.csv`: PSIS-LOO (`elpd_loo`, `p_loo`, Pareto k > 0.7, `elpd_diff` / `se_diff` against the best model), held-out 2020 (`elpd_test`, R², RMSE, MAE, 90% predictive coverage) and sampler diagnostics. The spatiotemporal model's LOO terms are conditional on the previous week (AR(1)); its 2020 score uses the marginal predictive, like the other two.
   - `zone_stats.csv`: mean peak week by zone, delays behind Southern, η² (share of between-site peak variance explained by zone)
   - `<model>/hyper.csv`, `<model>/beta.csv`
   - `peak_hist.png`, `season_curve.png`, `map_site_intercept.png`, `map_site_peak.png`

4. **Copy the results back [laptop].**
   ```bash
   rsync -av --exclude 'stan_csv' \
     <unityID>@login.hpc.ncsu.edu:cew_seasonal_bayes/results/step1_stan \
     ~/Desktop/cew_seasonal_bayes/results/
   ```

5. **Step 2 [Mac].** From the project directory, run `step2/step2_granger.R`. It takes about a minute up to tens of minutes, depending on the number of regions, and does four things:
   - **Clustering:** groups sites into regions within each zone, using space plus each site's mean weekly curve over 2008–2019.
   - **Zone level:** tests Southern lags of 1–8 weeks against the Transitional and Northern zones, with leads and the reverse direction as falsification checks.
   - **Region level:** finds which Southern regions predict which Transitional and Northern regions, and at which lag, then compares those lags with distance and the step-1 peak delays.
   - **DAGs:** tests South → North directly against South → Transitional → North, using conditional Granger tests (`S -> N | T` against `T -> N | S`) at zone level and for each Northern region.

   Every regression is fit to 2008–2019. Significance comes from a year-permutation null (Benjamini–Hochberg across region pairs), and 2020 gives the out-of-sample gain.

   Step 2 still expects `model_data.rds`, `predictions.rds` and `peak_draws.rds` in `results/step1/` from the archived brms comparison; the Stan summary does not write them yet.

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
