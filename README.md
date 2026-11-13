# cew_seasonal_bayes
Bayesian modeling of Corn earworm population timing: weekly spatial models of
H. zea (CEW) trap counts and the influence of Southern source populations.

```
cew_seasonal_bayes/
├── original scripts (unchanged; they source each other from this directory):
│     cew_fourier_covariates.R, cew_landscapemetrics.R, cew_source_sink_functions.R,
│     cew_source_sink_pipeline.R, cew_stage1_data.R, cew_stage1_hpc.R,
│     CEW_JIPM_Figures.R (publication figures)
├── cew_covariates_hpc.R        covariate extraction (Daymet, CDL) -> covariates/   [Hazel]
├── covariates/                 extracted covariates and the site table (cew_sites.csv)
├── notebooks/evi_extraction.ipynb   MODIS EVI per site (Earth Engine), keyed by location2
├── step1/step1_data.R          Step 1 modelling data from covariates/ (weekly log(CEW + 1),
│                               weeks 1-52 of every monitored site-year)
├── step1_stan/                 Step 1: three nested Stan models, 2020 held out
│   ├── fixed.stan                intercept, trap type, EVI, GDD, crop shares
│   ├── fourier.stan              + shared harmonic seasonal curve (<= 2 harmonics) + year intercept
│   ├── spatiotemporal.stan       + site intercept (spatial Matern + independent), spatial site
│   │                             slopes on cos1 / sin1 (raw sites), AR(1) residual in site-year
│   ├── stan_common.R             data, design matrices, priors, scoring and plots
│   └── step1_stan.R              modes dryrun / fit / summary -> results/step1_stan/
│                                 (Hazel: hpc/step1_stan.sbatch; hpc/README.md section 11)
├── step2/                      Step 2: regions, Granger causality S -> T / N, DAGs   [laptop]
├── hpc/                        Slurm scripts and the Hazel guide (hpc/README.md)
├── archive/                    kept for reference; paths inside are not updated
│   ├── step1_v1/                 first step 1 comparison (models, HPC driver, kriging test)
│   ├── step1_v2/                 brms comparison (fixed / seasonal / spatial GP) and its Slurm scripts
│   └── step1_explore/            R-INLA prototypes, earlier Stan / brms models (space-time knots,
│                                 zone harmonic, spacetime_simple with its brms check) and their results
└── gadm/, renv/, renv.lock     boundaries and the R package environment
```

Run order:
1. Covariates on Hazel: `./hpc/submit_covariates.sh` (details in hpc/README.md)
2. Step 1 on Hazel: `STAN_MODE=dryrun sbatch ... hpc/step1_stan.sbatch`, then the fit and
   summary (hpc/README.md section 11)
3. Step 2 from the step 1 results: `step2/step2_granger.R` (it still expects the files the
   archived brms comparison wrote: model_data.rds, predictions.rds, peak_draws.rds)

All scripts run from this directory.

The archived R-INLA prototypes (archive/step1_explore/step1_inla/) need INLA, which is not on
CRAN; for R 4.4 on Windows the newest binary is 24.12.11:
```r
utils::install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable",
                        type = "binary", lib = .libPaths()[1], dependencies = FALSE)
renv::install(c("fmesher", "sn"))
```
