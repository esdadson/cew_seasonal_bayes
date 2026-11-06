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
├── step1/                      Step 1: weekly Bayesian model of log(CEW + 1)
│   ├── step1_data.R              modelling data from covariates/
│   ├── step1_models.R            candidate models (fixed, seasonal, spatial GP), priors, sampler, metrics
│   └── step1_compare.R           fits, LOO + held-out 2020 scores -> results/step1/
├── step2/                      Step 2: regions, Granger causality S -> T / N, DAGs   [laptop]
├── hpc/                        Slurm scripts and the Hazel guide (hpc/README.md)
├── archive/step1_v1/           first step 1 comparison (models, HPC driver, kriging test);
│                               kept for reference, paths inside are not updated
└── gadm/, renv/, renv.lock     boundaries and the R package environment
```

Run order:
1. Covariates on Hazel: `./hpc/submit_covariates.sh` (details in hpc/README.md)
2. Step 1 on Hazel: `bash hpc/submit_step1_compare.sh` (hpc/README.md section 11)
3. Step 2 from the step 1 results: `step2/step2_granger.R`

All scripts run from this directory.
