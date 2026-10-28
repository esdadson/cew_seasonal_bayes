# cew_seasonal_bayes
Bayesian modeling of Corn earworm population timing: weekly spatial models of
H. zea (CEW) trap counts and the influence of Southern source populations.

```
cew_seasonal_bayes/
├── cew_covariates_hpc.R        covariate extraction (Daymet, CDL) -> covariates/   [Hazel]
├── step1/                      Step 1: weekly spatial Bayesian model of log(CEW + 1)
│   ├── step1_data.R              modelling data from covariates/ + weeks 1-52 prediction grid
│   ├── step1_models.R            fixed | spatial (GP) | fourier models, priors, predictions, peaks
│   └── step1_hpc.R               Slurm driver: dryrun | fit | compare            [Hazel]
├── step2/                      Step 2: spatial clustering + Granger causality, S -> T / N
│   ├── step2_functions.R
│   └── step2_granger.R                                                          [laptop]
├── results/step1/, results/step2/   outputs of steps 1 and 2 (read by the figures)
├── CEW_JIPM_Figures.R          publication figures; step 1-2 figures at the end
├── hpc/                        Slurm scripts and the Hazel guide (hpc/README.md)
└── original scripts (unchanged): cew_fourier_covariates.R, cew_landscapemetrics.R,
    cew_source_sink_functions.R, cew_source_sink_pipeline.R, cew_stage1_data.R, cew_stage1_hpc.R
```

Run order (details in hpc/README.md, steps 1-9 and 11):
1. Covariates on Hazel: `./hpc/submit_covariates.sh`
2. Step 1 on Hazel: dry runs (`sbatch hpc/step1_dryrun.sbatch`), then `./hpc/submit_step1.sh all`
   (fits the three models, compares them with `loo_compare`, predicts the selected one over
   weeks 1-52 of every monitored site-year and extracts peak weeks)
3. Copy `results/step1/` back; run `step2/step2_granger.R`
4. Figures: the "Step 1 and step 2 figures" section at the end of `CEW_JIPM_Figures.R`

All scripts run from this directory.
