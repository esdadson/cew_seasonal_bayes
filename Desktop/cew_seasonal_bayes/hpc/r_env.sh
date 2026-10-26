# R setup for Hazel, sourced by the job scripts and by hpc/install_r_packages.sh.
#   R_SETUP=module (default): Hazel's R module + packages in a personal library
#   R_SETUP=conda:            the conda environment in $CONDA_ENV (hpc/cew_env.yml)
R_SETUP=${R_SETUP:-module}
# Required. `module load R` is R 4.4.0 on Hazel (hpc.ncsu.edu/Software/Software.php).
R_MODULES=${R_MODULES:-"R"}
# Loaded only if they exist: terra needs GDAL, PROJ and GEOS to compile and run.
# Hazel doesn't document these names; if `module avail gdal` shows a versioned
# name (e.g. gdal/3.8.4), put it here.
R_OPTIONAL_MODULES=${R_OPTIONAL_MODULES:-"gdal proj geos"}
export R_LIBS_USER=${R_LIBS_USER:-/usr/local/usrapps/cmjone25/$USER/libs/R}
# CmdStan for brms/cmdstanr (hpc/install_stan.sh). cmdstanr uses the newest
# cmdstan-<version> folder inside this directory.
export CMDSTAN=${CMDSTAN:-/usr/local/usrapps/cmjone25/$USER/cmdstan}
# The project's .Rprofile starts renv (the laptop's package library); on Hazel
# packages come from $R_LIBS_USER, so keep renv off
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
CONDA_ENV=${CONDA_ENV:-/share/cmjone25/$USER/cew_env}

source ~/.bashrc
if [ "$R_SETUP" = conda ]; then
  conda activate "$CONDA_ENV"         # Hazel docs: don't also `module load R`
else
  module load $R_MODULES || { echo "r_env.sh: could not load: $R_MODULES" >&2; return 1; }
  R_LOADED="$R_MODULES"
  for m in $R_OPTIONAL_MODULES; do
    module load "$m" >/dev/null 2>&1 && R_LOADED="$R_LOADED $m"
  done
  echo "r_env.sh: loaded modules: $R_LOADED"
  command -v gdal-config >/dev/null || echo "r_env.sh: no gdal-config on PATH (terra can't compile without GDAL)" >&2
fi
