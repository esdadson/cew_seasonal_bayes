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
# Where your packages and CmdStan live (Hazel docs: /usr/local/usrapps/$GROUP/$USER).
# You need write access there; if your usrapps group isn't rkmeente, set
# CEW_APPS in ~/.bashrc, e.g. export CEW_APPS=/usr/local/usrapps/<group>/$USER
CEW_APPS=${CEW_APPS:-/usr/local/usrapps/rkmeente/$USER}
export R_LIBS_USER=${R_LIBS_USER:-$CEW_APPS/libs/R}
# CmdStan for brms/cmdstanr (hpc/install_stan.sh). cmdstanr uses the newest
# cmdstan-<version> folder inside this directory.
export CMDSTAN=${CMDSTAN:-$CEW_APPS/cmdstan}
# The project's .Rprofile starts renv (the laptop's package library); on Hazel
# packages come from $R_LIBS_USER, so keep renv off
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
CONDA_ENV=${CONDA_ENV:-/share/rkmeente/$USER/cew_env}

# ~/.bashrc defines `module`/conda; a failing line in it must not silently end
# a caller running under `set -e` (the install scripts), so pause errexit
echo "r_env.sh: setting up R ($R_SETUP)"
case $- in *e*) _errexit=1; set +e ;; *) _errexit= ;; esac
[ -f ~/.bashrc ] && source ~/.bashrc
[ -n "$_errexit" ] && set -e
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
