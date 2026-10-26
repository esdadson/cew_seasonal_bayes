# R setup for Hazel, sourced by the job scripts and the install scripts.
#   R_SETUP=conda (default):  the conda environment in $CONDA_ENV (hpc/cew_env.yml).
#                             conda-forge ships terra with GDAL/PROJ/GEOS and
#                             brms/cmdstanr/CmdStan prebuilt - nothing compiles.
#   R_SETUP=module:           Hazel's R module + packages built into $R_LIBS_USER
R_SETUP=${R_SETUP:-conda}
# Where your software lives (Hazel docs: /usr/local/usrapps/$GROUP/$USER).
# You need write access there; if your usrapps group isn't rkmeente, set
# CEW_APPS in ~/.bashrc, e.g. export CEW_APPS=/usr/local/usrapps/<group>/$USER
CEW_APPS=${CEW_APPS:-/usr/local/usrapps/rkmeente/$USER}
CONDA_ENV=${CONDA_ENV:-$CEW_APPS/cew_env}
# module setup only. `module load R` is R 4.4.0 on Hazel. The optional modules
# are loaded if they exist (terra needs GDAL, PROJ and GEOS to compile).
R_MODULES=${R_MODULES:-"R"}
R_OPTIONAL_MODULES=${R_OPTIONAL_MODULES:-"gdal proj geos"}
# The project's .Rprofile starts renv (the laptop's package library); on Hazel
# packages come from conda or $R_LIBS_USER, so keep renv off
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

# ~/.bashrc defines `module`; a failing line in it must not silently end
# a caller running under `set -e` (the install scripts), so pause errexit
echo "r_env.sh: setting up R ($R_SETUP)"
case $- in *e*) _errexit=1; set +e ;; *) _errexit= ;; esac
[ -f ~/.bashrc ] && source ~/.bashrc
if [ "$R_SETUP" = conda ]; then
  # Hazel docs: don't also `module load R`, and no personal R library
  module load conda >/dev/null 2>&1
  eval "$(conda shell.bash hook 2>/dev/null)"
  unset R_LIBS_USER R_LIBS
  if [ -d "$CONDA_ENV" ]; then
    conda activate "$CONDA_ENV" || { echo "r_env.sh: could not activate $CONDA_ENV" >&2; return 1; }
    echo "r_env.sh: conda env $CONDA_ENV; R: $(command -v Rscript); CmdStan: ${CMDSTAN:-none}"
  else
    echo "r_env.sh: no conda env at $CONDA_ENV - create it from hpc/cew_env.yml (hpc/README.md step 3)" >&2
  fi
  [ -n "$_errexit" ] && set -e
else
  [ -n "$_errexit" ] && set -e
  export R_LIBS_USER=${R_LIBS_USER:-$CEW_APPS/libs/R}
  # CmdStan from hpc/install_stan.sh; cmdstanr uses the newest cmdstan-<version> here
  export CMDSTAN=${CMDSTAN:-$CEW_APPS/cmdstan}
  module load $R_MODULES || { echo "r_env.sh: could not load: $R_MODULES" >&2; return 1; }
  R_LOADED="$R_MODULES"
  for m in $R_OPTIONAL_MODULES; do
    module load "$m" >/dev/null 2>&1 && R_LOADED="$R_LOADED $m"
  done
  echo "r_env.sh: loaded modules: $R_LOADED"
  command -v gdal-config >/dev/null || echo "r_env.sh: no gdal-config on PATH (terra can't compile without GDAL) - use conda instead: unset R_SETUP (hpc/README.md step 3)" >&2
fi
