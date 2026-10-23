# R setup for Hazel, sourced by the job scripts and by hpc/install_r_packages.sh.
#   R_SETUP=module (default): Hazel's R module + packages in a personal library
#   R_SETUP=conda:            the conda environment in $CONDA_ENV (hpc/cew_env.yml)
# Edit the three settings below to match what `module avail` shows on Hazel.
R_SETUP=${R_SETUP:-module}
R_MODULES=${R_MODULES:-"R"}           # e.g. "R gdal proj geos" if terra must be compiled
export R_LIBS_USER=${R_LIBS_USER:-/usr/local/usrapps/cmjone25/$USER/libs/R}
CONDA_ENV=${CONDA_ENV:-/share/cmjone25/$USER/cew_env}

source ~/.bashrc
if [ "$R_SETUP" = conda ]; then
  conda activate "$CONDA_ENV"         # Hazel docs: don't also `module load R`
else
  module load $R_MODULES
fi
