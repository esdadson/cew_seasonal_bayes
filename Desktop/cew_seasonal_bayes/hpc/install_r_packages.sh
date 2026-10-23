#!/bin/bash
# Install the R packages cew_covariates_hpc.R needs into $R_LIBS_USER (see hpc/r_env.sh).
# Run on a LOGIN node: Hazel's docs require package installs there because
# compute nodes have no internet. It only installs packages; it runs no analysis.
#   ./hpc/install_r_packages.sh
set -eo pipefail
cd "$(dirname "$0")/.."
source hpc/r_env.sh
[ "$R_SETUP" = conda ] && { echo "R_SETUP=conda: build the environment from hpc/cew_env.yml instead"; exit 1; }
mkdir -p "$R_LIBS_USER"
echo "R: $(command -v Rscript)   library: $R_LIBS_USER   modules: $R_MODULES"
Rscript -e '
pk   <- c("terra", "landscapemetrics", "dplyr", "tidyr", "geosphere")
have <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
cat("Already available:", if (any(have)) pk[have] else "none", "\n")
if (any(!have)) install.packages(pk[!have], lib = Sys.getenv("R_LIBS_USER"),
                                 repos = "https://cloud.r-project.org", Ncpus = 2)
ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Failed to install: ", paste(pk[!ok], collapse = ", "),
                   " - see the compiler errors above (terra needs GDAL/PROJ/GEOS modules)")
cat("All packages available. terra", as.character(packageVersion("terra")), "\n")'
