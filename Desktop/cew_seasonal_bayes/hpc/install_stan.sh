#!/bin/bash
# Install what the stage-1 fit needs on top of hpc/install_r_packages.sh:
# brms, cmdstanr, posterior, zoo, lubridate into $R_LIBS_USER and CmdStan into
# $CMDSTAN (see hpc/r_env.sh). Run on a LOGIN node: Hazel's docs require
# package installs there because compute nodes have no internet. It only
# installs software; it runs no analysis.
#   ./hpc/install_stan.sh
set -eo pipefail
cd "$(dirname "$0")/.."
source hpc/r_env.sh
[ "$R_SETUP" = conda ] && { echo "R_SETUP=conda: add r-brms and r-cmdstanr to the environment instead"; exit 1; }
mkdir -p "$R_LIBS_USER" "$CMDSTAN"
echo "R: $(command -v Rscript)   library: $R_LIBS_USER   CmdStan: $CMDSTAN"
Rscript -e '
lib  <- Sys.getenv("R_LIBS_USER")
pk   <- c("brms", "posterior", "zoo", "lubridate", "cmdstanr")
have <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
cat("Already available:", if (any(have)) pk[have] else "none", "\n")
# cmdstanr is not on CRAN; brms compiles rstan as a dependency (can take ~1 h)
if (any(!have)) install.packages(pk[!have], lib = lib, Ncpus = 2,
                                 repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) stop("Failed to install: ", paste(pk[!ok], collapse = ", "), " - see the compiler errors above")
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  cmdstanr::check_cmdstan_toolchain(fix = TRUE)
  cmdstanr::install_cmdstan(dir = Sys.getenv("CMDSTAN"), cores = 2)   # ~10-20 min
}
cat("All available. brms", as.character(packageVersion("brms")),
    "- CmdStan", cmdstanr::cmdstan_version(), "at", cmdstanr::cmdstan_path(), "\n")'
