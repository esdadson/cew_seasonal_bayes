# =============================================================================
# cew_covariates_hpc.R
# Heavy raster extraction for sections 2-3 of cew_source_sink_pipeline.R,
# split into one unit of work per call so it runs as a Slurm job array.
# On Hazel, never run this on a login node: submit it with hpc/submit_covariates.sh.
#
# Always run inside a Slurm job (the script refuses otherwise on Linux):
#   sbatch hpc/dryrun.sbatch                            # dry run
#   ./hpc/submit_covariates.sh run                      # everything else
# Modes, as called by hpc/covariates.sbatch and hpc/dryrun.sbatch:
#   Rscript cew_covariates_hpc.R sites                 # site + site-year tables (run first)
#   Rscript cew_covariates_hpc.R daymet <year> <var>   # var = tmax | tmin | prcp
#   Rscript cew_covariates_hpc.R cdl <year>            # landscape metrics, that year's sites
#   Rscript cew_covariates_hpc.R dryrun [n_cdl_sites]  # check inputs, time one unit of each
#   Rscript cew_covariates_hpc.R local [year ...]      # workstation: everything (or some
#                                                      # years) in one session, resumable
#
# Environment:
#   CEW_DATA  data root containing Original/ and Raster/
#             (default: /rs1/researchers/c/cmjone25/Data on Hazel, see default_data_root())
#   CEW_OUT   output directory                             (default covariates)
#   CEW_SITES which sites each year's extraction covers       (default active)
#             active: CDL year Y -> sites trapped in Y; Daymet file Y -> sites
#                     trapped in Y-1, Y or Y+1 (the file holds the December
#                     before season Y+1 and ISO weeks that cross New Year)
#             all:    every site in every year
#   SLURM_CPUS_PER_TASK  cores for the CDL step (set by sbatch --cpus-per-task)
#
# Why it is faster than cew_landscapemetrics.R:
#   - Daymet: every site and day per file in one extract(), not one call per DOY;
#     each (year, variable) file is its own array job (42 jobs in parallel).
#   - CDL: no crop/mask/classify of whole states. Each site reads only a
#     32 x 32 km window around itself, which is classified and measured; sites
#     are split across cores and each year is its own array job.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(terra); library(landscapemetrics)
})
source("cew_source_sink_functions.R")

data_root <- default_data_root()
out_dir   <- Sys.getenv("CEW_OUT", "covariates")
n_cores   <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
site_mode <- match.arg(Sys.getenv("CEW_SITES", "active"), c("active", "all"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_path <- file.path(data_root, "Original/pest-occurrence/corn_earworm/doi_10_5061_dryad_m0cfxpp5x__v20221206/Lawton_etal_data.csv")
daymet_path <- function(var, yr) sprintf("%s/Original/Daymet/%s/daymet_v4_daily_na_%s_%d.nc", data_root, var, var, yr)
cdl_path    <- function(yr) sprintf("%s/Raster/USA/landcover/%d_30m_cdls/%d_30m_cdls.tif", data_root, yr, yr)

yrs         <- 2008:2020                     # keep in sync with the pipeline
daymet_yrs  <- (min(yrs) - 1):max(yrs)       # previous year for winter severity
daymet_vars <- c("tmax", "tmin", "prcp")
radii       <- c(pland = 1000, pland16k = 16100)   # metres; 1.61 and 16.1 km used in the literature

# -----------------------------------------------------------------------------
# CDL reclassification. Fixes to cew_landscapemetrics.R:
#  - code 25 (other small grains) was mapped to BOTH soybean and other crops
#  - double-crop classes (225-254) and several crop codes (34, 38, 51, 52,
#    56, 66, 75, 210-250) fell through to Non_Crop
# Double crops are assigned to the H. zea host they contain (corn > cotton > soy).
# Verify codes against the CDL legend for each year before running.
# -----------------------------------------------------------------------------
cls <- list(
  `1`  = c(1, 12, 13, 225, 226, 228, 237, 241),           # Corn (incl. dbl-crop corn)
  `2`  = c(2, 232, 238, 239),                             # Cotton (incl. dbl-crop cotton)
  `5`  = c(5, 26, 240, 254),                              # Soybeans (incl. dbl-crop soy)
  `10` = 10, `11` = 11,                                   # Peanuts, Tobacco
  `44` = c(3, 4, 6, 14, 21:25, 27:36, 38, 39, 41:57, 59:61, 66:72, 74:77,
           204:224, 227, 229:231, 233:236, 242:250)        # Other crops
)
rcl <- do.call(rbind, lapply(names(cls), function(k) cbind(cls[[k]], as.integer(k))))
stopifnot(!any(duplicated(rcl[, 1])))
rcl_complete <- rbind(rcl, cbind(setdiff(1:255, rcl[, 1]), 63L))   # everything else Non_Crop
class_names  <- c(`1` = "Corn", `2` = "Cotton", `5` = "Soybeans", `10` = "Peanuts",
                  `11` = "Tobacco", `44` = "Other_Crops", `63` = "Non_Crop")

read_out_csv <- function(name) {
  f <- file.path(out_dir, name)
  if (!file.exists(f)) stop("Run `Rscript cew_covariates_hpc.R sites` first (", f, " missing)")
  read.csv(f)
}
read_sites      <- function() read_out_csv("cew_sites.csv")
read_site_years <- function() read_out_csv("cew_site_years.csv")

# -----------------------------------------------------------------------------
# Which sites were trapped in which (ISO) year: >= 1 non-missing count, the
# same year definition as the pipeline. location2 comes from the ONE site table,
# so ids are stable across years even though the set of sites changes.
# -----------------------------------------------------------------------------
make_site_years <- function(counts_raw, sites) {
  counts_raw %>%
    mutate(date = as.Date(substr(date, 1, 10)),
           year = as.integer(format(date, "%G")),
           CEW  = suppressWarnings(as.numeric(CEW_sum))) %>%
    filter(!is.na(CEW), year %in% yrs) %>%
    inner_join(sites %>% select(longitude, latitude, location2), by = c("longitude", "latitude")) %>%
    distinct(location2, year) %>%
    arrange(year, location2)
}

# Sites to extract for one CDL year or one Daymet calendar-year file
sites_for <- function(sites, site_years, yr, product = c("cdl", "daymet")) {
  if (site_mode == "all") return(sites)
  span <- if (match.arg(product) == "daymet") (yr - 1):(yr + 1) else yr
  sites[sites$location2 %in% site_years$location2[site_years$year %in% span], ]
}

# -----------------------------------------------------------------------------
# Daymet: one variable, one year, the given sites -> long data frame
# -----------------------------------------------------------------------------
extract_daymet <- function(sites, yr, var) {
  if (!nrow(sites)) return(setNames(data.frame(integer(), as.Date(character()), numeric()),
                                    c("location2", "date", var)))
  r <- rast(daymet_path(var, yr))
  pts <- project(vect(sites, geom = c("longitude", "latitude"), crs = "EPSG:4326"), crs(r))
  v <- as.matrix(terra::extract(r, pts, ID = FALSE))     # sites x days
  dates <- as.Date(time(r))
  if (all(is.na(dates))) dates <- as.Date(sprintf("%d-01-01", yr)) + 0:(nlyr(r) - 1)
  setNames(data.frame(rep(sites$location2, times = ncol(v)), rep(dates, each = nrow(v)), as.vector(v)),
           c("location2", "date", var))
}

# -----------------------------------------------------------------------------
# CDL: metrics for a set of sites in one year. Opens its own raster handle so
# it is safe to call from forked workers.
# -----------------------------------------------------------------------------
cdl_sites <- function(sites_xy, yr) {
  r  <- rast(cdl_path(yr))
  rr <- max(radii)
  bind_rows(lapply(seq_len(nrow(sites_xy)), function(i) {
    x <- sites_xy$x[i]; y <- sites_xy$y[i]; id <- sites_xy$location2[i]
    tryCatch({
      win <- classify(crop(r, ext(x - rr, x + rr, y - rr, y + rr)), rcl_complete)
      pt  <- cbind(x, y)
      # sample_lsm's circles carry no CRS, so terra warns on every crop; harmless
      quiet <- function(expr) withCallingHandlers(expr, warning = function(w)
        if (grepl("CRS do not match", conditionMessage(w))) invokeRestart("muffleWarning"))
      m_small <- quiet(sample_lsm(win, y = pt, plot_id = id, size = radii[["pland"]], shape = "circle",
                            what = c("lsm_c_pland", "lsm_c_ed"), verbose = FALSE, progress = FALSE))
      m_large <- quiet(sample_lsm(win, y = pt, plot_id = id, size = radii[["pland16k"]], shape = "circle",
                            what = "lsm_c_pland", verbose = FALSE, progress = FALSE)) %>%
        mutate(metric = "pland16k")
      bind_rows(m_small, m_large) %>% transmute(location2 = plot_id, class, metric, value)
    }, error = function(e) {
      message("  site ", id, " (", yr, "): ", conditionMessage(e))   # e.g. outside CDL extent
      NULL
    })
  }))
}

run_cdl <- function(sites, yr, cores = n_cores) {
  empty <- tibble(location2 = integer(), year = integer())
  if (!nrow(sites)) return(empty)
  r0 <- rast(cdl_path(yr))
  xy <- crds(project(vect(sites, geom = c("longitude", "latitude"), crs = "EPSG:4326"), crs(r0)))
  sites_xy <- data.frame(location2 = sites$location2, x = xy[, 1], y = xy[, 2])
  chunks <- split(sites_xy, (seq_len(nrow(sites_xy)) - 1) %% cores)   # interleaved: spreads slow regions
  long <- bind_rows(parallel::mclapply(chunks, cdl_sites, yr = yr, mc.cores = 1))
  if (!nrow(long)) return(empty)                      # e.g. every site outside the CDL
  long %>%
    mutate(class = class_names[as.character(class)]) %>%
    filter(!is.na(class)) %>%
    pivot_wider(id_cols = location2, names_from = c(class, metric), values_from = value,
                values_fill = 0) %>%
    mutate(year = yr, .after = location2)
}

peak_mem_gb <- function() {          # Linux only: peak resident memory of this process
  s <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
  kb <- as.numeric(gsub("\\D", "", grep("^VmHWM", s, value = TRUE)))
  if (length(kb)) kb / 1024^2 else NA_real_
}

# =============================================================================
# Hazel AUP: no code on login nodes. On Linux, refuse to run outside a Slurm job.
if (Sys.info()[["sysname"]] == "Linux" && !nzchar(Sys.getenv("SLURM_JOB_ID")))
  stop("Not inside a Slurm job. Submit with `sbatch hpc/dryrun.sbatch` or ",
       "`./hpc/submit_covariates.sh run` - never run R on a login node.")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "dryrun"

write_sites <- function() {
  counts_raw <- read.csv(counts_path) %>% filter(longitude >= -105)
  sites      <- make_site_table(counts_raw)
  site_years <- make_site_years(counts_raw, sites)
  write.csv(sites, file.path(out_dir, "cew_sites.csv"), row.names = FALSE)
  write.csv(site_years, file.path(out_dir, "cew_site_years.csv"), row.names = FALSE)
  message(nrow(sites), " sites; sites trapped per year (CEW_SITES=", site_mode, "):")
  print(table(site_years$year))
  list(sites = sites, site_years = site_years)
}

if (mode == "sites") {
  write_sites()

} else if (mode == "local") {
  # Same outputs as the HPC jobs, made one after another on this machine.
  # Finished files are skipped, so an interrupted run resumes where it stopped.
  sel   <- if (length(args) > 1) as.integer(args[-1]) else c(daymet_yrs, yrs)
  st    <- write_sites()
  cores <- max(1, parallel::detectCores() - 2)
  t_all <- Sys.time()
  for (yr in intersect(daymet_yrs, sel)) for (var in daymet_vars) {
    f <- file.path(out_dir, sprintf("daymet_%s_%d.rds", var, yr))
    if (file.exists(f)) next
    t0 <- Sys.time()
    s_yr <- sites_for(st$sites, st$site_years, yr, "daymet")
    saveRDS(extract_daymet(s_yr, yr, var), f)
    message(sprintf("daymet %s %d: %d sites, %.1f min", var, yr, nrow(s_yr),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  for (yr in intersect(yrs, sel)) {
    f <- file.path(out_dir, sprintf("cdl_%d.rds", yr))
    if (file.exists(f)) next
    t0 <- Sys.time()
    s_yr <- sites_for(st$sites, st$site_years, yr, "cdl")
    out  <- run_cdl(s_yr, yr, cores)
    saveRDS(out, f)
    message(sprintf("cdl %d: %d/%d sites returned metrics on %d cores, %.1f min", yr, nrow(out),
                    nrow(s_yr), cores, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  message(sprintf("Done in %.1f min. Outputs in %s", as.numeric(difftime(Sys.time(), t_all, units = "mins")), out_dir))

} else if (mode == "daymet") {
  yr <- as.integer(args[2]); var <- args[3]
  stopifnot(var %in% daymet_vars)
  t0 <- Sys.time()
  s_yr <- sites_for(read_sites(), read_site_years(), yr, "daymet")
  out <- extract_daymet(s_yr, yr, var)
  saveRDS(out, file.path(out_dir, sprintf("daymet_%s_%d.rds", var, yr)))
  message(sprintf("daymet %s %d: %d sites, %d rows, %.1f min, peak %.2f GB", var, yr, nrow(s_yr), nrow(out),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), peak_mem_gb()))

} else if (mode == "cdl") {
  yr <- as.integer(args[2])
  t0 <- Sys.time()
  s_yr <- sites_for(read_sites(), read_site_years(), yr, "cdl")
  out <- run_cdl(s_yr, yr)
  saveRDS(out, file.path(out_dir, sprintf("cdl_%d.rds", yr)))
  message(sprintf("cdl %d: %d/%d sites returned metrics on %d cores, %.1f min, peak %.2f GB",
                  yr, nrow(out), nrow(s_yr), n_cores,
                  as.numeric(difftime(Sys.time(), t0, units = "mins")), peak_mem_gb()))

} else if (mode == "dryrun") {
  # 1. Inputs: every file the arrays will open
  need <- c(counts_path,
            as.vector(outer(daymet_vars, daymet_yrs, function(v, y) mapply(daymet_path, v, y))),
            sapply(yrs, cdl_path))
  miss <- need[!file.exists(need)]
  message(length(need) - length(miss), "/", length(need), " input files found")
  if (length(miss)) { message("MISSING:\n  ", paste(miss, collapse = "\n  ")); quit(status = 1) }

  counts_raw <- read.csv(counts_path) %>% filter(longitude >= -105)
  sites      <- make_site_table(counts_raw)
  site_years <- make_site_years(counts_raw, sites)
  n_dm  <- sapply(daymet_yrs, function(y) nrow(sites_for(sites, site_years, y, "daymet")))
  n_cdl <- sapply(yrs, function(y) nrow(sites_for(sites, site_years, y, "cdl")))
  message(nrow(sites), " sites in total (CEW_SITES=", site_mode, "). Sites per extraction year:")
  print(data.frame(year = daymet_yrs, daymet = n_dm, cdl = c(NA, n_cdl)), row.names = FALSE)

  # 2. One full Daymet unit (= exactly one array element), the year with most sites
  y_dm <- daymet_yrs[which.max(n_dm)]
  t0 <- Sys.time()
  dm <- extract_daymet(sites_for(sites, site_years, y_dm, "daymet"), y_dm, "tmax")
  t_dm <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("Daymet unit (%d, %d sites): %.1f min, %d NA of %d values",
                  y_dm, max(n_dm), t_dm, sum(is.na(dm$tmax)), nrow(dm)))

  # 3. A sample of the busiest CDL year's sites on one core, scaled to that year / cores
  y_cdl <- yrs[which.max(n_cdl)]
  s_cdl <- sites_for(sites, site_years, y_cdl, "cdl")
  n_try <- min(nrow(s_cdl), if (length(args) > 1) as.integer(args[2]) else 20)
  t0 <- Sys.time()
  cd <- run_cdl(s_cdl[round(seq(1, nrow(s_cdl), length.out = n_try)), ], y_cdl, cores = 1)
  t_cdl <- as.numeric(difftime(Sys.time(), t0, units = "mins")) / n_try
  mem <- peak_mem_gb()
  message(sprintf("CDL: %.2f min per site, %d/%d sites returned metrics", t_cdl, nrow(cd), n_try))
  print(cd[, 1:min(6, ncol(cd))])

  # 4. Suggested sbatch resources: 1.5x measured time + 10 min, 1.5x peak memory, min 2 GB
  pad_w   <- function(m) { m <- ceiling(1.5 * m + 10); sprintf("%02d:%02d:00", m %/% 60, m %% 60) }
  pad_mem <- function(g, k = 1) if (is.na(g)) "<1.5 x seff 'Memory Utilized'>" else
    sprintf("%dG", max(2, ceiling(1.5 * g * k)))
  cdl_cores <- 4
  message("\n--- Suggested settings for hpc/submit_covariates.sh ---")
  message("Peak memory this job (Daymet + CDL sample): ", round(mem, 2), " GB")
  message("DAYMET_TIME=", pad_w(t_dm), "   DAYMET_MEM=", pad_mem(mem))
  message("CDL_CPUS=", cdl_cores, "   CDL_TIME=", pad_w(t_cdl * max(n_cdl) / cdl_cores),
          "   CDL_MEM=", pad_mem(mem, cdl_cores), "   # whole job: each forked worker holds a window")
  message("Array sizes: Daymet ", length(daymet_yrs) * length(daymet_vars), " jobs, CDL ", length(yrs), " jobs")

} else {
  stop("Unknown mode: ", mode)
}
