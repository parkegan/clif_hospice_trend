#!/usr/bin/env Rscript
# =============================================================================
# setup.R  —  one-time environment setup with renv.
# -----------------------------------------------------------------------------
# Run this ONCE per machine, from the project root:
#     Rscript setup.R
#
# Behaviour:
#   * If renv.lock EXISTS (a collaborator committed it), restore those exact
#     package versions — this is the reproducible path every site should use.
#   * If renv.lock is ABSENT (you are the first to set the project up), initialise
#     renv, install the dependency set, and write renv.lock for everyone else.
#
# After this finishes, run the pipeline with:
#     Rscript code/run_pipeline.R
#
# Updated 2026.10.16 (UMN 261015 postmortem):
#   * Added a mixed model FIT test at the end. The loadability check below
#     (requireNamespace) proves a namespace LOADS, not that it works: an lme4
#     binary paired with a Matrix it was not compiled against loads cleanly
#     and then errors on every glmer() call ("function 'cholmod_factor_ldetA'
#     not provided by package 'Matrix'" or similar). That is exactly the
#     failure mode a fresh renv library can produce on Linux when BINARY
#     packages are restored at locked versions (prebuilt binaries only exist
#     for current versions, so a restored older Matrix can meet an lme4 built
#     against a newer one), and it silently degrades every mixed model in 04
#     to fixed effects. The test below fits a toy GLMM and prints the repair
#     (reinstall Matrix and lme4 from SOURCE together, restart, snapshot).
#
# Updated 2026.10.09 (buddy site feedback, restricted install environments):
#   * The package list is now READ FROM code/_dependencies.R at run time and
#     unioned with the fallback list below, so the two manifests can never
#     drift apart. Sites building environments by hand (conda, offline
#     mirrors) can treat either file as the complete list.
#   * Added a system time zone database check. tzdata is an OS level
#     dependency, not an R package: minimal containers and bare conda
#     environments can lack it, and every datetime operation with a named
#     zone (America/Chicago, America/Los_Angeles, ...) then fails. R cannot
#     install it; the check below tells you what to install if it is missing
#     (conda: `conda install tzdata`; Debian/Ubuntu: `apt-get install tzdata`;
#     RHEL/Rocky: `dnf install tzdata`).
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"),
        renv.config.auto.snapshot = FALSE)

# an anchor file so here::here() resolves to this directory (the project root)
if (!file.exists(".here")) file.create(".here")

# --- system time zone database check (tzdata; see header note) ---------------
tz_ok <- tryCatch({
  tzs <- OlsonNames()
  length(tzs) > 0 && "America/Chicago" %in% tzs
}, error = function(e) FALSE)
if (!tz_ok) {
  message("WARNING: the system time zone database (tzdata) is missing or ",
          "incomplete.\n  Datetime handling with named zones will fail at ",
          "runtime. Install it at the SYSTEM level before running the ",
          "pipeline:\n    conda:         conda install tzdata\n",
          "    Debian/Ubuntu: apt-get install tzdata\n",
          "    RHEL/Rocky:    dnf install tzdata")
}

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv ...")
  install.packages("renv")
}

# Fallback package list. The authoritative manifest is code/_dependencies.R;
# the block below scans it and unions the result with this list, so adding a
# package there is enough — this vector only matters if that file is absent.
pkgs <- c(
  "arrow", "cmprsk", "collapse", "comorbidity", "data.table", "dplyr", "fst",
  "ggplot2", "here", "lme4", "lubridate", "mclogit", "patchwork", "ps", "rlang",
  "rprojroot", "scales", "stringr", "survival", "yaml"
)

dep_file <- file.path("code", "_dependencies.R")
if (file.exists(dep_file)) {
  dep_lines <- readLines(dep_file, warn = FALSE)
  hits      <- regmatches(dep_lines,
                          regexpr("^\\s*library\\([A-Za-z0-9._]+\\)", dep_lines))
  scanned   <- gsub("^\\s*library\\(|\\)$", "", hits)
  extra     <- setdiff(scanned, pkgs)
  if (length(extra))
    message("Packages found in code/_dependencies.R beyond the fallback list: ",
            paste(extra, collapse = ", "))
  pkgs <- sort(union(pkgs, scanned))
} else {
  message("NOTE: code/_dependencies.R not found — using the fallback package ",
          "list in this script. Run setup.R from the project root.")
}
message("Package set (", length(pkgs), "): ", paste(pkgs, collapse = ", "))

if (file.exists("renv.lock")) {
  message("renv.lock found -> restoring recorded package versions ...")
  renv::restore(prompt = FALSE)
} else {
  message("No renv.lock -> initialising a fresh renv library ...")
  renv::init(bare = TRUE, restart = FALSE)
  message("Installing ", length(pkgs), " packages (this can take a while the first time) ...")
  renv::install(pkgs)
  renv::snapshot(prompt = FALSE)
  message("Wrote renv.lock — commit it so every site installs identical versions.")
}

# sanity check: are all packages loadable?
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("\nWARNING: these packages did not install cleanly: ",
          paste(missing, collapse = ", "),
          "\nOn Linux they usually need system libraries — see the README ",
          "(libcurl, libssl, libxml2, and Arrow's dependencies).")
} else {
  message("\nAll ", length(pkgs), " packages available.")
}

# --- mixed model FIT test (2026.10.16; see header note) -----------------------
# Loadable is not workable: fit one toy GLMM so an lme4/Matrix binary mismatch
# is caught HERE, during setup, instead of silently degrading every mixed model
# in 04 to fixed effects during a multi hour site run.
glmm_ok <- tryCatch({
  set.seed(1)
  toy <- data.frame(y = rbinom(400L, 1L, 0.3),
                    g = rep(sprintf("h%02d", 1:8), each = 50L),
                    x = rnorm(400L))
  invisible(suppressWarnings(
    lme4::glmer(y ~ x + (1 | g), data = toy,
                family = binomial(link = "logit"), nAGQ = 1)))
  TRUE
}, error = function(e) conditionMessage(e))

if (isTRUE(glmm_ok)) {
  message("Mixed model fit test passed (lme4::glmer works in this library).")
} else {
  message("\nWARNING: lme4 loads but CANNOT FIT a mixed model in this library:\n    ",
          glmm_ok,
          "\n  This is usually a Matrix / lme4 binary mismatch (common after a ",
          "renv restore\n  with prebuilt Linux binaries, or after any update of ",
          "Matrix alone). Repair it\n  BEFORE running the pipeline — 04 will ",
          "refuse to run until this passes:\n",
          "    install.packages(c(\"Matrix\", \"lme4\"), type = \"source\")\n",
          "  then RESTART R, re-run this setup, and renv::snapshot() the fix.")
}

message("\nSetup complete. Run the pipeline with:  Rscript code/run_pipeline.R")
