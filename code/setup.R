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
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"),
        renv.config.auto.snapshot = FALSE)

# an anchor file so here::here() resolves to this directory (the project root)
if (!file.exists(".here")) file.create(".here")

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing renv ...")
  install.packages("renv")
}

pkgs <- c(
  "arrow", "cmprsk", "collapse", "comorbidity", "data.table", "dplyr", "fst",
  "ggplot2", "here", "lme4", "lubridate", "mclogit", "patchwork", "ps", "rlang",
  "rprojroot", "scales", "stringr", "survival", "yaml"
)

if (file.exists("renv.lock")) {

  # Report the R version the lockfile targets vs the one running here. A mismatch
  # is the usual reason restore() struggles: CRAN may have no BINARY of a pinned
  # package version for a newer R, so renv tries to compile from source and can
  # fail. We therefore attempt restore(), then FILL ANY GAPS with versions that
  # are compatible with the R actually running -- so a newer (or older) R still
  # ends up with a working library. The shared renv.lock is left untouched.
  lock_r <- tryCatch(renv::lockfile_read("renv.lock")$R$Version, error = function(e) NA_character_)
  if (!is.na(lock_r) && lock_r != as.character(getRversion()))
    message(sprintf("Note: renv.lock targets R %s; you are running R %s. Any packages\n",
                    lock_r, getRversion()),
            "  whose pinned build is unavailable for your R will be replaced with a\n",
            "  compatible version (results are not expected to change).")

  message("Restoring recorded package versions from renv.lock ...")
  tryCatch(renv::restore(prompt = FALSE),
           error = function(e) message("  restore() reported: ", conditionMessage(e)))

  # anything the lock could not supply for THIS R -> install a compatible version
  gaps <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(gaps)) {
    message("Filling ", length(gaps), " package(s) the lockfile could not supply on this R:\n  ",
            paste(gaps, collapse = ", "))
    renv::install(gaps)     # newest version that builds on the current R
  }

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
message("\nSetup complete. Run the pipeline with:  Rscript code/run_pipeline.R")
