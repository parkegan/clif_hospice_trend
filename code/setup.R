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
message("\nSetup complete. Run the pipeline with:  Rscript code/run_pipeline.R")
